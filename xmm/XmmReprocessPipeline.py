import subprocess
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.esa.xmm_newton import XMMNewton
import os, re
import shutil, tarfile
from pathlib import Path
import pandas as pd
from glob import glob
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from regions import CircleSkyRegion, CircleAnnulusSkyRegion
from matplotlib.patches import Circle
from matplotlib.colors import LogNorm
from astropy.wcs.utils import proj_plane_pixel_scales
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[0]

class XMMCrossMatchPipeline:
    """
        This class downloads xmm data, and reprocess the raw data: generate science files, source detection, generate spectra etc.
    """

    def __init__(self,
                 project_root="/Users/juliaahlvind/Documents/projekt_3/"):      

        # Paths
        self.project_root = project_root

        # Pipeline och data directories relativt project_root
        self.pipeline_dir = ROOT
        self.data_directory = os.path.join(project_root, "data", "xmm")
        self.sample_dir =  os.path.join(project_root, "sample")
        self.input_file = os.path.join(project_root, "sample", "test_xmm.txt")
        #self.input_file = os.path.join(project_root, "sample", "cross_matches_optical_xmm_epoch_cuts_type_filtered.txt")


    @staticmethod
    def fk5_to_physical(fits_file:Path, ra:str, dec:str):
        """ Translate celestial coordinates to physical coordinates in image.
            This is necessary for eregionanalysis to run more efficiently. 
            
            output: SN coord in physical image coordinates (float)
            """

        coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg), frame='icrs')
        # Convert to degree decimal
        ra_deg = coord.ra.degree
        dec_deg = coord.dec.degree

        # read the image file
        hdul = fits.open(fits_file)
        hd = hdul[0].header
        
        physical_coord_ref_x = hd['REFXCRPX']  # phys
        wcs_coord_ref_deg_x = hd['REFXCRVL']  # deg
        wcs_pixel_size_x = hd['REFXCDLT']  # deg/pix

        physical_coord_ref_y = hd['REFYCRPX']  # phys
        wcs_coord_ref_deg_y = hd['REFYCRVL'] # deg
        wcs_pixel_size_y = hd['REFYCDLT']  # deg/pix

        # translate coordinates
        x = (ra_deg - wcs_coord_ref_deg_x)*np.cos(np.deg2rad(dec_deg))/wcs_pixel_size_x + physical_coord_ref_x
        y = (dec_deg - wcs_coord_ref_deg_y)/wcs_pixel_size_y + physical_coord_ref_y
        
        return x, y

    def physical_to_fk5(fits_file:Path, x, y):
        """ 
            """

        # read the image file
        hdul = fits.open(fits_file)
        hd = hdul[0].header
        
        physical_coord_ref_x = hd['REFXCRPX']  # phys
        wcs_coord_ref_deg_x = hd['REFXCRVL']  # deg
        wcs_pixel_size_x = hd['REFXCDLT']  # deg/pix

        physical_coord_ref_y = hd['REFYCRPX']  # phys
        wcs_coord_ref_deg_y = hd['REFYCRVL'] # deg
        wcs_pixel_size_y = hd['REFYCDLT']  # deg/pix

        # translate coordinates
        dec_deg = (y - physical_coord_ref_y)*wcs_pixel_size_y + wcs_coord_ref_deg_y
        
        ra_deg = (x - physical_coord_ref_x)*wcs_pixel_size_x/np.cos(np.deg2rad(dec_deg)) + wcs_coord_ref_deg_x

        return ra_deg, dec_deg

    @staticmethod
    def read_eregion_src_region_output(eregion_file: str):
        """
        Read a DS9 FK5 circle region file and extract center coordinates and radius.

        Returns
        -------
        center : astropy.coordinates.SkyCoord
            
            Sky coordinate of region center (FK5/J2000)
        radius : astropy.units.Quantity
            Radius in arcseconds
        """
        # Läs hela filen som text
        with open(eregion_file, "r") as f:
            text = f.read()

        # Regex för att hitta SASCIRCLE: (X,Y) in CIRCLE(X,Y,R)
        pattern = r"SASCIRCLE:\s*\(X,Y\)\s*in\s*CIRCLE\s*\(\s*([0-9.eE+-]+)\s*,\s*([0-9.eE+-]+)\s*,\s*([0-9.eE+-]+)\s*\)"


        match = re.search(pattern, text)
        if match:
            x, y, r = match.groups()
            print(f"CIRCLE coordinates: X={x}, Y={y}, R={r}")
            return x, y ,r
        else:
            print("Ingen CIRCLE hittades")

    @staticmethod
    def read_ds9_circle_region(reg_file: str):
        """
        Read a DS9 FK5 circle region file and extract center coordinates and radius.

        Returns
        -------
        center : astropy.coordinates.SkyCoord
            
            Sky coordinate of region center (FK5/J2000)
        radius : astropy.units.Quantity
            Radius in arcseconds
        """
        with open(reg_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith("circle"):
                    match = re.search(r"circle\(([^,]+),([^,]+),([^)]+)\)", line)
                    if not match:
                        raise ValueError("Could not parse circle region")

                    ra_str, dec_str, radius_str = match.groups()

                    return ra_str, dec_str, radius_str

        raise ValueError("No circle region found")

    def extract_tar_file(self, tar_file: str, path: Path):
        """Extract a tar-fil to a specific place."""

        print(f'path:{path}')
        print(f'file: {tar_file}')
        os.chdir(path)
        with tarfile.open(tar_file, 'r') as tar:
            tar.extractall()  # Se till att filerna extraheras till rätt plats
        

    def move_file(self, file_to_move: Path, new_path: Path):
        """ moves files from one place to another. """
        shutil.move(file_to_move,new_path)
        

    def find_file(self, path: Path, name_key: str) -> str:
        """ Find file at path that ends with name_key"""
        files = list(Path(path).glob(f'*{name_key}'))
        print(files)
        return str(files[0])

    # --- Step 1: Download raw data ---
    def download_xmm_data(self, obsid: str):
        """ Downloads the ODF files from webpage """

        print(f"Downloading data for: {obsid}")

        os.mkdir(f'{self.data_directory}/{obsid}')

        # go to folder 
        os.chdir(f'{self.data_directory}/{obsid}')
        
        subprocess.run(
            #f" curl -o files.tar https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?{obsid}&level=ODF",
            f' curl -o files.tar "https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno={obsid}&instname=PN&level=ODF"',
            shell=True, # Python runs the command in shell-process eg bash
            check=True # Python flags if bash command fails
        )
    
        os.chdir(self.pipeline_dir)

        print("done downloading data!")

    # --- Step 2: Restructure files within folder ---
    def xmm_data_is_fantastic_so_this_is_what_i_have_to_do(self, obsid: str):
        """ This is a merge of random things like moving and untaring files due to the wonderful xmm structure. """
        print(f"Moving around files and un tar for {obsid}")


        # Opens the tar file: this creates a folder with the name {obsid} which has two dolders within: odf and pps
        try:
            self.extract_tar_file(f'files.tar', f"{self.data_directory}/{obsid}")

        except Exception as e:
            print(f"Could not open tar file. Error: {e}")

            # save failed obsid to output file
            ff = open(f"{self.data_directory}/failed_downloads.txt", 'a')
            ff.write(f"{obsid}\n")
            ff.close()
            print("Continue with next obsid")
            os.chdir(self.pipeline_dir)

            return

        # Remove the tar file that was just opened
        os.remove(f'{self.data_directory}/{obsid}/files.tar')

        # move files from subfolder to parent folder for the reprocessing to work
        ##self.move_file(f'{self.data_directory}/{obsid}/odf/{obsid}.tar.gz', f'{self.data_directory}/{obsid}/{obsid}.tar.gz')

        # open the new tar file that just got moved
        ##self.extract_tar_file(f'{self.data_directory}/{obsid}/{obsid}.tar.gz', f'{self.data_directory}/{obsid}/')

        # remove the tar file that just got opened
        ##os.remove(f'{self.data_directory}/{obsid}/{obsid}.tar.gz')

        # identify the new TAR file that got created
        file = self.find_file(f'{self.data_directory}/{obsid}/', '.TAR')

        # open the TAR file
        self.extract_tar_file(file,  os.path.dirname(file))

        # remove the tar file
        os.remove(f'{file}')

        # remove the empthy odf folder and pps(this is not needed for data processing)
        ##shutil.rmtree(f'{self.data_directory}/{obsid}/odf')
        ##shutil.rmtree(f'{self.data_directory}/{obsid}/pps')

        os.chdir(self.pipeline_dir)

        print("done moving files!")

    # --- Step 3: Pre process raw data to generate event files etc. ---
    def reprocess_xmm_EPIC_data(self, obsid: str, sn_name: str, preprocess: bool):
        """ runs the bash SASReduction: first reprocessing. Only reprocessing EPIC-pn data 
        
            outputs: 1. science files: image(obs_id_pn_img_filt_gti_05_10.fits), eventf(obs_id_pn_filt_gti_05_10.fits)
                     2. spectral files: name_obsid_pn_spec_grp_bin1.fits (name_obsid_pn_spec_bkg_bin1.fits, name_obsid_pn_rmf_bin1.fits, name_obsid_pn_arf_bin1.fits)
                     3. ewavelet files: source_list05_10_pn_3σ_name_obs_id.fits
                     """
        print(f"Reprocessing first step for {sn_name} {obsid}...")
        subprocess.run(
            f"{self.pipeline_dir}/SASReduction.sh {obsid} {sn_name} {preprocess}",
            shell=True, # Python runs the command in shell-process eg bash
            check=True # Python flags if bash command fails
        )

        os.chdir(self.pipeline_dir)
        print("Done reprocessing first steps")

    # --- Step 4: Create source region file by invoking eregionanalysis ---
    def generte_src_region_file(self, sn_name:str, obsid:str, ra:str, dec:str):
        """ Create src region file by providing eregionanalys with physical 
            coordinates (image-physical coord) for the SN position.

            The default src region size is 10".

            output:
            src_sn_name_obsid.reg
            large_src_flagged_src.txt
            """
        
        # path fits file
        fits_path = glob(f"{self.data_directory}/{obsid}/{obsid}_pn_img_filt_gti_05_10.fits")[0]

        # caclulate physical coordinates
        x_physical, y_physical = self.fk5_to_physical(fits_path, ra, dec)

        self.call_eregionanalyse(obsid, sn_name, x_physical, y_physical)

        # --- create source region based on eregion output ---
        # NOTE!! Flag if src reg > 200/15" could be that the SN is positioned near a brighter source
        # these are flagged and saved to a txt file to check later!
        eregion_file = glob(f"{self.data_directory}/{obsid}/{sn_name}_{obsid}_eregion_pn_05_10.txt")[0]
        x_physical_eregion, y_physical_eregion, radi_physical_eregion = self.read_eregion_src_region_output(eregion_file)

        # print src region file
        ff = open(f"{self.data_directory}/{obsid}/src_{sn_name}_{obsid}.reg", 'w+')
        ff.write(f"circle({x_physical_eregion},{y_physical_eregion},{radi_physical_eregion})")
        ff.close()

        print(float(radi_physical_eregion))
        # save flag if large region
        if float(radi_physical_eregion) > 300: # > 15"
            ff = open(f"{self.sample_dir}/large_src_flagged_src.txt", 'a')
            ff.write(f"{sn_name}\t{obsid}\t{radi_physical_eregion}\n")
            ff.close()


    # --- Call eregionanalyse ---
    def call_eregionanalyse(self, obsid: str, sn_name: str, x_physical:float, y_physical:float):
        """ runs the SAS task eregionanalyse to extract aoutomatic source regions based on src and bkg regions 
        
            output:
            sn_name_obsid_eregion_pn_05_10.txt
            """
            
        subprocess.run(
            f"{self.pipeline_dir}/eregionanalyse.sh {obsid} {sn_name} {x_physical} {y_physical} ",
            shell=True, # Python runs the command in shell-process eg bash
            check=True # Python flags if bash command fails
        )
        
    # --- Step 6: Generate background region file, annulus ---
    def generate_bkg_region(self, sn_name:str, obsid:str):
        """ Generate a background region file.
            shape annulus: iner radius = src+0.5 arcsec
                           outer radius = innder_rad*2 
                           
        output:
        bkg_annulus_sn_name_obsid.reg"
        
        """

        print(f"Creating bkg annulus file for {sn_name} {obsid}")
        # source region
        src_file=f"{self.data_directory}/{obsid}/src_{sn_name}_{obsid}.reg"
        x_physical_src, y_physical_src, radius_str = self.read_ds9_circle_region(src_file)

        radius_float = float(radius_str)
        inner_rad = f'{radius_float+10}'
        outer_rad = f'{(radius_float+10)*2}'

        ff = open(f"{self.data_directory}/{obsid}/bkg_annulus_{sn_name}_{obsid}.reg", '+w')
        ff.write(f"annulus({x_physical_src},{y_physical_src},{inner_rad},{outer_rad})\n")

        ff.close()

        print("finnished creating bkg file!")
        print(" ")


    def open_image_with_ds9(self, sn_name: str, obsid: str, ra:str, dec:str):

        image = f"{self.data_directory}/{obsid}/{obsid}_pn_img_filt_gti_05_10.fits"
        src_reg_file = f"{self.data_directory}/{obsid}/src_{sn_name}_{obsid}.reg"
        
        print(f"Pleas check src region for {sn_name} {obsid}...")
        print(f"ra:{ra}, dec:{dec}")

        # run ds9
        subprocess.run(
            f". ds9 {image} -region {src_reg_file} -scale log ",
            shell=True, # Python runs the command in shell-process eg bash
            check=True # Python flags if bash command fails
        )
        
        # Wait for the user to interact with DS9
        input(f"Press Enter after inspecting the image in DS9...{sn_name}, {obsid}")
        
        # close the ds9 window
        #subprocess.run('pkill ds9', shell=True, check=True)
        



    def plot_save_image_with_src_bkg_reg(self, sn_name:str, obsid:str):
        """ This function plots the image with source and background region and 
            saves the figure output without appearing on screen """

        src_file = f"{self.data_directory}/{obsid}/src_{sn_name}_{obsid}.reg"
        x_physical, y_physical, r_phys = self.read_ds9_circle_region(src_file)

        x_physical= float(x_physical)
        y_physical= float(y_physical)
        r_phys = float(r_phys)
        # open image
        image_file = f"{self.data_directory}/{obsid}/{obsid}_pn_img_filt_gti_05_10.fits"
        with fits.open(image_file) as hdul:
            data = hdul[0].data
            header = hdul[0].header

            wcs_phys = WCS(header, key='L')

        x_img, y_img = wcs_phys.wcs_world2pix(x_physical, y_physical,0)
        x2_img, _ = wcs_phys.wcs_world2pix(x_physical + r_phys, y_physical, 0)
        r1_img = abs(x2_img - x_img)

        # plot image
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=wcs_phys)
        ax.imshow(data, origin='lower', cmap='viridis', norm=LogNorm())

        # plot src region
        circle = Circle((x_img, y_img), r1_img, edgecolor='green', facecolor='none', lw=1, label='src reg')
        ax.add_patch(circle)

        # retrive bkg data
        outer_rad = (r_phys+10)*2
        x_img, y_img = wcs_phys.wcs_world2pix(x_physical, y_physical,0)
        x2_img, _ = wcs_phys.wcs_world2pix(x_physical + outer_rad, y_physical, 0)
        r2_img = abs(x2_img - x_img)

        # plot bkg region
        circle = Circle((x_img, y_img), r2_img, edgecolor='red', facecolor='none', lw=1, label='bkg reg')
        ax.add_patch(circle)

        # zoom: 10 x outer radius
        zoom_rad_pix = (10.0 * r2_img) 

        # zoom
        ax.set_xlim(x_img - zoom_rad_pix, x_img + zoom_rad_pix)
        ax.set_ylim(y_img - zoom_rad_pix, y_img + zoom_rad_pix)

        plt.title(f'{sn_name} {obsid}: inner={round(r1_img,2)}, outer={round(r2_img,2)}"')

        plt.savefig(f"{self.data_directory}/ds9_saves/{sn_name}_{obsid}.png")
        plt.clf()


    def process_input_file(self):
        """Read input file containing obsids to download and reprocess."""

        # load Chandra data
        xmm = os.path.join(self.input_file)
        cols_with_epoch_xmm = ["index", "name", "type", "dist_Mpc", "ra", "dec", "optical_date",
                                "obsid", "sep_arcmin", "expt_ks", "X_ray_obs_date", "epoch_days",
                                "PIname", "target"]
        df_xmm = pd.read_csv(xmm, sep="\t", skiprows=1, names=cols_with_epoch_xmm, dtype={"obsid": str})


        for i,(sn_name,obsid) in enumerate(zip(df_xmm['name'],df_xmm['obsid'])):
            ra, dec = df_xmm['ra'][i], df_xmm['dec'][i]

            print(f"Running for index:{df_xmm['index'][i]}")
            # if folder already downloaded, skip this step
            obsid_folder = os.path.join(self.project_root, "data", "xmm")
            if obsid not in os.listdir(obsid_folder): # proceed

                # 1. Download data
                self.download_xmm_data(obsid)
                
                # 2. move and open all files
                self.xmm_data_is_fantastic_so_this_is_what_i_have_to_do(obsid)
                    
                # 3. reprocess the EPIC-pn/mos data
                self.reprocess_xmm_EPIC_data(obsid, sn_name, preprocess=True)

            #4. Create source regions
            self.generte_src_region_file(sn_name, obsid, ra, dec)

            # 5. STOP!!
            # Check the sus src regions in large_src_flagged_src and adjust accordingly 
            # before moving on to the next step!
            #with open(f"{self.sample_dir}/large_src_flagged_src.txt", "r") as ff:
            #    text = ff.read()

            #pattern = f"{sn_name}\t{obsid}"
            #match = re.search(pattern, text)
            #if match:
            self.open_image_with_ds9(sn_name, obsid, ra, dec)

            input_continue = input("Input 'y' to continue with the next steps for this source, or 'n' to skip to the next source: ")

            if input_continue.lower() == 'y':

                # 6. Create background region annulus
                self.generate_bkg_region(sn_name, obsid)

                # plot ds9 image with src and bkg regions witout vewing it live -> store image
                self.plot_save_image_with_src_bkg_reg(sn_name,obsid)

                # 7. Generate spectra
                self.reprocess_xmm_EPIC_data(obsid, sn_name, preprocess=False)

            # Optional. Run edetect for image source detection
            #self.reprocess_xmm_EPIC_data(obsid, sn_name, preprocess='false', run_ewavelet=True)







pipeline = XMMCrossMatchPipeline()
pipeline.process_input_file()

            
os.chdir("/Users/juliaahlvind/Documents/projekt_3/pipelines")



