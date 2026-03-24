import subprocess
import os, time, re
from astropy.io import ascii, fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
import pandas as pd
from regions import CircleSkyRegion, CircleAnnulusSkyRegion
from matplotlib.colors import LogNorm
from astropy.wcs.utils import proj_plane_pixel_scales
from pathlib import Path

ROOT = Path(__file__).resolve().parents[0]


class ChandraDataProcessor:
    """A class to handle chandra download of raw data, reprocessing and extract source region."""

    def __init__(self, 
                project_root: str = "/Users/juliaahlvind/Documents/projekt_3/"):

        self.project_root = project_root

        # Pipeline och data directories relativt project_root
        self.pipeline_dir = ROOT
        self.data_directory = os.path.join(project_root, "data", "chandra")
        self.sample_dir = os.path.join(project_root, "sample")
        #self.input_file = os.path.join(project_root, "sample", "test_chandra.txt")
        self.input_file = os.path.join(project_root, "sample", "cross_matches_optical_chandra_epoch_cuts_type_filtered_data_quality_cuts.txt")

    @staticmethod
    def fk5_to_physical(fits_file:Path, ra:str, dec:str, r_arcsec:str):
        """ Translate celestial coordinates to physical coordinates in image.
            This is necessary to more efficiently create a smaller image centred on src reg to run wavdetect on. 
            
            ra: src centre in hour angle
            dec: src centre in deg
            r_arcsec: src redius in arcsec
            output: SN coord in physical image coordinates (float: (x,y))
            """

        coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg), frame='icrs')
        # Convert to degree decimal
        ra_deg = coord.ra.degree
        dec_deg = coord.dec.degree

        # read the image file
        hdul = fits.open(fits_file)
        hd = hdul[0].header

        crpix_x = hd['CRPIX1'] # ref pixel for x (physical)
        crval_x = hd['CRVAL1'] # ref ra (deg)
        cdelt_x = hd['CDELT1'] # deg/pixel in x

        crpix_y = hd['CRPIX2'] # ref pixel for y (physical)
        crval_y = hd['CRVAL2'] # ref dec (deg)
        cdelt_y = hd['CDELT2'] # deg/pixel in y

        # translate coordinates
        x = (ra_deg - crval_x) * np.cos(np.deg2rad(dec_deg)) / cdelt_x + crpix_x
        y = (dec_deg - crval_y) / cdelt_y + crpix_y
        
        arcsec_per_pixel = abs(cdelt_x) * 3600.0
        rad_pixl = float(r_arcsec) / arcsec_per_pixel

        return x, y, rad_pixl

    # --- Step 1: Downlaod data ---
    def download_obsid(self, obsid: str):
        """ Download raw chandra data for specific obsID. """
        
        print(f"Running CIAO tool download_chandra_data for {obsid}...")
        
        if obsid not in os.listdir(self.data_directory):
            os.chdir(self.data_directory)
            subprocess.run(
                f'conda run -n ciao download_chandra_obsid {obsid}',
                shell=True,
                check=True
            )
        print("finnished downloading Chandra data!")
        print(" ")


    # --- Step 2: Reprocess the raw chandra data ---
    def process_obsid(self, obsid: str):
        """ Reprocess the data. If folder repro excists, reprecessing already done
        
            output:
            images (image_058_bin1.fits)
            event files (acisfobsid_repro_evt2_03_10_clean.fits)
            """
        
        print(f"Reprocessing the data of {obsid}...")

        if "repro" not in os.listdir(f"{self.data_directory}/{obsid}"):

            full_obsid = f"{int(obsid):0>5}"
            print(full_obsid)

            subprocess.run(
                f"{self.pipeline_dir}/repro.sh {obsid} {full_obsid}",
                shell=True,
                check=True
            )
        print("finnished reprocessing data!")
        print(" ")

  

    # --- Optional step: Create small image around SN position ---
    def run_create_small_image(self, sn_name:str, obsid:str):
        """ Create a small images around the SN position to run wavdetect more quickly. 
            
            output:
            small_image_058_bin1_sn_name.fits
            """

        print(f"Creating a smaller image for: {sn_name} {obsid}:")

        src_file = f"{self.data_directory}/{obsid}/repro/psfsize_src_{sn_name}_{obsid}.reg"
        ra_centre, dec_centre, radius_str = self.read_ds9_circle_region(src_file)
        radius_value = radius_str.replace('"', '')

        fits_file_image = f"{self.data_directory}/{obsid}/repro/image_058_bin1.fits"
        xx, yy, rad_pixl = self.fk5_to_physical(fits_file_image,ra_centre, dec_centre, radius_value)

        subprocess.run(
            f"{self.pipeline_dir}/create_small_images.sh {sn_name} {obsid} {xx} {yy} {rad_pixl*20}", ## OBS!! small image =20 gr radie av source reg
            shell=True,
            check=True
        )    
        print(f"finnished creating small image!")
        print(" ")

    # --- Optional step: Run wavdetect ---
    def run_wavdetect(self, sn_name: str, obsid: str):
        """Run wavdetect for a specific obsid.
        
            outfile:
            sources_058_wavdetect_obsid.fits
            """

        print(f"Running wavdetect for: {sn_name} {obsid}:")
        subprocess.run(
            f'{self.pipeline_dir}/wavdetect.sh {sn_name} {obsid}',
            shell=True,
            check=True
        )    
        print(f"finnished running wavetect!")
        print(" ")



    def plot_wavdetect_results(self, sn_name:str, obsid:str):
        # FITS file with counts from wavdetect of found
        fits_file = f'{self.data_directory}/{obsid}/repro/sources_058_wavdetect_{sn_name}.fits'#outfile_{name}.fits'

        hdulist = fits.open(fits_file)
        data = hdulist[1].data


        # region file to check if it is indeed one of the detected
        src_file = f"{self.data_directory}/{obsid}/repro/psfsize_src_{sn_name}_{obsid}.reg"
        ra_centre, dec_centre, radius_str = self.read_ds9_circle_region(src_file)

        # Parse the region string to get the region coordinates
        region_coords = SkyCoord(ra=ra_centre, dec=dec_centre, frame='fk5', unit=(u.hourangle, u.deg)) 
        RA, DEC = region_coords.ra.deg, region_coords.dec.deg 

        # list all detected sources coordinates
        sources_coords = SkyCoord(data['RA'],data['DEC'], frame='fk5', unit=(u.deg, u.deg))

        if len(sources_coords)==0:
            return float('NaN'),float('NaN'),float('NaN'),float('NaN')
        else:
            # find if there is a match with coordinates of detected sources and src_reg, by a separation limit
            index_best, separation_best, _ = region_coords.match_to_catalog_sky(sources_coords)
            closest_det_src_cord_RA, closest_det_src_cord_DEC = sources_coords[index_best].ra.deg, sources_coords[index_best].dec.deg
            sep = separation_best.deg[0]

            fig, ax = plt.subplots(1,1)
            fig.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.92)
            ax.set_xlabel('RA (degrees)')
            ax.set_ylabel('DEC (degrees)')
            p = ax.scatter(data['RA'], data['DEC'],s=data['SRC_SIGNIFICANCE']*10, c=data['SRC_SIGNIFICANCE'], cmap='viridis', marker='.')
            p.set_clip_path(ax.patch)
            for i in range(len(data['RA'])):
                t = ax.text(data['RA'][i], data['DEC'][i], s=data['SRC_SIGNIFICANCE'][i])
                t.set_clip_on(True)

            ax.plot(float(RA),float(DEC),'+',color='red', markersize=8, label='SN pos.')
            circle1 = plt.Circle((float(RA),float(DEC)),10/3600,color='red', fill=False, label='10" radi ')
            ax.add_patch(circle1)

        

            ax.set_title(f'{sn_name}, {obsid}, min sep. {round(sep*3600,3)}", 3σ')

            ax.plot([], [], ' ', label='FOV = 20" × 20"')
            ax.legend()
            ax.set_xlim([float(RA)-20/3600, float(RA)+20/3600])
            ax.set_ylim([float(DEC)-20/3600, float(DEC)+20/3600])

            #plt.show()
            plt.savefig(f"{self.data_directory}/wavdetect_plots/{sn_name}_{obsid}.png")
            plt.clf()

            min_sep_arcsec = round(sep*3600,3)
            min_sep_ra, min_sep_dec = closest_det_src_cord_RA, closest_det_src_cord_DEC
            min_sep_significance = data['SRC_SIGNIFICANCE'][index_best]
            
            return min_sep_arcsec, min_sep_ra, min_sep_dec, min_sep_significance

    # --- Step 3 : Create src region via psfsize ---
    def extract_source_region(self, sn_name: str, obsid: str, ra:str, dec:str):
        """ Extract source region for specified obsID. 
        
            output:
            psfsize_src_sn_name_obsid.reg
            """

        print(f"Running psfsize for {sn_name} {obsid}...")
        subprocess.run(
            f'{self.pipeline_dir}/psfsize_srcs.sh {sn_name} {obsid} {ra} {dec}',
            shell=True,
            check=True
        )
        print("finnished running psfsize!")
        print(" ")



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

    # --- Step 4: Create background regions ---
    def generate_bkg_region(self, sn_name:str, obsid:str):
        """ Generate a background region file.
            shape annulus: iner radius = psfsize+0.5 arcsec
                           outer radius = innder_rad*2 
                           
        output:
        bkg_annulus_sn_name_obsid.reg
        """

        print(f"Creating bkg annulus file for {sn_name} {obsid}")
        # source region
        src_file = f"{self.data_directory}/{obsid}/repro/psfsize_src_{sn_name}_{obsid}.reg"
        ra_centre, dec_centre, radius_str = self.read_ds9_circle_region(src_file)
        radius_value = radius_str.replace('"', '')

        radius_arcsec = float(radius_value)
        inner_rad = f'{radius_arcsec+0.5}"'

        if (radius_arcsec+0.5)*2 < 5: #if the outer radius is smaller than 5, put deafult to 5 since a too small bkg region can have zero counts and thus no spectra can be extracted
            outer_rad = '5"'
        else:
            outer_rad = f'{(radius_arcsec+0.5)*2}"'

        ff = open(f"{self.data_directory}/{obsid}/repro/bkg_annulus_{sn_name}_{obsid}.reg", '+w')
        ff.write(f"annulus({ra_centre},{dec_centre},{inner_rad},{outer_rad})\n")

        ff.close()

        print("finnished creating bkg file!")
        print(" ")


    # --- Optional step: Open image in ds9 to inspect image ---
    def open_image_with_ds9(self, sn_name: str, obsid: str):

        image = f"{self.data_directory}/{obsid}/repro/image_058_bin1.fits"
        psf_reg_file = f"{self.data_directory}/{obsid}/repro/psfsize_src_{sn_name}_{obsid}.reg"
        bkg_reg_file = f"{self.data_directory}/{obsid}/repro/bkg_annulus_{sn_name}_{obsid}.reg"
        
        # run ds9
        subprocess.run(
            f". ds9 {image} -region {psf_reg_file} -region {bkg_reg_file} -scale log",
            shell=True, # Python runs the command in shell-process eg bash
            check=True # Python flags if bash command fails
        )
        
        # Wait for the user to interact with DS9
        input(f"Press Enter after inspecting the image in DS9...{sn_name}, {obsid}")
        
        # close the ds9 window
        #subprocess.run('pkill ds9', shell=True, check=True)
        


    # --- Step 5: Extract spectrum ---
    def spacextract_bash(self, index, sn_name: str, obsid: str):
        """ calls bash script spacextract.sh which creates grouped spectra
        
            output: 
                grouped spectr: spectra_sn_name_obsid_grp.pi 
                           bkg: spectra_sn_name_obsid_bkg.pi
                           rmf: spectra_sn_name_obsid.rmf
                           arf: spectra_sn_name_obsid.corr.arf
        
        """

        print(f"Creating spectra for {sn_name} {obsid}...")

        try:
            result = subprocess.run(
                f'{self.pipeline_dir}/specextract.sh {sn_name} {obsid}',
                check=True,
                shell=True,
            )
            print("finnished creating spectra!")
            print(" ")

        except subprocess.CalledProcessError as e:
            if e.returncode == 42:
                print("Could not create spectra due to zero counts in src reg.")
                ff = open(f"{self.data_directory}/zero_count_src_no_spectra.txt",'a')
                ff.write(f"{index}\t{sn_name}\t{obsid}\n")
                ff.close()
                return True


    # --- Step 6: Calculate flux limits with srcflux ---
    def srcflux_bash(self, sn_name:str, obsid: str):
        """ calls bash script srcflux.sh which calculates flux limits of given source region. 
        
            output:
            srcflux_output_05_8_sn_name_obsid.txt
        """
        
        print(f"Running srcflux: aka image source deteciton on {sn_name} {obsid}...")
        result = subprocess.run(
            f"{self.pipeline_dir}/srcflux.sh {sn_name} {obsid}",
            check=True,
            shell=True,
        )

        print("finnished running image source detection!")
        print(" ")


    def plot_save_image_with_src_bkg_reg(self, sn_name:str, obsid:str):
        """ This function plots the image with source and background region and 
            saves the figure output without appearing on screen """

        # open image
        image_file = f"{self.data_directory}/{obsid}/repro/image_058_bin1.fits"
        hdu = fits.open(image_file)[0]
        wcs = WCS(hdu.header)

        # plot image
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=wcs)
        ax.imshow(hdu.data, origin='lower', cmap='viridis', norm=LogNorm())

        # retrive src and bkg info
        src_file = f"{self.data_directory}/{obsid}/repro/psfsize_src_{sn_name}_{obsid}.reg"
        ra_centre, dec_centre, radius_str = self.read_ds9_circle_region(src_file)
        radius_value = float(radius_str.replace('"', ''))

        inner_rad = radius_value+0.5
        outer_rad = (radius_value+0.5)*2

        source_pos = SkyCoord(ra=ra_centre, dec=dec_centre, unit=(u.hourangle, u.deg), frame='icrs')

        # add src and bkg
        source_reg = CircleSkyRegion(center=source_pos, radius=radius_value*u.arcsec)
        bg_reg = CircleAnnulusSkyRegion(center=source_pos, inner_radius=inner_rad*u.arcsec, outer_radius=outer_rad*u.arcsec)

        bg_reg.to_pixel(wcs).plot(ax=ax, edgecolor='red', lw=1, label='Background')
        source_reg.to_pixel(wcs).plot(ax=ax, edgecolor='green', lw=1, label='Source')

        # centre pixel coord
        x_cen, y_cen = wcs.world_to_pixel(source_pos)

        # pixelscale (arcsec / pixel)
        pixscale = proj_plane_pixel_scales(wcs)[0] * 3600.0

        # zoom: 10 x outer radius
        zoom_rad_pix = (10.0 * outer_rad) / pixscale

        # zoom
        ax.set_xlim(x_cen - zoom_rad_pix, x_cen + zoom_rad_pix)
        ax.set_ylim(y_cen - zoom_rad_pix, y_cen + zoom_rad_pix)

        plt.title(f'{sn_name} {obsid}: inner={radius_str}, outer={outer_rad}"')
        plt.savefig(f"{self.data_directory}/ds9_saves/{sn_name}_{obsid}.png")
        plt.clf()


    def process_spectra_from_input_file(self):
        """Read input file containing obsids to download and reprocess."""

        # load Chandra data
        chandra = self.input_file
        cols_with_epoch_chandra = [
                    "index", "name", "type", "dist_Mpc", "ra", "dec", "optical_date",
                    "obsid", "sep_arcmin", "instr", "grating", "expt_ks",
                    "X_ray_obs_date", "epoch_days", "PIname", "target"
                ]
        df_chandra = pd.read_csv(chandra, sep="\t", skiprows=1, names=cols_with_epoch_chandra)


        for i,(sn_name,obsid) in enumerate(zip(df_chandra['name'],df_chandra['obsid'])):
            ra, dec = df_chandra['ra'][i], df_chandra['dec'][i]

            df_row = df_chandra.iloc[i]
            print(f"Running for index:{df_chandra['index'][i]}")

            # if folder already downloaded, skip this step
            obsid_folder = os.path.join(self.data_directory)

            """if f"{obsid}" not in os.listdir(obsid_folder): # proceed
                # 1. Download data
                self.download_obsid(obsid)
                
                # 2. Reprocess data
                self.process_obsid(obsid)

            # 3. Create src region via psfsize
            self.extract_source_region(sn_name, obsid, ra, dec)
            
            self.open_image_with_ds9(sn_name, obsid) # live time viewing image

            # 4. Generate background region
            self.generate_bkg_region(sn_name, obsid)"""

            # optional
            self.run_create_small_image(sn_name, obsid)
            self.run_wavdetect(sn_name, obsid)
            min_sep_arcsec, min_sep_ra, min_sep_dec, min_sep_significance = self.plot_wavdetect_results(sn_name, obsid)
            # save ewavelet output in separat txt file
            with open(f"{self.sample_dir}/wavdetect_output.txt", "a") as ff:    
                ff.write("\t".join(str(v) for v in df_row.values) + "\t")
                ff.write(f"{min_sep_arcsec}\t{min_sep_ra}\t{min_sep_dec}\t{min_sep_significance}\n")
            ff.close()


            # Open image with ds9
            """self.plot_save_image_with_src_bkg_reg(sn_name, obsid) # save output image

            # 5. Extract spectrum
            src_zero_counts = self.spacextract_bash(df_chandra['index'][i],sn_name, obsid)

            if not src_zero_counts:
                # 6. Calculate flux limits with srcflux
                self.srcflux_bash(sn_name, obsid)"""


pipeline = ChandraDataProcessor()
pipeline.process_spectra_from_input_file()

            
os.chdir(ROOT)




