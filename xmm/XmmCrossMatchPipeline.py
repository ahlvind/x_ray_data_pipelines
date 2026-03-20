
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.esa.xmm_newton import XMMNewton
import os
import pandas as pd


class XMMCrossMatchPipeline:
    """
    This class runs the cross matching cleanup for XMM, creates txt files with different cuts to the sample:
        1. Creates a clean output from the raw output(from inputting a list of optical coordinates in https://heasarc.gsfc.nasa.gov/cgi-bin/W3Browse/w3browse.pl) : clean_cross_match_output()
        2. Calculates and adds epochs to the output list: add_epoch_to_cross_matches()
        3. Removes the obsid with invalid epochs, a.k.a epochs <0days: removes_invalid_epochs()
        4. Removes the SN types that are not type II: removes_incorrect_SNtypes()
    """

    def __init__(self,
                 project_root="/Users/juliaahlvind/Documents/projekt_3/sample/",
                 osnc_path="/Users/juliaahlvind/Documents/projekt_1/OSNC/OSNC_60MPc.txt"):      

        # Paths
        self.project_root = project_root
        self.osnc_path = osnc_path

        # output file names
        self.raw_cross = os.path.join(project_root, "raw_cross_match_heasarc_xmm_optical.txt")
        self.cleaned = os.path.join(project_root, "cross_matches_optical_xmm.txt")
        self.epoch_filtered = os.path.join(project_root, "cross_matches_optical_xmm_epoch_cuts.txt")
        self.type_filtered = os.path.join(project_root, "cross_matches_optical_xmm_epoch_cuts_type_filtered.txt")

        # Allowed SN types for filtering
        self.allowed_types = ['II_P', 'II', 'II_P?', 'II_L', '.IIP', 'II_Pec?', 'II_Pec']

        self.cols_in = ["trash", "obsid", "status", "target_name", "ra_pnt", "dec_pnt", "obs_date", "expt_s",
                        "pi_lname", "pi_fname", "public_date", "data_in_heasarc", "search_offset_arcmin", "trash2"]

        self.cols_clean = ["index", "name", "type", "dist_Mpc", "ra", "dec", "optical_date",
                            "obsid", "sep_arcmin", "expt_ks",
                            "X_ray_obs_date", "PIname", "target"]

        self.cols_with_epoch = ["index", "name", "type", "dist_Mpc", "ra", "dec", "optical_date",
                            "obsid", "sep_arcmin", "expt_ks",
                            "X_ray_obs_date", "epoch_days", "PIname", "target"]


    @staticmethod
    def parse_offset_column(
        df: pd.DataFrame,
        col: str = "search_offset_arcmin",
        sep_col: str = "sep_arcmin",
        ra_col: str = "ra",
        dec_col: str = "dec",
        coord_col: str = "coord",
        ) -> pd.DataFrame:
        """
        Parse separation (arcmin), RA, DEC from a line like:
        '1.062 (23:04:54.96 +12:18:20.1)'

        Returns a copy of df with new columns added.
        """

        pattern = r"([\d.]+)\s*\(([^ ]+)\s+([^)]+)\)"

        out = df.copy()

        out[[sep_col, ra_col, dec_col]] = (
            out[col]
            .astype(str)
            .str.extract(pattern)
        )

        out[sep_col] = out[sep_col].astype(float)

        out[coord_col] = SkyCoord(
            ra=out[ra_col],
            dec=out[dec_col],
            unit=(u.hourangle, u.deg)
        )

        return out


    def clean_cross_match_output(self):

        # Read OSNC filewith panda: SN_name ra, dec, optical_date
        # input format: SNname	type	dist_Mpc	ra	dec	obsdate
        sn_df = pd.read_csv(self.osnc_path, delim_whitespace=True)



        # Read input file
        df = pd.read_csv(
            self.raw_cross,
            sep="|",
            skiprows=3,           # skipp info data
            names=self.cols_in,
            skipinitialspace=True,
            dtype={"obsid": str}, 
            engine="python"       # recomended?
        )

        # remove potential empty kolumns
        df = df.loc[:, df.columns.notnull()]
        # Remove empty columnes that are created due to separation with | and each row starts and ends with it
        df = df.drop(columns=["trash", "trash2"])

        # add three columns with ra, dec and sep_arcmin based on the column search_offset_arcmin
        df = self.parse_offset_column(df) # "sep_arcmin", "ra", "dec"

        # remove potentil spaces that would interfere with string matching between files
        # by creating two new columns
        sn_df["ra_str"]  = sn_df["ra"].str.strip()
        sn_df["dec_str"] = sn_df["dec"].str.strip()

        df["ra_str"]  = df["ra"].str.strip()
        df["dec_str"] = df["dec"].str.strip()

        # match OSNC line with coorinates from df
        final_df = df.merge(
            sn_df,
            on=["ra_str", "dec_str"],
            how="left",
            suffixes=("", "_sn")
        )

        # ensure keeping 0 in obsid
        final_df["obsid"] = final_df["obsid"].astype(str)
        
        # remove SN that were not in OSNC
        out_df = final_df[final_df["SNname"].notna()].copy()

        # create a copy of the columns we are interested in printing
        cols_out = ["SNname", "type", "dist_Mpc", "ra", "dec", "obsdate", "obsid", "sep_arcmin", 
                    "expt_s", "obs_date", "pi_lname", "target_name"]

        # create a copy of the panda frame
        out_df = out_df[cols_out].copy()

        # add index in the first column
        out_df.insert(0, "index", range(len(out_df)))

        # remove empty spaces in string name
        out_df['pi_lname']=out_df['pi_lname'].str.strip()
        out_df['target_name']=out_df['target_name'].str.strip()

        # strip the time stamp from X-ray obs date
        out_df['obs_date']=(pd.to_datetime(out_df['obs_date']).dt.date)

        # remove "SN" before SN name for each row
        out_df["SNname"] = out_df["SNname"].str.replace("SN", "", regex=False)

        # ensure keeping 0 in obsid
        out_df["obsid"] = out_df["obsid"].astype(str)
        print(out_df["obsid"].head())
        print(out_df["obsid"].apply(type).unique())

        # write the clean output file
        out_df.to_csv(
            self.cleaned,
            sep=" ",
            index=False,
            header=self.cols_clean
        )


# --- Step 2: Add epochs(days between optical and X-ray observations) to the output ---
    def add_epoch_to_cross_matches(self):
        """ Calculate epochs and update the cross match file """
        
        df = pd.read_csv(self.cleaned, sep=" ", skiprows=1, names=self.cols_clean, dtype={"obsid": str})


        # reformat strings to dates
        df['optical_date'] = pd.to_datetime(df['optical_date'])
        df['X_ray_obs_date'] = pd.to_datetime(df['X_ray_obs_date'])

        # Insert epoch_days after X_ray_obs_date
        df.insert(df.columns.get_loc("X_ray_obs_date") + 1,
            "epoch_days",    
            (df['X_ray_obs_date'] - df['optical_date']).dt.days # calculates epochs
            )

        # rewrite index list
        df = df.reset_index(drop=True)
        df['index'] = df.index

        df["obsid"] = df["obsid"].astype(str)

        # Save updated file
        df.to_csv(self.cleaned, sep="\t", index=False)
        print("Added epochs to cross match file")


    # --- Step 3: Remove invalid epochs where X-ray obs is before optical ---
    def removes_invalid_epochs(self):
        """ Creates the list of all valid X-ray observations, in practise this means removing observations with epochs < 0days"""
        # read the input file
        df = pd.read_csv(self.cleaned, sep="\t", skiprows=1, names=self.cols_with_epoch, dtype={"obsid": str})

        # Filter: keep only epochs_days >= 0
        df_clean = df[df['epoch_days'] >= 0].copy()

        # revwrite index list
        df_clean = df_clean.reset_index(drop=True)
        df_clean['index'] = df_clean.index

        df_clean["obsid"] = df_clean["obsid"].astype(str)

        # Save to new file
        df_clean.to_csv(
            self.epoch_filtered,
            sep="\t",
            index=False
        )

        print("Removed invalid observations based on epochs")


    # --- Stepr 4: Filter output to correct SN type ---
    def removes_incorrect_SNtypes(self):
        """ Creates the final list of obsid that should be downloaded, in practise this means removing observations for SN with incorrect type.

        All types (from  df['type'].unique() ): 'Ia/c', 'Ic', 'Ib', 'Ib/c?', 'IIb', 'Ic_BL', 'IIb/Ib/Ic(Ca_rich)',
        'Ib/c', 'Ibn', 'Ib_Pec', 'IIb/Ib/Ic', 'Ic_Pec', 'II_P', 'II',
        'II_P?', 'II_L', '.IIP', 'II_Pec?', 'II_Pec'

        All correct types: 'II_P', 'II', 'II_P?', 'II_L', '.IIP', 'II_Pec?', 'II_Pec'

        """
        df = pd.read_csv(self.epoch_filtered, sep="\t", skiprows=1, names=self.cols_with_epoch, dtype={"obsid": str})

        # filter panda frame
        df_type_filtered = df[df["type"].isin(self.allowed_types)].copy()

        # rewrite index list
        df_type_filtered = df_type_filtered.reset_index(drop=True)
        df_type_filtered["index"] = df_type_filtered.index

        df_type_filtered["obsid"] = df_type_filtered["obsid"].astype(str)

        # save to new output
        df_type_filtered.to_csv(self.type_filtered, sep="\t", index=False)
        print(f"Filtered by SN type, svaed in output: {self.type_filtered}")

        print(f"Final sample size, unique SN:{len(df_type_filtered['name'].unique())}, obs{len(df_type_filtered)}")

    def run_all(self):
        self.clean_cross_match_output()
        self.add_epoch_to_cross_matches()
        self.removes_invalid_epochs()
        self.removes_incorrect_SNtypes()
        print("Full Chandra cross match pipeline completed!")


pipeline = XMMCrossMatchPipeline()
pipeline.run_all()


