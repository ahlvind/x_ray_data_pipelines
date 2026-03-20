
import pandas as pd
import os
import re
import time
import subprocess
from pathlib import Path


class ChandraCrossMatchPipeline:
    """ 
    This class runs the cross matching for Chandra, creates txt files with different cuts to the sample:
        1. Cross matches with coordinates of optical SNe, this is calling a  bash script: run_find_chandra_obsid()
        2. Creates a clean output from the previous raw output: clean_cross_match_output()
        3. Calculates and adds epochs to the output list: add_epoch_to_cross_matches()
        4. Removes the obsid with invalid epochs, a.k.a epochs <0days: removes_invalid_epochs()
        5. Removes the SN types that are not type II: removes_incorrect_SNtypes()
    """
    def __init__(self,
                 project_root="/Users/juliaahlvind/Documents/projekt_3/sample/",
                 osnc_path="/Users/juliaahlvind/Documents/projekt_1/OSNC/OSNC_60MPc.txt"):
        
        # Paths
        self.project_root = project_root
        self.osnc_path = osnc_path
        
        # output file names
        self.raw_cross = os.path.join(project_root, "raw_cross_matches_optical_chandra.txt")
        self.cleaned = os.path.join(project_root, "cross_matches_optical_chandra.txt")
        self.epoch_filtered = os.path.join(project_root, "cross_matches_optical_chandra_epoch_cuts.txt")
        self.type_filtered = os.path.join(project_root, "cross_matches_optical_chandra_epoch_cuts_type_filtered.txt")

        # Allowed SN types for filtering
        self.allowed_types = ['II_P', 'II', 'II_P?', 'II_L', '.IIP', 'II_Pec?', 'II_Pec']

        # Column definitions
        self.cols_no_epoch = [
            "index", "name", "type", "dist_Mpc", "ra", "dec", "optical_date",
            "obsid", "sep_arcmin", "instr", "grating", "expt_ks",
            "X_ray_obs_date", "PIname", "target"
        ]

        self.cols_with_epoch = [
            "index", "name", "type", "dist_Mpc", "ra", "dec", "optical_date",
            "obsid", "sep_arcmin", "instr", "grating", "expt_ks",
            "X_ray_obs_date", "epoch_days", "PIname", "target"
        ]

    @staticmethod
    def parse_observation_line(line):
        """ Read Chandra cross match file and extract data.
        obsid, separation, instrument, grating, exposure time, obs date X-ray, PI, target of observation. """
        
        # this function is used by clean_cross_match_output()

        parts = line.split()
        if len(parts) < 8:
            return None

        obsid = parts[0]
        sepn = parts[1]
        inst = parts[2]
        grat = parts[3]
        time = parts[4]
        obsdate = parts[5]
        piname = parts[6]
        target = " ".join(parts[7:]).strip('"')

        return obsid, sepn, inst, grat, time, obsdate, piname, target


    # --- Step 0: Run the external bash script to find cross matches with coordinates ---
    def run_find_chandra_obsid(self):
        """ Running the find_chandra_obsid. Cross match optical SN from OSNC_60Mpc with Chandra archive. """
        start = time.time()
        cmd = f"./find_chandra_obsid.sh"
        try: 
            subprocess.run(
                cmd,
                shell=True,
                check=True
            )
            print(f"Pipeline finished successfully in {time.time() - start:.2f} s")

        except subprocess.CalledProcessError as e:
            print("find_chandra_obsid pipeline failed:", e)
            raise


    # --- Step 1: Clean the raw output from cross match ---
    def clean_cross_match_output(self):
        """ This function reads the output of the find_chandra_obsid and generates a clean txt"""

        print("Cleaning raw cross-match output...")
        

        # Read OSNC file to create dictionary with: SN → {ra, dec, optical_date}
        # Format:SN  ra  dec  optical_date
        # input format: SNname	type	dist_Mpc	ra	dec	obsdate
        # Load OSNC file
        extra = {}
        with open(self.osnc_path, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 6:
                    sn = parts[0].replace("SN", "")
                    extra[sn] = {
                        "SNtype": parts[1],
                        "dist_Mpc": parts[2],
                        "ra": parts[3],
                        "dec": parts[4],
                        "optical_date": parts[5]
                    }

        output_rows = []
        current_sn = None
        reading_block = False
        index = 0

        # reading the input file from cross match of chandra obsid.
        with open(self.raw_cross, "r") as f:
            for line in f:
                line = line.rstrip()

                # Find new SN
                if re.match(r"^\d{4}[A-Za-z]+$", line): #Ankare för början av raden. Säger: “mönstret måste börja här”.\d{4} Exakt fyra siffror (årtalet), t.ex. 1991. [A-Za-z]+En eller flera latinska bokstäver (A–Z eller a–z). Det täcker både stora och små bokstäver, och är skiftlägeskänsligt (dvs A ≠ a i jämförelse, men båda accepteras). $ Ankare för slutet av raden. Säger: “mönstret måste sluta här”. Det gör att vi inte accepterar extra tecken efter.
                    current_sn = line.strip()
                    reading_block = False
                    continue

                # Header-linje, mark that block with observations begins here
                if line.startswith("# obsid"):
                    reading_block = True
                    continue

                # if entering a block of obsid: save output data per obs
                if reading_block and line.strip() and not line.startswith("#"): # skip first line with data column headers like obsid etc...
                    parsed = self.parse_observation_line(line) # read obsid line
                    if parsed:
                        obsid, sepn, inst, grat, time, obsdate, piname, target = parsed

                        # extract info from OSNC that matches current SN 
                        ra = extra.get(current_sn, {}).get("ra", "NA")
                        dec = extra.get(current_sn, {}).get("dec", "NA")
                        optical_date = extra.get(current_sn, {}).get("optical_date", "NA")
                        sn_type = extra.get(current_sn, {}).get("SNtype", "NA")
                        dist_Mpc = extra.get(current_sn, {}).get("dist_Mpc", "NA")

                        row = (
                            f"{index}\t{current_sn}\t{sn_type}\t{dist_Mpc}\t"
                            f"{ra}\t{dec}\t{optical_date}\t"
                            f"{obsid}\t{sepn}\t{inst}\t{grat}\t{time}\t{obsdate}\t"
                            f"{piname}\t{target}"
                        )
                        output_rows.append(row)
                        index += 1


        # Write formatted output
        with open(self.cleaned, "w") as f:
            header = " ".join(self.cols_no_epoch) + "\n"
            f.write(header)
            for row in output_rows:
                f.write(row + "\n")

        print(f" Cleaned cross match output finnished. Saved to: {self.cleaned}")

        # OBS! If you have an observation from Santos Lleo, manually change the txt after running clean_cross_match_output from "" to none ! :) beacus this is the only one that saved their name with "" in the database
        # ooor sometimes it works


    # --- Step 2: Add epochs(days between optical and X-ray observations) to the output ---
    def add_epoch_to_cross_matches(self):
        """ Calculate epochs and update the cross match file """
        
        df = pd.read_csv(self.cleaned, sep="\t", skiprows=1, names=self.cols_no_epoch)

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

        # Save updated file
        df.to_csv(self.cleaned, sep="\t", index=False)
        print("Added epochs to cross match file")


    # --- Step 3: Remove invalid epochs where X-ray obs is before optical ---
    def removes_invalid_epochs(self):
        """ Creates the list of all valid X-ray observations, in practise this means removing observations with epochs < 0days"""
        # read the input file
        df = pd.read_csv(self.cleaned, sep="\t", skiprows=1, names=self.cols_with_epoch)

        # Filter: keep only epochs_days >= 0
        df_clean = df[df['epoch_days'] >= 0].copy()

        # revwrite index list
        df_clean = df_clean.reset_index(drop=True)
        df_clean['index'] = df_clean.index

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
        df = pd.read_csv(self.epoch_filtered, sep="\t", skiprows=1, names=self.cols_with_epoch)

        # filter panda frame
        df_type_filtered = df[df["type"].isin(self.allowed_types)].copy()

        # rewrite index list
        df_type_filtered = df_type_filtered.reset_index(drop=True)
        df_type_filtered["index"] = df_type_filtered.index

        # save to new output
        df_type_filtered.to_csv(self.type_filtered, sep="\t", index=False)
        print(f"Filtered by SN type, svaed in output: {self.type_filtered}")


    def run_all(self):
        self.clean_cross_match_output()
        self.add_epoch_to_cross_matches()
        self.removes_invalid_epochs()
        self.removes_incorrect_SNtypes()
        print("Full Chandra cross match pipeline completed!")


    # --- Plott the final list of SNe ---
    def plot_output(self):
        df = pd.read_csv(self.type_filtered, sep="\t", skiprows=1, names=self.cols_with_epoch)

        print(len(df['name']))
        """sns.lineplot(
            data=df,
            x="epoch_days",
            y="name",
            hue="name",
            marker="o"
        )

        plt.show()"""


pipeline = ChandraCrossMatchPipeline()
pipeline.run_all()
pipeline.plot_output()



