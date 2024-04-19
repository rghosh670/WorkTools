"""
Class definition for a BED object.
"""

import os
from pprint import pprint

import datatable as dt
import numpy as np
import pandas as pd
import pysam
from icecream import ic

import worktools.api as api
import worktools.api.utils as utils


class BEDObject:
    def __init__(
        self,
        bed_file_path=None,
        bed_cols=None,
        bed_df=None,
        add_rt_style_index=False,
        ref_alt_cols=None,
    ):
        self.bed_df = bed_df
        ic(self.bed_df)
        self.bed_file_path = bed_file_path
        self.bed_cols = bed_cols

        if bed_file_path and bed_df is None:
            self.read_bed(bed_file_path)

        if add_rt_style_index:
            self.add_rt_style_index(ref_alt_cols)

        ic(self.bed_df)
        return

    def __str__(self):
        if self.bed_df is None:
            return "Empty BED object"

        return self.bed_df.__str__()

    def __repr__(self):
        if self.bed_df is None:
            return "Empty BED object"

        return self.bed_df.__str__()

    def has_chr_prefix(self):
        """
        Checks if the '#CHROM' column of the VCF file has the 'chr' prefix.

        Returns:
        - True if the 'chr' prefix is present, False otherwise.
        """
        return self.bed_df["#CHROM"][0].startswith("chr")

    def sort_bed_df(self):
        chrom_mapping = {f"chr{i}": i for i in range(1, 23)}  # Numeric chromosomes
        chrom_mapping.update(
            {"chrX": 23, "chrY": 24, "chrM": 25}
        )  # Non-numeric chromosomes

        # Identify all unique chromosomes in the DataFrame
        all_chroms = set(self.bed_df["#CHROM"].unique())

        # Find chromosomes not already in the mapping
        unmapped_chroms = [chrom for chrom in all_chroms if chrom not in chrom_mapping]

        # Assign these unmapped chromosomes a numeric key starting from 26
        next_chrom_key = 26
        for chrom in sorted(unmapped_chroms):  # Sorting ensures consistent mapping
            chrom_mapping[chrom] = next_chrom_key
            next_chrom_key += 1

        # Step 1: Create a sortable numeric key based on the updated mapping
        self.bed_df["chrom_key"] = self.bed_df["#CHROM"].map(chrom_mapping)

        # Step 2: Sort by this numeric key and by position
        self.bed_df = self.bed_df.sort_values(by=["chrom_key", "chromStart"])

        # Optionally, remove the temporary 'chrom_key' if you don't need it anymore
        self.bed_df.drop("chrom_key", axis=1, inplace=True)
        return

    def read_bed(self, bed_file_path):
        """
        Reads a BED file into a pandas DataFrame.

        Parameters:
        - bed_file_path: Path to the VCF file.

        Returns:
        - A pandas DataFrame containing the VCF data.
        """
        self.bed_file_path = bed_file_path
        first_line = pd.read_csv(bed_file_path, sep="\t", nrows=1, header=None).iloc[0]
        has_header = type(first_line.iloc[1]) != np.int64

        # Example: {'#CHROM': 0, 'chromStart': 1, 'chromEnd': 1}
        if self.bed_cols:
            self.bed_df = dt.fread(
                bed_file_path, sep="\t", header=has_header, columns=dt.str32
            ).to_pandas()
            bed_df_cols = self.bed_df.columns.to_list()

            if 0 not in self.bed_cols.values():
                smallest_col_num = min(self.bed_cols.values())
                cols_to_move = bed_df_cols[:smallest_col_num]
                bed_df_cols = bed_df_cols[smallest_col_num:]
                bed_df_cols = bed_df_cols + cols_to_move
                self.bed_df = self.bed_df[bed_df_cols]
                self.bed_cols = {
                    key: val - smallest_col_num for key, val in self.bed_cols.items()
                }

            for key, val in self.bed_cols.items():
                bed_df_cols[val] = key

            self.bed_df.columns = bed_df_cols

            if "chromEnd" not in bed_df_cols:
                self.bed_df["chromEnd"] = self.bed_df.iloc[
                    :, self.bed_cols["chromStart"]
                ]
                self.bed_df = self.bed_df.drop("chromStart", axis=1)

            self.bed_df["chromStart"] = (
                self.bed_df.iloc[:, self.bed_cols["chromEnd"]].astype(int) - 1
            )
            columns = list(self.bed_df.columns)
            columns.remove("chromStart")
            columns.insert(1, "chromStart")
            self.bed_df = self.bed_df[columns]

            if not self.has_chr_prefix():
                self.bed_df["#CHROM"] = "chr" + self.bed_df["#CHROM"]

            return
        else:
            self.bed_df = dt.fread(
                bed_file_path, sep="\t", header=has_header
            ).to_pandas()
            self.bed_df.columns = ["#CHROM", "chromStart", "chromEnd"] + [
                i for i in self.bed_df.columns.to_list()[3:]
            ]
        return

    def write_to_file(self, file_path, include_header=False, tabix=False):
        """
        Writes the VCF data to a file.

        Parameters:
        - file_path: Path to the output file.
        """
        utils.run_shell_command(f"> {file_path}")
        self.bed_df.to_csv(
            file_path, sep="\t", index=False, header=include_header, mode="a"
        )

        if tabix:
            pysam.tabix_index(file_path, force=True, preset="bed")
        return

    def has_chr_prefix(self):
        """
        Checks if the '#CHROM' column of the BED file has the 'chr' prefix.

        Returns:
        - True if the 'chr' prefix is present, False otherwise.
        """
        return self.bed_df["#CHROM"][0].startswith("chr")

    def is_original_file_valid(self):
        pass

    def lift_over_bed(self, genome_build):
        return api.lift_over_bed(self, genome_build)

    # TODO - Make sure that if the file is modified, and not written to file, this returns False (this should be an atrribute, not a function)
    def is_original_bed_file_valid(self):
        if self.bed_file_path is None:
            return False

        return True

    def add_rt_style_index(self, ref_alt_cols=None, force=False):
        if (
            "REF" not in self.bed_df.columns or "ALT" not in self.bed_df.columns
        ) and ref_alt_cols is None:
            raise ValueError(
                "BED file must have 'REF' and 'ALT' columns or their integer locations must be passed in to ref_alt_cols to add an 'rt-style' index."
            )
        if ref_alt_cols:
            ref_col_num = ref_alt_cols[0]
            alt_col_num = ref_alt_cols[1]
        else:
            ref_col_num = self.bed_df.columns.get_loc("REF")
            alt_col_num = self.bed_df.columns.get_loc("ALT")

        self.bed_df.index = (
            self.bed_df["#CHROM"]
            + ";"
            + self.bed_df["chromEnd"].astype(str)
            + ";"
            + self.bed_df.iloc[:, ref_col_num].astype(str)
            + ";"
            + self.bed_df.iloc[:, alt_col_num].astype(str)
        )
        return
