"""
Class definition for a VCF object.
Provides utilities for reading, manipulating, and analyzing VCF files.
"""

import os
import re
import subprocess
import tempfile
from pprint import pprint

import datatable as dt
import pandas as pd
import pysam
from Bio import Entrez
from cyvcf2 import VCF
from icecream import ic
from loguru import logger
from pyliftover import LiftOver

import worktools
import worktools.api as api
import worktools.api.utils as utils

GRCh38_CHROM_LENGTH_DICT = {
    "chr1": "##contig=<ID=chr1,length=248956422,assembly=gnomAD_GRCh38>",
    "chr2": "##contig=<ID=chr2,length=242193529,assembly=gnomAD_GRCh38>",
    "chr3": "##contig=<ID=chr3,length=198295559,assembly=gnomAD_GRCh38>",
    "chr4": "##contig=<ID=chr4,length=190214555,assembly=gnomAD_GRCh38>",
    "chr5": "##contig=<ID=chr5,length=181538259,assembly=gnomAD_GRCh38>",
    "chr6": "##contig=<ID=chr6,length=170805979,assembly=gnomAD_GRCh38>",
    "chr7": "##contig=<ID=chr7,length=159345973,assembly=gnomAD_GRCh38>",
    "chr8": "##contig=<ID=chr8,length=145138636,assembly=gnomAD_GRCh38>",
    "chr9": "##contig=<ID=chr9,length=138394717,assembly=gnomAD_GRCh38>",
    "chr10": "##contig=<ID=chr10,length=133797422,assembly=gnomAD_GRCh38>",
    "chr11": "##contig=<ID=chr11,length=135086622,assembly=gnomAD_GRCh38>",
    "chr12": "##contig=<ID=chr12,length=133275309,assembly=gnomAD_GRCh38>",
    "chr13": "##contig=<ID=chr13,length=114364328,assembly=gnomAD_GRCh38>",
    "chr14": "##contig=<ID=chr14,length=107043718,assembly=gnomAD_GRCh38>",
    "chr15": "##contig=<ID=chr15,length=101991189,assembly=gnomAD_GRCh38>",
    "chr16": "##contig=<ID=chr16,length=90338345,assembly=gnomAD_GRCh38>",
    "chr17": "##contig=<ID=chr17,length=83257441,assembly=gnomAD_GRCh38>",
    "chr18": "##contig=<ID=chr18,length=80373285,assembly=gnomAD_GRCh38>",
    "chr19": "##contig=<ID=chr19,length=58617616,assembly=gnomAD_GRCh38>",
    "chr20": "##contig=<ID=chr20,length=64444167,assembly=gnomAD_GRCh38>",
    "chr21": "##contig=<ID=chr21,length=46709983,assembly=gnomAD_GRCh38>",
    "chr22": "##contig=<ID=chr22>",
    "chrX": "##contig=<ID=chrX,length=156040895,assembly=gnomAD_GRCh38>",
    "chrY": "##contig=<ID=chrY,length=57227415,assembly=gnomAD_GRCh38>",
    "chrM": "##contig=<ID=chrM,length=16569,assembly=gnomAD_GRCh38>",
}


class VCFObject:
    """
    A class for handling operations on VCF files including reading, lifting over, and annotating VCF files.

    Attributes:
        vcf_df (DataFrame): Pandas DataFrame containing VCF data.
        pysam_vcf (VariantFile): pysam object to interact with VCF data.
        vcf_file_path (str): Path to the VCF file.
        lift_over_obj (LiftOver): pyliftover object for genome coordinate conversions.
        has_rt_style_index (bool): Indicates if the VCF DataFrame has an RT-style index.
        metadata (str): Metadata from the VCF file.
        info_col_df (DataFrame): DataFrame derived from INFO column of VCF data.
    """
    def __init__(self, vcf_file_path=None, vcf_df=None, read_file=True):
        self.vcf_df = vcf_df
        self.pysam_vcf = None
        self.vcf_file_path = vcf_file_path
        self.lift_over_obj = None
        self.has_rt_style_index = False
        self.metadata = None
        self.info_col_df = None

        if vcf_file_path and read_file:
            self.read_vcf(vcf_file_path)

    def __str__(self):
        if self.vcf_df is None:
            return f"Empty VCF object :("

        metadata_str = ""
        if self.metadata:
            num_metadata_lines = self.metadata.count("\n")

            max_metadata_lines = 10
            half_max_metadata_lines = max_metadata_lines // 2

            if num_metadata_lines > max_metadata_lines:
                metadata_lines = self.metadata.split("\n")
                metadata_lines = (
                    metadata_lines[:half_max_metadata_lines]
                    + ["..."]
                    + metadata_lines[-half_max_metadata_lines:]
                )
                metadata_str = "\n".join(metadata_lines).strip()

        return_str = f"{metadata_str}\n{self.vcf_df.__str__()}"

        return return_str

    def __repr__(self):
        return self.__str__()

    def read_vcf(self, vcf_file_path):
        """
        Reads a VCF file into a pandas DataFrame.

        Parameters:
        - vcf_file_path: Path to the VCF file.

        Returns:
        - A pandas DataFrame containing the VCF data.
        """
        self.vcf_file_path = vcf_file_path
        column_types = {
            "#CHROM": dt.str32,
            "POS": dt.int32,
            "ID": dt.str32,
            "REF": dt.str32,
            "ALT": dt.str32,
            "QUAL": dt.str32,
            "FILTER": dt.str32,
            "INFO": dt.str32,
        }

        fread_cmd = f"bcftools view {self.vcf_file_path} | grep -v '##'"
        if vcf_file_path.endswith(".vcf.gz"):
            if not utils.is_bgzipped(vcf_file_path):
                raise ValueError("VCF file is not bgzipped but ends with .gz")
            fread_cmd = f"zgrep -v '^##' {vcf_file_path}"
        elif vcf_file_path.endswith(".vcf"):
            fread_cmd = f"grep -v '^##' {vcf_file_path}"
        elif vcf_file_path.endswith(".bcf"):
            fread_cmd = f"bcftools view {self.vcf_file_path} | grep -v '##'"

        self.vcf_df = dt.fread(
            cmd=fread_cmd, sep="\t", header=True, columns=column_types
        ).to_pandas()

        logger.info("Read VCF file %s", vcf_file_path)

        if not self.has_chr_prefix():
            self.adjust_chr_prefix_in_vcf("add", True)
            logger.info(
                "'chr' prefix was not found, and was added (AND WRITTEN) to the file. Run rt.adjust_chr_prefix_in_vcf(VCFObject, 'remove') to remove it."
            )
        self.pysam_vcf = pysam.VariantFile(vcf_file_path, threads=os.cpu_count())
        if not self.has_index() or self.is_index_outdated():
            self.index_variant_file()

        self.metadata = self.read_metadata()
        self.vcf_df = self.vcf_df.fillna(".")

        return

    def is_index_outdated(self):
        if self.has_index():
            data_file_mod_time = os.path.getmtime(self.vcf_file_path)
            suffix = ".csi" if self.vcf_file_path.endswith(".bcf") else ".tbi"
            index_file_mod_time = os.path.getmtime(f"{self.vcf_file_path}{suffix}")
            return index_file_mod_time < data_file_mod_time
        return True

    def has_chr_prefix(self):
        """
        Checks if the '#CHROM' column of the VCF file has the 'chr' prefix.

        Returns:
        - True if the 'chr' prefix is present, False otherwise.
        """
        ret = True
        try:
            ret = self.vcf_df["#CHROM"].iloc[0].startswith("chr")
        except:
            pass
        return ret

    def write_metadata_to_file(self, file_path, do_not_generate_metadata=False):
        if not do_not_generate_metadata:
            self.generate_metadata_from_vcf_df()

        with open(file_path, "w") as file:
            if self.metadata:
                file.write(self.metadata.strip() + "\n")

    def write_to_file(
        self,
        file_path,
        write_metadata=True,
        bgzip=False,
        tabix=False,
        do_not_generate_metadata=False,
        bcf=False,
    ):
        """
        Writes the VCF data to a file.

        Parameters:
        - file_path: Path to the output file.
        """
        utils.run_shell_command(f"> {file_path}")

        if write_metadata:
            if not do_not_generate_metadata:
                self.generate_metadata_from_vcf_df()

            with open(file_path, "w") as file:
                if self.metadata:
                    file.write(self.metadata.strip() + "\n")

        self.vcf_df = self.vcf_df.fillna(".")
        self.sort_vcf_df()
        self.vcf_df.to_csv(file_path, sep="\t", index=False, mode="a", na_rep=".")
        self.vcf_file_path = file_path

        if tabix:
            os.system(f"mv {self.vcf_file_path} {self.vcf_file_path.rstrip('.gz')}")
            self.vcf_file_path = self.vcf_file_path.rstrip(".gz")
            self.index_variant_file()
            return

        if bgzip:
            if file_path.endswith(".gz"):
                file_path = file_path.rstrip(".gz")
            bgzip_cmd = f"bgzip -f {file_path}"
            utils.run_shell_command(bgzip_cmd)

        if bcf:
            with tempfile.NamedTemporaryFile(
                mode="w+", dir="/tmp", delete=False, suffix=".bcf"
            ) as bcf_temp_file:
                bcf_cmd = f"bcftools view -Ob -o {bcf_temp_file.name} {file_path}"
                utils.run_shell_command(bcf_cmd)
                bcf_temp_file.flush()

                mv_cmd = f"mv {bcf_temp_file.name} {file_path}"
                utils.run_shell_command(mv_cmd)

                index_cmd = f"bcftools index -f {file_path}"
                utils.run_shell_command(index_cmd)
        return

    def index_variant_file(self):
        if self.vcf_file_path.endswith(".vcf"):
            self.vcf_file_path = pysam.tabix_index(
                self.vcf_file_path, force=True, preset="vcf"
            )
            logger.info(
                "No tabix index was found for VCF file. File was bgzipped (%s) and indexed (%s).",
                self.vcf_file_path,
                self.vcf_file_path + ".tbi",
            )
        elif self.vcf_file_path.endswith(".bcf"):
            index_cmd = f"bcftools index -f {self.vcf_file_path}"
            logger.info(
                "No csi index was found for BCF file. File (%s) was indexed (%s).",
                self.vcf_file_path,
                self.vcf_file_path + ".csi",
            )
            utils.run_shell_command(index_cmd)
        elif self.vcf_file_path.endswith(".vcf.gz"):
            index_cmd = f"tabix -p vcf -f {self.vcf_file_path}"
            utils.run_shell_command(index_cmd)
            logger.info(
                "No tabix index was found for VCF file. File (%s) was indexed (%s).",
                self.vcf_file_path,
                self.vcf_file_path + ".tbi",
            )
        return

    def adjust_chr_prefix_in_vcf(self, action, write_to_file=False):
        """
        Adjusts the 'chr' prefix in the '#CHROM' column of a VCF dataframe.

        Parameters:
        - vcf_df: A pandas DataFrame containing VCF data.
        - action: A string specifying the action to perform: 'add' to add 'chr' prefix, 'remove' to remove it.
        """

        match action:
            case "add":
                if not self.has_chr_prefix():
                    self.vcf_df["#CHROM"] = "chr" + self.vcf_df["#CHROM"]

                    if write_to_file:
                        if utils.is_bgzipped(self.vcf_file_path):
                            temp_file = tempfile.NamedTemporaryFile(
                                delete=False, suffix=".vcf"
                            )
                            add_chr_command = f"zcat {self.vcf_file_path} | sed -i '/^[^#]/s/^/chr/' {self.vcf_file_path} | bgzip > {temp_file.name} && mv {temp_file.name} {self.vcf_file_path}"
                            temp_file.close()
                        else:
                            add_chr_command = (
                                f"sed -i '/^[^#]/s/^/chr/' {self.vcf_file_path}"
                            )
                        utils.run_shell_command(add_chr_command)
            case "remove":
                if self.has_chr_prefix():
                    self.vcf_df["#CHROM"] = self.vcf_df["#CHROM"].str.lstrip("chr")

                    if write_to_file:
                        if utils.is_bgzipped(self.vcf_file_path):
                            temp_file = tempfile.NamedTemporaryFile(
                                delete=False, suffix=".vcf"
                            )
                            add_chr_command = f"zcat {self.vcf_file_path} | sed -i '/^[^#]/s/^chr//' {self.vcf_file_path} | bgzip > {temp_file.name} && mv {temp_file.name} {self.vcf_file_path}"
                            temp_file.close()
                        else:
                            add_chr_command = (
                                f"sed -i '/^[^#]/s/^chr//' {self.vcf_file_path}"
                            )
                        utils.run_shell_comand(add_chr_command)
            case _:
                raise ValueError("Invalid action specified. Choose 'add' or 'remove'.")

        return

    def get_genome_build_vcf(self):
        """
        Returns the genome build of the VCF file.

        Returns:
        - A string containing the genome build (GRCh38 or GRCh37).
        """
        variant_ids = self.vcf_df["ID"]
        if not variant_ids.str.startswith("rs").any():
            raise NotImplementedError(
                "Getting the genome build without any variant with rsIDs is not yet supported."
            )

        query_variant_index = variant_ids[
            variant_ids.str.startswith("rs")
        ].first_valid_index()
        query_variant_id = variant_ids[query_variant_index]

        try:
            Entrez.email = "rohit.ghosh@yale.edu"
            handle = Entrez.efetch(db="snp", id=query_variant_id, retmode="text")
            snp_data = handle.read()

            hg38_match = re.search(r"<CHRPOS>(.*?)</CHRPOS>", snp_data)
            hg19_match = re.search(
                r"<CHRPOS_PREV_ASSM>(.*?)</CHRPOS_PREV_ASSM>", snp_data
            )

            if hg38_match:
                grch38_match = hg38_match.group(1)
                grch38_chrom = "chr" + grch38_match.split(":")[0]
                grch38_pos = grch38_match.split(":")[1]
            else:
                logger.error("GRCh38 position not found")
                return

            if hg19_match:
                grch37_match = hg19_match.group(1)
                grch37_chrom = "chr" + grch37_match.split(":")[0]
                grch37_pos = grch37_match.split(":")[1]
            else:
                logger.error("GRCh37 position not found")
                return

        except Exception as e:
            logger.error("Error fetching SNP data: %s", e)
            return
        finally:
            handle.close()

        query_variant_chrom = self.vcf_df["#CHROM"][query_variant_index]
        query_variant_pos = str(self.vcf_df["POS"][query_variant_index])

        if query_variant_chrom != grch38_chrom and query_variant_chrom != grch37_chrom:
            logger.error("Chromosome doesn't match rsID")
            return

        if query_variant_pos != grch38_pos and query_variant_pos != grch37_pos:
            logger.error("Position doesn't match hg38 or hg19")
            return

        if query_variant_pos == grch38_pos:
            return "hg38"

        return "hg19"

    def lift_over_vcf(self, genome_build):
        return api.lift_over_vcf(self, genome_build)

    def get_malinois_scores(scores="all"):
        return api.get_malinois_scores(self, scores)

    def has_index(self):
        if self.vcf_file_path.endswith(".vcf.gz"):
            return os.path.exists(self.vcf_file_path + ".tbi")
        elif self.vcf_file_path.endswith(".bcf"):
            return os.path.exists(self.vcf_file_path + ".csi")
        return False

    def fetch(self, chrom, start, stop):
        rows = []

        for record in self.pysam_vcf.fetch(chrom, start, stop):
            rows.append(record.__str__().split("\t"))

        fetch_result_df = pd.DataFrame(rows)
        fetch_result_df.columns = self.vcf_df.columns
        return fetch_result_df

    def get_liftover_object(self, original_genome_build, new_genome_build):
        if self.lift_over_obj is not None and self.lift_over_obj[0] == (
            f"{original_genome_build}_{new_genome_build}"
        ):
            return self.lift_over_obj[1]

        # Dictionary to map genome build aliases
        build_aliases = {
            "GRCh37": "hg19",
            "GRCh38": "hg38",
        }

        # Replace genome build names if they are in aliases
        original_genome_build = build_aliases.get(
            original_genome_build, original_genome_build
        )
        new_genome_build = build_aliases.get(new_genome_build, new_genome_build)

        # Validate genome builds
        valid_builds = ["hg19", "hg38"]
        if (
            original_genome_build not in valid_builds
            or new_genome_build not in valid_builds
        ):
            raise ValueError(
                "Invalid genome build. Only 'hg19', 'hg38', 'GRCh37', and 'GRCh38' are supported."
            )

        # Create the LiftOver object
        self.lift_over_obj = (
            f"{original_genome_build}_{new_genome_build}",
            LiftOver(original_genome_build, new_genome_build),
        )
        return self.lift_over_obj[1]

    def get_info_col_df(self, tags=None, cache=False):
        if cache and self.info_col_df is not None:
            return self.info_col_df

        vcf_df = self.vcf_df

        if tags is not None:
            vcf_df = self.vcf_df.copy(deep=True)
            condition = vcf_df["INFO"].str.contains(tags[0])
            for tag in tags[1:]:
                condition = condition & vcf_df["INFO"].str.contains(tag)

            # Filter the DataFrame based on the condition
            vcf_df = vcf_df.loc[condition]

        # Extract all key-value pairs, keeping values as strings
        pattern = r"(?P<Key>\w+)=([^;]+);?"

        # Extracting the matches
        matches = vcf_df["INFO"].str.extractall(pattern)

        # Aggregate each key into a list within a dictionary comprehension
        data_dict = {key: group[1].tolist() for key, group in matches.groupby("Key")}

        # Construct the DataFrame from the aggregated dictionary
        df = pd.DataFrame(data_dict)

        self.info_col_df = df

        return df

    def sort_vcf_df(self):
        chrom_mapping = {f"chr{i}": i for i in range(1, 23)}  # Numeric chromosomes
        chrom_mapping.update(
            {"chrX": 23, "chrY": 24, "chrM": 25}
        )  # Non-numeric chromosomes

        # Identify all unique chromosomes in the DataFrame
        all_chroms = set(self.vcf_df["#CHROM"].unique())

        # Find chromosomes not already in the mapping
        unmapped_chroms = [chrom for chrom in all_chroms if chrom not in chrom_mapping]

        # Assign these unmapped chromosomes a numeric key starting from 26
        next_chrom_key = 26
        for chrom in sorted(unmapped_chroms):  # Sorting ensures consistent mapping
            chrom_mapping[chrom] = next_chrom_key
            next_chrom_key += 1

        # Step 1: Create a sortable numeric key based on the updated mapping
        self.vcf_df["chrom_key"] = self.vcf_df["#CHROM"].map(chrom_mapping)

        # Step 2: Sort by this numeric key and by position
        self.vcf_df = self.vcf_df.sort_values(by=["chrom_key", "POS"])

        # Optionally, remove the temporary 'chrom_key' if you don't need it anymore
        self.vcf_df.drop("chrom_key", axis=1, inplace=True)
        self.vcf_df.reset_index(drop=True, inplace=True)
        return

    # TODO - Make sure that if the file is modified, and not written to file, this returns False (this should be an atrribute, not a function)
    def is_original_vcf_file_valid(self):
        if self.bed_file_path is None:
            return False

        return True

    def add_rt_style_index(self, force=False):
        if force:
            self.vcf_df.index = (
                self.vcf_df["#CHROM"]
                + ";"
                + self.vcf_df["POS"].astype(str)
                + ";"
                + self.vcf_df["REF"]
                + ";"
                + self.vcf_df["ALT"]
            )
            self.has_rt_style_index = True
            return

        if not self.has_rt_style_index:
            self.vcf_df.index = (
                self.vcf_df["#CHROM"]
                + ";"
                + self.vcf_df["POS"].astype(str)
                + ";"
                + self.vcf_df["REF"]
                + ";"
                + self.vcf_df["ALT"]
            )
            self.has_rt_style_index = True
        return

    def drop_column(self, column_name, replace_chr="."):
        if column_name not in self.vcf_df.columns:
            raise ValueError(f"Column {column_name} not found in VCF data.")

        self.vcf_df[column_name] = [replace_chr] * len(self.vcf_df)

    def keep_in_info_column(self, values_to_keep):
        # self.vcf_df['INFO'] = self.vcf_df['INFO'].apply(lambda x: ';'.join([i for i in x.split(';') if i.split('=')[0] in values_to_keep]))
        self.vcf_df["INFO"] = self.vcf_df["INFO"].apply(
            lambda x: ";".join(
                [
                    i
                    for i in x.split(";")
                    if "=" in i and i.split("=")[0] in values_to_keep
                ]
            )
        )
        return

    def add_to_info_column(self, col):
        self.vcf_df["INFO"] = self.vcf_df["INFO"] + ";" + col
        return

    def remove_from_info_column(self, values_to_remove):
        # self.vcf_df['INFO'] = self.vcf_df['INFO'].apply(lambda x: ';'.join([i for i in x.split(';') if i.split('=')[0] not in values_to_remove]))
        self.vcf_df["INFO"] = self.vcf_df["INFO"].apply(
            lambda x: ";".join(
                [
                    i
                    for i in x.split(";")
                    if "=" in i and i.split("=")[0] not in values_to_remove
                ]
            )
        )
        return

    def drop_rows_without_info_tag(self, tags):
        self.vcf_df = self.vcf_df.loc[self.vcf_df["INFO"].str.contains("&".join(tags))]
        return

    def add_prefix_to_info_values(self, prefix):
        self.vcf_df["INFO"] = self.vcf_df["INFO"].apply(
            lambda x: ";".join([f"{prefix}_{i}" for i in x.split(";")])
        )
        return

    def remove_columns_from_file(
        self,
        cols_to_remove,
        update_file=True,
        update_vcf_df=True,
        output_file_format="z",
        output_file_path=None,
    ):
        cols_to_remove = ",".join(cols_to_remove)

        suffix = None

        match output_file_format:
            case "z":
                suffix = ".vcf.gz"
            case "b":
                suffix = ".bcf"
            case "v":
                suffix = ".vcf"

        with tempfile.NamedTemporaryFile(
            mode="w+", dir="/tmp", delete=False, suffix=suffix
        ) as vcf_temp_file:
            annotate_cmd = f"bcftools annotate --remove {cols_to_remove} {self.vcf_file_path} -O{output_file_format} -o {vcf_temp_file.name}"

            vcf_temp_file.flush()
            annotate_cmd_output = utils.run_shell_command(annotate_cmd)

            if "Error" in annotate_cmd_output:
                raise RuntimeError(f"Error annotating VCF file: {annotate_cmd_output}")

            logger.info("%s columns were removed the source VCF file.", cols_to_remove)

        if output_file_path:
            mv_cmd = f"mv {vcf_temp_file.name} {output_file_path}"
            utils.run_shell_command(mv_cmd)

        if update_file:
            mv_cmd = f"mv {vcf_temp_file.name} {self.vcf_file_path}"
            utils.run_shell_command(mv_cmd)

        if update_vcf_df:
            self.read_vcf(self.vcf_file_path)
        return

    def generate_metadata_from_vcf_df(self):
        chroms = self.vcf_df["#CHROM"].unique()
        ic(self.vcf_df)
        ic()
        pprint(chroms)
        metadata = [
            GRCh38_CHROM_LENGTH_DICT.get(
                chrom,
                f"##contig=<ID=chr{chrom}, length=248956422,assembly=gnomAD_GRCh38>",
            )
            for chrom in chroms
        ]
        metadata = ["##fileformat=VCFv4.2"] + metadata
        info_strings = self.vcf_df["INFO"]
        ic(info_strings.head())

        # Split by semicolon to get all key-value pairs in a list for each row
        key_value_pairs = info_strings.str.split(";").explode()

        # Split by equals sign to separate keys and values, then take the keys
        keys = key_value_pairs[key_value_pairs.str.contains("=")].str.split("=").str[0]

        # Find unique keys
        unique_keys = [re.sub(r"\W", "_", key) for key in keys.unique()]

        unique_keys = [f"_{key}" if key[0].isdigit() else key for key in keys.unique()]
        ic(unique_keys)

        info_metadata = [
            f'##INFO=<ID={tag},Number=.,Type=String,Description="Poop scoop">'
            if tag != "."
            else ""
            for tag in unique_keys
        ]
        metadata.extend(info_metadata)
        self.metadata = "\n".join(metadata)
        return

    def filter_gnomAD_variants(
        self,
        output_file_path=None,
        output_format="z",
        update_file=True,
        update_vcf_df=True,
    ):
        if output_file_path is None:
            with tempfile.NamedTemporaryFile(
                mode="w+", dir="/tmp", delete=False, suffix=".vcf.gz"
            ) as temp_output_file:
                output_file_path = temp_output_file.name
                filter_cmd = f"bcftools view -i 'FILTER=\"PASS\" && INFO/AF > 0 && INFO/AF < 1 && INFO/AN >= 76156' {self.vcf_file_path} -O{output_format} -o {output_file_path} --write-index"
                utils.run_shell_command(filter_cmd)
                temp_output_file.flush()
        else:
            filter_cmd = f"bcftools view -i 'FILTER=\"PASS\" && INFO/AF > 0 && INFO/AF < 1 && INFO/AN >= 76156' {self.vcf_file_path} -O{output_format} -o {output_file_path} --write-index"
            utils.run_shell_command(filter_cmd)

        if update_file:
            mv_cmd = f"mv {output_file_path} {self.vcf_file_path}"
            utils.run_shell_command(mv_cmd)

        if update_vcf_df:
            self.read_vcf(self.vcf_file_path)

        return

    # TODO: make this sequentially check against every chr file and then concat the ones that pass for each chr
    def filter_against_gnomAD_QC(
        self,
        output_file_path=None,
        output_format="z",
        update_file=True,
        update_vcf_df=True,
    ):
        # if output_file_path is None:
        gnomAD_QC_bcf_obj = worktools.core.large_vcf_object.LargeVCFObject(
            "/gpfs/gibbs/pi/reilly/VariantEffects/results/rohit_results/gnomad_snps_background_nonexonic_pass_QC.bcf"
        )
        # with tempfile.NamedTemporaryFile(mode='w+', dir='/tmp', delete=False, suffix='.vcf.gz') as temp_output_file:
        api.get_variants_in_common(
            self,
            gnomAD_QC_bcf_obj,
            output_file_path=output_file_path,
            output_format=output_format,
        )

        # if update_file:
        #     mv_cmd = f"mv {output_file_path} {self.vcf_file_path}"
        #     utils.run_shell_command(mv_cmd)

        # if update_vcf_df:
        #     self.read_vcf(self.vcf_file_path)
        return

    def is_tag_in_file(self, tag):
        query_cmd = f"bcftools query -i 'INFO/{tag} != \".\"' -f '%INFO/{tag}\n' {self.vcf_file_path} | wc -l"
        num_lines = utils.run_shell_command(query_cmd)
        return num_lines > 0

    def transform_filter_expression_with_abs(self, input_expression):
        # Replace 'and' and 'or' with their bcftools equivalents '&&' and '||'
        transformed_expression = input_expression.replace(" and ", " && ").replace(
            " or ", " || "
        )

        def repl(match):
            tag = match.group(1)  # Capture the tag name correctly
            operator = match.group(2)
            # Adjust for the correct detection of 'abs' usage
            abs_prefix = "ABS(" if "abs" in match.group(0) else ""
            abs_suffix = ")" if "abs" in match.group(0) else ""
            # Prefix INFO tags correctly
            return f"{abs_prefix}INFO/{tag}{abs_suffix} {operator} "

        # Adjusted regex to match 'abs' function usage and INFO tags correctly
        transformed_expression = re.sub(
            r"abs\('([a-zA-Z0-9_]+)'\)\s*(>=|<=|!=|==|>|<)",
            repl,
            transformed_expression,
        )
        # Transform non-abs prefixed tags
        transformed_expression = re.sub(
            r"'([a-zA-Z0-9_]+)'\s*(>=|<=|!=|==|>|<)",
            lambda m: f"INFO/{m.group(1)} {m.group(2)} ",
            transformed_expression,
        )

        # Extract all unique tags mentioned in the transformed_expression
        all_tags = set(re.findall(r"INFO/([a-zA-Z0-9_]+)", transformed_expression))

        if all_tags:
            self.modify_header_for_tags(all_tags)

        return transformed_expression

    def modify_header_for_tags(self, tags):
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_header:
            # Export the current header to the temporary file
            subprocess.run(
                f"bcftools view -h {self.vcf_file_path} > {temp_header.name}",
                shell=True,
                check=True,
            )
            temp_header.flush()
            # os.system(f"cat {temp_header.name}")
            # Read the header from the start
            temp_header.seek(0)
            lines = temp_header.readlines()

            with tempfile.NamedTemporaryFile(
                mode="w+", delete=False
            ) as real_temp_header:
                for line in lines:
                    modified_line = line
                    for tag in tags:
                        if (
                            line.startswith(f"##INFO=<ID={tag},")
                            and "Type=String" in line
                        ):
                            modified_line = line.replace("Type=String", "Type=Float")
                    real_temp_header.write(modified_line)
                    real_temp_header.flush()

                os.system(f"mv {real_temp_header.name} {temp_header.name}")

            # Apply the modified header
            with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_vcf:
                subprocess.run(
                    f"bcftools reheader -h {temp_header.name} -o {temp_vcf.name} {self.vcf_file_path}",
                    shell=True,
                    check=True,
                )

                # Update the VCF file path to the modified file
                os.system(f"mv {temp_vcf.name} {self.vcf_file_path}")
                self.index_variant_file()
        return

    # Example usage: vcf_obj.filter_by_info_tag_condition(['AF > 0.01' or 'AN > 1000' and 'AC < 100'])
    # Example usage: vcf_obj.filter_by_info_tag_condition("(abs('K562__ref') > 0.5 or abs('K562__alt') > 0.5) and 'K562__skew' > 0.1")
    def filter_by_info_tag_condition(
        self,
        filter_string,
        output_file_path=None,
        output_format="z",
        update_file=False,
        update_vcf_df=False,
    ):
        condition_string = self.transform_filter_expression_with_abs(filter_string)
        ic(condition_string)
        if output_file_path is None:
            with tempfile.NamedTemporaryFile(
                mode="w+", dir="/tmp", delete=False, suffix=".vcf.gz"
            ) as temp_output_file:
                output_file_path = temp_output_file.name
                filter_cmd = f"bcftools view -i '{condition_string}' {self.vcf_file_path} -O{output_format} -o {output_file_path} --write-index"
                utils.run_shell_command(filter_cmd)
                self.vcf_file_path = output_file_path
                temp_output_file.flush()
        else:
            filter_cmd = f"bcftools view -i '{condition_string}' {self.vcf_file_path} -O{output_format} -o {output_file_path}  --write-index"
            utils.run_shell_command(filter_cmd)
            # if output_format == "z":
            #     pysam.tabix_index(output_file_path, force=True, preset="vcf")
            # elif self.vcf_file_path.endswith(".bcf"):
            #     index_cmd = f"bcftools index -f {output_file_path}"
            #     utils.run_shell_command(index_cmd)

        if update_file:
            mv_cmd = f"mv {output_file_path} {self.vcf_file_path}"
            utils.run_shell_command(mv_cmd)

        if update_vcf_df:
            self.read_vcf(self.vcf_file_path)
        pass

    def get_promoter_variants(
        self, output_file_path=None, update_file=True, update_vcf_df=True
    ):
        return api.get_promoter_variants(
            self, output_file_path, update_file, update_vcf_df
        )

    def get_dels_variants(
        self, output_file_path=None, update_file=True, update_vcf_df=True
    ):
        return api.get_dels_variants(self, output_file_path, update_file, update_vcf_df)

    def get_pels_variants(
        self, output_file_path=None, update_file=True, update_vcf_df=True
    ):
        return api.get_pels_variants(self, output_file_path, update_file, update_vcf_df)

    def get_ccre_variants(
        self, output_file_path=None, update_file=True, update_vcf_df=True
    ):
        return api.get_ccre_variants(self, output_file_path, update_file, update_vcf_df)

    def get_variants_within_250_bp_of_tss(
        self, output_file_path=None, update_file=True, update_vcf_df=True
    ):
        return api.get_variants_within_250_bp_of_tss(
            self, output_file_path, update_file, update_vcf_df
        )

    def get_variants_within_500_bp_of_tss(
        self, output_file_path=None, update_file=True, update_vcf_df=True
    ):
        return api.get_variants_within_500_bp_of_tss(
            self, output_file_path, update_file, update_vcf_df
        )

    def get_variants_within_750_bp_of_tss(
        self, output_file_path=None, update_file=True, update_vcf_df=True
    ):
        return api.get_variants_within_750_bp_of_tss(
            self, output_file_path, update_file, update_vcf_df
        )

    def get_variants_within_1kb_of_tss(
        self, output_file_path=None, update_file=True, update_vcf_df=True
    ):
        return api.get_variants_within_1kb_of_tss(
            self, output_file_path, update_file, update_vcf_df
        )

    def get_exonic_variants(
        self, output_file_path=None, update_file=True, update_vcf_df=True
    ):
        return api.get_exonic_variants(
            self, output_file_path, update_file, update_vcf_df
        )

    def get_non_exonic_variants(
        self, output_file_path=None, output_file_format="z", update_file=True, update_vcf_df=True
    ):
        return api.get_non_exonic_variants(
            self, output_file_path, output_file_format, update_file, update_vcf_df
        )

    def filter_by_ccre(
        self,
        ccre,
        output_file_path=None,
        output_file_format="z",
        update_file=True,
        update_vcf_df=True,
    ):
        return api.filter_by_ccre(
            vcf_obj=self,
            ccre=ccre,
            output_file_path=output_file_path,
            output_file_format=output_file_format,
            update_file=update_file,
            update_vcf_df=update_vcf_df,
        )

    def remove_non_snps(self, output_file_path=None):
        self.vcf_df = self.vcf_df.loc[
            (self.vcf_df["REF"].str.len() == 1) & (self.vcf_df["ALT"].str.len() == 1)
        ]
        return

    def update_vcf(self, file_path):
        self.vcf_file_path = file_path
        self.read_vcf(file_path)
        return

    def flush_to_temp_file(self, output_file_format="b"):
        suffix = None
        bcf = False
        tabix = False

        match output_file_format:
            case "z":
                suffix = ".vcf.gz"
                tabix = True
            case "b":
                suffix = ".bcf"
                bcf = True
            case "v":
                suffix = ".vcf"

        with tempfile.NamedTemporaryFile(
            mode="w+", dir="/tmp", delete=False, suffix=suffix
        ) as temp_file:
            self.write_to_file(
                temp_file.name,
                write_metadata=True,
                bgzip=False,
                tabix=tabix,
                do_not_generate_metadata=False,
                bcf=bcf,
            )
            temp_file.flush()

        # utils.run_shell_command(f'bcftools view {temp_file.name}')
        self.vcf_file_path = temp_file.name
        # self.read_vcf(temp_file.name)
        # ic(self.metadata)
        # ic(self.vcf_df)
        # self.__init__(vcf_file_path = temp_file.name, vcf_df = None, read_file = True)

        return temp_file.name

    def read_metadata(self):
        grep_zgrep = "zgrep" if utils.is_bgzipped(self.vcf_file_path) else "grep"

        if self.vcf_file_path.endswith(".bcf"):
            metadata_cmd = f"bcftools view -h {self.vcf_file_path} | grep '^##'"
            self.metadata = utils.run_shell_command(metadata_cmd)
            return

        metadata_cmd = f"{grep_zgrep} '^##' {self.vcf_file_path}"
        self.metadata = utils.run_shell_command(metadata_cmd)

        return

    def get_num_variants(self, from_df=True):
        if from_df:
            return len(self.vcf_df)

        if self.vcf_file_path.endswith(".bcf"):
            num_variants = utils.run_shell_command(
                f"bcftools view {self.vcf_file_path} | grep -v '^#' | wc -l"
            )
        elif self.vcf_file_path.endswith(".vcf.gz"):
            num_variants = utils.run_shell_command(
                f"zcat {self.vcf_file_path} | grep -v '^#' | wc -l"
            )
        else:
            num_variants = utils.run_shell_command(
                f"grep -v '^#' {self.vcf_file_path} | wc -l"
            )
        return int(num_variants)

    def get_duplicate_identical_pos_df(self):
        return self.vcf_df[self.vcf_df.duplicated(subset=["#CHROM", "POS"], keep=False)]
