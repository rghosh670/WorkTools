import os
import re
import tempfile

import datatable as dt
import kipoiseq
import numpy as np
import pandas as pd
import pysam
from icecream import ic
from kipoiseq import Interval, Variant

from worktools.core import BEDObject, FASTAObject, LargeVCFObject, VCFObject
from worktools.core.fasta_object import FastaStringExtractor

HG38_REFERENCE_GENOME_FILE_PATH = "/home/rg972/project/Data/hg38.fa"


def convert_bed_to_vcf(bed_obj, rt_style_index=None):
    bed_df = bed_obj.bed_df

    if rt_style_index is None:
        if bed_df.index.dtype == np.int64:
            raise NotImplementedError(
                "BED file must have an 'rt-style' index formatted as such: 'chr;pos;ref;alt'"
            )

    vcf_df = pd.Series(rt_style_index.to_list()).str.split(";", expand=True)
    vcf_df.columns = ["#CHROM", "chromEnd", "REF", "ALT"]
    vcf_df.index = rt_style_index
    ic(vcf_df)

    merge_df = pd.merge(vcf_df, bed_df, on="chromEnd", how="right")
    ic(merge_df)
    columns_to_drop = ["#CHROM_y", "chromStart", "C3", "C4", "C5"]

    # Only drop columns that exist in the DataFrame

    ic(bed_df)
    ic(len(rt_style_index))
    bed_df.index = (
        bed_df["#CHROM"]
        + ";"
        + bed_df["chromEnd"].astype(str)
        + ";"
        + bed_df.iloc[:, 2].astype(str)
        + ";"
        + bed_df.iloc[:, 3].astype(str)
    )
    ic(bed_df)
    base_columns = ["#CHROM", "chromStart", "chromEnd"]

    # Find any additional columns that contain the specified strings
    additional_columns = [
        col
        for col in bed_df.columns
        if not any(base_col in col for base_col in base_columns)
    ]

    # Combine the base columns with the additional columns, ensuring uniqueness
    selected_columns = base_columns + list(set(additional_columns) - set(base_columns))
    ic(selected_columns)

    # Use the combined list of columns to select from the DataFrame
    bed_df = bed_df[selected_columns]
    ic(bed_df)

    # bed_df.index = pd.Series(rt_style_index).iloc[matching_indices]

    bed_df.iloc[:, 3:] = bed_df.iloc[:, 3:].astype(str)

    # Construct the column-value pairs in a vectorized manner
    col_names = bed_df.columns[3:]  # Adjust the index as necessary
    col_names = pd.Index([re.sub(r"\W", "_", col) for col in col_names])
    col_names = pd.Index([f"_{col}" if col[0].isdigit() else col for col in col_names])
    col_prefixes = col_names + "="
    col_values = bed_df.iloc[:, 3:]

    # Use numpy broadcasting to prepend column names to their values
    combined = col_prefixes + col_values.to_numpy()

    # Join the column-value pairs for each row and create a Series instead of a list
    info_col_series = [";".join(row) for row in combined]
    ic(info_col_series[:5])

    bed_df = merge_df.drop(
        columns=merge_df.columns.intersection(columns_to_drop), axis=1
    )
    bed_df = bed_df.iloc[:, :4]
    ic(bed_df)

    bed_df.columns = vcf_df.columns

    vcf_df["chromEnd"] = vcf_df["chromEnd"].astype(int)
    vcf_df["merge_key"] = (
        vcf_df["#CHROM"].astype(str) + "_" + vcf_df["chromEnd"].astype(str)
    )
    bed_df["merge_key"] = (
        bed_df["#CHROM"].astype(str) + "_" + bed_df["chromEnd"].astype(str)
    )
    vcf_df = vcf_df[vcf_df["merge_key"].isin(bed_df["merge_key"])]
    vcf_df.drop("merge_key", axis=1, inplace=True)
    bed_df.drop("merge_key", axis=1, inplace=True)

    vcf_df.columns = ["#CHROM", "POS", "REF", "ALT"]

    possible_id_columns = ["name", "id", "ID", "Name", "rsid"]

    # Find the first column that exists in bed_df from the possible_id_columns list
    id_column = next(
        (col for col in possible_id_columns if col in bed_df.columns), None
    )
    # If any of the specified columns are found, set vcf_df['ID'] to it
    if id_column is not None:
        vcf_df["ID"] = bed_df[id_column].values
    else:
        vcf_df["ID"] = ["."] * len(vcf_df)

    vcf_df["QUAL"] = ["."] * len(vcf_df)
    vcf_df["FILTER"] = ["."] * len(vcf_df)
    vcf_df["INFO"] = info_col_series  # ['.'] * len(vcf_df)
    ic(vcf_df)
    vcf_df = vcf_df.dropna()
    ic(vcf_df)
    vcf_df = vcf_df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]]

    ret_vcf = VCFObject(vcf_df=vcf_df)
    ret_vcf.remove_from_info_column([])  # remove invalid INFO values
    ret_vcf.sort_vcf_df()
    ic(ret_vcf)
    return ret_vcf


def convert_vcf_to_bed(vcf_obj, expand_info=False):
    if not isinstance(vcf_obj, LargeVCFObject):
        vcf_obj.add_rt_style_index()
        vcf_df = vcf_obj.vcf_df
        bed_df = vcf_df[["#CHROM", "POS"]].copy()
        bed_df.columns = ["#CHROM", "chromEnd"]
        bed_df["chromStart"] = bed_df["chromEnd"] - 1
        bed_df["name"] = vcf_df["ID"]
        bed_df["score"] = [0] * len(bed_df)
        bed_df["strand"] = ["."] * len(bed_df)
        bed_df = bed_df[["#CHROM", "chromStart", "chromEnd", "name", "score", "strand"]]

        bed_df = bed_df.astype(
            {
                "#CHROM": "str",
                "chromStart": "int",
                "chromEnd": "int",
                "name": "str",
                "score": "float",  # Assuming 'score' should be a floating point number
                "strand": "str",
            }
        )

        bed_df.index = vcf_df.index

        if expand_info:
            pattern = r"(?P<key>\w+__[a-z]+)=(?P<value>[-+]?\d*\.\d+|\d+)"
            info_expanded = vcf_df["INFO"].str.extractall(pattern)

            # Step 2: Pivot the extracted DataFrame to wide format
            info_expanded.reset_index(inplace=True)
            info_wide = info_expanded.pivot(
                index="level_0", columns="key", values="value"
            )

            # Step 3: Convert the 'value' column to numeric
            info_wide = info_wide.apply(pd.to_numeric)

            # Step 4: Join the expanded info back to the original DataFrame
            bed_df = bed_df.join(info_wide)

        return BEDObject(bed_df=bed_df)
    else:
        cmd = f"zgrep -v '^##' {vcf_obj.vcf_file_path}"

        if vcf_obj.vcf_file_path.endswith(".bcf"):
            cmd = f"bcftools view {vcf_obj.vcf_file_path} | grep -v '##'"

        vcf_columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

        columns_dict = dict(
            zip(vcf_columns, (([dt.str32, dt.int32, dt.str32]) + ([None] * 5)))
        )

        bed_df = dt.fread(
            cmd=cmd, columns=columns_dict, sep="\t", header=True
        ).to_pandas()

        ic(bed_df)
        bed_df.columns = ["#CHROM", "chromEnd", "name"]
        bed_df["chromStart"] = (bed_df["chromEnd"] - 1).astype(str)
        ic(bed_df)
        bed_df["score"] = [0] * len(bed_df)
        bed_df["strand"] = ["."] * len(bed_df)
        bed_df = bed_df[["#CHROM", "chromStart", "chromEnd", "name", "score", "strand"]]
        ic(bed_df)

        bed_df = bed_df.astype(
            {
                "#CHROM": "str",
                "chromStart": "int",
                "chromEnd": "int",
                "name": "str",
                "score": "float",  # Assuming 'score' should be a floating point number
                "strand": "str",
            }
        )

        bed_df.index = vcf_obj.get_rt_style_index()
        ic(bed_df)

        if expand_info:
            pattern = r"(?P<key>\w+__[a-z]+)=(?P<value>[-+]?\d*\.\d+|\d+)"
            info_expanded = vcf_df["INFO"].str.extractall(pattern)

            # Step 2: Pivot the extracted DataFrame to wide format
            info_expanded.reset_index(inplace=True)
            info_wide = info_expanded.pivot(
                index="level_0", columns="key", values="value"
            )

            # Step 3: Convert the 'value' column to numeric
            info_wide = info_wide.apply(pd.to_numeric)

            # Step 4: Join the expanded info back to the original DataFrame
            bed_df = bed_df.join(info_wide)

        ic(bed_df)
        return BEDObject(bed_df=bed_df)


def get_fasta_from_vcf(
    vcf_obj, genome_fasta_file=HG38_REFERENCE_GENOME_FILE_PATH, sequence_length=30
):
    vcf_bcf_file_path = vcf_obj.vcf_file_path
    if os.path.getsize(vcf_bcf_file_path) == 0:
        raise ValueError(f"VCF/BCF file {vcf_bcf_file_path} is empty.")

    seq_extractor = kipoiseq.extractors.VariantSeqExtractor(
        reference_sequence=FastaStringExtractor(genome_fasta_file)
    )

    reference_sequences = {}
    alternate_sequences = {}

    # Using pysam to open both VCF and BCF files
    with pysam.VariantFile(vcf_bcf_file_path) as vcf_bcf:
        for record in vcf_bcf:
            chrom = record.chrom
            pos = record.pos
            id = record.id
            ref = record.ref
            alts = record.alts  # alts is a tuple of alternate alleles

            for alt in alts:
                interval = Interval(chrom, pos, pos)
                interval = interval.resize(sequence_length)
                center = interval.center() - interval.start

                variant = Variant(chrom=chrom, pos=pos, ref=ref, alt=alt, id=id)
                reference_sequences[f">{interval.chrom}_{interval.end}"] = (
                    seq_extractor.extract(interval, [], anchor=center)
                )
                alternate_sequences[f">{interval.chrom}_{interval.end}"] = (
                    seq_extractor.extract(interval, [variant], anchor=center)
                )

    # Temporary file handling for reference and alternate sequences remains the same

    with tempfile.NamedTemporaryFile(
        mode="w+", dir="/tmp", delete=False, suffix=".fasta"
    ) as ref_temp_file:
        for id, seq in reference_sequences.items():
            ref_temp_file.write(f"{id}\n{seq}\n")
        ref_temp_file.flush()
        ref_fasta_obj = FASTAObject(fasta_file_path=ref_temp_file.name)

    with tempfile.NamedTemporaryFile(
        mode="w+", dir="/tmp", delete=False, suffix=".fasta"
    ) as alt_temp_file:
        for id, seq in alternate_sequences.items():
            alt_temp_file.write(f"{id}\n{seq}\n")

        alt_temp_file.flush()
        alt_fasta_obj = FASTAObject(fasta_file_path=alt_temp_file.name)

    return ref_fasta_obj, alt_fasta_obj
