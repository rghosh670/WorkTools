import logging
import tempfile

import datatable as dt
import numpy as np
import pandas as pd
from icecream import ic

import worktools.api.utils as utils

LIFTOVER_SCRIPT_FILE_PATH = "/home/rg972/project/Downloaded_Tools/ucsc-tools/liftOver"
HG19_TO_HG38_CHAIN_FILE_PATH = (
    "/home/rg972/project/Downloaded_Tools/ucsc-tools/hg19ToHg38.over.chain.gz"
)

logger = logging.getLogger("worktools")


def adjust_chr_prefix_in_bed(bed_obj, action):
    """
    Adjusts the 'chr' prefix in the '#CHROM' column of a BED dataframe.

    Parameters:
    - action: A string specifying the action to perform: 'add' to add 'chr' prefix, 'remove' to remove it.
    """

    match action:
        case "add":
            if not bed_obj.has_chr_prefix():
                bed_obj.bed_df["#CHROM"] = "chr" + bed_obj.bed_df["#CHROM"]
        case "remove":
            if bed_obj.has_chr_prefix():
                bed_obj.bed_df["#CHROM"] = bed_obj.bed_df["#CHROM"].str.lstrip("chr")
        case _:
            raise ValueError("Invalid action specified. Choose 'add' or 'remove'.")

    return


def get_genome_build_bed(bed_obj):
    if not variant_ids.str.startswith("rs").any():
        raise NotImplementedError(
            "Getting the genome build without any variant with rsIDs is not yet supported."
        )

    variant_ids = bed_obj.bed_df.iloc[:, 5]
    query_variant_index = variant_ids[
        variant_ids.str.startswith("rs")
    ].first_valid_index()
    query_variant_id = variant_ids[query_variant_index]

    try:
        handle = Entrez.efetch(db="snp", id="query_variant_id", retmode="text")
        snp_data = handle.read()

        hg38_match = re.search(r"<CHRPOS>(.*?)</CHRPOS>", snp_data)
        hg19_match = re.search(r"<CHRPOS_PREV_ASSM>(.*?)</CHRPOS_PREV_ASSM>", snp_data)

        if hg38_match:
            grch38_match = hg38_match.group(1)
            grch38_chrom = grch38_match.split(":")[0]
            grch38_pos = grch38_match.split(":")[1]
        else:
            logger.error("GRCh38 position not found")
            return

        if hg19_match:
            grch37_match = hg19_match.group(1)
            grch37_chrom = grch37_match.split(":")[0]
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
    query_variant_pos = self.vcf_df["#CHROM"][query_variant_index]

    if query_variant_chrom != grch38_chrom and query_variant_chrom != grch37_chrom:
        logger.error("Chromosome doesn't match rsID")
        return

    if query_variant_pos != grch38_pos and query_variant_pos != grch37_pos:
        logger.error("Position doesn't match hg38 or hg19")
        return

    if query_variant_pos == grch38_pos:
        return "GRCh38"

    return "GRCh37"


def lift_over_bed(bed_obj, genome_build):
    raise NotImplementedError(
        "This function is not yet implemented. Refer to lift_over_vcf for an example of how to implement this."
    )
