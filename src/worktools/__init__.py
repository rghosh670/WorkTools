"""Top-level package for worktools."""

__author__ = """Rohit Ghosh"""  # Indicates the author of the package
__email__ = "rohit.ghosh@yale.edu"  # Contact email for the author
__version__ = "0.1.0"  # Current version of the package

import os

# IceCream is used for debugging: it prints variable names and their values automatically
from icecream import ic, install

install()

# Importing specific functions from the worktools API module
from worktools.api import (adjust_chr_prefix_in_bed, adjust_chr_prefix_in_vcf,
                           annotate_vcf_with_info_from_other,
                           annotate_vcf_with_info_from_other_any_allele,
                           convert_bed_to_vcf, convert_vcf_to_bed,
                           fishers_test_decile_overlap_enrichment,
                           fishers_test_fimo_motif_enrichment,
                           fishers_test_overlap_enrichment,
                           fishers_test_trait_enrichment, get_fasta_from_vcf,
                           get_fimo_motif_changes, get_genome_build_bed,
                           get_hocomoco_motif_lengths, get_malinois_scores,
                           get_variants_in_common,
                           index_all_variant_files_in_dir,
                           keep_variants_in_bed, lift_over_bed, lift_over_vcf,
                           move_columns_from_one_vcf_to_another,
                           plot_prc_curves, remove_variants_in_bed,
                           run_fimo_with_hocomoco)
# Importing data structures used across the worktools package
from worktools.core import BEDObject, FASTAObject, LargeVCFObject, VCFObject

# Configure IceCream to always include context in its output, aiding in debugging
ic.configureOutput(includeContext=True)

# Defines a list of all publicly accessible objects in the module
__all__ = [
    "BEDObject",
    "FASTAObject",
    "VCFObject",
    "LargeVCFObject",
    "adjust_chr_prefix_in_bed",
    "adjust_chr_prefix_in_vcf",
    "convert_bed_to_vcf",
    "convert_vcf_to_bed",
    "get_fasta_from_vcf",
    "run_fimo_with_hocomoco",
    "get_fimo_motif_changes",
    "move_columns_from_one_vcf_to_another",
    "get_variants_in_common",
    "annotate_vcf_with_info_from_other",
    "annotate_vcf_with_info_from_other_any_allele",
    "get_malinois_scores",
    "get_hocomoco_motif_lengths",
    "keep_variants_in_bed",
    "remove_variants_in_bed",
    "index_all_variant_files_in_dir",
    "fishers_test_trait_enrichment",
    "fishers_test_overlap_enrichment",
    "fishers_test_fimo_motif_enrichment",
    "fishers_test_decile_overlap_enrichment",
    "plot_prc_curves",
]
