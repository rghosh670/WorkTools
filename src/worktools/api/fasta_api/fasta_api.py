import logging
import os
import tempfile
from itertools import chain
from pprint import pprint

import datatable as dt
import pandas as pd
from icecream import ic

import worktools.api.utils as utils

FIMO_SCRIPT_FILE_PATH = (
    "/home/rg972/project/Worktools/src/worktools/lib/motif_liquidator.sh"
)
HOCOMOCO_FILE_PATH = "/home/rg972/project/Worktools/src/worktools/data/HOCOMOCOv11_core_HUMAN_mono_meme_format.txt"

logger = logging.getLogger("worktools")


def run_fimo_with_hocomoco(
    fasta_obj, filter_by_pval=True, pval_filter=1e-5, output_file_path=None
):
    if not fasta_obj.fasta_file_path.endswith(".fasta"):
        raise ValueError(
            f"Invalid file extension for fasta file: {fasta_obj.fasta_file_path}"
        )

    if output_file_path is None:
        # Create a temporary file and set output_file_path to its path
        temp_output_file = tempfile.NamedTemporaryFile(
            delete=False, suffix=".txt", mode="w"
        )
        output_file_path = temp_output_file.name
        temp_output_file.close()  # It's important to close the file to flush and allow other processes to access it
    output_file_path = os.path.abspath(output_file_path)
    logger.info(
        f"Running FIMO with HOCOMOCO on {fasta_obj.fasta_file_path} and writing results to {output_file_path}"
    )
    # Construct and run the command to execute the motif liquidator script
    motif_liq_command = (
        f"{FIMO_SCRIPT_FILE_PATH} {fasta_obj.fasta_file_path} {output_file_path}"
    )
    utils.run_shell_command(motif_liq_command)
    logger.debug(f"Running command: {motif_liq_command}")

    if filter_by_pval:
        with tempfile.NamedTemporaryFile(
            mode="w+", dir="/tmp", delete=False, suffix=".txt"
        ) as temp_file:
            filter_cmd = f"mawk 'BEGIN{{FS=OFS=\"\\t\"}} $7 < {pval_filter} || NR == 1' {output_file_path} > {temp_file.name}"
            logger.debug(f"Running command: {filter_cmd}")
            utils.run_shell_command(filter_cmd)
            temp_file.flush()
            logger.info(
                f"Removed {utils.get_num_lines_in_file(output_file_path) - utils.get_num_lines_in_file(temp_file.name)} FIMO results with p-value >= {pval_filter}"
            )
            mv_cmd = f"mv {temp_file.name} {output_file_path}"
            logger.debug(f"Running command: {mv_cmd}")
            utils.run_shell_command(mv_cmd)

    # Load the results into a pandas DataFrame
    fimo_result_df = dt.fread(output_file_path, sep="\t", header=True).to_pandas()
    logger.debug(
        f"Loaded FIMO results into DataFrame with shape {fimo_result_df.shape}"
    )
    fimo_result_df.columns = fimo_result_df.columns.str.replace(" ", "_")

    # Pre-compute fasta sequences to avoid repetitive access
    fasta_sequences = {
        name: fasta_obj.get_sequence_from_fasta_id(name)
        for name in fimo_result_df["sequence_name"].unique()
    }
    fimo_result_df["fasta_sequence"] = fimo_result_df["sequence_name"].map(
        fasta_sequences
    )
    fimo_result_df["start"] = fimo_result_df.apply(
        lambda row: str(row["fasta_sequence"]).find(str(row["matched_sequence"])),
        axis=1,
    )
    fimo_result_df["end"] = (
        fimo_result_df["start"] + fimo_result_df["matched_sequence"].str.len()
    )
    fimo_result_df["middle_pos"] = fimo_result_df["fasta_sequence"].str.len() // 2
    fimo_result_df = fimo_result_df[
        (fimo_result_df["start"] <= fimo_result_df["middle_pos"])
        & (fimo_result_df["middle_pos"] <= fimo_result_df["end"])
    ]

    fimo_result_df = fimo_result_df.drop(
        [
            "q-value",
            "fasta_sequence",
            "start",
            "end",
            "middle_pos",
        ],
        axis="columns",
    )
    # Extract the gene name from the motif pattern name
    fimo_result_df["#pattern_name"] = (
        fimo_result_df["#pattern_name"].str.split("_HUMAN").str[0]
    )
    # Write the possibly filtered and processed DataFrame back to a file
    fimo_result_df.to_csv(output_file_path, sep="\t", index=False)
    logger.info(f"Wrote FIMO results to {output_file_path}")
    fasta_obj.fimo_result_df = fimo_result_df
    return fimo_result_df


def get_fimo_motif_changes(ref_fasta_obj, alt_fasta_obj):
    if ref_fasta_obj.fimo_result_df is None:
        run_fimo_with_hocomoco(ref_fasta_obj)

    if alt_fasta_obj.fimo_result_df is None:
        run_fimo_with_hocomoco(alt_fasta_obj)

    ref_fimo_df = ref_fasta_obj.fimo_result_df
    alt_fimo_df = alt_fasta_obj.fimo_result_df

    variant_rows = []
    for variant, ref_variant_df in ref_fimo_df.groupby(by="sequence_name"):
        alt_variant_df = alt_fimo_df[alt_fimo_df["sequence_name"] == variant]
        ref_motifs = set(ref_variant_df["#pattern_name"])
        alt_motifs = set(alt_variant_df["#pattern_name"])
        variant_row = pd.Series(
            data={
                "ref_motifs": ref_motifs,
                "alt_motifs": alt_motifs,
                "unchanged": ref_motifs.intersection(alt_motifs),
                "disrupted": ref_motifs - alt_motifs,
                "gained": alt_motifs - ref_motifs,
            },
            name=variant,
        )
        variant_rows.append(variant_row)

    summary_df = pd.concat(variant_rows, axis=1).T
    disrupted_motifs = pd.Series(
        chain(*[list(motif_set) for motif_set in summary_df["disrupted"]])
    )
    disrupted_motif_proportions = disrupted_motifs.value_counts() / len(
        disrupted_motifs
    )

    gained_motifs = pd.Series(
        chain(*[list(motif_set) for motif_set in summary_df["gained"]])
    )
    gained_motif_proportions = gained_motifs.value_counts() / len(gained_motifs)
    return summary_df, disrupted_motif_proportions, gained_motif_proportions
