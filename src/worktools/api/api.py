import ast
import os
import re
import tempfile
from copy import deepcopy

import datatable as dt
import numpy as np
import pandas as pd
from icecream import ic
from scipy.stats import fisher_exact

from worktools.api.conversion import (convert_bed_to_vcf, convert_vcf_to_bed)
from worktools.api.vcf_api import (get_variants_in_common,
                           move_columns_from_one_vcf_to_another)
import worktools.api.utils as utils
from worktools.core import BEDObject, VCFObject

HOCOMOCO_FILE_PATH = "/home/rg972/project/WorkTools/worktools/data/HOCOMOCOv11_core_HUMAN_mono_meme_format.txt"
GENOME_FILE_PATH = "/home/rg972/project/Data/hg38.genome"
BACKGROUND_VARS_TOTAL = 536377520


def get_hocomoco_motif_lengths():
    grep_cmd = f"grep -B 1 'w=' {HOCOMOCO_FILE_PATH}"
    grep_cmd_output = utils.run_shell_command(grep_cmd)
    # Define a pattern to capture the motif name and the 'w' value
    pattern = re.compile(r"MOTIF (\S+).*?w= (\d+)", re.DOTALL)

    # Find all matches and create the dictionary
    motif_lengths = pd.Series(
        {
            match.group(1): int(match.group(2))
            for match in pattern.finditer(grep_cmd_output)
        }
    )
    motif_lengths = motif_lengths.sort_values(ascending=False)
    return motif_lengths


def remove_keep_helper(
    vcf_obj, bed_obj, keep=True, output_file_path=None, output_file_format="z"
):
    ic(output_file_path)
    vcf_file_path = vcf_obj.vcf_file_path
    bed_file_path = bed_obj.bed_file_path
    keep_remove_string = "-v" if not keep else ""

    bed_obj.sort_bed_df()
    ic(vcf_obj)
    ic(bed_obj)

    with tempfile.NamedTemporaryFile(
        mode="w+", dir="/tmp", delete=False, suffix=".bed"
    ) as temp_bed_file:
        bed_obj.write_to_file(temp_bed_file.name)
        temp_bed_file.flush()
        bed_file_path = temp_bed_file.name

    vcf_to_bed_obj = convert_vcf_to_bed(vcf_obj)
    vcf_to_bed_obj.sort_bed_df()

    with tempfile.NamedTemporaryFile(
        mode="w+", dir="/tmp", delete=False, suffix=".bed"
    ) as temp_vcf_bed_file:
        vcf_to_bed_obj.write_to_file(temp_vcf_bed_file.name)
        temp_vcf_bed_file.flush()
        vcf_to_bed_file_path = temp_vcf_bed_file.name

    # TODO - MAKE THIS SORTED AGAIN
    with tempfile.NamedTemporaryFile(
        mode="w+", dir="/tmp", delete=False, suffix=".bed"
    ) as temp_output_file:
        if keep:
            intersect_cmd = f"bedtools intersect -a {vcf_to_bed_file_path} -b {bed_file_path} -g {GENOME_FILE_PATH} > {temp_output_file.name}"
        else:
            intersect_cmd = f"bedtools intersect -a {vcf_to_bed_file_path} -b {bed_file_path} -g {GENOME_FILE_PATH} -v > {temp_output_file.name}"
        utils.run_shell_command(intersect_cmd)

        # os.system(f"cp {vcf_to_bed_file_path} /gpfs/gibbs/pi/reilly/VariantEffects/vcf_to_bed.bed")
        # os.system(f"cp {bed_file_path} /gpfs/gibbs/pi/reilly/VariantEffects/bed.bed")
        # os.system(f"cp {temp_output_file.name} /gpfs/gibbs/pi/reilly/VariantEffects/output.bed")
        # os.system(intersect_cmd)
        temp_output_file.flush()
        bed_output_file_path = temp_output_file.name

    ic(vcf_to_bed_obj)
    ic(bed_obj)
    ic(bed_output_file_path)
    if os.path.getsize(bed_output_file_path) == 0:
        os.system(f"> {output_file_path}")
        vcf_obj.write_metadata_to_file(output_file_path, do_not_generate_metadata=True)
        with open(output_file_path, "w") as file:
            file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        return
    output_bed_obj = BEDObject(
        bed_file_path=bed_output_file_path,
        bed_cols={"#CHROM": 0, "chromStart": 1, "chromEnd": 2},
    )
    output_vcf_obj = convert_bed_to_vcf(
        output_bed_obj, rt_style_index=vcf_to_bed_obj.bed_df.index
    )

    ic(output_vcf_obj)
    move_columns_from_one_vcf_to_another(
        vcf_obj,
        output_vcf_obj,
        output_file_format="b",
        output_file_path=(bed_output_file_path.rstrip(".bed") + ".bcf"),
    )
    ic(output_vcf_obj)

    if output_file_path is None:
        with tempfile.NamedTemporaryFile(
            mode="w+", dir="/tmp", delete=False, suffix=".vcf"
        ) as temp_output_file:
            output_file_path = temp_output_file.name

    ic(output_vcf_obj)
    ic(output_file_path)
    # bcf = output_file_path.endswith('.bcf')
    # tabix = output_file_path.endswith('.vcf.gz')
    # ic(bcf,tabix)
    
    bcf = False
    tabix = True

    if output_file_format == "z":
        tabix = True
        bcf = False
    elif output_file_format == "b":
        tabix = False
        bcf = True

    output_vcf_obj.write_to_file(output_file_path, bcf=bcf, tabix=tabix)

    # output_file_path SHOULD BE BCF FILE!

    return output_file_path


def keep_variants_in_bed(
    vcf_obj, bed_obj, output_file_path=None, output_file_format="z"
):
    return remove_keep_helper(
        vcf_obj=vcf_obj,
        bed_obj=bed_obj,
        keep=True,
        output_file_path=output_file_path,
        output_file_format=output_file_format,
    )


def remove_variants_in_bed(
    vcf_obj, bed_obj, output_file_path=None, output_file_format="z"
):
    return remove_keep_helper(
        vcf_obj=vcf_obj,
        bed_obj=bed_obj,
        keep=False,
        output_file_path=output_file_path,
        output_file_format=output_file_format,
    )


def index_all_variant_files_in_dir(dir_path):
    vcf_files = [
        os.path.join(dir_path, file)
        for file in os.listdir(dir_path)
        if file.endswith((".vcf", ".vcf.gz"))
    ]
    bcf_files = [
        os.path.join(dir_path, file)
        for file in os.listdir(dir_path)
        if file.endswith(".bcf")
    ]

    for vcf_file in vcf_files:
        vcf_file_path = os.path.join(dir_path, vcf_file)
        utils.run_shell_command(f"tabix -f -p vcf {vcf_file_path}")

    for bcf_file in bcf_files:
        vcf_file_path = os.path.join(dir_path, vcf_file)
        utils.run_shell_command(f"bcftools index -f {bcf_file_path}")
    return


def compute_trait_enrichment(set_anno_counts, background_anno_counts, trait):
    set_trait_count = set_anno_counts.get(trait, 0)
    background_trait_count = background_anno_counts.get(trait, 0)
    set_not_trait_count = sum(set_anno_counts) - set_trait_count
    background_not_trait_count = sum(background_anno_counts) - background_trait_count

    # ic(trait, set_trait_count, background_trait_count, set_not_trait_count, background_not_trait_count)
    if any(
        x == 0
        for x in [
            set_trait_count,
            background_trait_count,
            set_not_trait_count,
            background_not_trait_count,
        ]
    ):
        return np.nan, np.nan, np.nan, np.nan, np.nan, False

    oddsratio, pvalue = fisher_exact(
        [
            [set_trait_count, set_not_trait_count],
            [background_trait_count, background_not_trait_count],
        ]
    )
    lower_ci_bound = np.exp(
        np.log(oddsratio)
        - 1.96
        * np.sqrt(
            (1 / set_trait_count)
            + (1 / background_trait_count)
            + (1 / set_not_trait_count)
            + (1 / background_not_trait_count)
        )
    )
    upper_ci_bound = np.exp(
        np.log(oddsratio)
        + 1.96
        * np.sqrt(
            (1 / set_trait_count)
            + (1 / background_trait_count)
            + (1 / set_not_trait_count)
            + (1 / background_not_trait_count)
        )
    )
    confidence_interval = (
        f"({'{:.3f}'.format(lower_ci_bound)}, {'{:.3f}'.format(upper_ci_bound)})"
    )
    pvalue = format(pvalue, ".2e")
    significant = lower_ci_bound > 1
    # ic(oddsratio, pvalue, lower_ci_bound, upper_ci_bound, confidence_interval, significant)
    return (
        oddsratio,
        pvalue,
        lower_ci_bound,
        upper_ci_bound,
        confidence_interval,
        significant,
    )


def fishers_test_trait_enrichment(
    vcf, background_vcf, col_name, pip_cutoff=None, pip_col_name=None
):
    info_df = vcf.get_info_col_df()
    background_info_df = background_vcf.get_info_col_df(cache=True)

    if pip_cutoff is not None:
        info_df = info_df[info_df[pip_col_name].astype(float) > pip_cutoff]
        background_info_df = background_info_df[
            background_info_df[pip_col_name].astype(float) > pip_cutoff
        ]

    set_anno_counts = info_df[col_name].value_counts()
    background_anno_counts = background_info_df[col_name].value_counts()

    unique_traits = set_anno_counts.index.union(background_anno_counts.index)

    fishers_test_res = [
        compute_trait_enrichment(set_anno_counts, background_anno_counts, trait)
        for trait in unique_traits
    ]

    df_cols = [
        "oddsratio",
        "pvalue",
        "lower_ci_bound",
        "upper_ci_bound",
        "confidence_interval",
        "significant",
    ]
    fishers_test_df = pd.DataFrame(
        fishers_test_res, columns=df_cols, index=unique_traits
    )
    # fishers_test_df["pvalue"] = fishers_test_df["pvalue"].round(5)
    fishers_test_df = fishers_test_df.sort_values("pvalue")
    return fishers_test_df


def get_counts(info_df_col):
    info_df_col = info_df_col.apply(ast.literal_eval)
    return pd.Series([item for subset in info_df_col for item in subset]).value_counts()


def fishers_test_fimo_motif_enrichment(vcf, background_vcf):
    info_df = vcf.get_info_col_df(tags=["disrupted", "gained"])

    if info_df.empty:
        ic("INFO DF EMPTY!")
        return pd.DataFrame(), pd.DataFrame()

    background_info_df = background_vcf.get_info_col_df(tags=["disrupted", "gained"])

    ic(info_df)
    ic(background_info_df)
    # disrupted
    set_tf_counts = get_counts(info_df["disrupted"])
    background_tf_counts = get_counts(background_info_df["disrupted"])

    unique_tfs = set_tf_counts.index.union(background_tf_counts.index)
    ic(unique_tfs)
    ic(len(unique_tfs))

    fishers_test_res = [
        compute_trait_enrichment(set_tf_counts, background_tf_counts, tf)
        for tf in unique_tfs
    ]

    df_cols = [
        "oddsratio",
        "pvalue",
        "lower_ci_bound",
        "upper_ci_bound",
        "confidence_interval",
        "significant",
    ]
    fishers_test_df = pd.DataFrame(fishers_test_res, columns=df_cols, index=unique_tfs)
    # fishers_test_df["pvalue"] = fishers_test_df["pvalue"].round(5)
    disrupted_fishers_test_df = fishers_test_df.sort_values("pvalue")

    # gained
    set_tf_counts = get_counts(info_df["gained"])
    background_tf_counts = get_counts(background_info_df["gained"])
    unique_tfs = set_tf_counts.index.union(background_tf_counts.index)

    fishers_test_res = [
        compute_trait_enrichment(set_tf_counts, background_tf_counts, tf)
        for tf in unique_tfs
    ]

    df_cols = [
        "oddsratio",
        "pvalue",
        "lower_ci_bound",
        "upper_ci_bound",
        "confidence_interval",
        "significant",
    ]
    fishers_test_df = pd.DataFrame(fishers_test_res, columns=df_cols, index=unique_tfs)
    # fishers_test_df["pvalue"] = fishers_test_df["pvalue"].round(5)
    gained_fishers_test_df = fishers_test_df.sort_values("pvalue")
    return disrupted_fishers_test_df, gained_fishers_test_df


def fishers_test_overlap_enrichment(
    vcf, background_vcf, enriched_for_vcf, pip_cutoff=None, pip_col_name=None
):
    if pip_cutoff is not None:
        with tempfile.NamedTemporaryFile(
            mode="w+", suffix=".bcf", delete=False
        ) as temp_vcf:
            filter_string = f"{pip_col_name} > {pip_cutoff}"
            enriched_for_vcf.filter_by_info_tag_condition(
                self,
                filter_string,
                output_file_path=temp_vcf.name,
                output_format="b",
                update_file=False,
                update_vcf_df=False,
            )
            enriched_for_vcf = VCFObject(temp_vcf.name)

    in_vcf_in_set_oi = get_variants_in_common(
        vcf, enriched_for_vcf, ret_num_variants=True
    )
    in_vcf_not_in_set_oi = vcf.get_num_variants() - in_vcf_in_set_oi

    in_background_in_set_oi = get_variants_in_common(
        background_vcf, enriched_for_vcf, ret_num_variants=True
    )
    in_background_not_in_set_oi = BACKGROUND_VARS_TOTAL - in_background_in_set_oi

    if any(
        x == 0
        for x in [
            in_vcf_in_set_oi,
            in_vcf_not_in_set_oi,
            in_background_in_set_oi,
            in_background_not_in_set_oi,
        ]
    ):
        return np.nan, np.nan, np.nan, np.nan, np.nan, False

    oddsratio, pvalue = fisher_exact(
        [
            [in_vcf_in_set_oi, in_vcf_not_in_set_oi],
            [in_background_in_set_oi, in_background_not_in_set_oi],
        ]
    )
    lower_ci_bound = np.exp(
        np.log(oddsratio)
        - 1.96
        * np.sqrt(
            (1 / in_vcf_in_set_oi)
            + (1 / in_background_in_set_oi)
            + (1 / in_vcf_not_in_set_oi)
            + (1 / in_background_not_in_set_oi)
        )
    )
    upper_ci_bound = np.exp(
        np.log(oddsratio)
        + 1.96
        * np.sqrt(
            (1 / in_vcf_in_set_oi)
            + (1 / in_background_in_set_oi)
            + (1 / in_vcf_not_in_set_oi)
            + (1 / in_background_not_in_set_oi)
        )
    )
    confidence_interval = (
        f"({'{:.3f}'.format(lower_ci_bound)}, {'{:.3f}'.format(upper_ci_bound)})"
    )
    pvalue = format(pvalue, ".2e")
    significant = lower_ci_bound > 1

    result_index = [
        "oddsratio",
        "pvalue",
        "lower_ci_bound",
        "upper_ci_bound",
        "confidence_interval",
        "significant",
    ]
    result_vals = [
        oddsratio,
        pvalue,
        lower_ci_bound,
        upper_ci_bound,
        confidence_interval,
        significant,
    ]
    result_series = pd.Series(result_vals, index=result_index)
    # result_series["pvalue"] = f"{float(result_series['pvalue']):.5f}"
    return result_series


# example: vcf = background_vcf, background_vcf = background_vcf, enriched_for_vcf = ccre.bcf, tag_to_decile = "K562__skew"
def fishers_test_decile_overlap_enrichment(
    vcf, background_vcf, enriched_for_vcf, tag_to_decile, abs=True
):
    tag_values_cmd = f"bcftools query -f '%INFO/{tag_to_decile}\n' {vcf.vcf_file_path}"
    tag_values_df = dt.fread(cmd=tag_values_cmd).to_pandas()
    # ic(tag_values_df)
    if tag_values_df.empty:
        return pd.DataFrame()

    tag_values = tag_values_df.iloc[:, 0].astype(float)
    tag_values = tag_values.abs()
    tag_values.name = tag_to_decile

    tag_deciles, bin_edges = pd.qcut(tag_values, q=10, retbins=True)

    result_rows = []

    in_background_in_set_oi = get_variants_in_common(
        background_vcf, enriched_for_vcf, ret_num_variants=True
    )
    in_background_not_in_set_oi = BACKGROUND_VARS_TOTAL - in_background_in_set_oi

    for idx, val in enumerate(bin_edges):
        if idx == len(bin_edges) - 1:
            break

        lower_bound = val
        upper_bound = bin_edges[idx + 1]

        filter_string = (
            f"{tag_to_decile} > {lower_bound} and {tag_to_decile} <= {upper_bound}"
        )
        with tempfile.NamedTemporaryFile(
            mode="w+", dir="/tmp", delete=False, suffix=".bcf"
        ) as bcf_filtered:
            vcf.filter_by_info_tag_condition(
                filter_string=filter_string,
                output_format="b",
                output_file_path=bcf_filtered.name,
            )

        filtered_vcf = VCFObject(bcf_filtered.name)
        # ic(filtered_vcf)
        in_vcf_in_set_oi = get_variants_in_common(
            filtered_vcf, enriched_for_vcf, ret_num_variants=True
        )
        in_vcf_not_in_set_oi = filtered_vcf.get_num_variants() - in_vcf_in_set_oi

        if any(
            x == 0
            for x in [
                in_vcf_in_set_oi,
                in_vcf_not_in_set_oi,
                in_background_in_set_oi,
                in_background_not_in_set_oi,
            ]
        ):
            result = pd.Series(
                [np.nan, np.nan, np.nan, np.nan, np.nan, False],
                index=[
                    "oddsratio",
                    "pvalue",
                    "lower_ci_bound",
                    "upper_ci_bound",
                    "confidence_interval",
                    "significant",
                ],
            )
            result.name = f"{idx}:({lower_bound}, {upper_bound}]"
            result_rows.append(result)
            continue

        oddsratio, pvalue = fisher_exact(
            [
                [in_vcf_in_set_oi, in_vcf_not_in_set_oi],
                [in_background_in_set_oi, in_background_not_in_set_oi],
            ]
        )
        lower_ci_bound = np.exp(
            np.log(oddsratio)
            - 1.96
            * np.sqrt(
                (1 / in_vcf_in_set_oi)
                + (1 / in_background_in_set_oi)
                + (1 / in_vcf_not_in_set_oi)
                + (1 / in_background_not_in_set_oi)
            )
        )
        upper_ci_bound = np.exp(
            np.log(oddsratio)
            + 1.96
            * np.sqrt(
                (1 / in_vcf_in_set_oi)
                + (1 / in_background_in_set_oi)
                + (1 / in_vcf_not_in_set_oi)
                + (1 / in_background_not_in_set_oi)
            )
        )
        confidence_interval = (
            f"({'{:.3f}'.format(lower_ci_bound)}, {'{:.3f}'.format(upper_ci_bound)})"
        )
        pvalue = format(pvalue, ".2e")
        significant = lower_ci_bound > 1

        result_index = [
            "oddsratio",
            "pvalue",
            "lower_ci_bound",
            "upper_ci_bound",
            "confidence_interval",
            "significant",
        ]
        result_vals = [
            oddsratio,
            pvalue,
            lower_ci_bound,
            upper_ci_bound,
            confidence_interval,
            significant,
        ]
        result_series = pd.Series(result_vals, index=result_index)
        # result_series["pvalue"] = f"{float(result_series['pvalue']):.5f}"
        result_series.name = f"{idx}:({lower_bound}, {upper_bound}]"
        result_rows.append(result_series)
        ic(result_series)

    ret_df = pd.DataFrame(result_rows)
    ic(ret_df)
    return ret_df
