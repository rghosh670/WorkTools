import concurrent.futures
import logging
import os
import subprocess
import tempfile
from functools import partial
from pprint import pprint

import datatable as dt
import numpy as np
import pandas as pd
from icecream import ic
from pyliftover import LiftOver

import worktools.api as api
import worktools.api.utils as utils
from worktools.api.conversion import convert_bed_to_vcf, convert_vcf_to_bed
from worktools.core import BEDObject, LargeVCFObject, VCFObject

logger = logging.getLogger("worktools")

MALINOIS_BCF_DIR_PATH = (
    "/home/rg972/project/WorkData/gnomAD_variants/malinois_predictions/bcf_files"
)
MALINOIS_VCF_DIR_PATH = (
    "/home/rg972/project/WorkData/gnomAD_variants/malinois_predictions/vcf_files"
)
MALINOIS_BED_DIR_PATH = (
    "/home/rg972/project/WorkData/gnomAD_variants/malinois_predictions/bed_files"
)
MALINOIS_BCF_FILE_PATH = "/home/rg972/project/WorkData/gnomAD_variants/malinois_predictions/gnomad_variants_combined.bcf"

CCRE_BED_FILE_PATH = (
    "/home/rg972/project/WorkData/genome_annotations/GRCh38-cCREs.V4.bed.gz"
)
DELS_BED_FILE_PATH = (
    "/home/rg972/project/WorkData/genome_annotations/GRCh38-dELS.V4.bed.gz"
)
PELS_BED_FILE_PATH = (
    "/home/rg972/project/WorkData/genome_annotations/GRCh38-pELS.V4.bed.gz"
)
PLS_BED_FILE_PATH = (
    "/home/rg972/project/WorkData/genome_annotations/GRCh38-PLS.V4.bed.gz"
)

EXONS_BED_FILE_PATH = "/home/rg972/project/WorkData/genome_annotations/gencode.v44.protein.coding.exons.autosomes.v2.bed"
TWO_FIFTY_BP_BED_FILE_PATH = "/home/rg972/project/WorkData/genome_annotations/gencode.v44.protein.coding.250bp.promoters.autosomes.v2.bed"
FIVE_HUNDRED_BP_BED_FILE_PATH = "/home/rg972/project/WorkData/genome_annotations/gencode.v44.protein.coding.500bp.promoters.autosomes.v2.bed"
SEVEN_FIFTY_BP_BED_FILE_PATH = "/home/rg972/project/WorkData/genome_annotations/gencode.v44.protein.coding.750bp.promoters.autosomes.v2.bed"
THOUSAND_BP_BED_FILE_PATH = "/home/rg972/project/WorkData/genome_annotations/gencode.v44.protein.coding.1kb.promoters.autosomes.v2.bed"


def adjust_chr_prefix_in_vcf(vcf_obj, action):
    """
    Adjusts the 'chr' prefix in the '#CHROM' column of a VCF dataframe.

    Parameters:
    - vcf_df: A pandas DataFrame containing VCF data.
    - action: A string specifying the action to perform: 'add' to add 'chr' prefix, 'remove' to remove it.
    """

    match action:
        case "add":
            if not vcf_obj.has_chr_prefix():
                vcf_obj.vcf_df["#CHROM"] = "chr" + vcf_obj.vcf_df["#CHROM"]
        case "remove":
            if vcf_obj.has_chr_prefix():
                vcf_obj.vcf_df["#CHROM"] = vcf_obj.vcf_df["#CHROM"].str.lstrip("chr")
        case _:
            raise ValueError("Invalid action specified. Choose 'add' or 'remove'.")

    return


def lift_over_vcf(vcf_obj, genome_build, chunksize=10000):
    """
    Lifts over the VCF data to a new genome build using a chain file, processing in chunks.

    Parameters:
    - genome_build: The genome build to lift over to (GRCh38 or GRCh37).
    - chunksize: The number of rows to process in each chunk.
    """
    if genome_build not in ["hg19", "hg38", "GRCh38", "GRCh37"]:
        raise ValueError(
            "Invalid genome build specified. Choose from 'hg19', 'hg38', 'GRCh38', or 'GRCh37'"
        )

    build_aliases = {
        "GRCh37": "hg19",
        "GRCh38": "hg38",
    }

    genome_build = build_aliases.get(genome_build, genome_build)
    current_genome_build = vcf_obj.get_genome_build_vcf()

    try:
        if genome_build == current_genome_build:
            logger.info("VCF file is already in the specified genome build.")
            return
    except NotImplementedError:
        pass

    lo = vcf_obj.get_liftover_object(current_genome_build, genome_build)
    vcf_df = vcf_obj.vcf_df

    num_unmapped_variants = 0
    num_variants_with_multiple_mappings = 0

    def process_chunk(chunk):
        nonlocal num_unmapped_variants
        nonlocal num_variants_with_multiple_mappings
        lifted_chunk = []
        for _, row in chunk.iterrows():
            converted_coord = lo.convert_coordinate(row["#CHROM"], row["POS"])
            if not converted_coord:
                num_unmapped_variants += 1
                continue
            # Handle both single and multiple conversion results
            if len(converted_coord) > 1:
                num_variants_with_multiple_mappings += 1

            for coord in converted_coord:
                new_row = [coord[0], coord[1]] + row.tolist()[
                    2:
                ]  # Adjust based on the structure of your VCF
                lifted_chunk.append(new_row)
        return lifted_chunk

    # Process in chunks
    chunks = [
        vcf_df.iloc[i : i + chunksize] for i in range(0, vcf_df.shape[0], chunksize)
    ]
    lifted_df_rows = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_to_chunk = {
            executor.submit(process_chunk, chunk): chunk for chunk in chunks
        }
        for future in concurrent.futures.as_completed(future_to_chunk):
            lifted_df_rows.extend(future.result())

    lifted_df = pd.DataFrame(lifted_df_rows, columns=vcf_df.columns)
    vcf_obj.vcf_df = lifted_df
    vcf_obj.sort_vcf_df()
    logging.info(
        "Lifted over VCF data to %s. %d variants had no mappings. %d variants had multiple mappings",
        genome_build,
        num_unmapped_variants,
        num_variants_with_multiple_mappings,
    )
    return


# TODO: Accomodate largevcfobjects
def move_columns_from_one_vcf_to_another(
    source_vcf_obj,
    dest_vcf_obj,
    columns=["ID", "QUAL", "FILTER", "INFO"],
    transfer_metadata=True,
    update_file=True,
    update_vcf_df=True,
    output_file_format="z",
    output_file_path=None,
):
    # if isinstance(source_vcf_obj, LargeVCFObject):
    #     source_vcf_obj.read_vcf(source_vcf_obj.vcf_file_path, call_super = True)

    # if isinstance(dest_vcf_obj, LargeVCFObject):
    #     dest_vcf_obj.read_vcf(dest_vcf_obj.vcf_file_path, call_super = True)

    # if isinstance(source_vcf_obj, VCFObject):
    #     source_vcf_obj.add_rt_style_index(force = True)
    #     dest_vcf_obj.add_rt_style_index(force = True)

    #     source_df = source_vcf_obj.vcf_df
    #     dest_df = dest_vcf_obj.vcf_df

    #     temp_df = source_df[columns].loc[dest_df.index.intersection(source_df.index)]
    #     dest_vcf_obj.vcf_df.update(temp_df)

    # elif isinstance(source_vcf, str):
    #     pass
    cols_string = ",".join(columns)
    source_vcf_obj_fp = source_vcf_obj.vcf_file_path
    dest_vcf_obj_fp = dest_vcf_obj.vcf_file_path

    if dest_vcf_obj_fp is None:
        dest_vcf_obj.flush_to_temp_file()
        dest_vcf_obj_fp = dest_vcf_obj.vcf_file_path
        dest_vcf_obj.read_vcf(dest_vcf_obj_fp)
        ic(dest_vcf_obj)

    annotate_cmd = f"bcftools annotate {dest_vcf_obj_fp} -a {source_vcf_obj_fp} -c {cols_string} -O{output_file_format} -o {output_file_path} --threads {os.cpu_count()}"
    utils.run_shell_command(annotate_cmd)

    dest_vcf_obj.vcf_file_path = output_file_path
    dest_vcf_obj.index_variant_file()
    dest_vcf_obj.read_vcf(output_file_path)

    # os.system(f"bcftools view {output_file_path}")
    if transfer_metadata:
        if source_vcf_obj.metadata is None:
            source_vcf_obj.read_metadata()

        dest_vcf_obj.metadata = source_vcf_obj.metadata.strip()

    return


def modify_alleles_to_n(vcf_path, output_path):
    cmd = f'bcftools view {vcf_path} | awk \'{{if($0 ~ /^#/) print; else {{$4="N"; $5="N"; print}}}}\' > {output_path}'
    subprocess.run(cmd, shell=True, check=True)
    return


def get_variants_in_common(
    vcf_obj1,
    vcf_obj2,
    output_file_path=None,
    ret_num_variants=False,
    output_format="z",
    ignore_alleles=False,
):
    vcf1_file_path = vcf_obj1.vcf_file_path
    vcf2_file_path = vcf_obj2.vcf_file_path

    # if ignore_alleles:
    # # Create temporary files for modified VCFs
    #     _, temp_vcf1_path = tempfile.mkstemp(suffix='.vcf')
    #     _, temp_vcf2_path = tempfile.mkstemp(suffix='.vcf')

    #     # Modify REF and ALT alleles to 'N' for both VCFs
    #     modify_alleles_to_n(vcf1_file_path, temp_vcf1_path)
    #     modify_alleles_to_n(vcf2_file_path, temp_vcf2_path)

    #     # Update paths to use modified VCFs
    #     vcf1_file_path = temp_vcf1_path
    #     vcf2_file_path = temp_vcf2_path

    if ret_num_variants:
        ignore_alleles_string = "-c both" if ignore_alleles else ""
        isec_cmd = f'bcftools isec -n=2 -w1 {vcf1_file_path} {vcf2_file_path} {ignore_alleles_string} | grep -v "^#" | wc -l'
        num_variants = utils.run_shell_command(isec_cmd)
        return int(num_variants)
    else:
        isec_cmd = f"bcftools isec -n=2 -w1 {vcf1_file_path} {vcf2_file_path} -O{output_format} -o {output_file_path}"
        utils.run_shell_command(isec_cmd)
        temp_vcf = LargeVCFObject(output_file_path)
        temp_vcf.index_variant_file()
    return


def get_info_tags(vcf_obj):
    # Command to get the header of the VCF/BCF file
    command = ["bcftools", "view", "-h", vcf_obj.vcf_file_path]

    # Execute the command
    result = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    # Check for errors
    if result.returncode != 0:
        raise Exception(f"Error in bcftools command: {result.stderr}")

    # Extract INFO tags from the header
    info_tags = []
    for line in result.stdout.splitlines():
        if line.startswith("##INFO="):
            # Extract the ID attribute from the INFO definition
            info_id = line.split("ID=")[1].split(",")[0]
            info_tags.append(info_id)

    return info_tags


def annotate_vcf_with_info_from_other_any_allele(
    source_vcf,
    dest_vcf,
    columns=None,
    output_file_path=None,
    output_format="z",
    update_file=False,
    update_vcf_df=False,
):
    df = dest_vcf.vcf_df

    if df.empty:
        logger.warn("VCF file is empty. No variants to annotate.")
        return

    df["GROUP"] = df.groupby(["#CHROM", "POS"]).cumcount()

    max_duplicates = df["GROUP"].max() + 1

    concat_file_paths = []

    for i in range(max_duplicates):
        # with tempfile.NamedTemporaryFile(mode='w+', dir='/tmp', delete=False, suffix='.vcf') as vcf_temp_file:
        df_filtered = df[df["GROUP"] == i]
        df_filtered = df_filtered.drop(columns=["GROUP"])
        current_vcf = VCFObject(vcf_df=df_filtered)
        current_vcf.flush_to_temp_file(output_file_format="v")
        with tempfile.NamedTemporaryFile(
            mode="w+", dir="/tmp", delete=False, suffix=".bcf"
        ) as anno_vcf:
            errno = annotate_vcf_with_info_from_other(
                source_vcf=source_vcf,
                dest_vcf=current_vcf,
                columns=columns,
                output_file_path=anno_vcf.name,
                output_format="b",
                update_file=update_file,
                update_vcf_df=update_vcf_df,
                ignore_alleles=True,
            )
            if errno == -1:
                continue
            annotated_vcf = VCFObject(anno_vcf.name)
            concat_file_paths.append(
                annotated_vcf.flush_to_temp_file(output_file_format="b")
            )

    if len(concat_file_paths) == 0:
        logger.warn("No variants were annotated.")
        return

    if len(concat_file_paths) == 1:
        mv_cmd = f"mv {concat_file_paths[0]} {output_file_path}"
        utils.run_shell_command(mv_cmd)
        return

    concat_cmd = f"bcftools concat -a {' '.join(concat_file_paths)} -O{output_format} -o {output_file_path}"
    return utils.run_shell_command(concat_cmd)


def index_and_bgzip_file(file_path):
    os.system(f"bgzip -f {file_path}")
    os.system(f"tabix -p vcf {file_path}.gz")
    return f"{file_path}.gz"


def annotate_vcf_with_info_from_other(
    source_vcf,
    dest_vcf,
    columns=None,
    output_file_path=None,
    output_format="z",
    update_file=False,
    update_vcf_df=False,
    ignore_alleles=False,
):
    source_vcf_file_path = (
        source_vcf.vcf_file_path
    )  # if isinstance(source_vcf, VCFObject) else source_vcf
    dest_vcf_file_path = (
        dest_vcf.vcf_file_path
    )  # if isinstance(dest_vcf, VCFObject) else source_vcf

    if source_vcf_file_path.endswith(".vcf"):
        source_vcf_file_path = index_and_bgzip_file(source_vcf_file_path)

    if dest_vcf_file_path.endswith(".vcf"):
        dest_vcf_file_path = index_and_bgzip_file(dest_vcf_file_path)

    ic(source_vcf)
    ic(dest_vcf)
    # if isinstance(source_vcf, VCFObject):
    #     num_variants_in_common = get_variants_in_common(source_vcf, dest_vcf, ret_num_variants = True, ignore_alleles = True)
    #     if num_variants_in_common == 0:
    #         logger.warn("The two VCF files have no variants in common.")
    #         return -1

    if columns is None:
        columns = get_info_tags(source_vcf)

    info_cols = [f"INFO/{col}" for col in columns]
    info_cols = ",".join(info_cols)

    if output_file_path is None:
        with tempfile.NamedTemporaryFile(
            mode="w+", dir="/tmp", delete=False, suffix=".vcf.gz"
        ) as vcf_temp_file:
            annotate_cmd = f"bcftools annotate {dest_vcf_file_path} -a {source_vcf_file_path} -c {info_cols} -O{output_format} -o {vcf_temp_file.name} --threads {os.cpu_count()} --write-index"
            if ignore_alleles:
                annotate_cmd += " --pair-logic both"

            vcf_temp_file.flush()
            annotate_cmd_output = utils.run_shell_command(annotate_cmd)

            if "Error" in annotate_cmd_output:
                raise RuntimeError(f"Error annotating VCF file: {annotate_cmd_output}")

            mv_cmd = f"mv {vcf_temp_file.name} {dest_vcf_file_path}"

            if update_file:
                utils.run_shell_command(mv_cmd)

                if isinstance(source_vcf, VCFObject):
                    logger.info(
                        "%d variants were annotated with INFO from the source VCF file.",
                        num_variants_in_common,
                    )
                else:
                    logger.info(
                        "VCF file was annotated with INFO from the source VCF file."
                    )
    else:
        annotate_cmd = f"bcftools annotate {dest_vcf_file_path} -a {source_vcf_file_path} -c {info_cols} -O{output_format} -o {output_file_path} --threads {os.cpu_count()} --write-index"
        if ignore_alleles:
            annotate_cmd += " --pair-logic both"
        annotate_cmd_output = utils.run_shell_command(annotate_cmd)
        return

    if update_vcf_df and isinstance(dest_vcf, VCFObject) and not update_file:
        dest_vcf.read_vcf(vcf_temp_file.name)

    return vcf_temp_file.name


def get_malinois_scores(vcf_obj, scores="all"):
    cell_type = ["K562", "HepG2", "SKNSH"]

    if scores == "all":
        scores = ["ref", "alt", "skew"]

    info_tags = [f"{a}__{b}" for a in cell_type for b in scores]

    # Get a list of all .bcf files in the directory
    bcf_files = [
        os.path.join(MALINOIS_BCF_DIR_PATH, f)
        for f in os.listdir(MALINOIS_BCF_DIR_PATH)
        if f.endswith(".bcf")
    ]

    for bcf_file in bcf_files:
        annotate_vcf_with_info_from_other(bcf_file, vcf_obj, info_tags)

    return


def get_promoter_variants(
    vcf_obj,
    output_file_path=None,
    output_file_format="z",
    update_file=True,
    update_vcf_df=True,
):
    promoter_bed = BEDObject(bed_file_path=PLS_BED_FILE_PATH)
    filtered_vcf_output_file_path = api.keep_variants_in_bed(
        vcf_obj,
        promoter_bed,
        output_file_path=output_file_path,
        output_file_format=output_file_format,
    )

    if update_file:
        mv_cmd = f"mv {filtered_vcf_output_file_path} {vcf_obj.vcf_file_path}"
        utils.run_shell_command(mv_cmd)

    if update_vcf_df:
        vcf_obj.read_vcf(vcf_obj.vcf_file_path)
    return


def get_dels_variants(
    vcf_obj,
    output_file_path=None,
    output_file_format="z",
    update_file=True,
    update_vcf_df=True,
):
    dels_bed = BEDObject(bed_file_path=DELS_BED_FILE_PATH)
    filtered_vcf_output_file_path = api.keep_variants_in_bed(
        vcf_obj,
        dels_bed,
        output_file_path=output_file_path,
        output_file_format=output_file_format,
    )

    if update_file:
        mv_cmd = f"mv {filtered_vcf_output_file_path} {vcf_obj.vcf_file_path}"
        utils.run_shell_command(mv_cmd)

    if update_vcf_df:
        vcf_obj.read_vcf(vcf_obj.vcf_file_path)
    return


def get_pels_variants(
    vcf_obj,
    output_file_path=None,
    output_file_format="z",
    update_file=True,
    update_vcf_df=True,
):
    pels_bed = BEDObject(bed_file_path=PELS_BED_FILE_PATH)
    filtered_vcf_output_file_path = api.keep_variants_in_bed(
        vcf_obj,
        pels_bed,
        output_file_path=output_file_path,
        output_file_format=output_file_format,
    )

    if update_file:
        mv_cmd = f"mv {filtered_vcf_output_file_path} {vcf_obj.vcf_file_path}"
        utils.run_shell_command(mv_cmd)

    if update_vcf_df:
        vcf_obj.read_vcf(vcf_obj.vcf_file_path)
    return


def get_ccre_variants(
    vcf_obj,
    output_file_path=None,
    output_file_format="z",
    update_file=True,
    update_vcf_df=True,
):
    ccre_bed = BEDObject(bed_file_path=CCRE_BED_FILE_PATH)
    filtered_vcf_output_file_path = api.keep_variants_in_bed(
        vcf_obj,
        ccre_bed,
        output_file_path=output_file_path,
        output_file_format=output_file_format,
    )

    if update_file:
        mv_cmd = f"mv {filtered_vcf_output_file_path} {vcf_obj.vcf_file_path}"
        utils.run_shell_command(mv_cmd)

    if update_vcf_df:
        vcf_obj.read_vcf(vcf_obj.vcf_file_path)
    return


def get_variants_within_250_bp_of_tss(
    vcf_obj,
    output_file_path=None,
    output_file_format="z",
    update_file=True,
    update_vcf_df=True,
):
    two_fifty_bed = BEDObject(bed_file_path=TWO_FIFTY_BP_BED_FILE_PATH)
    filtered_vcf_output_file_path = api.keep_variants_in_bed(
        vcf_obj,
        two_fifty_bed,
        output_file_path=output_file_path,
        output_file_format=output_file_format,
    )

    if update_file:
        mv_cmd = f"mv {filtered_vcf_output_file_path} {vcf_obj.vcf_file_path}"
        utils.run_shell_command(mv_cmd)

    if update_vcf_df:
        vcf_obj.read_vcf(vcf_obj.vcf_file_path)
    return


def get_variants_within_500_bp_of_tss(
    vcf_obj,
    output_file_path=None,
    output_file_format="z",
    update_file=True,
    update_vcf_df=True,
):
    five_fifty_bed = BEDObject(bed_file_path=FIVE_HUNDRED_BP_BED_FILE_PATH)
    filtered_vcf_output_file_path = api.keep_variants_in_bed(
        vcf_obj,
        five_fifty_bed,
        output_file_path=output_file_path,
        output_file_format=output_file_format,
    )

    if update_file:
        mv_cmd = f"mv {filtered_vcf_output_file_path} {vcf_obj.vcf_file_path}"
        utils.run_shell_command(mv_cmd)

    if update_vcf_df:
        vcf_obj.read_vcf(vcf_obj.vcf_file_path)
    return


def get_variants_within_750_bp_of_tss(
    vcf_obj,
    output_file_path=None,
    output_file_format="z",
    update_file=True,
    update_vcf_df=True,
):
    seven_fifty_bed = BEDObject(bed_file_path=SEVEN_FIFTY_BP_BED_FILE_PATH)
    filtered_vcf_output_file_path = api.keep_variants_in_bed(
        vcf_obj,
        seven_fifty_bed,
        output_file_path=output_file_path,
        output_file_format=output_file_format,
    )

    if update_file:
        mv_cmd = f"mv {filtered_vcf_output_file_path} {vcf_obj.vcf_file_path}"
        utils.run_shell_command(mv_cmd)

    if update_vcf_df:
        vcf_obj.read_vcf(vcf_obj.vcf_file_path)
    return


def get_variants_within_1kb_of_tss(
    vcf_obj,
    output_file_path=None,
    output_file_format="z",
    update_file=True,
    update_vcf_df=True,
):
    thousand_bed = BEDObject(bed_file_path=THOUSAND_BP_BED_FILE_PATH)
    filtered_vcf_output_file_path = api.keep_variants_in_bed(
        vcf_obj,
        thousand_bed,
        output_file_path=output_file_path,
        output_file_format=output_file_format,
    )

    if update_file:
        mv_cmd = f"mv {filtered_vcf_output_file_path} {vcf_obj.vcf_file_path}"
        utils.run_shell_command(mv_cmd)

    if update_vcf_df:
        vcf_obj.read_vcf(vcf_obj.vcf_file_path)
    return


def get_exonic_variants(
    vcf_obj,
    output_file_path=None,
    output_file_format="z",
    update_file=True,
    update_vcf_df=True,
):
    exon_bed = BEDObject(bed_file_path=EXONS_BED_FILE_PATH)
    filtered_vcf_output_file_path = api.keep_variants_in_bed(
        vcf_obj,
        exon_bed,
        output_file_path=output_file_path,
        output_file_format=output_file_format,
    )

    if update_file:
        mv_cmd = f"mv {filtered_vcf_output_file_path} {vcf_obj.vcf_file_path}"
        utils.run_shell_command(mv_cmd)

    if update_vcf_df:
        vcf_obj.read_vcf(vcf_obj.vcf_file_path)
    return


def get_non_exonic_variants(
    vcf_obj,
    output_file_path=None,
    output_file_format="z",
    update_file=True,
    update_vcf_df=True,
):
    exon_bed = BEDObject(bed_file_path=EXONS_BED_FILE_PATH)
    filtered_vcf_output_file_path = api.remove_variants_in_bed(
        vcf_obj,
        exon_bed,
        output_file_path=output_file_path,
        output_file_format=output_file_format,
    )
    ic(filtered_vcf_output_file_path)

    if update_file:
        mv_cmd = f"mv {filtered_vcf_output_file_path} {vcf_obj.vcf_file_path}"
        utils.run_shell_command(mv_cmd)

    if update_vcf_df:
        vcf_obj.read_vcf(vcf_obj.vcf_file_path)
    return


def filter_by_ccre(
    vcf_obj,
    ccre,
    output_file_format="z",
    output_file_path=None,
    update_file=True,
    update_vcf_df=True,
):
    match ccre.lower():
        case "dels":
            get_dels_variants(
                vcf_obj,
                output_file_path,
                output_file_format,
                update_file,
                update_vcf_df,
            )
        case "pels":
            get_pels_variants(
                vcf_obj,
                output_file_path,
                output_file_format,
                update_file,
                update_vcf_df,
            )
        case "pls" | "promoter":
            get_promoter_variants(
                vcf_obj,
                output_file_path,
                output_file_format,
                update_file,
                update_vcf_df,
            )
        case "ccre":
            get_ccre_variants(
                vcf_obj,
                output_file_path,
                output_file_format,
                update_file,
                update_vcf_df,
            )

        case "250_bp_of_tss" | "250_bp":
            get_variants_within_250_bp_of_tss(
                vcf_obj,
                output_file_path,
                output_file_format,
                update_file,
                update_vcf_df,
            )
        case "500_bp_of_tss" | "500_bp":
            get_variants_within_500_bp_of_tss(
                vcf_obj,
                output_file_path,
                output_file_format,
                update_file,
                update_vcf_df,
            )
        case "750_bp_of_tss" | "750_bp":
            get_variants_within_750_bp_of_tss(
                vcf_obj,
                output_file_path,
                output_file_format,
                update_file,
                update_vcf_df,
            )
        case "1_kb_of_tss" | "1000_bp":
            get_variants_within_1kb_of_tss(
                vcf_obj,
                output_file_path,
                output_file_format,
                update_file,
                update_vcf_df,
            )
        case _:
            raise ValueError(
                f"{ccre} - Invalid CCRE specified. Choose from 'dels', 'pels', 'pls', 'promoter', 'ccre', '250_bp_of_tss', '500_bp_of_tss', '750_bp_of_tss', or '1_kb_of_tss'"
            )
