import logging
import os
import subprocess
import tempfile

import datatable as dt
import pandas as pd
import pysam
from icecream import ic

import worktools.api.utils as utils
from worktools.core.vcf_object import VCFObject

logger = logging.getLogger("worktools")


class LargeVCFObject(VCFObject):
    def __init__(self, vcf_file_path, vcf_df=None):
        self.vcf_file_path = vcf_file_path  # Save vcf_file_path as an instance variable

        # if not super().has_index() or super().is_index_outdated():
        #     super().index_variant_file()

        self.pysam_vcf = pysam.VariantFile(vcf_file_path, threads=os.cpu_count())

        if vcf_file_path:
            super().read_metadata()

        super().__init__(
            vcf_file_path, vcf_df=None, read_file=False
        )  # Corrected super syntax

    def __str__(self):
        tail_cmd = None
        if utils.is_bgzipped(self.vcf_file_path):
            head_cmd = f"zcat {self.vcf_file_path} | head -n 10"
            # tail_cmd = f"zcat {self.vcf_file_path} | tail -n 10"
        else:
            head_cmd = f"head -n 50 {self.vcf_file_path}"
            tail_cmd = f"tail -n 10 {self.vcf_file_path}"

        head_str = utils.run_shell_command(head_cmd)
        if tail_cmd:
            tail_str = utils.run_shell_command(tail_cmd)
        else:
            tail_str = "..."

        return_str = f"{head_str}\n...\n{tail_str}"
        return return_str

    def read_vcf(self, vcf_file_path, call_super=False):
        if call_super:
            super().read_vcf(vcf_file_path)
            return

        super().read_metadata()

        ic(self.metadata)
        return

    def remove_non_snps(self, output_file_path=None):
        remove_cmd = (
            f"bcftools view -V indels,mnps,other -Oz -o {output_file_path} {self.vcf_file_path}"
        )
        utils.run_shell_command(remove_cmd)

        tabix_cmd = f"tabix -p vcf {output_file_path}"
        utils.run_shell_command(tabix_cmd)
        return

    def get_rt_style_index(self):
        cmd = f"zgrep -v '^##' {self.vcf_file_path}"

        if self.vcf_file_path.endswith(".bcf"):
            cmd = f"bcftools view {self.vcf_file_path} | grep -v '##'"

        vcf_columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

        columns_dict = dict(
            zip(
                vcf_columns,
                (([dt.str32, dt.str32, None, dt.str32, dt.str32]) + ([None] * 3)),
            )
        )

        dt_df = dt.fread(cmd=cmd, columns=columns_dict, sep="\t", header=True)
        ic(dt_df)

        dt_df[
            :, dt.update(index=dt.f[0] + ";" + dt.f[1] + ";" + dt.f[2] + ";" + dt.f[3])
        ]
        # Convert the 'index' column to a pandas Series
        index_series = dt_df[:, "index"].to_pandas().iloc[:, 0]

        return pd.Index(index_series)

    def get_num_variants(self):
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
