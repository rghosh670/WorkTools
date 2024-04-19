""" "
Class definition for a FASTA object.
"""

import logging
import os

import datatable as dt
import pandas as pd
import pyfaidx
from icecream import ic
from kipoiseq import Interval

import worktools.api as api
import worktools.api.utils as utils

logger = logging.getLogger("worktools")


class FastaStringExtractor:
    def __init__(self, fasta_file):
        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

    def extract(self, interval: Interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = Interval(
            interval.chrom,
            max(interval.start, 0),
            min(interval.end, chromosome_length),
        )
        # pyfaidx wants a 1-based interval
        sequence = str(
            self.fasta.get_seq(
                trimmed_interval.chrom,
                trimmed_interval.start + 1,
                trimmed_interval.stop,
            ).seq
        ).upper()
        # Fill truncated values with N's.
        pad_upstream = "N" * max(-interval.start, 0)
        pad_downstream = "N" * max(interval.end - chromosome_length, 0)
        return pad_upstream + sequence + pad_downstream

    def close(self):
        return self.fasta.close()


class FASTAObject:
    def __init__(self, fasta_file_path=None):
        self.fasta = None
        self.fasta_file_path = (
            os.path.abspath(fasta_file_path) if fasta_file_path is not None else None
        )
        self.fimo_result_df = None
        if fasta_file_path is not None:
            self.read_fasta(fasta_file_path)

    def __str__(self):
        if self.fasta is None:
            return "Empty FASTA object"

        seq_list = []
        for name in self.fasta.keys():
            seq = self.fasta[name].__str__()
            seq_list.extend([f">{name}", seq])

        max_lines = 200
        half_max_lines = max_lines // 2

        if len(seq_list) > max_lines:
            output = (
                "\n".join(seq_list[:half_max_lines])
                + "\n...\n"
                + "\n".join(seq_list[-half_max_lines:])
            )
            return output

        return "\n".join(seq_list)

    def __repr__(self):
        return self.__str__()

    def __del__(self):
        # if self.fasta_file_path.startswith('/tmp'):
        #     self.fasta_file_path.close()

        if self.fasta:
            self.fasta.close()

        return

    def read_fasta(self, fasta_file_path):
        """
        Reads a FASTA file into a python dictt.

        Parameters:
        - fasta_file_path: Path to the FASTA file.

        Returns:
        - A python dict containing the FASTA data.
        """
        self.fasta_file_path = fasta_file_path
        logger.info(f"Reading FASTA file {fasta_file_path}")
        self.fasta = pyfaidx.Fasta(fasta_file_path)
        logger.info(f"Read FASTA file {fasta_file_path}")
        return

    def write_to_file(self, file_path):
        """
        Writes the FASTA data to a file.

        Parameters:
        - file_path: Path to the output file.
        """
        utils.run_shell_command(f"> {file_path}")

        with open(file_path, "w") as output_file:
            output_file.write(self.__repr__())

        return

    def run_fimo_with_hocomoco(self, output_file_path):
        self.fimo_result_df = api.run_fimo_with_hocomoco(
            self, output_file_path=output_file_path
        )
        return

    def get_sequence_from_fasta_id(self, fasta_id):
        if self.fasta is None:
            raise ValueError("No FASTA data loaded.")

        if fasta_id not in self.fasta:
            raise ValueError(f"FASTA ID {fasta_id} not found in FASTA data.")

        return self.fasta[fasta_id].__str__()
