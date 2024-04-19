import logging
import os
import subprocess
import tempfile
from pprint import pprint

import datatable as dt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from icecream import ic
from sklearn.metrics import (auc, precision_recall_curve, roc_auc_score,
                             roc_curve)

import worktools
import worktools.api as api
import worktools.api.utils as utils

logger = logging.getLogger("worktools")


def plot_prc_curves(
    vcf_obj,
    y_class_col,
    y_pred_cols,
    only_average=True,
    only_skew=True,
    output_file_path=None,
    output_file_dir=None,
    qc_col=None,
    abs=True,
    title=None,
    cell_type=None,
    show=False,
):
    if isinstance(y_pred_cols, str):
        y_pred_cols = [y_pred_cols]

    info_df = vcf_obj.get_info_col_df()

    if y_pred_cols == ["all"]:
        y_pred_cols = info_df.columns
        if qc_col is not None:
            y_pred_cols = [i for i in y_pred_cols if i != qc_col]

        if y_class_col in y_pred_cols:
            y_pred_cols = [i for i in y_pred_cols if i != y_class_col]

    if cell_type:
        cell_type = "K562" if "k562" == cell_type.lower() else cell_type
        cell_type = "HepG2" if "hepg2" == cell_type.lower() else cell_type
        cell_type = "SKNSH" if "sknsh" == cell_type.lower() else cell_type
        y_pred_cols = [
            i for i in y_pred_cols if cell_type in i and not i.startswith("CAGE")
        ]

    if only_skew:
        y_pred_cols = [i for i in y_pred_cols if not ("__" in i and "skew" not in i)]

    if qc_col is not None:
        info_df[qc_col] = [i == "True" for i in info_df[qc_col]]
        info_df = info_df.loc[info_df[qc_col]]

    info_df[y_class_col] = [i == "True" for i in info_df[y_class_col]]
    # y_class_col = np.array(info_df[y_class_col])
    y_class_col = np.array(info_df[y_class_col].astype(int))
    # print(y_class_col)

    for y_pred in y_pred_cols:
        if y_pred.startswith("_") and only_average:
            continue

        y_pred_col = info_df[y_pred].astype(float)
        y_pred_col = np.array(y_pred_col)
        y_pred_col = np.abs(y_pred_col)

        precision, recall, thresholds = precision_recall_curve(
            y_class_col, y_pred_col, pos_label=1
        )
        auc_score = auc(recall, precision)
        plt.plot(
            recall, precision, label=f"{y_pred.strip('_')} (AUC = {auc_score:.5f})"
        )
        plt.grid(False)
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        title = title if title is not None else f"Precision-Recall Curve"
        plt.title(title)
        plt.legend()
        # plt.grid()

    if title is None:
        output_file_name = ".".join(y_pred_cols) + "_prc.png"
    else:
        output_file_name = title.replace(" ", "_") + "_prc.png"

    if output_file_path is None:
        if not (output_file_dir is None or output_file_name is None):
            output_file_path = os.path.join(output_file_dir, output_file_name)

    if not show:
        plt.savefig(output_file_path, bbox_inches="tight")
    else:
        plt.show()
    plt.clf()
