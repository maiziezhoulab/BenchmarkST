import argparse
import scanpy as sc
import numpy as np
import torch
import random
import copy
import os


def average_impute_(unimputed_m):

    row_num, col_num = unimputed_m.shape
    unimputed_m_transpose = unimputed_m.T
    imputed_m = copy.deepcopy(unimputed_m_transpose)
    for col_idx in range(col_num):
        sum_ = np.sum(unimputed_m_transpose[col_idx])
        denom_ = np.count_nonzero(unimputed_m_transpose[col_idx])
        if denom_ == 0:
            pass
        else:
            # print(sum_, denom_)
            avg_ = sum_ / denom_
            imputed_m[col_idx] = np.where(imputed_m[col_idx] == 0, avg_, imputed_m[col_idx])

    return imputed_m.T


def impute_(cluster_num, exp_m, pred_label_list, barcode_list):
    dict_ = {}

    """
    doc num of spots for each cluster
    """
    spots_num_sep = []
    for cluster_idx in range(cluster_num):
        spots_num_sep.append(0)

    imputed_exp_m = np.zeros_like(exp_m)
    for idx_ in range(len(barcode_list)):
        # e1: idx, e2: pred_label, e3: expression
        dict_[idx_] = [int(pred_label_list[idx_]), exp_m[idx_]]
        # print(int(pred_label_list[idx_]))
        spots_num_sep[int(pred_label_list[idx_])] += 1

    # print(spots_num_sep)
    """
    initialize empty matrix for each cluster
    """
    matrix_list = []
    for cluster_idx in range(cluster_num):
        # print((spots_num_sep[cluster_idx], exp_m.shape[1]))
        matrix_list.append(np.zeros((spots_num_sep[cluster_idx], exp_m.shape[1])))
        spots_num_sep[cluster_idx] = 0

    """
    assign spots within same cluster to each matrix
    """

    # for cluster_idx in range(cluster_num):
    for k in range(len(barcode_list)):
        # if int(dict_[k][0]) == cluster_idx:
        # print(k)
        # print(matrix_list)
        # print(k)
        # print(dict_[k][0])
        # print(dict_[k][0], spots_num_sep[dict_[k][0]])
        matrix_list[dict_[k][0]][spots_num_sep[dict_[k][0]]] = dict_[k][1]
        # print()
        spots_num_sep[dict_[k][0]] += 1
    print(spots_num_sep)
    """
    imputation within each smaller matrix
    """
    for cluster_idx in range(cluster_num):
        matrix_list[cluster_idx] = average_impute_(matrix_list[cluster_idx])
        spots_num_sep[cluster_idx] = 0

    """
    assign imputed values in each smaller matrix back to the large one
    """
    for k in range(len(barcode_list)):

        imputed_exp_m[k] = matrix_list[dict_[k][0]][spots_num_sep[dict_[k][0]]]
        spots_num_sep[dict_[k][0]] += 1
    # imputed_matrices.append(imputed_exp_m)
    #
    # total_m = 6
    # avg_m = np.zeros_like(imputed_matrices[0])
    # for m in imputed_matrices:
    #     avg_m += m
    # # final inputed matrix
    # avg_m /= total_m
    # # avg_m is the final output
    # print(avg_m)
    return imputed_exp_m


if __name__ == '__main__':

    section_id = '151673'
    """
    toy input
    """
    input_dir = os.path.join('/home/yunfei/Spatial/SphereLoss/withAnnotations', section_id)
    adata_ = sc.read_visium(path=input_dir, count_file=section_id + '_filtered_feature_bc_matrix.h5')

    barcode_list_ = adata_.obs.index.values.tolist()
    print(barcode_list_)
    # pred_label_list = random.sample(range(0, cluster_num), len(barcode_list))
    # values = [5, 6, 7, 8, 9, 10]
    values = [0, 1, 2, 3, 4, 5, 6]
    pred_label_list_ = random.choices(values, k=len(barcode_list_))
    exp_m_ = adata_.X.toarray()  # spots by genes
    m_ = impute_(7, exp_m_, pred_label_list_, barcode_list_)
