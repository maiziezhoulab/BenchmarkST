"""
main function is currently ongoing debugging
"""
import warnings
import pandas as pd
import numpy as np
import scanpy as sc
from scipy import sparse
import os
from imputation.impute import impute_
import GAAE
import argparse
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import time
import seaborn as sns 
from GAAE.utils import impute, DE_num_calc, initialize, filter_num_calc, downstream_analyses
warnings.filterwarnings("ignore")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str, default="./", help="root dir for input data")
    parser.add_argument('--gt_dir', type=str, default="./", help="root dir for data ground truth")
    parser.add_argument('--input_data', type=str, default="151673", help="input data section id")
    parser.add_argument('--impute_cluster_num', type=str, default="7", help="diff cluster numbers for imputation")
    parser.add_argument('--cluster_num', type=int, default=7, help="input data cluster number")
    parser.add_argument('--radius', type=int, default=150, help="input data radius")
    parser.add_argument("--de_candidates", type=str, default="200", help="candidate de list for imputation, separated by comma")
    parser.add_argument('--no_de', type=int, default=0, help='switch on/off DEG selection module')
    parser.add_argument("--use_mean", type=int, default=0, help="use mean value in de list or not")
    parser.add_argument("--impute_runs", type=int, default=1, help="time of runs for imputation")
    parser.add_argument("--runs", type=int, default=2, help="total runs for the data")
    parser.add_argument('--gt', type=int, default=1, help="ground truth for the input data")
    parser.add_argument('--use_hvgs', type=int, default=3000, help="select highly variable genes before training")
    parser.add_argument('--use_preprocessing', type=int, default=1, help='use preprocessed input or raw input')
    parser.add_argument('--save_fig', type=int, default=1, help='saving output visualization')
    parser.add_argument('--filter_nzr', type=float, default=0.15, help='suggested nzr threshold for filtering')
    parser.add_argument('--filter_num', type=int, default=None, help='suggested gene threshold for filtering')
    parser.add_argument('--de_nzr_min', type=float, default=0.299, help='suggested min nzr threshold after de selection')
    parser.add_argument('--de_nzr_max', type=float, default=0.399, help='suggested max nzr threshold after de selection')
    parser.add_argument('--use_gpu_id', type=str, default='0', help='use which GPU, only applies when you have multiple gpu')
    args = parser.parse_args()
    args.impute_cluster_num = args.impute_cluster_num.split(",")  # ["5", "6", "7"]
    root_d = args.data_dir

    if args.input_data not in ['20180417_BZ5_control', '20180419_BZ9_control', '20180424_BZ14_control', 'STARmap_20180505_BY3_1k.h5ad']:
        filter_num = filter_num_calc(args, args.filter_num)
        print("optimized filter number = ", filter_num)
    else:
        filter_num = 0
    adata, adata_ori = initialize(args, filter_num)

    if args.de_candidates == "None":
        if os.path.exists('./cache/' + args.data_dir.split('/')[-2] + dataset + '.txt'):
            with open('./cache/' + args.data_dir.split('/')[-2] + '.txt', 'r') as fp:
                line = fp.readlines()[0]
                split_ = line.strip().split(",")
                de_top_k_list = []
                for e in split_:
                    de_top_k_list.append(int(e))
            print("previously cached de list = ", de_top_k_list)
        else:
            de_top_k_list = DE_num_calc(args, adata)
            print("optimized de list = ", de_top_k_list)
            with open('./cache/' + args.data_dir.split('/')[-2] + '.txt', 'a+') as fp:
                # fp.write('de list: ')
                fp.write(','.join([str(i) for i in de_top_k_list]))
                # fp.write('\n')
    else:
        split_ = args.de_candidates.strip().split(",")
        de_top_k_list = []
        for e in split_:
            de_top_k_list.append(int(e))
        print("manually defined de list = ", de_top_k_list)
    
    de_list_epoch = []
    if de_top_k_list != []:
        print("performing DEGs selection")
        adata_list = []
        for de_ in de_top_k_list:
            for cluster_n in args.impute_cluster_num:
                print("cluster_n = ", cluster_n)
                GAAE.get_kNN(adata, rad_cutoff=args.radius)

                ari_ini, ari_final, de_list, adata_out = GAAE.train_ADEPT_use_DE(adata, n_epochs=1000,
                                                                            num_cluster=int(cluster_n),
                                                                            dif_k=de_, device_id=args.use_gpu_id)
                de_list_epoch.append(de_list)
                adata_list.append(adata_out)
        g_union = set.union(*de_list_epoch)
        imputed_ad = impute(args, adata_list, g_union, de_top_k_list)
    else:
        print("skip performing DEGs selection")
        imputed_ad = adata

    """result of imputed data"""
    GAAE.get_kNN(imputed_ad, rad_cutoff=args.radius)
    ari_ini, ARI, de_list, adata_out = GAAE.train_ADEPT_use_DE(imputed_ad, n_epochs=1000, num_cluster=args.cluster_num, device_id=args.use_gpu_id)

    print('Dataset:', dataset)
    print('ARI:', ARI)
    aris.append(ARI)
    print(aris)

    if args.save_fig:
        if args.input_data in ['20180417_BZ5_control', '20180419_BZ9_control', '20180424_BZ14_control', 'STARmap_20180505_BY3_1k.h5ad']:
            timestr = time.strftime("%Y%m%d-%H%M%S")
            plt.rcParams["figure.figsize"] = (3, 3)
            sc.pl.spatial(adata_out, color=["mclust_impute", "Ground Truth"],
                            title=['ADEPT (ARI=%.2f)' % ari_ini, "Ground Truth"], spot_size=95)
            plt.savefig(os.path.join(root_d, args.input_data + '_viz', "_" + timestr + "_" + str(i) + ".pdf"))
            downstream_analyses(args.input_data, adata_out, ari_ini, root_d, args.input_data + "_" + timestr, imputed_=1)
        if args.input_data in ['151507', '151508', '151509', '151510', '151673', '151674', '151675', '151676']:
            timestr = time.strftime("%Y%m%d-%H%M%S")
            plt.rcParams["figure.figsize"] = (3, 3)
            sc.pl.spatial(adata_out, color=["mclust_impute", "Ground Truth"],
                            title=['ADEPT (ARI=%.2f)' % ari_ini, "Ground Truth"], spot_size=55)
            plt.savefig(os.path.join(root_d, args.input_data + '_viz', "_" + timestr + "_" + str(i) + ".pdf"))
            downstream_analyses(args.input_data, adata_out, ari_ini, root_d, args.input_data + "_" + timestr, imputed_=1)
        if args.input_data == 'section1':
            timestr = time.strftime("%Y%m%d-%H%M%S")
            plt.rcParams["figure.figsize"] = (3, 3)
            sc.pl.spatial(adata_out, color=["mclust_impute", "Ground Truth"],
                            title=['ADEPT (ARI=%.2f)' % ari_ini, "Ground Truth"], spot_size=150, color_map='viridis')
            print(adata_out.uns)
            print(adata_out.uns['mclust_impute_colors'])
            adata_out.uns['mclust_impute_colors'] = ['#440154', '#481467', '#482576', '#453781', '#404688',
                                                    '#39558c', '#33638d', '#2d718e', '#287d8e', '#238a8d',
                                                    '#1f968b', '#20a386', '#29af7f', '#3dbc74', '#56c667',
                                                    '#75d054', '#95d840', '#bade28', '#dde318', '#fde725']
            print(adata_out.uns['mclust_impute_colors'])
            sc.pl.spatial(adata_out, color=["mclust_impute", "Ground Truth"],
                            title=['ADEPT (ARI=%.2f)' % ari_ini, "Ground Truth"], spot_size=150, color_map='viridis')
            plt.savefig(os.path.join(root_d, args.input_data + '_viz', "_" + timestr + "_" + str(i) + ".pdf"))
            downstream_analyses(args.input_data, adata_out, ari_ini, root_d, args.input_data + "_" + timestr, imputed_=1)
