import scanpy as sc
import os
import pandas as pd
import numpy as np
import anndata
import matplotlib.pyplot as plt


# for loading DLPFC12 data
def load_DLPFC(root_dir='../benchmarking_data/DLPFC12', section_id='151507'):
    # 151507, ..., 151676 12 in total
    ad = sc.read_visium(path=os.path.join(root_dir, section_id), count_file=section_id+'_filtered_feature_bc_matrix.h5', load_images=True)
    ad.var_names_make_unique()

    gt_dir = os.path.join(root_dir, section_id, 'gt')
    gt_df = pd.read_csv(os.path.join(gt_dir, 'tissue_positions_list_GTs.txt'), sep=',', header=None, index_col=0)
    ad.obs['original_clusters'] = gt_df.loc[:, 6]
    keep_bcs = ad.obs.dropna().index
    ad = ad[keep_bcs].copy()
    ad.obs['original_clusters'] = ad.obs['original_clusters'].astype(int).astype(str)
    # print(ad.obs)
    return ad


# for loading BC data
# cluster = 20
def load_BC(root_dir='../benchmarking_data/BC', section_id='section1'):
    # section1
    ad = sc.read_visium(path=os.path.join(root_dir, section_id), count_file=section_id+'_filtered_feature_bc_matrix.h5', load_images=True)
    ad.var_names_make_unique()
    
    gt_dir = os.path.join(root_dir, section_id, 'gt')
    gt_df = pd.read_csv(os.path.join(gt_dir, 'tissue_positions_list_GTs.txt'), sep=',', header=None, index_col=0)
    ad.obs['original_clusters'] = gt_df.loc[:, 6].astype(int)
    ad.obs['original_clusters'] += 1
    keep_bcs = ad.obs.dropna().index
    ad = ad[keep_bcs].copy()
    ad.obs['original_clusters'] = ad.obs['original_clusters'].astype(int).astype(str)
    return ad


# for loading mouse_PFC data
def load_mPFC(root_dir = '/home/yunfei/spatial_benchmarking/benchmarking_data/STARmap_mouse_PFC', section_id='20180417_BZ5_control'):  
    # section id = '20180417_BZ5_control', '20180419_BZ9_control', '20180424_BZ14_control' 3 in total
    # cluster       4                       4                       4
    info_file = os.path.join(root_dir, 'starmap_mpfc_starmap_info.xlsx')
    cnts_file = os.path.join(root_dir, 'starmap_mpfc_starmap_cnts.xlsx')
    
    xls_cnts = pd.ExcelFile(cnts_file)
    df_cnts = pd.read_excel(xls_cnts, section_id)
    
    xls_info = pd.ExcelFile(info_file)
    df_info = pd.read_excel(xls_info, section_id)

    spatial_X = df_info.to_numpy()
    obs_ = df_info
    obs_.columns = ['psuedo_barcodes', 'x', 'y', 'gt', 'original_clusters']
    obs_.index = obs_['psuedo_barcodes'].tolist()
    
    var_ = df_cnts.iloc[:, 0]
    var_ = pd.DataFrame(var_)
    
    ad = anndata.AnnData(X=df_cnts.iloc[:,1:].T, obs=obs_, var=var_)
    ad.obs['original_clusters'] = ad.obs['original_clusters'].astype(int).astype(str)
    spatial = np.vstack((ad.obs['x'].to_numpy(), ad.obs['y'].to_numpy()))
    ad.obsm['spatial'] = spatial.T
    return ad


# for loading mouse_visual_cortex data
# cluster = 7
def load_mVC(root_dir='../benchmarking_data/STARmap_mouse_visual_cortex', section_id='STARmap_20180505_BY3_1k.h5ad'):
    ad = sc.read(os.path.join(root_dir, section_id))
    ad.var_names_make_unique()
    ad.obs.columns = ['Total_counts', 'imagerow', 'imagecol', 'original_clusters']
    return ad


# for loading mHypothalamus data
# already preprocessed? Xs are floats
def load_mHypothalamus(root_dir='/home/yunfei/spatial_benchmarking/benchmarking_data/mHypothalamus', section_id='0.26'):
    # section id = '0.26', '0.21', '0.16', '0.11', '0.06', '0.01', '-0.04', '-0.09', '-0.14', '-0.19', '-0.24', '-0.29' 12 in total
    # cluster =     15      15      14      15      15      15      14       15       15       15      16        15
    info_file = os.path.join(root_dir, 'MERFISH_Animal1_info.xlsx')
    cnts_file = os.path.join(root_dir, 'MERFISH_Animal1_cnts.xlsx')
    xls_cnts = pd.ExcelFile(cnts_file)
    # print(xls_cnts.sheet_names)
    df_cnts = pd.read_excel(xls_cnts, section_id)
    
    xls_info = pd.ExcelFile(info_file)
    df_info = pd.read_excel(xls_info, section_id)
    # print(df_cnts, df_info)
    spatial_X = df_info.to_numpy()
    obs_ = df_info
    if len(df_info.columns) == 5:
        obs_.columns = ['psuedo_barcodes', 'x', 'y', 'original_clusters', 'Neuron_cluster_ID']
    elif len(df_info.columns) == 6:
        obs_.columns = ['psuedo_barcodes', 'x', 'y', 'cell_types', 'Neuron_cluster_ID', 'original_clusters']
        # print(section_id)
        # print(obs_['z'].nunique())
    obs_.index = obs_['psuedo_barcodes'].tolist()
    # print(obs_)

    var_ = df_cnts.iloc[:, 0]
    var_ = pd.DataFrame(var_)
    # print(var_)
    
    ad = anndata.AnnData(X=df_cnts.iloc[:,1:].T, obs=obs_, var=var_)
    spatial = np.vstack((ad.obs['x'].to_numpy(), ad.obs['y'].to_numpy()))
    ad.obsm['spatial'] = spatial.T
    return ad

# for loading her2_tumor data
def load_her2_tumor(root_dir='../benchmarking_data/Her2_tumor', section_id='A1'):
    # section id = A1(348) B1(295) C1(177) D1(309) E1(587) F1(692) G2(475) H1(613) ~J1(254), 8 in total
    # clusters =   6       5       4       4       4       4       7       7
    cnts_dir = os.path.join(root_dir, 'ST-cnts')
    gt_dir = os.path.join(root_dir, 'ST-pat/lbl')
    gt_file_name = section_id+'_labeled_coordinates.tsv'
    cnt_file_name = section_id+'.tsv'
    cnt_df = pd.read_csv(os.path.join(cnts_dir, cnt_file_name), sep='\t', header=0)
    # print(cnt_file_name)
    # print(cnt_df)
    gt_df = pd.read_csv(os.path.join(gt_dir, gt_file_name), sep='\t', header=0)
    # print(gt_file_name)
    # print(gt_df)
    keep_bcs = gt_df.dropna().index
    gt_df = gt_df.iloc[keep_bcs]
    xs = gt_df['x'].tolist()
    ys = gt_df['y'].tolist()
    # print(xs)
    # print(ys)
    rounded_xs = [round(elem) for elem in xs]
    # print(rounded_xs)
    rounded_ys = [round(elem) for elem in ys]
    # print(rounded_ys)

    res = [str(i) + 'x' + str(j) for i, j in zip(rounded_xs, rounded_ys)]
    # print(len(set(res)))
    gt_df['Row.names'] = res
    # print(gt_df)

    spatial_X = cnt_df.to_numpy()
    obs_ = gt_df
    obs_ = obs_.sort_values(by=['Row.names'])
    obs_ = obs_.loc[obs_['Row.names'].isin(cnt_df['Unnamed: 0'])]
    obs_ = obs_.reset_index(drop=True)
    # print(obs_)

    var_ = cnt_df.iloc[0, 1:]
    var_ = pd.DataFrame(var_)

    ad = anndata.AnnData(X=cnt_df.iloc[:,1:], obs=obs_, var=var_, dtype=np.int64)
    ad.obs['original_clusters'] = ad.obs['label']
    spatial = np.vstack((ad.obs['pixel_x'].to_numpy(), ad.obs['pixel_y'].to_numpy()))
    ad.obsm['spatial'] = spatial.T
    return ad


# 
# cluster = 52
def load_mMAMP(root_dir='/home/yunfei/spatial_benchmarking/benchmarking_data/mMAMP', section_id='MA'):
    ad = sc.read_visium(path=os.path.join(root_dir, section_id), count_file=section_id+'_filtered_feature_bc_matrix.h5', load_images=True)
    ad.var_names_make_unique()

    gt_dir = os.path.join(root_dir, section_id, 'gt')
    gt_df = pd.read_csv(os.path.join(gt_dir, 'tissue_positions_list_GTs.txt'), sep='\t', header=0, index_col=0)
    ad.obs = gt_df
    ad.obs['original_clusters'] = ad.obs['ground_truth']
    keep_bcs = ad.obs.dropna().index
    ad = ad[keep_bcs].copy()
    return ad


# visualize everything for sanity check
def anndata_visualization(ad, fname, save_folder='/home/yunfei/spatial_benchmarking/benchmarking_data/gt_visualization', col_name='Ground Truth', spot_size=150):
    sc.pl.spatial(ad,
                  color=[col_name],
                  title=[col_name],
                  show=True, spot_size=spot_size)
    plt.savefig(os.path.join(save_folder, fname + ".pdf"))


if __name__ == '__main__':
    # for sec in ['151507', '151508', '151509', '151510', '151669', '151670', '151671', '151672', '151673', '151674', '151675', '151676']:
    #     ad = load_DLPFC(root_dir='/home/yunfei/spatial_benchmarking/benchmarking_data/DLPFC12', section_id=sec)
    #     # print(ad.obs)
    #     # exit(-1)
    #     # spatial = np.vstack((ad.obs['pixel_x'].to_numpy(), ad.obs['pixel_y'].to_numpy()))
    #     # ad.obsm['spatial'] = spatial.T
    #     anndata_visualization(ad, 'dlpfc_' + sec, col_name='Ground Truth', spot_size=75)
    #     print(sec)
    #     print(ad.obs['Ground Truth'].nunique())
    # print("dlpfc test passed")

    # ad = load_BC(root_dir='/home/yunfei/spatial_benchmarking/benchmarking_data/BC', section_id='section1')
    # anndata_visualization(ad, 'bc_section1', col_name='Ground Truth', spot_size=155)
    # # print(sec)
    # print(ad.obs['Ground Truth'].nunique())
    # print("bc test passed")

    # for sec in ['20180417_BZ5_control', '20180419_BZ9_control', '20180424_BZ14_control']:
    #     ad = load_mPFC(root_dir = '/home/yunfei/spatial_benchmarking/benchmarking_data/STARmap_mouse_PFC', section_id=sec)
    #     spatial = np.vstack((ad.obs['x'].to_numpy(), ad.obs['y'].to_numpy()))
    #     ad.obsm['spatial'] = spatial.T
    #     # print(ad.obs)
    #     anndata_visualization(ad, 'mpfc_' + sec, col_name='z', spot_size=175)
    #     print(sec)
    #     print(ad.obs['GT'].nunique())
    #     print(ad.obs['z'].nunique())
    # print("mPFC test passed")

    for sec in ['0.26', '0.21', '0.16', '0.11', '0.06', '0.01', '-0.04', '-0.09', '-0.14', '-0.19', '-0.24', '-0.29']:
        ad = load_mHypothalamus(root_dir='/home/yunfei/spatial_benchmarking/benchmarking_data/mHypothalamus', section_id=sec)
        print(ad.obs)
        spatial = np.vstack((ad.obs['x'].to_numpy(), ad.obs['y'].to_numpy()))
        ad.obsm['spatial'] = spatial.T
        # exit(-1)
        # anndata_visualization(ad, 'mHypo_' + sec, col_name='GT', spot_size=30)
        # print(sec)
        # print(ad.obs['GT'].nunique())
    print("mH test passed")

    # ad = load_mVC(root_dir='/home/yunfei/spatial_benchmarking/benchmarking_data/STARmap_mouse_visual_cortex', section_id='STARmap_20180505_BY3_1k.h5ad')
    # print(ad)
    # anndata_visualization(ad, 'mvc_' + 'STARmap_20180505_BY3_1k.h5ad', col_name='Ground Truth', spot_size=175)
    # print(sec)
    # print(ad.obs['Ground Truth'].nunique())
    # print("mVC test passed")

    # for sec in ['A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G2', 'H1']:
    #     ad = load_her2_tumor(root_dir='/home/yunfei/spatial_benchmarking/benchmarking_data/Her2_tumor', section_id=sec)
    #     spatial = np.vstack((ad.obs['pixel_x'].to_numpy(), ad.obs['pixel_y'].to_numpy()))
    #     ad.obsm['spatial'] = spatial.T
    #     anndata_visualization(ad, 'her2tumor_' + sec, col_name='label', spot_size=75)
    #     print(sec)
    #     print(ad.obs['label'].nunique())
    # print("her2tumor test passed")

    # ad = load_mMAMP()
    # print(ad.obs)
    # anndata_visualization(ad, 'mMAMP_anterior', col_name='ground_truth', spot_size=80)
    # print(ad.obs['ground_truth'].nunique())
