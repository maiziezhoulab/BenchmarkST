import torch
import argparse
import random
import numpy as np
import pandas as pd
from src.graph_func import graph_construction
from src.utils_func import mk_dir, adata_preprocess, load_ST_file, res_search_fixed_clus, plot_clustering
from src.training import conST_training

import anndata
from sklearn import metrics
import matplotlib.pyplot as plt
import scanpy as sc
import os
import warnings
warnings.filterwarnings('ignore')
from st_loading_utils import load_DLPFC, load_BC, load_mVC, load_mPFC, load_mHypothalamus, load_her2_tumor, load_mMAMP


parser = argparse.ArgumentParser()
parser.add_argument('--k', type=int, default=10, help='parameter k in spatial graph')
parser.add_argument('--knn_distanceType', type=str, default='euclidean',
                    help='graph distance type: euclidean/cosine/correlation')
parser.add_argument('--epochs', type=int, default=200, help='Number of epochs to train.')
parser.add_argument('--cell_feat_dim', type=int, default=300, help='Dim of PCA')
parser.add_argument('--feat_hidden1', type=int, default=100, help='Dim of DNN hidden 1-layer.')
parser.add_argument('--feat_hidden2', type=int, default=20, help='Dim of DNN hidden 2-layer.')
parser.add_argument('--gcn_hidden1', type=int, default=32, help='Dim of GCN hidden 1-layer.')
parser.add_argument('--gcn_hidden2', type=int, default=8, help='Dim of GCN hidden 2-layer.')
parser.add_argument('--p_drop', type=float, default=0.2, help='Dropout rate.')
parser.add_argument('--use_img', type=bool, default=False, help='Use histology images.')
parser.add_argument('--img_w', type=float, default=0.1, help='Weight of image features.')
parser.add_argument('--use_pretrained', type=bool, default=True, help='Use pretrained weights.')
parser.add_argument('--using_mask', type=bool, default=False, help='Using mask for multi-dataset.')
parser.add_argument('--feat_w', type=float, default=10, help='Weight of DNN loss.')
parser.add_argument('--gcn_w', type=float, default=0.1, help='Weight of GCN loss.')
parser.add_argument('--dec_kl_w', type=float, default=10, help='Weight of DEC loss.')
parser.add_argument('--gcn_lr', type=float, default=0.01, help='Initial GNN learning rate.')
parser.add_argument('--gcn_decay', type=float, default=0.01, help='Initial decay rate.')
parser.add_argument('--dec_cluster_n', type=int, default=10, help='DEC cluster number.')
parser.add_argument('--dec_interval', type=int, default=20, help='DEC interval nnumber.')
parser.add_argument('--dec_tol', type=float, default=0.00, help='DEC tol.')

parser.add_argument('--seed', type=int, default=0, help='random seed')
parser.add_argument('--beta', type=float, default=100, help='beta value for l2c')
parser.add_argument('--cont_l2l', type=float, default=0.3, help='Weight of local contrastive learning loss.')
parser.add_argument('--cont_l2c', type=float, default= 0.1, help='Weight of context contrastive learning loss.')
parser.add_argument('--cont_l2g', type=float, default= 0.1, help='Weight of global contrastive learning loss.')

parser.add_argument('--edge_drop_p1', type=float, default=0.1, help='drop rate of adjacent matrix of the first view')
parser.add_argument('--edge_drop_p2', type=float, default=0.1, help='drop rate of adjacent matrix of the second view')
parser.add_argument('--node_drop_p1', type=float, default=0.2, help='drop rate of node features of the first view')
parser.add_argument('--node_drop_p2', type=float, default=0.3, help='drop rate of node features of the second view')

# ______________ Eval clustering Setting ______________
parser.add_argument('--eval_resolution', type=int, default=1, help='Eval cluster number.')
parser.add_argument('--eval_graph_n', type=int, default=20, help='Eval graph kN tol.') 

params =  parser.parse_args(args=['--k', '20', '--knn_distanceType', 'euclidean', '--epochs', '200'])

np.random.seed(params.seed)
torch.manual_seed(params.seed)
torch.cuda.manual_seed(params.seed)
device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
print('Using device: ' + device)
params.device = device

def seed_torch(seed):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

iters=20
# """DLPFC"""
# # the number of clusters
# setting_combinations = [[7, '151507'], [7, '151508'], [7, '151509'], [7, '151510'], [5, '151669'], [5, '151670'], [5, '151671'], [5, '151672'], [7, '151673'], [7, '151674'], [7, '151675'], [7, '151676']]
# # setting_combinations = [[7, '151674'], [7, '151675'], [7, '151676']]
# for setting_combi in setting_combinations:
#     path = '/home/yunfei/spatial_benchmarking/benchmarking_data/DLPFC12'
#     adata_h5 = load_DLPFC(root_dir=path, section_id=setting_combi[1])
#     adata_X = adata_preprocess(adata_h5, min_cells=5, pca_n_comps=params.cell_feat_dim)
#     graph_dict = graph_construction(adata_h5.obsm['spatial'], adata_h5.shape[0], params)
    
#     dataset = data_name = setting_combi[1]
#     n_clusters = setting_combi[0]
#     aris = []
#     save_root = './output/spatialLIBD/'
#     # data_root = '../spatialLIBD'
#     params.save_path = mk_dir(f'{save_root}/{data_name}/conST')

#     params.cell_num = adata_h5.shape[0]

#     for iter_ in range(iters):
#         seed_torch(params.seed)
        
#         if params.use_img:
#             img_transformed = np.load('./MAE-pytorch/extracted_feature.npy')
#             img_transformed = (img_transformed - img_transformed.mean()) / img_transformed.std() * adata_X.std() + adata_X.mean()
#             conST_net = conST_training(adata_X, graph_dict, params, n_clusters, img_transformed)
#         else:
#             conST_net = conST_training(adata_X, graph_dict, params, n_clusters)

#         conST_net.pretraining()
#         conST_net.major_training()

#         conST_embedding = conST_net.get_embedding()

#         # np.save(f'{params.save_path}/conST_result.npy', conST_embedding)
#         # clustering
#         adata_conST = anndata.AnnData(conST_embedding, obs=adata_h5.obs)
#         adata_conST.uns['spatial'] = adata_h5.uns['spatial']
#         adata_conST.obs['original_clusters'] = adata_h5.obs['original_clusters']
#         adata_conST.obsm['spatial'] = adata_h5.obsm['spatial']

#         sc.pp.neighbors(adata_conST, n_neighbors=params.eval_graph_n)

#         eval_resolution = res_search_fixed_clus(adata_conST, n_clusters)
#         print(eval_resolution)
#         cluster_key = "conST_leiden"
#         sc.tl.leiden(adata_conST, key_added=cluster_key, resolution=eval_resolution)

#         keep_bcs = adata_conST.obs.dropna().index
#         adata_conST = adata_conST[keep_bcs].copy()
#         ARI = metrics.adjusted_rand_score(adata_conST.obs[cluster_key], adata_conST.obs['original_clusters'])

#         print('Dataset:', dataset)
#         print('ARI:', ARI)
#         aris.append(ARI)
#     print('Dataset:', dataset)
#     print(aris)
#     print(np.mean(aris))
#     with open('const_aris.txt', 'a+') as fp:
#         fp.write('DLPFC' + dataset + ' ')
#         fp.write(' '.join([str(i) for i in aris]))
#         fp.write('\n')


# """BC"""
# # the number of clusters
# setting_combinations = [[20, 'section1']]
# for setting_combi in setting_combinations:
#     path = '/home/yunfei/spatial_benchmarking/benchmarking_data/BC'
#     adata_h5 = load_DLPFC(root_dir=path, section_id=setting_combi[1])
#     adata_X = adata_preprocess(adata_h5, min_cells=5, pca_n_comps=params.cell_feat_dim)
#     graph_dict = graph_construction(adata_h5.obsm['spatial'], adata_h5.shape[0], params)

#     dataset = data_name = setting_combi[1]
#     n_clusters = setting_combi[0]
#     aris = []
#     save_root = './output/BC/'
#     # data_root = '../BC'
#     params.save_path = mk_dir(f'{save_root}/{data_name}/conST')

#     params.cell_num = adata_h5.shape[0]

#     for iter_ in range(iters):
#         seed_torch(params.seed)
        
#         if params.use_img:
#             img_transformed = np.load('./MAE-pytorch/extracted_feature.npy')
#             img_transformed = (img_transformed - img_transformed.mean()) / img_transformed.std() * adata_X.std() + adata_X.mean()
#             conST_net = conST_training(adata_X, graph_dict, params, n_clusters, img_transformed)
#         else:
#             conST_net = conST_training(adata_X, graph_dict, params, n_clusters)

#         conST_net.pretraining()
#         conST_net.major_training()

#         conST_embedding = conST_net.get_embedding()

#         # np.save(f'{params.save_path}/conST_result.npy', conST_embedding)
#         # clustering
#         adata_conST = anndata.AnnData(conST_embedding, obs=adata_h5.obs)
#         adata_conST.uns['spatial'] = adata_h5.uns['spatial']
#         adata_conST.obs['original_clusters'] = adata_h5.obs['original_clusters']
#         adata_conST.obsm['spatial'] = adata_h5.obsm['spatial']

#         sc.pp.neighbors(adata_conST, n_neighbors=params.eval_graph_n)

#         eval_resolution = res_search_fixed_clus(adata_conST, n_clusters)
#         print(eval_resolution)
#         cluster_key = "conST_leiden"
#         sc.tl.leiden(adata_conST, key_added=cluster_key, resolution=eval_resolution)

#         keep_bcs = adata_conST.obs.dropna().index
#         adata_conST = adata_conST[keep_bcs].copy()
#         ARI = metrics.adjusted_rand_score(adata_conST.obs[cluster_key], adata_conST.obs['original_clusters'])

#         print('Dataset:', dataset)
#         print('ARI:', ARI)
#         aris.append(ARI)
#     print('Dataset:', dataset)
#     print(aris)
#     print(np.mean(aris))
#     with open('const_aris.txt', 'a+') as fp:
#         fp.write('BC' + dataset + ' ')
#         fp.write(' '.join([str(i) for i in aris]))
#         fp.write('\n')


# """MA"""
# # the number of clusters
# setting_combinations = [[52, 'MA']]
# for setting_combi in setting_combinations:
#     path = '/home/yunfei/spatial_benchmarking/benchmarking_data/mMAMP'
#     adata_h5 = load_mMAMP(root_dir=path, section_id=setting_combi[1])
#     adata_X = adata_preprocess(adata_h5, min_cells=5, pca_n_comps=params.cell_feat_dim)
#     graph_dict = graph_construction(adata_h5.obsm['spatial'], adata_h5.shape[0], params)

#     dataset = data_name = setting_combi[1]
#     n_clusters = setting_combi[0]
#     aris = []
#     save_root = './output/MA/'
#     # data_root = '../BC'
#     params.save_path = mk_dir(f'{save_root}/{data_name}/conST')

#     params.cell_num = adata_h5.shape[0]

#     for iter_ in range(iters):
#         seed_torch(params.seed)
        
#         if params.use_img:
#             img_transformed = np.load('./MAE-pytorch/extracted_feature.npy')
#             img_transformed = (img_transformed - img_transformed.mean()) / img_transformed.std() * adata_X.std() + adata_X.mean()
#             conST_net = conST_training(adata_X, graph_dict, params, n_clusters, img_transformed)
#         else:
#             conST_net = conST_training(adata_X, graph_dict, params, n_clusters)

#         conST_net.pretraining()
#         conST_net.major_training()

#         conST_embedding = conST_net.get_embedding()

#         # np.save(f'{params.save_path}/conST_result.npy', conST_embedding)
#         # clustering
#         adata_conST = anndata.AnnData(conST_embedding, obs=adata_h5.obs)
#         adata_conST.uns['spatial'] = adata_h5.uns['spatial']
#         adata_conST.obs['original_clusters'] = adata_h5.obs['original_clusters']
#         adata_conST.obsm['spatial'] = adata_h5.obsm['spatial']

#         sc.pp.neighbors(adata_conST, n_neighbors=params.eval_graph_n)

#         eval_resolution = res_search_fixed_clus(adata_conST, n_clusters)
#         print(eval_resolution)
#         cluster_key = "conST_leiden"
#         sc.tl.leiden(adata_conST, key_added=cluster_key, resolution=eval_resolution)

#         keep_bcs = adata_conST.obs.dropna().index
#         adata_conST = adata_conST[keep_bcs].copy()
#         ARI = metrics.adjusted_rand_score(adata_conST.obs[cluster_key], adata_conST.obs['original_clusters'])

#         print('Dataset:', dataset)
#         print('ARI:', ARI)
#         aris.append(ARI)
#     print('Dataset:', dataset)
#     print(aris)
#     print(np.mean(aris))
#     with open('const_aris.txt', 'a+') as fp:
#         fp.write('mAB ')
#         fp.write(' '.join([str(i) for i in aris]))
#         fp.write('\n')


# """Her2st"""
# # the number of clusters
# setting_combinations = [[5, 'B1'], [4, 'C1'], [4, 'D1'], [4, 'E1'], [4, 'F1'], [7, 'G2'], [7, 'H1']]
# # setting_combinations = [[7, '151674'], [7, '151675'], [7, '151676']] [6, 'A1'], 
# for setting_combi in setting_combinations:
#     args = parser.parse_args()
#     # seed
#     seed_torch(1)

#     path = args.path = '/home/yunfei/spatial_benchmarking/benchmarking_data/Her2_tumor'
#     adata_h5 = load_her2_tumor(root_dir=path, section_id=setting_combi[1])
#     if params.cell_feat_dim > len(adata_h5.obs.index):
#         params.cell_feat_dim = len(adata_h5.obs.index)-1
#         print(params.cell_feat_dim)
#     adata_X = adata_preprocess(adata_h5, min_cells=5, pca_n_comps=params.cell_feat_dim)
#     graph_dict = graph_construction(adata_h5.obsm['spatial'], adata_h5.shape[0], params)

#     dataset = data_name = setting_combi[1]
#     n_clusters = setting_combi[0]
#     aris = []
#     save_root = './output/her2tumor/'
#     # data_root = '../BC'
#     params.save_path = mk_dir(f'{save_root}/{data_name}/conST')

#     params.cell_num = adata_h5.shape[0]

#     for iter_ in range(iters):
#         seed_torch(params.seed)
        
#         if params.use_img:
#             img_transformed = np.load('./MAE-pytorch/extracted_feature.npy')
#             img_transformed = (img_transformed - img_transformed.mean()) / img_transformed.std() * adata_X.std() + adata_X.mean()
#             conST_net = conST_training(adata_X, graph_dict, params, n_clusters, img_transformed)
#         else:
#             conST_net = conST_training(adata_X, graph_dict, params, n_clusters)

#         conST_net.pretraining()
#         conST_net.major_training()

#         conST_embedding = conST_net.get_embedding()

#         # np.save(f'{params.save_path}/conST_result.npy', conST_embedding)
#         # clustering
#         adata_conST = anndata.AnnData(conST_embedding, obs=adata_h5.obs)
#         # adata_conST.uns['spatial'] = adata_h5.uns['spatial']
#         adata_conST.obs['original_clusters'] = adata_h5.obs['original_clusters']
#         adata_conST.obsm['spatial'] = adata_h5.obsm['spatial']

#         sc.pp.neighbors(adata_conST, n_neighbors=params.eval_graph_n)

#         eval_resolution = res_search_fixed_clus(adata_conST, n_clusters)
#         print(eval_resolution)
#         cluster_key = "conST_leiden"
#         sc.tl.leiden(adata_conST, key_added=cluster_key, resolution=eval_resolution)

#         keep_bcs = adata_conST.obs.dropna().index
#         adata_conST = adata_conST[keep_bcs].copy()
#         ARI = metrics.adjusted_rand_score(adata_conST.obs[cluster_key], adata_conST.obs['original_clusters'])

#         print('Dataset:', dataset)
#         print('ARI:', ARI)
#         aris.append(ARI)
#     print('Dataset:', dataset)
#     print(aris)
#     print(np.mean(aris))
#     with open('const_aris.txt', 'a+') as fp:
#         fp.write('Her2tumor' + dataset + ' ')
#         fp.write(' '.join([str(i) for i in aris]))
#         fp.write('\n')


# """mVC"""
# # the number of clusters
# setting_combinations = [[7, 'STARmap_20180505_BY3_1k.h5ad']]
# for setting_combi in setting_combinations:
#     args = parser.parse_args()
#     # seed
#     seed_torch(1)

#     path = args.path = '/home/yunfei/spatial_benchmarking/benchmarking_data/STARmap_mouse_visual_cortex'
#     adata_h5 = load_mVC(root_dir=path, section_id=setting_combi[1])
#     adata_X = adata_preprocess(adata_h5, min_cells=5, pca_n_comps=params.cell_feat_dim)
#     graph_dict = graph_construction(adata_h5.obsm['spatial'], adata_h5.shape[0], params)

#     dataset = data_name = setting_combi[1]
#     n_clusters = setting_combi[0]
#     aris = []
#     save_root = './output/her2tumor/'
#     # data_root = '../BC'
#     params.save_path = mk_dir(f'{save_root}/{data_name}/conST')

#     params.cell_num = adata_h5.shape[0]

#     for iter_ in range(iters):
#         seed_torch(params.seed)
        
#         if params.use_img:
#             img_transformed = np.load('./MAE-pytorch/extracted_feature.npy')
#             img_transformed = (img_transformed - img_transformed.mean()) / img_transformed.std() * adata_X.std() + adata_X.mean()
#             conST_net = conST_training(adata_X, graph_dict, params, n_clusters, img_transformed)
#         else:
#             conST_net = conST_training(adata_X, graph_dict, params, n_clusters)

#         conST_net.pretraining()
#         conST_net.major_training()

#         conST_embedding = conST_net.get_embedding()

#         # np.save(f'{params.save_path}/conST_result.npy', conST_embedding)
#         # clustering
#         adata_conST = anndata.AnnData(conST_embedding, obs=adata_h5.obs)
#         # adata_conST.uns['spatial'] = adata_h5.uns['spatial']
#         adata_conST.obs['original_clusters'] = adata_h5.obs['original_clusters']
#         adata_conST.obsm['spatial'] = adata_h5.obsm['spatial']

#         sc.pp.neighbors(adata_conST, n_neighbors=params.eval_graph_n)

#         eval_resolution = res_search_fixed_clus(adata_conST, n_clusters)
#         print(eval_resolution)
#         cluster_key = "conST_leiden"
#         sc.tl.leiden(adata_conST, key_added=cluster_key, resolution=eval_resolution)

#         keep_bcs = adata_conST.obs.dropna().index
#         adata_conST = adata_conST[keep_bcs].copy()
#         ARI = metrics.adjusted_rand_score(adata_conST.obs[cluster_key], adata_conST.obs['original_clusters'])

#         print('Dataset:', dataset)
#         print('ARI:', ARI)
#         aris.append(ARI)
#     print('Dataset:', dataset)
#     print(aris)
#     print(np.mean(aris))
#     with open('const_aris.txt', 'a+') as fp:
#         fp.write('mVC ')
#         fp.write(' '.join([str(i) for i in aris]))
#         fp.write('\n')


# """mPFC"""
# # the number of clusters
# setting_combinations = [[4, '20180417_BZ5_control'], [4, '20180419_BZ9_control'], [4, '20180424_BZ14_control']]
# for setting_combi in setting_combinations:
#     args = parser.parse_args()
#     # seed
#     seed_torch(1)
#     path = args.path = '/home/yunfei/spatial_benchmarking/benchmarking_data/STARmap_mouse_PFC'
#     adata_h5 = load_mPFC(root_dir=path, section_id=setting_combi[1])
#     if params.cell_feat_dim > len(adata_h5.var.index):
#         params.cell_feat_dim = len(adata_h5.var.index)-1
#         print(params.cell_feat_dim)
#     adata_X = adata_preprocess(adata_h5, min_cells=5, pca_n_comps=params.cell_feat_dim)
#     graph_dict = graph_construction(adata_h5.obsm['spatial'], adata_h5.shape[0], params)

#     dataset = data_name = setting_combi[1]
#     n_clusters = setting_combi[0]
#     aris = []
#     save_root = './output/her2tumor/'
#     # data_root = '../BC'
#     params.save_path = mk_dir(f'{save_root}/{data_name}/conST')

#     params.cell_num = adata_h5.shape[0]

#     for iter_ in range(iters):
#         seed_torch(params.seed)
        
#         if params.use_img:
#             img_transformed = np.load('./MAE-pytorch/extracted_feature.npy')
#             img_transformed = (img_transformed - img_transformed.mean()) / img_transformed.std() * adata_X.std() + adata_X.mean()
#             conST_net = conST_training(adata_X, graph_dict, params, n_clusters, img_transformed)
#         else:
#             conST_net = conST_training(adata_X, graph_dict, params, n_clusters)

#         conST_net.pretraining()
#         conST_net.major_training()

#         conST_embedding = conST_net.get_embedding()

#         # np.save(f'{params.save_path}/conST_result.npy', conST_embedding)
#         # clustering
#         adata_conST = anndata.AnnData(conST_embedding, obs=adata_h5.obs)
#         # adata_conST.uns['spatial'] = adata_h5.uns['spatial']
#         adata_conST.obs['original_clusters'] = adata_h5.obs['original_clusters']
#         adata_conST.obsm['spatial'] = adata_h5.obsm['spatial']

#         sc.pp.neighbors(adata_conST, n_neighbors=params.eval_graph_n)

#         eval_resolution = res_search_fixed_clus(adata_conST, n_clusters)
#         print(eval_resolution)
#         cluster_key = "conST_leiden"
#         sc.tl.leiden(adata_conST, key_added=cluster_key, resolution=eval_resolution)

#         keep_bcs = adata_conST.obs.dropna().index
#         adata_conST = adata_conST[keep_bcs].copy()
#         ARI = metrics.adjusted_rand_score(adata_conST.obs[cluster_key], adata_conST.obs['original_clusters'])

#         print('Dataset:', dataset)
#         print('ARI:', ARI)
#         aris.append(ARI)
#     print('Dataset:', dataset)
#     print(aris)
#     print(np.mean(aris))
#     with open('const_aris.txt', 'a+') as fp:
#         fp.write('mPFC' + dataset + ' ')
#         fp.write(' '.join([str(i) for i in aris]))
#         fp.write('\n')


"""mHypo"""
setting_combinations = [[8, '-0.04'], [8, '-0.09'], [8, '-0.14'], [8, '-0.19'], [8, '-0.24'], [8, '-0.29']]
for setting_combi in setting_combinations:
    args = parser.parse_args()
    # seed
    seed_torch(1)
    path = args.path = '/home/yunfei/spatial_benchmarking/benchmarking_data/mHypothalamus'
    adata_h5 = load_mHypothalamus(root_dir=path, section_id=setting_combi[1])
    if params.cell_feat_dim > len(adata_h5.var.index):
        params.cell_feat_dim = len(adata_h5.var.index)-1
        # print(params.cell_feat_dim)
    if params.cell_feat_dim > len(adata_h5.obs.index):
        params.cell_feat_dim = len(adata_h5.obs.index)-1
        # print(params.cell_feat_dim)
    adata_X = adata_preprocess(adata_h5, min_cells=5, pca_n_comps=params.cell_feat_dim)
    graph_dict = graph_construction(adata_h5.obsm['spatial'], adata_h5.shape[0], params)

    dataset = data_name = setting_combi[1]
    n_clusters = setting_combi[0]
    aris = []
    save_root = './output/her2tumor/'
    # data_root = '../BC'
    params.save_path = mk_dir(f'{save_root}/{data_name}/conST')

    params.cell_num = adata_h5.shape[0]

    for iter_ in range(iters):
        seed_torch(params.seed)
        
        if params.use_img:
            img_transformed = np.load('./MAE-pytorch/extracted_feature.npy')
            img_transformed = (img_transformed - img_transformed.mean()) / img_transformed.std() * adata_X.std() + adata_X.mean()
            conST_net = conST_training(adata_X, graph_dict, params, n_clusters, img_transformed)
        else:
            conST_net = conST_training(adata_X, graph_dict, params, n_clusters)

        conST_net.pretraining()
        conST_net.major_training()

        conST_embedding = conST_net.get_embedding()

        # np.save(f'{params.save_path}/conST_result.npy', conST_embedding)
        # clustering
        adata_conST = anndata.AnnData(conST_embedding, obs=adata_h5.obs)
        # adata_conST.uns['spatial'] = adata_h5.uns['spatial']
        adata_conST.obs['original_clusters'] = adata_h5.obs['original_clusters']
        adata_conST.obsm['spatial'] = adata_h5.obsm['spatial']

        sc.pp.neighbors(adata_conST, n_neighbors=params.eval_graph_n)

        eval_resolution = res_search_fixed_clus(adata_conST, n_clusters)
        print(eval_resolution)
        cluster_key = "conST_leiden"
        sc.tl.leiden(adata_conST, key_added=cluster_key, resolution=eval_resolution)

        keep_bcs = adata_conST.obs.dropna().index
        adata_conST = adata_conST[keep_bcs].copy()
        ARI = metrics.adjusted_rand_score(adata_conST.obs[cluster_key], adata_conST.obs['original_clusters'])

        print('Dataset:', dataset)
        print('ARI:', ARI)
        aris.append(ARI)
    print('Dataset:', dataset)
    print(aris)
    print(np.mean(aris))
    with open('const_aris.txt', 'a+') as fp:
        fp.write('mHypothalamus' + dataset + ' ')
        fp.write(' '.join([str(i) for i in aris]))
        fp.write('\n')