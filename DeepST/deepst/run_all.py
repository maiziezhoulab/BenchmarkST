import os 
from DeepST import run
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.metrics import adjusted_rand_score
import numpy as np
from st_loading_utils import load_DLPFC, load_BC, load_mVC, load_mPFC, load_mHypothalamus, load_her2_tumor, load_mMAMP, load_spacelhBC

iters = 1 # for boxplotting

save_dir_r = './clustering_result/'


if not os.path.exists(save_dir_r):
    os.makedirs(save_dir_r)

# """DLPFC"""
# # the number of clusters
# # [7, '151507'], [7, '151508'], [7, '151509'], [7, '151510'], [5, '151669']
# setting_combinations = [[7, '151507'], [7, '151508'], [7, '151509'], [7, '151510'], [5, '151669'], [5, '151670'], [5, '151671'], [5, '151672'], [7, '151673'], [7, '151674'], [7, '151675'], [7, '151676']]
# # setting_combinations = [[7, '151674'], [7, '151675'], [7, '151676']]
# for setting_combi in setting_combinations:
#    n_clusters = setting_combi[0]  # 7

#    dataset = setting_combi[1]  # '151673'
#    aris = []
#    dir_ = '/data/maiziezhou_lab/Datasets/ST_datasets/DLPFC12'

#    save_path = '../results/' + dataset + '/'
#    deepen = run(save_path = save_path, 
#                   platform = "Visium",
#                   pca_n_comps = 200,
#                   pre_epochs = 800, #### According to your own hardware, choose the number of training
#                   epochs = 1000, #### According to your own hardware, choose the number of training
#                   Conv_type="GCNConv", #### you can choose GNN types. 
#                   )
#    adata_ = deepen._get_adata(dir_, dataset)
   
#    for iter_ in range(iters):
#       adata = deepen._get_augment(adata_, adjacent_weight = 0.3, neighbour_k = 4)
#       graph_dict = deepen._get_graph(adata.obsm["spatial"], distType="BallTree", k=12)
#       adata = deepen._fit(adata, graph_dict, pretrain = False)
#       adata = deepen._get_cluster_data(adata, n_domains = n_clusters, priori=True) ###### without using prior knowledge, setting priori = False.
#       print(adata.obs)

#       ARI = adjusted_rand_score(adata.obs["DeepST_refine_domain"], adata.obs["original_clusters"])
#       save_df = adata.obs[['original_clusters', 'DeepST_refine_domain']]
#       save_df.to_csv(os.path.join(save_dir_r, 'deepst_' + dataset + '_iter_' + str(iter_) + '_out.csv'), index=True)
#       aris.append(ARI)
#       print(iter_)
#       print('Dataset:', dataset)
#       print(ARI)
#    print('Dataset:', dataset)
#    print(aris)
#    print(np.mean(aris))
#    with open('deepst_aris.txt', 'a+') as fp:
#       fp.write('DLPFC' + dataset + ' ')
#       fp.write(' '.join([str(i) for i in aris]))
#       fp.write('\n')
# # exit(-1)

# """BC"""
# # the number of clusters
# setting_combinations = [[20, 'section1']]
# for setting_combi in setting_combinations:
#    n_clusters = setting_combi[0]  # 7

#    dataset = setting_combi[1]  #
#    aris = []
#    dir_ = '/data/maiziezhou_lab/Datasets/ST_datasets/BC'

#    save_path = '../results/' + dataset + '/'
#    deepen = run(save_path = save_path, 
#                    platform = "Visium",
#                    pca_n_comps = 200,
#                    pre_epochs = 800, #### According to your own hardware, choose the number of training
#                    epochs = 1000, #### According to your own hardware, choose the number of training
#                    Conv_type="GCNConv", #### you can choose GNN types. 
#                    )
#    adata_ = deepen._get_adata(dir_, dataset)
   
#    for iter_ in range(iters):
#       adata = deepen._get_augment(adata_, adjacent_weight = 0.3, neighbour_k = 4)
#       graph_dict = deepen._get_graph(adata.obsm["spatial"], distType="BallTree", k=12)
#       adata = deepen._fit(adata, graph_dict, pretrain = False)
#       adata = deepen._get_cluster_data(adata, n_domains = n_clusters, priori=True) ###### without using prior knowledge, setting priori = False.
#       print(adata.obs)

#       ARI = adjusted_rand_score(adata.obs["DeepST_refine_domain"], adata.obs["original_clusters"])
#       save_df = adata.obs[['original_clusters', 'DeepST_refine_domain']]
#       save_df.to_csv(os.path.join(save_dir_r, 'deepst_' + dataset + '_iter_' + str(iter_) + '_out.csv'), index=True)
#       aris.append(ARI)
#       print(iter_)
#       print('Dataset:', dataset)
#       print(ARI)
#    print('Dataset:', dataset)
#    print(aris)
#    print(np.mean(aris))
#    with open('deepst_aris.txt', 'a+') as fp:
#       fp.write('HBRC1 ')
#       fp.write(' '.join([str(i) for i in aris]))
#       fp.write('\n')


# """load mMAMP ma section"""
# # the number of clusters
# # 6       5       4       4       4       4       7       7
# setting_combinations = [[52, 'MA']]
# for setting_combi in setting_combinations:
#    n_clusters = setting_combi[0]  # 7

#    dataset = setting_combi[1]  # '151673'
#    aris = []
#    dir_ = '/data/maiziezhou_lab/Datasets/ST_datasets/mMAMP'

#    save_path = '../results/' + dataset + '/'
#    deepen = run(save_path = save_path, 
#                    platform = "Visium",
#                    pca_n_comps = 200,
#                    pre_epochs = 800, #### According to your own hardware, choose the number of training
#                    epochs = 1000, #### According to your own hardware, choose the number of training
#                    Conv_type="GCNConv", #### you can choose GNN types. 
#                    )
#    adata_ = deepen._get_adata(dir_, dataset)
   
#    for iter_ in range(iters):
#       adata = deepen._get_augment(adata_, adjacent_weight = 0.3, neighbour_k = 4)
#       graph_dict = deepen._get_graph(adata.obsm["spatial"], distType="BallTree", k=12)
#       adata = deepen._fit(adata, graph_dict, pretrain = False)
#       adata = deepen._get_cluster_data(adata, n_domains = n_clusters, priori=True) ###### without using prior knowledge, setting priori = False.
#       print(adata.obs)

#       ARI = adjusted_rand_score(adata.obs["DeepST_refine_domain"], adata.obs["original_clusters"])
#       save_df = adata.obs[['original_clusters', 'DeepST_refine_domain']]
#       save_df.to_csv(os.path.join(save_dir_r, 'deepst_' + dataset + '_iter_' + str(iter_) + '_out.csv'), index=True)
#       aris.append(ARI)
#       print(iter_)
#       print('Dataset:', dataset)
#       print(ARI)
#    print('Dataset:', dataset)
#    print(aris)
#    print(np.mean(aris))
#    with open('deepst_aris.txt', 'a+') as fp:
#       fp.write('mABC ')
#       fp.write(' '.join([str(i) for i in aris]))
#       fp.write('\n')


# """mVC"""
# # the number of clusters
# setting_combinations = [[7, 'STARmap_20180505_BY3_1k.h5ad']]
# for setting_combi in setting_combinations:
#    n_clusters = setting_combi[0]

#    dataset = setting_combi[1]
#    aris = []
#    dir_ = '/data/maiziezhou_lab/Datasets/ST_datasets/STARmap_mouse_visual_cortex'
   
#    save_path = '../results/' + dataset + '/'
#    deepen = run(save_path = save_path, 
#                    platform = "benchmark_test",
#                    pca_n_comps = 200,
#                    pre_epochs = 800, #### According to your own hardware, choose the number of training
#                    epochs = 1000, #### According to your own hardware, choose the number of training
#                    )
#    # adata_ = deepen._get_adata(dir_, dataset)
#    adata_, graph_dict = deepen._get_single_adata(dir_, dataset, weights="weights_matrix_nomd") #### Augmentation without using morphological information
#    # adata_ = deepen._get_augment(adata_, adjacent_weight = 0.3, neighbour_k = 4)
#    # graph_dict = deepen._get_graph(adata_.obsm["spatial"], distType="BallTree", k=12)
#    for iter_ in range(iters):
      
      
#       adata = deepen._fit(adata_, graph_dict, pretrain = False)
#       adata = deepen._get_cluster_data(adata, n_domains = n_clusters, priori=True) ###### without using prior knowledge, setting priori = False.
#       print(adata.obs)

#       ARI = adjusted_rand_score(adata.obs["DeepST_refine_domain"], adata.obs["original_clusters"])
#       save_df = adata.obs[['original_clusters', 'DeepST_refine_domain']]
#       save_df.to_csv(os.path.join(save_dir_r, 'deepst_' + dataset + '_iter_' + str(iter_) + '_out.csv'), index=True)
#       aris.append(ARI)
#       print(iter_)
#       print('Dataset:', dataset)
#       print(ARI)
#    print('Dataset:', dataset)
#    print(aris)
#    print(np.mean(aris))
#    with open('deepst_aris.txt', 'a+') as fp:
#       fp.write('mVC ')
#       fp.write(' '.join([str(i) for i in aris]))
#       fp.write('\n')


# """mPFC"""
# # the number of clusters
# setting_combinations = [[4, '20180417_BZ5_control'], [4, '20180419_BZ9_control'], [4, '20180424_BZ14_control']]
# for setting_combi in setting_combinations:
#    n_clusters = setting_combi[0]

#    dataset = setting_combi[1]
#    aris = []
#    dir_ = '/data/maiziezhou_lab/Datasets/ST_datasets/STARmap_mouse_PFC'
   
#    save_path = '../results/' + dataset + '/'
#    deepen = run(save_path = save_path, 
#                    platform = "benchmark_test",
#                    pca_n_comps = 200,
#                    pre_epochs = 800, #### According to your own hardware, choose the number of training
#                    epochs = 1000, #### According to your own hardware, choose the number of training
#                    )
#    # adata_ = deepen._get_adata(dir_, dataset)
#    adata_, graph_dict = deepen._get_single_adata(dir_, dataset, weights="weights_matrix_nomd") #### Augmentation without using morphological information
#    # adata_ = deepen._get_augment(adata_, adjacent_weight = 0.3, neighbour_k = 4)
#    # graph_dict = deepen._get_graph(adata_.obsm["spatial"], distType="BallTree", k=12)
#    for iter_ in range(iters):
      
      
#       adata = deepen._fit(adata_, graph_dict, pretrain = False)
#       adata = deepen._get_cluster_data(adata, n_domains = n_clusters, priori=True) ###### without using prior knowledge, setting priori = False.
#       print(adata.obs)

#       ARI = adjusted_rand_score(adata.obs["DeepST_refine_domain"], adata.obs["original_clusters"])
#       save_df = adata.obs[['original_clusters', 'DeepST_refine_domain']]
#       save_df.to_csv(os.path.join(save_dir_r, 'deepst_' + dataset + '_iter_' + str(iter_) + '_out.csv'), index=True)
#       aris.append(ARI)
#       print(iter_)
#       print('Dataset:', dataset)
#       print(ARI)
#    print('Dataset:', dataset)
#    print(aris)
#    print(np.mean(aris))
#    with open('spagcn_aris.txt', 'a+') as fp:
#       fp.write('mPFC' + dataset + ' ')
#       fp.write(' '.join([str(i) for i in aris]))
#       fp.write('\n')


# """mHypo"""
# # the number of clusters
# # 15      15      14      15      15      15      14       15       15       15      16        15
# # setting_combinations = [[15, '0.26'], [15, '0.21'], [14, '0.16'], [15, '0.11'], [15, '0.06'], [15, '0.01'], [14, '-0.04'], [15, '-0.09'], [15, '-0.14'], [15, '-0.19'], [16, '-0.24'], [15, '-0.29']]
# setting_combinations = [[8, '-0.04'], [8, '-0.09'], [8, '-0.14'], [8, '-0.19'], [8, '-0.24'], [8, '-0.29']]
# # setting_combinations = [[8, '-0.24'], [15, '-0.29']]
# for setting_combi in setting_combinations:
#    n_clusters = setting_combi[0]  # 7

#    dataset = setting_combi[1]  #
#    aris = []
#    dir_ = '/data/maiziezhou_lab/Datasets/ST_datasets/mHypothalamus'
   
#    save_path = '../results/' + dataset + '/'
#    deepen = run(save_path = save_path, 
#                    platform = "benchmark_test",
#                    pca_n_comps = 200,
#                    pre_epochs = 800, #### According to your own hardware, choose the number of training
#                    epochs = 1000, #### According to your own hardware, choose the number of training
#                    )
#    # adata_ = deepen._get_adata(dir_, dataset)
#    adata_, graph_dict = deepen._get_single_adata(dir_, dataset, weights="weights_matrix_nomd") #### Augmentation without using morphological information
#    # adata_ = deepen._get_augment(adata_, adjacent_weight = 0.3, neighbour_k = 4)
#    # graph_dict = deepen._get_graph(adata_.obsm["spatial"], distType="BallTree", k=12)
#    for iter_ in range(iters):
      
      
#       adata = deepen._fit(adata_, graph_dict, pretrain = False)
#       adata = deepen._get_cluster_data(adata, n_domains = n_clusters, priori=True) ###### without using prior knowledge, setting priori = False.
#       print(adata.obs)

#       ARI = adjusted_rand_score(adata.obs["DeepST_refine_domain"], adata.obs["original_clusters"])
#       save_df = adata.obs[['original_clusters', 'DeepST_refine_domain']]
#       save_df.to_csv(os.path.join(save_dir_r, 'deepst_' + dataset + '_iter_' + str(iter_) + '_out.csv'), index=True)
#       aris.append(ARI)
#       print(iter_)
#       print('Dataset:', dataset)
#       print(ARI)
#    print('Dataset:', dataset)
#    print(aris)
#    print(np.mean(aris))
#    with open('deepst_aris.txt', 'a+') as fp:
#       fp.write('mHypothalamus' + dataset + ' ')
#       fp.write(' '.join([str(i) for i in aris]))
#       fp.write('\n')


"""spacel hBC"""
# the number of clusters
# 10 for all
setting_combinations = [[10, 'human_bc_spatial_1142243F'], 
                        [10, 'human_bc_spatial_1160920F'], 
                        [10, 'human_bc_spatial_CID4290'], 
                        [10, 'human_bc_spatial_CID4465'], 
                        [10, 'human_bc_spatial_CID4535'], 
                        [10, 'human_bc_spatial_CID44971'], 
                        [10, 'human_bc_spatial_Parent_Visium_Human_BreastCancer'], 
                        [10, 'human_bc_spatial_V1_Breast_Cancer_Block_A_Section_1'],
                        [10, 'human_bc_spatial_V1_Breast_Cancer_Block_A_Section_2'],
                        [10, 'human_bc_spatial_V1_Human_Invasive_Ductal_Carcinoma'],
                        [10, 'human_bc_spatial_Visium_FFPE_Human_Breast_Cancer']]
for setting_combi in setting_combinations:
   n_clusters = setting_combi[0]  # 7

   dataset = setting_combi[1]  # '151673'
   aris = []
   # dir_ = '/data/maiziezhou_lab/Datasets/ST_datasets/Her2_tumor'
   dir_ = '/home/yunfei/spatial_benchmarking/benchmarking_data/visium_human_breast_cancer'
   
   save_path = '../results/' + dataset + '/'
   deepen = run(save_path = save_path, 
                   platform = "benchmark_test",
                   pca_n_comps = 200,
                   pre_epochs = 800, #### According to your own hardware, choose the number of training
                   epochs = 1000, #### According to your own hardware, choose the number of training
                   )
   # adata_ = deepen._get_adata(dir_, dataset)
   adata_, graph_dict = deepen._get_single_adata(dir_, dataset, weights="weights_matrix_nomd") #### Augmentation without using morphological information
   # adata_ = deepen._get_augment(adata_, adjacent_weight = 0.3, neighbour_k = 4)
   # graph_dict = deepen._get_graph(adata_.obsm["spatial"], distType="BallTree", k=12)
   for iter_ in range(iters):
      
      
      adata = deepen._fit(adata_, graph_dict, pretrain = False)
      adata = deepen._get_cluster_data(adata, n_domains = n_clusters, priori=True) ###### without using prior knowledge, setting priori = False.
      print(adata.obs)

      ARI = adjusted_rand_score(adata.obs["DeepST_refine_domain"], adata.obs["original_clusters"])
      save_df = adata.obs[['original_clusters', 'DeepST_refine_domain']]
      save_df.to_csv(os.path.join(save_dir_r, 'deepst_' + dataset + '_iter_' + str(iter_) + '_out.csv'), index=True)
      aris.append(ARI)
      print(iter_)
      print('Dataset:', dataset)
      print(ARI)
   print('Dataset:', dataset)
   print(aris)
   print(np.mean(aris))
   with open('deepst_aris.txt', 'a+') as fp:
      fp.write('spacel_hBC' + dataset + ' ')
      fp.write(' '.join([str(i) for i in aris]))
      fp.write('\n')