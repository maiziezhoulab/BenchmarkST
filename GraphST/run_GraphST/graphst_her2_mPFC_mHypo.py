import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
import numpy as np
from GraphST import GraphST
from GraphST.utils import clustering
from st_loading_utils import load_mPFC, load_mHypothalamus, load_her2_tumor, load_mMAMP

# Run device, by default, the package is implemented on 'cpu'. We recommend using GPU.
device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

# the location of R, which is necessary for mclust algorithm. Please replace the path below with local R installation path
# os.environ['R_HOME'] = '/home/yunfei/anaconda3/envs/GraphST/lib/R'
# os.environ['R_USER'] = '/home/yunfei/anaconda3/envs/GraphST/lib/python3.8/site-packages'

# """mPFC"""
# # the number of clusters
# setting_combinations = [[4, '20180417_BZ5_control'], [4, '20180419_BZ9_control'], [4, '20180424_BZ14_control']]
# for setting_combi in setting_combinations:
#    n_clusters = setting_combi[0]

#    dataset = setting_combi[1]
   
#    dir_ = '/home/yunfei/spatial_benchmarking/benchmarking_data/STARmap_mouse_PFC'
#    ad = load_mPFC(root_dir=dir_, section_id=dataset)

#    aris = []
#    for iter in range(5):

      
#       # print(ad)

#       # define model
#       model = GraphST.GraphST(ad, device=device)

#       # train model
#       ad = model.train()

#       # print(ad)

#       # set radius to specify the number of neighbors considered during refinement
#       radius = 50

#       tool = 'mclust' # mclust, leiden, and louvain

#       # clustering
#       if tool == 'mclust':
#          clustering(ad, n_clusters, radius=radius, method=tool, refinement=True) # For DLPFC dataset, we use optional refinement step.
#       elif tool in ['leiden', 'louvain']:
#          clustering(ad, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)

#       # filter out NA nodes
#       ad = ad[~pd.isnull(ad.obs['original_clusters'])]

#       # calculate metric ARI
#       ARI = metrics.adjusted_rand_score(ad.obs['domain'], ad.obs['original_clusters'])
#       ad.uns['ARI'] = ARI

#       # print('Dataset:', dataset)
#       # print('ARI:', ARI)
#       aris.append(ARI)
#    print('Dataset:', dataset)
#    print(aris)
#    print(np.mean(aris))
#    with open('graphst_aris.txt', 'a+') as fp:
#       fp.write('mPFC' + dataset + ' ')
#       fp.write(' '.join([str(i) for i in aris]))
#       fp.write('\n')


"""mHypo"""
# the number of clusters
# 15      15      14      15      15      15      14       15       15       15      16        15
# setting_combinations = [[15, '0.26'], [15, '0.21'], [14, '0.16'], [15, '0.11'], [15, '0.06'], [15, '0.01'], [14, '-0.04'], [15, '-0.09'], [15, '-0.14'], [15, '-0.19'], [16, '-0.24'], [15, '-0.29']]
setting_combinations = [[8, '-0.04'], [8, '-0.09'], [8, '-0.14'], [8, '-0.19'], [8, '-0.24'], [8, '-0.29']]
for setting_combi in setting_combinations:
   n_clusters = setting_combi[0]

   dataset = setting_combi[1]  #
   
   dir_ = '/home/yunfei/spatial_benchmarking/benchmarking_data/mHypothalamus'
   ad = load_mHypothalamus(root_dir=dir_, section_id=dataset)

   aris = []
   for iter in range(5):

      
      # print(ad)

      # define model
      model = GraphST.GraphST(ad, device=device)

      # train model
      ad = model.train()

      # print(ad)

      # set radius to specify the number of neighbors considered during refinement
      radius = 50

      tool = 'mclust' # mclust, leiden, and louvain

      # clustering
      if tool == 'mclust':
         clustering(ad, n_clusters, radius=radius, method=tool, refinement=True) # For DLPFC dataset, we use optional refinement step.
      elif tool in ['leiden', 'louvain']:
         clustering(ad, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)

      # filter out NA nodes
      ad = ad[~pd.isnull(ad.obs['original_clusters'])]

      # calculate metric ARI
      ARI = metrics.adjusted_rand_score(ad.obs['domain'], ad.obs['original_clusters'])
      ad.uns['ARI'] = ARI

      # print('Dataset:', dataset)
      # print('ARI:', ARI)
      aris.append(ARI)
   print('Dataset:', dataset)
   print(aris)
   print(np.mean(aris))
   with open('graphst_aris.txt', 'a+') as fp:
      fp.write('mHypothalamus' + dataset + ' ')
      fp.write(' '.join([str(i) for i in aris]))
      fp.write('\n')


# """Her2"""
# # the number of clusters
# # 6       5       4       4       4       4       7       7
# setting_combinations = [[6, 'A1'], [5, 'B1'], [4, 'C1'], [4, 'D1'], [4, 'E1'], [4, 'F1'], [7, 'G2'], [7, 'H1']]
# for setting_combi in setting_combinations:
#    n_clusters = setting_combi[0]  # 7

#    dataset = setting_combi[1]  # '151673'
   
#    dir_ = '/home/yunfei/spatial_benchmarking/benchmarking_data/Her2_tumor'
#    ad = load_her2_tumor(root_dir=dir_, section_id=dataset)

#    aris = []
#    for iter in range(5):

      
#       # print(ad)

#       # define model
#       model = GraphST.GraphST(ad, device=device)

#       # train model
#       ad = model.train()

#       # print(ad)

#       # set radius to specify the number of neighbors considered during refinement
#       radius = 50

#       tool = 'mclust' # mclust, leiden, and louvain

#       # clustering
#       if tool == 'mclust':
#          clustering(ad, n_clusters, radius=radius, method=tool, refinement=True) # For DLPFC dataset, we use optional refinement step.
#       elif tool in ['leiden', 'louvain']:
#          clustering(ad, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)

#       # filter out NA nodes
#       ad = ad[~pd.isnull(ad.obs['original_clusters'])]

#       # calculate metric ARI
#       ARI = metrics.adjusted_rand_score(ad.obs['domain'], ad.obs['original_clusters'])
#       ad.uns['ARI'] = ARI

#       # print('Dataset:', dataset)
#       # print('ARI:', ARI)
#       aris.append(ARI)
#    print('Dataset:', dataset)
#    print(aris)
#    print(np.mean(aris))
#    with open('graphst_aris.txt', 'a+') as fp:
#       fp.write('Her2tumor' + dataset + ' ')
#       fp.write(' '.join([str(i) for i in aris]))
#       fp.write('\n')


# """load mMAMP ma section"""
# # the number of clusters
# # 6       5       4       4       4       4       7       7
# setting_combinations = [[52, 'MA']]
# for setting_combi in setting_combinations:
#    n_clusters = setting_combi[0]  # 7

#    dataset = setting_combi[1]  # '151673'
   
#    dir_ = '/home/yunfei/spatial_benchmarking/benchmarking_data/mMAMP'
#    ad = load_mMAMP(root_dir=dir_, section_id=dataset)

#    aris = []
#    for iter in range(5):

      
#       # print(ad)

#       # define model
#       model = GraphST.GraphST(ad, device=device)

#       # train model
#       ad = model.train()

#       # print(ad)

#       # set radius to specify the number of neighbors considered during refinement
#       radius = 50

#       tool = 'mclust' # mclust, leiden, and louvain

#       # clustering
#       if tool == 'mclust':
#          clustering(ad, n_clusters, radius=radius, method=tool, refinement=True) # For DLPFC dataset, we use optional refinement step.
#       elif tool in ['leiden', 'louvain']:
#          clustering(ad, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)

#       # filter out NA nodes
#       ad = ad[~pd.isnull(ad.obs['original_clusters'])]

#       # calculate metric ARI
#       ARI = metrics.adjusted_rand_score(ad.obs['domain'], ad.obs['original_clusters'])
#       ad.uns['ARI'] = ARI

#       # print('Dataset:', dataset)
#       # print('ARI:', ARI)
#       aris.append(ARI)
#    print('Dataset:', dataset)
#    print(aris)
#    print(np.mean(aris))
#    with open('graphst_aris.txt', 'a+') as fp:
#       fp.write('mABC' + dataset + ' ')
#       fp.write(' '.join([str(i) for i in aris]))
#       fp.write('\n')
