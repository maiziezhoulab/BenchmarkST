import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
import numpy as np
from GraphST import GraphST
from GraphST.utils import clustering
from st_loading_utils import load_DLPFC, load_BC, load_mVC

# Run device, by default, the package is implemented on 'cpu'. We recommend using GPU.
device = torch.device('cuda:3' if torch.cuda.is_available() else 'cpu')

# the location of R, which is necessary for mclust algorithm. Please replace the path below with local R installation path
# os.environ['R_HOME'] = '/home/yunfei/anaconda3/envs/GraphST/lib/R'
# os.environ['R_USER'] = '/home/yunfei/anaconda3/envs/GraphST/lib/python3.8/site-packages'

# """DLPFC"""
# # the number of clusters
# # setting_combinations = [[7, '151507'], [7, '151508'], [7, '151509'], [7, '151510'], [5, '151669'], [5, '151670'], [5, '151671'], [5, '151672'], [7, '151673'], [7, '151674'], [7, '151675'], [7, '151676']]
# setting_combinations = [[7, '151674'], [7, '151675'], [7, '151676']]
# for setting_combi in setting_combinations:
#    n_clusters = setting_combi[0]  # 7

#    dataset = setting_combi[1]  # '151673'
   
#    dir_ = '/home/yunfei/spatial_benchmarking/benchmarking_data/DLPFC12'
#    ad = load_DLPFC(root_dir=dir_, section_id=dataset)

#    aris = []
#    for iter in range(20):

      
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

#       print('Dataset:', dataset)
#       print('ARI:', ARI)
#       aris.append(ARI)
#    print('Dataset:', dataset)
#    print(aris)
#    print(np.mean(aris))
#    with open('graphst_aris.txt', 'a+') as fp:
#       fp.write('DLPFC' + dataset + ' ')
#       fp.write(' '.join([str(i) for i in aris]))
#       fp.write('\n')


"""BC"""
# the number of clusters
setting_combinations = [[20, 'section1']]
for setting_combi in setting_combinations:
   n_clusters = setting_combi[0]  # 7

   dataset = setting_combi[1]  #
   
   dir_ = '/home/yunfei/spatial_benchmarking/benchmarking_data/BC'
   ad = load_BC(root_dir=dir_, section_id=dataset)

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
      fp.write('HBRC1 ')
      fp.write(' '.join([str(i) for i in aris]))
      fp.write('\n')


"""mVC"""
# the number of clusters
setting_combinations = [[7, 'STARmap_20180505_BY3_1k.h5ad']]
for setting_combi in setting_combinations:
   n_clusters = setting_combi[0]

   dataset = setting_combi[1]
   
   dir_ = '/home/yunfei/spatial_benchmarking/benchmarking_data/STARmap_mouse_visual_cortex'
   ad = load_mVC(root_dir=dir_, section_id=dataset)

   aris = []
   for iter in range(5):

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
      fp.write('mVC ')
      fp.write(' '.join([str(i) for i in aris]))
      fp.write('\n')
