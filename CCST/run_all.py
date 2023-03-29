import os
import sys
import matplotlib
matplotlib.use('Agg')
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
#import pylab as pl
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import time
from sklearn import metrics
from sklearn.metrics import adjusted_rand_score, fowlkes_mallows_score
from scipy import sparse
#from sklearn.metrics import roc_curve, auc, roc_auc_score
from st_loading_utils import load_mPFC, load_mHypothalamus, load_her2_tumor, load_mMAMP, load_DLPFC, load_BC, load_mVC
import numpy as np
import pickle
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, ChebConv, GATConv, DeepGraphInfomax, global_mean_pool, global_max_pool  # noqa
from torch_geometric.data import Data, DataLoader
from datetime import datetime
import argparse
parser = argparse.ArgumentParser()
# ================Specify data type firstly===============
parser.add_argument( '--data_type', default='nsc', help='"sc" or "nsc", \
   refers to single cell resolution datasets(e.g. MERFISH) and \
   non single cell resolution data(e.g. ST) respectively') 
# =========================== args ===============================
parser.add_argument( '--data_name', type=str, default='V1_Breast_Cancer_Block_A_Section_1', help="'MERFISH' or 'V1_Breast_Cancer_Block_A_Section_1") 
parser.add_argument( '--lambda_I', type=float, default=0.3) #0.8 on MERFISH, 0.3 on ST
parser.add_argument( '--data_path', type=str, default='generated_data/', help='data path')
parser.add_argument( '--model_path', type=str, default='model') 
parser.add_argument( '--embedding_data_path', type=str, default='Embedding_data') 
parser.add_argument( '--result_path', type=str, default='results') 
parser.add_argument( '--DGI', type=int, default=1, help='run Deep Graph Infomax(DGI) model, otherwise direct load embeddings')
parser.add_argument( '--load', type=int, default=0, help='Load pretrained DGI model')
parser.add_argument( '--num_epoch', type=int, default=5000, help='numebr of epoch in training DGI')
parser.add_argument( '--hidden', type=int, default=256, help='hidden channels in DGI') 
parser.add_argument( '--PCA', type=int, default=1, help='run PCA or not')
parser.add_argument( '--cluster', type=int, default=1, help='run cluster or not')
parser.add_argument( '--n_clusters', type=int, default=5, help='number of clusters in Kmeans, when ground truth label is not avalible.') #5 on MERFISH, 20 on Breast
parser.add_argument( '--draw_map', type=int, default=1, help='run drawing map')
parser.add_argument( '--diff_gene', type=int, default=0, help='Run differential gene expression analysis')
parser.add_argument( '--batch_size', type=int, default=512, help='training batch size')
parser.add_argument( '--gpu_id', type=str, default="0", help='default gpu id')
args = parser.parse_args()
# iters=2 # for script testing
iters = 5 # for boxplotting
# the location of R, which is necessary for mclust algorithm. Please replace the path below with local R installation path
# os.environ['R_HOME'] = '/home/yunfei/anaconda3/envs/GraphST/lib/R'
# os.environ['R_USER'] = '/home/yunfei/anaconda3/envs/GraphST/lib/python3.8/site-packages'
args.embedding_data_path = './generated_data'
# save_dir_v = './integration_visualization/'
save_dir_r = './clustering_results/'

# if not os.path.exists(save_dir_v):
#     os.makedirs(save_dir_v)

if not os.path.exists(save_dir_r):
    os.makedirs(save_dir_r)


def purity_score(y_true, y_pred):
    # compute contingency matrix (also called confusion matrix)
    contingency_matrix = metrics.cluster.contingency_matrix(y_true, y_pred)
    # return purity
    return np.sum(np.amax(contingency_matrix, axis=0)) / np.sum(contingency_matrix)


"""DLPFC"""
# the number of clusters
setting_combinations = [[7, '151507'], [7, '151508'], [7, '151509'], [7, '151510'], [5, '151669'], [5, '151670'], [5, '151671'], [5, '151672'], [7, '151673'], [7, '151674'], [7, '151675'], [7, '151676']]
# setting_combinations = [[7, '151674'], [7, '151675'], [7, '151676']]

a = time.time()
for setting_combi in setting_combinations:
   args.n_clusters = setting_combi[0]  # 7

   args.data_name = setting_combi[1]  # '151673'
   dataset = setting_combi[1]
   args.data_type = 'nsc'
   dir_ = '/data/maiziezhou_lab/Datasets/ST_datasets/DLPFC12'
   save_dir_r = './'
   ad = load_DLPFC(root_dir=dir_, section_id=args.data_name)
   aris = []
   pss = []
   fmss = []
   args.embedding_data_path = args.embedding_data_path +'/'+ args.data_name +'/'
   args.model_path = args.model_path +'/'+ args.data_name +'/'
   args.result_path = args.result_path +'/'+ args.data_name +'/'
   if not os.path.exists(args.embedding_data_path):
      os.makedirs(args.embedding_data_path) 
   if not os.path.exists(args.model_path):
      os.makedirs(args.model_path) 
   args.result_path = args.result_path+'lambdaI'+str(args.lambda_I) +'/'
   if not os.path.exists(args.result_path):
      os.makedirs(args.result_path) 
   print ('------------------------Model and Training Details--------------------------')
   print(args) 
   
   for iter_ in range(iters):

      
      if args.data_type == 'sc': # should input a single cell resolution dataset, e.g. MERFISH
         from CCST_merfish_utils import CCST_on_MERFISH
         CCST_on_MERFISH(args)
      elif args.data_type == 'nsc': # should input a non-single cell resolution dataset, e.g. V1_Breast_Cancer_Block_A_Section_1
         from CCST_ST_utils import CCST_on_ST
         preds = CCST_on_ST(args)
      else:
         print('Data type not specified')

      # calculate metric ARI
      obs_df = ad.obs.dropna()
      obs_df['ccst_preds'] = np.array(preds)[:, 1]
      obs_df.to_csv(os.path.join(save_dir_r, 'ccst_' + args.data_name + '_iter_' + str(iter_) + '_out.csv'), index=True)
      # print(preds)
      # print(obs_df['original_clusters'].to_list())
      ARI = adjusted_rand_score(np.array(preds)[:, 1], obs_df['original_clusters'].to_list())
      FMS = fowlkes_mallows_score(np.array(preds)[:, 1], obs_df['original_clusters'].to_list())
      PS = purity_score(np.array(preds)[:, 1], obs_df['original_clusters'].to_list())
      print('Dataset:', dataset)
      print('ARI:', ARI)
      aris.append(ARI)
      print('Dataset:', dataset)
      print('PS:', PS)
      pss.append(PS)
      print('Dataset:', dataset)
      print('FMS:', FMS)
      fmss.append(FMS)

   print('Dataset:', dataset)
   print(aris)
   print(np.mean(aris))
   print('Dataset:', dataset)
   print(pss)
   print(np.mean(pss))
   print('Dataset:', dataset)
   print(fmss)
   print(np.mean(fmss))
   with open('ccst_aris.txt', 'a+') as fp:
      fp.write('DLPFC' + dataset + ' ')
      fp.write(' '.join([str(i) for i in aris]))
      fp.write('\n')
   with open('ccst_pss.txt', 'a+') as fp:
      fp.write('DLPFC' + dataset + ' ')
      fp.write(' '.join([str(i) for i in pss]))
      fp.write('\n')
   with open('ccst_fmss.txt', 'a+') as fp:
      fp.write('DLPFC' + dataset + ' ')
      fp.write(' '.join([str(i) for i in fmss]))
      fp.write('\n')
b = time.time()
print("running time = ", b-a, " s")
# exit(-1)

"""BC"""
# the number of clusters
setting_combinations = [[20, 'section1']]
a = time.time()
for setting_combi in setting_combinations:
   args.n_clusters = setting_combi[0]  # 7

   args.data_name = setting_combi[1]  # '151673'
   dataset = setting_combi[1]
   args.data_type = 'nsc'
   dir_ = '/data/maiziezhou_lab/Datasets/ST_datasets/BC'
   ad = load_BC(root_dir=dir_, section_id=args.data_name)
   aris = []
   chis = []
   sss = []
   args.embedding_data_path = args.embedding_data_path +'/'+ args.data_name +'/'
   args.model_path = args.model_path +'/'+ args.data_name +'/'
   args.result_path = args.result_path +'/'+ args.data_name +'/'
   if not os.path.exists(args.embedding_data_path):
      os.makedirs(args.embedding_data_path) 
   if not os.path.exists(args.model_path):
      os.makedirs(args.model_path) 
   args.result_path = args.result_path+'lambdaI'+str(args.lambda_I) +'/'
   if not os.path.exists(args.result_path):
      os.makedirs(args.result_path) 
   print ('------------------------Model and Training Details--------------------------')
   print(args) 
   
   for iter_ in range(iters):

      
      if args.data_type == 'sc': # should input a single cell resolution dataset, e.g. MERFISH
         from CCST_merfish_utils import CCST_on_MERFISH
         CCST_on_MERFISH(args)
      elif args.data_type == 'nsc': # should input a non-single cell resolution dataset, e.g. V1_Breast_Cancer_Block_A_Section_1
         from CCST_ST_utils import CCST_on_ST
         preds = CCST_on_ST(args)
      else:
         print('Data type not specified')

      # calculate metric ARI
      obs_df = ad.obs.dropna()
      obs_df['ccst_preds'] = np.array(preds)[:, 1]
      obs_df.to_csv(os.path.join(save_dir_r, 'ccst_' + args.data_name + '_iter_' + str(iter_) + '_out.csv'), index=True)
      ARI = adjusted_rand_score(np.array(preds)[:, 1], obs_df['original_clusters'].to_list())
      FMS = fowlkes_mallows_score(np.array(preds)[:, 1], obs_df['original_clusters'].to_list())
      PS = purity_score(np.array(preds)[:, 1], obs_df['original_clusters'].to_list())
      print('Dataset:', dataset)
      print('ARI:', ARI)
      aris.append(ARI)
      print('Dataset:', dataset)
      print('PS:', PS)
      pss.append(PS)
      print('Dataset:', dataset)
      print('FMS:', FMS)
      fmss.append(FMS)

   print('Dataset:', dataset)
   print(aris)
   print(np.mean(aris))
   print('Dataset:', dataset)
   print(pss)
   print(np.mean(pss))
   print('Dataset:', dataset)
   print(fmss)
   print(np.mean(fmss))
   with open('ccst_aris.txt', 'a+') as fp:
      fp.write('HRBC1' + dataset + ' ')
      fp.write(' '.join([str(i) for i in aris]))
      fp.write('\n')
   with open('ccst_pss.txt', 'a+') as fp:
      fp.write('HRBC1' + dataset + ' ')
      fp.write(' '.join([str(i) for i in pss]))
      fp.write('\n')
   with open('ccst_fmss.txt', 'a+') as fp:
      fp.write('HRBC1' + dataset + ' ')
      fp.write(' '.join([str(i) for i in fmss]))
      fp.write('\n')
b = time.time()
print("running time = ", b-a, " s")

"""mVC"""
# the number of clusters
a = time.time()
setting_combinations = [[7, 'STARmap_20180505_BY3_1k.h5ad']]
for setting_combi in setting_combinations:
   args.n_clusters = setting_combi[0]  # 7

   args.data_name = setting_combi[1]  # '151673'
   dataset = setting_combi[1]
   args.data_type = 'nsc'
   dir_ = '/data/maiziezhou_lab/Datasets/ST_datasets/STARmap_mouse_visual_cortex'
   ad = load_mVC(root_dir=dir_, section_id=args.data_name)
   aris = []
   sss = []
   chis = []
   args.embedding_data_path = args.embedding_data_path +'/'+ args.data_name +'/'
   args.model_path = args.model_path +'/'+ args.data_name +'/'
   args.result_path = args.result_path +'/'+ args.data_name +'/'
   if not os.path.exists(args.embedding_data_path):
      os.makedirs(args.embedding_data_path) 
   if not os.path.exists(args.model_path):
      os.makedirs(args.model_path) 
   args.result_path = args.result_path+'lambdaI'+str(args.lambda_I) +'/'
   if not os.path.exists(args.result_path):
      os.makedirs(args.result_path) 
   print ('------------------------Model and Training Details--------------------------')
   print(args) 
   
   for iter_ in range(iters):

      
      if args.data_type == 'sc': # should input a single cell resolution dataset, e.g. MERFISH
         from CCST_merfish_utils import CCST_on_MERFISH
         CCST_on_MERFISH(args)
      elif args.data_type == 'nsc': # should input a non-single cell resolution dataset, e.g. V1_Breast_Cancer_Block_A_Section_1
         from CCST_ST_utils import CCST_on_ST
         preds = CCST_on_ST(args)
      else:
         print('Data type not specified')

      # calculate metric ARI
      obs_df = ad.obs.dropna()
      obs_df['ccst_preds'] = np.array(preds)[:, 1]
      obs_df.to_csv(os.path.join(save_dir_r, 'ccst_' + args.data_name + '_iter_' + str(iter_) + '_out.csv'), index=True)
      ARI = adjusted_rand_score(np.array(preds)[:, 1], obs_df['original_clusters'].to_list())
      FMS = fowlkes_mallows_score(np.array(preds)[:, 1], obs_df['original_clusters'].to_list())
      PS = purity_score(np.array(preds)[:, 1], obs_df['original_clusters'].to_list())
      print('Dataset:', dataset)
      print('ARI:', ARI)
      aris.append(ARI)
      print('Dataset:', dataset)
      print('PS:', PS)
      pss.append(PS)
      print('Dataset:', dataset)
      print('FMS:', FMS)
      fmss.append(FMS)

   print('Dataset:', dataset)
   print(aris)
   print(np.mean(aris))
   print('Dataset:', dataset)
   print(pss)
   print(np.mean(pss))
   print('Dataset:', dataset)
   print(fmss)
   print(np.mean(fmss))
   with open('ccst_aris.txt', 'a+') as fp:
      fp.write('mVC ')
      fp.write(' '.join([str(i) for i in aris]))
      fp.write('\n')
   with open('ccst_pss.txt', 'a+') as fp:
      fp.write('mVC ')
      fp.write(' '.join([str(i) for i in pss]))
      fp.write('\n')
   with open('ccst_fmss.txt', 'a+') as fp:
      fp.write('mVC ')
      fp.write(' '.join([str(i) for i in fmss]))
      fp.write('\n')
b = time.time()
print("running time = ", b-a, " s")


"""mPFC"""
# the number of clusters
setting_combinations = [[4, '20180417_BZ5_control'], [4, '20180419_BZ9_control'], [4, '20180424_BZ14_control']]
for setting_combi in setting_combinations:
   args.n_clusters = setting_combi[0]  # 7

   args.data_name = setting_combi[1]  # '151673'
   dataset = setting_combi[1]
   args.data_type = 'nsc'
   dir_ = '/data/maiziezhou_lab/Datasets/ST_datasets/STARmap_mouse_PFC'
   ad = load_mPFC(root_dir=dir_, section_id=args.data_name)
   aris = []
   args.embedding_data_path = args.embedding_data_path +'/'+ args.data_name +'/'
   args.model_path = args.model_path +'/'+ args.data_name +'/'
   args.result_path = args.result_path +'/'+ args.data_name +'/'
   if not os.path.exists(args.embedding_data_path):
      os.makedirs(args.embedding_data_path) 
   if not os.path.exists(args.model_path):
      os.makedirs(args.model_path) 
   args.result_path = args.result_path+'lambdaI'+str(args.lambda_I) +'/'
   if not os.path.exists(args.result_path):
      os.makedirs(args.result_path) 
   print ('------------------------Model and Training Details--------------------------')
   print(args) 
   
   for iter_ in range(iters):

      
      if args.data_type == 'sc': # should input a single cell resolution dataset, e.g. MERFISH
         from CCST_merfish_utils import CCST_on_MERFISH
         CCST_on_MERFISH(args)
      elif args.data_type == 'nsc': # should input a non-single cell resolution dataset, e.g. V1_Breast_Cancer_Block_A_Section_1
         from CCST_ST_utils import CCST_on_ST
         preds = CCST_on_ST(args)
      else:
         print('Data type not specified')

      # calculate metric ARI
      obs_df = ad.obs.dropna()
      obs_df['ccst_preds'] = np.array(preds)[:, 1]
      obs_df.to_csv(os.path.join(save_dir_r, 'ccst_' + args.data_name + '_iter_' + str(iter_) + '_out.csv'), index=True)
      ARI = adjusted_rand_score(np.array(preds)[:, 1], obs_df['original_clusters'].to_list())

      print('Dataset:', dataset)
      print('ARI:', ARI)
      aris.append(ARI)
   print('Dataset:', dataset)
   print(aris)
   print(np.mean(aris))
   with open('ccst_aris.txt', 'a+') as fp:
      fp.write('mPFC' + dataset + ' ')
      fp.write(' '.join([str(i) for i in aris]))
      fp.write('\n')


"""mHypo"""
# the number of clusters
# 15      15      14      15      15      15      14       15       15       15      16        15
# setting_combinations = [[15, '0.26'], [15, '0.21'], [14, '0.16'], [15, '0.11'], [15, '0.06'], [15, '0.01'], [14, '-0.04'], [15, '-0.09'], [15, '-0.14'], [15, '-0.19'], [16, '-0.24'], [15, '-0.29']]
setting_combinations = [[8, '-0.14'], [8, '-0.19'], [8, '-0.24'], [8, '-0.29']]
for setting_combi in setting_combinations:
   args.n_clusters = setting_combi[0]  # 7

   args.data_name = setting_combi[1]  # '151673'
   dataset = setting_combi[1]
   args.data_type = 'nsc'
   dir_ = '/data/maiziezhou_lab/Datasets/ST_datasets/mHypothalamus'
   ad = load_mHypothalamus(root_dir=dir_, section_id=args.data_name)
   aris = []
   args.embedding_data_path = args.embedding_data_path +'/'+ args.data_name +'/'
   args.model_path = args.model_path +'/'+ args.data_name +'/'
   args.result_path = args.result_path +'/'+ args.data_name +'/'
   if not os.path.exists(args.embedding_data_path):
      os.makedirs(args.embedding_data_path) 
   if not os.path.exists(args.model_path):
      os.makedirs(args.model_path) 
   args.result_path = args.result_path+'lambdaI'+str(args.lambda_I) +'/'
   if not os.path.exists(args.result_path):
      os.makedirs(args.result_path) 
   print ('------------------------Model and Training Details--------------------------')
   print(args) 
   
   for iter_ in range(iters):

      
      if args.data_type == 'sc': # should input a single cell resolution dataset, e.g. MERFISH
         from CCST_merfish_utils import CCST_on_MERFISH
         CCST_on_MERFISH(args)
      elif args.data_type == 'nsc': # should input a non-single cell resolution dataset, e.g. V1_Breast_Cancer_Block_A_Section_1
         from CCST_ST_utils import CCST_on_ST
         preds = CCST_on_ST(args)
      else:
         print('Data type not specified')

      # calculate metric ARI
      obs_df = ad.obs.dropna()
      obs_df['ccst_preds'] = np.array(preds)[:, 1]
      obs_df.to_csv(os.path.join(save_dir_r, 'ccst_' + args.data_name + '_iter_' + str(iter_) + '_out.csv'), index=True)
      # print(obs_df)
      # print(np.array(preds).shape)
      ARI = adjusted_rand_score(np.array(preds)[:, 1], ad.obs['original_clusters'].to_list())
      # exit(-1)
      print('Dataset:', dataset)
      print('ARI:', ARI)
      aris.append(ARI)
   print('Dataset:', dataset)
   print(aris)
   print(np.mean(aris))
   with open('ccst_aris.txt', 'a+') as fp:
      fp.write('mHypothalamus' + dataset + ' ')
      fp.write(' '.join([str(i) for i in aris]))
      fp.write('\n')


"""Her2"""
# the number of clusters
# 6       5       4       4       4       4       7       7
setting_combinations = [[6, 'A1'], [5, 'B1'], [4, 'C1'], [4, 'D1'], [4, 'E1'], [4, 'F1'], [7, 'G2'], [7, 'H1']]
for setting_combi in setting_combinations:
   args.n_clusters = setting_combi[0]  # 7

   args.data_name = setting_combi[1]  # '151673'
   dataset = setting_combi[1]
   args.data_type = 'nsc'
   dir_ = '/data/maiziezhou_lab/Datasets/ST_datasets/Her2_tumor'
   ad = load_her2_tumor(root_dir=dir_, section_id=args.data_name)
   aris = []
   args.embedding_data_path = args.embedding_data_path +'/'+ args.data_name +'/'
   args.model_path = args.model_path +'/'+ args.data_name +'/'
   args.result_path = args.result_path +'/'+ args.data_name +'/'
   if not os.path.exists(args.embedding_data_path):
      os.makedirs(args.embedding_data_path) 
   if not os.path.exists(args.model_path):
      os.makedirs(args.model_path) 
   args.result_path = args.result_path+'lambdaI'+str(args.lambda_I) +'/'
   if not os.path.exists(args.result_path):
      os.makedirs(args.result_path) 
   print ('------------------------Model and Training Details--------------------------')
   print(args) 
   
   for iter_ in range(iters):

      
      if args.data_type == 'sc': # should input a single cell resolution dataset, e.g. MERFISH
         from CCST_merfish_utils import CCST_on_MERFISH
         CCST_on_MERFISH(args)
      elif args.data_type == 'nsc': # should input a non-single cell resolution dataset, e.g. V1_Breast_Cancer_Block_A_Section_1
         from CCST_ST_utils import CCST_on_ST
         preds = CCST_on_ST(args)
      else:
         print('Data type not specified')

      # calculate metric ARI
      obs_df = ad.obs.dropna()
      obs_df['ccst_preds'] = np.array(preds)[:, 1]
      obs_df.to_csv(os.path.join(save_dir_r, 'ccst_' + args.data_name + '_iter_' + str(iter_) + '_out.csv'), index=True)
      ARI = adjusted_rand_score(np.array(preds)[:, 1], obs_df['original_clusters'].to_list())

      print('Dataset:', dataset)
      print('ARI:', ARI)
      aris.append(ARI)
   print('Dataset:', dataset)
   print(aris)
   print(np.mean(aris))
   with open('ccst_aris.txt', 'a+') as fp:
      fp.write('Her2tumor' + dataset + ' ')
      fp.write(' '.join([str(i) for i in aris]))
      fp.write('\n')


"""load mMAMP ma section"""
# the number of clusters
# 6       5       4       4       4       4       7       7
setting_combinations = [[52, 'MA']]
for setting_combi in setting_combinations:
   args.n_clusters = setting_combi[0]  # 7

   args.data_name = setting_combi[1]  # '151673'
   dataset = setting_combi[1]
   args.data_type = 'nsc'
   dir_ = '/data/maiziezhou_lab/Datasets/ST_datasets/mMAMP'
   ad = load_mMAMP(root_dir=dir_, section_id=args.data_name)
   aris = []
   args.embedding_data_path = args.embedding_data_path +'/'+ args.data_name +'/'
   args.model_path = args.model_path +'/'+ args.data_name +'/'
   args.result_path = args.result_path +'/'+ args.data_name +'/'
   if not os.path.exists(args.embedding_data_path):
      os.makedirs(args.embedding_data_path) 
   if not os.path.exists(args.model_path):
      os.makedirs(args.model_path) 
   args.result_path = args.result_path+'lambdaI'+str(args.lambda_I) +'/'
   if not os.path.exists(args.result_path):
      os.makedirs(args.result_path) 
   print ('------------------------Model and Training Details--------------------------')
   print(args) 
   
   for iter_ in range(iters):

      
      if args.data_type == 'sc': # should input a single cell resolution dataset, e.g. MERFISH
         from CCST_merfish_utils import CCST_on_MERFISH
         CCST_on_MERFISH(args)
      elif args.data_type == 'nsc': # should input a non-single cell resolution dataset, e.g. V1_Breast_Cancer_Block_A_Section_1
         from CCST_ST_utils import CCST_on_ST
         preds = CCST_on_ST(args)
      else:
         print('Data type not specified')

      # calculate metric ARI
      obs_df = ad.obs.dropna()
      
      ARI = adjusted_rand_score(np.array(preds)[:, 1], obs_df['original_clusters'].to_list())
      obs_df['ccst_preds'] = np.array(preds)[:, 1]
      obs_df.to_csv(os.path.join(save_dir_r, 'ccst_' + args.data_name + '_iter_' + str(iter_) + '_out.csv'), index=True)
      print('Dataset:', dataset)
      print('ARI:', ARI)
      aris.append(ARI)
   print('Dataset:', dataset)
   print(aris)
   print(np.mean(aris))
   with open('ccst_aris.txt', 'a+') as fp:
      fp.write('mABC ')
      fp.write(' '.join([str(i) for i in aris]))
      fp.write('\n')
   