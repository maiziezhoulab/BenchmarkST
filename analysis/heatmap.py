import os
import numpy as np
import pickle

# plot heatmap #
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import pandas as pd
import glob
from collections import defaultdict

from collections import OrderedDict
import math


def plot_heatmap(heatmap_df,save_dir,data_key,save_df=True):

    if save_df:
        with open(save_dir+'/'+data_key+'_across_datasets_heatmap_data.pkl','wb') as hmpk:
            pickle.dump(heatmap_df,hmpk)
    
    fig, ax = plt.subplots(figsize=(len(heatmap_df.columns)+1,len(heatmap_df.index)))
    #sns.heatmap(heatmap_df, annot=True, ax=ax,cmap='viridis',vmin=0, vmax=1,fmt='.1f',mask=heatmap_df.isnull())
    sns.heatmap(heatmap_df, annot=True, ax=ax,cmap='viridis',vmin=0, vmax=1,mask=heatmap_df.isnull())
    ax.xaxis.tick_top() # x axis on top
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    #plt.title(data_key)
    mpl.rcParams['pdf.fonttype'] = 42
    plt.savefig(save_dir+'/'+data_key+'_across_datasets_heatmap.pdf',bbox_inches='tight')
    plt.close()
    return 0


def read_txt_to_df(results_path='/home/ethan/BenchmarkST/ari_results/*.txt'):
    tools = []
    results = defaultdict(list)
    datasets = ['DLPFC151507', 'DLPFC151508', 'DLPFC151509', 'DLPFC151510', 
                'DLPFC151669', 'mVC', 'mPFC20180417_BZ5_control', 'mPFC20180419_BZ9_control', 
                'mPFC20180424_BZ14_control', 'mHypothalamus-0.04', 'mHypothalamus-0.09', 
                'mHypothalamus-0.14', 'mHypothalamus-0.19', 'mHypothalamus-0.24', 'mHypothalamus-0.29', 
                'Her2tumorA1', 'Her2tumorB1', 'Her2tumorC1', 'Her2tumorD1', 'Her2tumorE1', 
                'Her2tumorF1', 'Her2tumorG2', 'Her2tumorH1', 'DLPFC151670', 'DLPFC151671', 
                'DLPFC151672','DLPFC151673', 'DLPFC151674', 'DLPFC151675', 'DLPFC151676', 'HBRC1', 'mABC']
    for path in glob.glob(results_path):
        tool = os.path.basename(path)[:-4]
        tools.append(tool)
        # data = defaultdict(list)
        with open(path, 'r') as f:
            for line in f.readlines():
                res = line.strip().split()
                dataset = res[0]
                runs = []
                for i in range(1, len(res)):
                    runs.append(0.0 if float(res[i]) <= 0.0 else float(res[i]))
                if sum(runs) == 0:
                    results[dataset].append(float('nan'))
                else:
                    results[dataset].append(sum(runs) / len(runs))
    # sort_res = [[] ]
    # res_dict = dict(zip(tools, [results]))
    heatmap_df = pd.DataFrame.from_dict(results, orient='columns')
    heatmap_df.index = tools
    return heatmap_df


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input_dir', type=str,default='/home/ethan/BenchmarkST/ari_results')
    parser.add_argument('--data_type', type=str,default="ari")
    parser.add_argument('--save_dir', type=str,default="/home/ethan/BenchmarkST/data/")
    args = parser.parse_args()

    input_dir = args.input_dir
    data_type = args.data_type
    save_dir = args.save_dir


    heatmap_df = read_txt_to_df()
    plot_heatmap(heatmap_df,save_dir,data_type,save_df=True)