import os

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def read_data(path_):
    dict_ = {}
    with open(path_, 'r') as f:
        lines = f.readlines()

    for line in lines:
        key = line.split(' ')[0]
        if key in ['mHypothalamus0.01', 'mHypothalamus0.06', 'mHypothalamus0.11', 'mHypothalamus0.16', 'mHypothalamus0.21', 'mHypothalamus0.26']:
            continue
        values = line.rstrip('\n').split(' ')[1:]
        # print(values)
        if len(values) != 20:
            values = [values[0]] * 20
        if values == ['']:
            values = [0.5] * 20
        values = [float(i) for i in values]
        dict_[key] = values
        # split_ = line.split(' ')
        # ari.append(float(split_[-2]))
    return dict_


if __name__ == '__main__':
    dir_ = '/home/yunfei/spatial_benchmarking/BenchmarkST/ari_results'
    stagate_ = 'stagate_aris.txt'
    ccst_ = 'ccst_aris.txt'
    sedr_ = 'sedr_aris.txt'
    spagcn_ = 'spagcn_aris.txt'
    graphst_ = 'graphst_aris.txt'
    deepst_ = 'deepst_aris.txt'

    stagate_data = read_data(os.path.join(dir_, stagate_))
    name_list1 = ['STAGATE', 'STAGATE', 'STAGATE', 'STAGATE', 'STAGATE',
                  'STAGATE', 'STAGATE', 'STAGATE', 'STAGATE', 'STAGATE',
                  'STAGATE', 'STAGATE', 'STAGATE', 'STAGATE', 'STAGATE',
                  'STAGATE', 'STAGATE', 'STAGATE', 'STAGATE', 'STAGATE']
    print(stagate_data.keys())
    # exit(-1)
    ccst_data = read_data(os.path.join(dir_, ccst_))
    name_list2 = ['CCST', 'CCST', 'CCST', 'CCST', 'CCST',
                  'CCST', 'CCST', 'CCST', 'CCST', 'CCST',
                  'CCST', 'CCST', 'CCST', 'CCST', 'CCST',
                  'CCST', 'CCST', 'CCST', 'CCST', 'CCST'
                  ]
    print(ccst_data.keys())

    sedr_data = read_data(os.path.join(dir_, sedr_))
    name_list3 = ['SEDR', 'SEDR', 'SEDR', 'SEDR', 'SEDR',
                  'SEDR', 'SEDR', 'SEDR', 'SEDR', 'SEDR',
                  'SEDR', 'SEDR', 'SEDR', 'SEDR', 'SEDR',
                  'SEDR', 'SEDR', 'SEDR', 'SEDR', 'SEDR']
    print(sedr_data.keys())

    spagcn_data = read_data(os.path.join(dir_, spagcn_))
    name_list4 = ['SpaGCN', 'SpaGCN', 'SpaGCN', 'SpaGCN', 'SpaGCN',
                  'SpaGCN', 'SpaGCN', 'SpaGCN', 'SpaGCN', 'SpaGCN',
                  'SpaGCN', 'SpaGCN', 'SpaGCN', 'SpaGCN', 'SpaGCN',
                  'SpaGCN', 'SpaGCN', 'SpaGCN', 'SpaGCN', 'SpaGCN']
    print(spagcn_data.keys())

    graphst_data = read_data(os.path.join(dir_, graphst_))
    name_list5 = ['GraphST', 'GraphST', 'GraphST', 'GraphST', 'GraphST',
                  'GraphST', 'GraphST', 'GraphST', 'GraphST', 'GraphST',
                  'GraphST', 'GraphST', 'GraphST', 'GraphST', 'GraphST',
                  'GraphST', 'GraphST', 'GraphST', 'GraphST', 'GraphST']
    print(graphst_data.keys())

    deepst_data = read_data(os.path.join(dir_, deepst_))
    name_list6 = ['DeepST', 'DeepST', 'DeepST', 'DeepST', 'DeepST',
                  'DeepST', 'DeepST', 'DeepST', 'DeepST', 'DeepST',
                  'DeepST', 'DeepST', 'DeepST', 'DeepST', 'DeepST',
                  'DeepST', 'DeepST', 'DeepST', 'DeepST', 'DeepST']
    print(deepst_data.keys())

    # exit(-2)
    r_dict = {}
    for k in spagcn_data.keys():
        r_dict[k] = []
    for k in spagcn_data.keys():
        # print(k)
        # print("**********stagate")
        # print(np.mean(stagate_data[k]))
        # print(np.std(stagate_data[k]))
        r_dict[k].append(np.mean(stagate_data[k]))

        # print("**********ccst")
        # print(np.mean(ccst_data[k]))
        # print(np.std(ccst_data[k]))
        r_dict[k].append(np.mean(ccst_data[k]))

        # print("**********sedr")
        # print(np.mean(sedr_data[k]))
        # print(np.std(sedr_data[k]))
        r_dict[k].append(np.mean(sedr_data[k]))

        # print("**********spagcn")
        # print(np.mean(spagcn_data[k]))
        # print(np.std(spagcn_data[k]))
        r_dict[k].append(np.mean(spagcn_data[k]))

        # print("**********deepst")
        # print(np.mean(deepst_data[k]))
        # print(np.std(deepst_data[k]))
        r_dict[k].append(np.mean(deepst_data[k]))

        # print("**********graphst")
        # print(np.mean(graphst_data[k]))
        # print(np.std(graphst_data[k]))
        r_dict[k].append(np.mean(graphst_data[k]))

        # print("**********")
        d = {'ARI': deepst_data[k] + stagate_data[k] + ccst_data[k] + sedr_data[k] + spagcn_data[k] + graphst_data[k],
             k: name_list6 + name_list1 + name_list2 + name_list3 + name_list4 + name_list5}
        df = pd.DataFrame(data=d)

        fig, ax = plt.subplots()
        # df = sns.load_dataset('iris')
        print(df)
        sns.boxplot(x=df[k], y=df["ARI"])
        ax.set_ylim([0, 0.85])
        # plt.savefig(os.path.join(dir_, k + '_benchmarking_boxplot0224.pdf'))
    print(r_dict)
    # print(ccst_data)
