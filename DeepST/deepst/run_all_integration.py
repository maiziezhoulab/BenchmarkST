import os 
from DeepST import run
import matplotlib.pyplot as plt
from pathlib import Path
from st_loading_utils import load_BC, load_DLPFC
import scanpy as sc
import numpy as np
from sklearn import metrics
from sklearn.metrics import adjusted_rand_score, fowlkes_mallows_score
import json

save_dir_v = './integration_visualization/'
save_dir_r = './integration_result/'

if not os.path.exists(save_dir_v):
    os.makedirs(save_dir_v)

if not os.path.exists(save_dir_r):
    os.makedirs(save_dir_r)


def purity_score(y_true, y_pred):
    # compute contingency matrix (also called confusion matrix)
    contingency_matrix = metrics.cluster.contingency_matrix(y_true, y_pred)
    # return purity
    return np.sum(np.amax(contingency_matrix, axis=0)) / np.sum(contingency_matrix)


iter_=20
dict_out = {}
"""pair-wise result"""
sets_ = [['151507', '151508', '151509', '151510'], ['151669', '151670', '151671', '151672'], ['151673', '151674', '151675', '151676']]
# section_list = ['151673', '151674', '151675', '151676']
for section_list in sets_:
    for index_ in range(1, len(section_list)):
        sec1 = section_list[index_-1]
        sec2 = section_list[index_]
        temp_section_list = [sec1, sec2]

        if sec1 + '_' + sec2 not in dict_out.keys():
            dict_out[sec1 + '_' + sec2] = {}
        

        for e in range(iter_):
            data_path = "/home/yunfei/spatial_benchmarking/benchmarking_data/DLPFC12" 
            data_name_list = temp_section_list
            save_path = "../Results" 
            n_domains = 7 
            deepen = run(save_path = save_path, 
                pca_n_comps = 200,
                pre_epochs = 800, #### According to your own hardware, choose the number of training
                epochs = 1000, #### According to your own hardware, choose the number of training
                platform = "Visium",
                )
            adata, graph_dict, domains = deepen._get_multiple_adata(data_path, data_name_list)
            adata = deepen._fit(adata, graph_dict, domains, pretrain = True)
            adata = deepen._get_cluster_data(adata, n_domains = n_domains, priori=True)


            np.save(os.path.join(save_dir_r, 'deepst_' + sec1 + '_' + sec2 + '_iter_' + str(e) + '_deepst_embedding.npy'), adata.obsm['DeepST_embed'])
                        
            # get/save ARI, purity, FMS
            for section_id in temp_section_list:
                if section_id not in dict_out[sec1 + '_' + sec2].keys():
                    dict_out[sec1 + '_' + sec2][section_id] = {}
                temp_adata = adata.obs[adata.obs['batch_name'] == section_id]
                temp_obs = temp_adata.dropna()
                save_df = temp_obs[['original_clusters', 'DeepST_refine_domain']]
                save_df.to_csv(os.path.join(save_dir_r, 'deepst_' + section_id + '_in_' + sec1 + '_' + sec2 + '_iter_' + str(e) + '_deepst_out.csv'), index=True)

                temp_ARI = adjusted_rand_score(temp_obs['DeepST_refine_domain'], temp_obs['original_clusters'])
                FMS = fowlkes_mallows_score(temp_obs['DeepST_refine_domain'], temp_obs['original_clusters'])
                PS = purity_score(temp_obs['DeepST_refine_domain'], temp_obs['original_clusters'])
                print('ARI of section ID %s: %.3f' %(section_id, temp_ARI))
                print('FMS of section ID %s: %.3f' %(section_id, FMS))
                print('PS of section ID %s: %.3f' %(section_id, PS))
                if 'ARI' not in dict_out[sec1 + '_' + sec2][section_id].keys():
                    dict_out[sec1 + '_' + sec2][section_id]['ARI'] = [temp_ARI]
                else:
                    dict_out[sec1 + '_' + sec2][section_id]['ARI'].append(temp_ARI)
                
                if 'FMS' not in dict_out[sec1 + '_' + sec2][section_id].keys():
                    dict_out[sec1 + '_' + sec2][section_id]['FMS'] = [FMS]
                else:
                    dict_out[sec1 + '_' + sec2][section_id]['FMS'].append(FMS)
                
                if 'PS' not in dict_out[sec1 + '_' + sec2][section_id].keys():
                    dict_out[sec1 + '_' + sec2][section_id]['PS'] = [PS]
                else:
                    dict_out[sec1 + '_' + sec2][section_id]['PS'].append(PS)


with open('deepst_result_dlpfc_pairwise.json', 'w', encoding='utf-8') as f:
    json.dump(dict_out, f, ensure_ascii=False, indent=4)


"""all slices result"""
sets_ = [['151507', '151508', '151509', '151510'], ['151669', '151670', '151671', '151672'], ['151673', '151674', '151675', '151676']]
# section_list = ['151673', '151674', '151675', '151676']
for section_list in sets_:
    temp_section_list = section_list
    sec1, sec2, sec3, sec4 = temp_section_list

    if sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 not in dict_out.keys():
        dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4] = {}

    for e in range(iter_):
        data_path = "/home/yunfei/spatial_benchmarking/benchmarking_data/DLPFC12" 
        data_name_list = temp_section_list
        save_path = "../Results" 
        n_domains = 7 
        deepen = run(save_path = save_path, 
            pca_n_comps = 200,
            pre_epochs = 800, #### According to your own hardware, choose the number of training
            epochs = 1000, #### According to your own hardware, choose the number of training
            platform = "Visium",
            )
        adata, graph_dict, domains = deepen._get_multiple_adata(data_path, data_name_list)
        adata = deepen._fit(adata, graph_dict, domains, pretrain = True)
        adata = deepen._get_cluster_data(adata, n_domains = n_domains, priori=True)

        np.save(os.path.join(save_dir_r, 'deepst_' + sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 + '_iter_' + str(e) + '_deepst_embedding.npy'), adata.obsm['DeepST_embed'])
                    
        # get/save ARI, purity, FMS
        for section_id in temp_section_list:
            if section_id not in dict_out[sec1 + '_' + sec2].keys():
                dict_out[sec1 + '_' + sec2][section_id] = {}
            temp_adata = adata.obs[adata.obs['batch_name'] == section_id]
            temp_obs = temp_adata.dropna()
            save_df = temp_obs[['original_clusters', 'DeepST_refine_domain']]
            save_df.to_csv(os.path.join(save_dir_r, 'deepst_' + sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 + '_iter_' + str(e) + '_deepst_out.csv'), index=True)

            temp_ARI = adjusted_rand_score(temp_obs['DeepST_refine_domain'], temp_obs['original_clusters'])
            FMS = fowlkes_mallows_score(temp_obs['DeepST_refine_domain'], temp_obs['original_clusters'])
            PS = purity_score(temp_obs['DeepST_refine_domain'], temp_obs['original_clusters'])
            print('ARI of section ID %s: %.3f' %(section_id, temp_ARI))
            print('FMS of section ID %s: %.3f' %(section_id, FMS))
            print('PS of section ID %s: %.3f' %(section_id, PS))
            if 'ARI' not in dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4][section_id].keys():
                dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4][section_id]['ARI'] = [temp_ARI]
            else:
                dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4][section_id]['ARI'].append(temp_ARI)
            
            if 'FMS' not in dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4][section_id].keys():
                dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4][section_id]['FMS'] = [FMS]
            else:
                dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4][section_id]['FMS'].append(FMS)
            
            if 'PS' not in dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4][section_id].keys():
                dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4][section_id]['PS'] = [PS]
            else:
                dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4][section_id]['PS'].append(PS)

with open('deepst_result_dlpfc_all.json', 'w', encoding='utf-8') as f:
    json.dump(dict_out, f, ensure_ascii=False, indent=4)


dict_out = {}
"""pair-wise result BC"""
section_list = ['section1', 'section2']
for index_ in range(1, len(section_list)):
    sec1 = section_list[index_-1]
    sec2 = section_list[index_]
    temp_section_list = [sec1, sec2]

    if sec1 + '_' + sec2 not in dict_out.keys():
        dict_out[sec1 + '_' + sec2] = {}

    for section_id in temp_section_list:
        temp_section_list = [sec1, sec2]

        if sec1 + '_' + sec2 not in dict_out.keys():
            dict_out[sec1 + '_' + sec2] = {}
        

        for e in range(iter_):
            data_path = "/home/yunfei/spatial_benchmarking/benchmarking_data/BC" 
            data_name_list = temp_section_list
            save_path = "../Results" 
            n_domains = 20 
            deepen = run(save_path = save_path, 
                pca_n_comps = 200,
                pre_epochs = 800, #### According to your own hardware, choose the number of training
                epochs = 1000, #### According to your own hardware, choose the number of training
                platform = "Visium",
                )
            adata, graph_dict, domains = deepen._get_multiple_adata(data_path, data_name_list)
            adata = deepen._fit(adata, graph_dict, domains, pretrain = True)
            adata = deepen._get_cluster_data(adata, n_domains = n_domains, priori=True)


            np.save(os.path.join(save_dir_r, 'deepst_' + sec1 + '_' + sec2 + '_iter_' + str(e) + '_deepst_embedding.npy'), adata.obsm['DeepST_embed'])
                        
            # get/save ARI, purity, FMS
            for section_id in temp_section_list:
                if section_id not in dict_out[sec1 + '_' + sec2].keys():
                    dict_out[sec1 + '_' + sec2][section_id] = {}
                temp_adata = adata.obs[adata.obs['batch_name'] == section_id]
                temp_obs = temp_adata.dropna()
                save_df = temp_obs[['original_clusters', 'DeepST_refine_domain']]
                save_df.to_csv(os.path.join(save_dir_r, 'deepst_' + section_id + '_in_' + sec1 + '_' + sec2 + '_iter_' + str(e) + '_deepst_out.csv'), index=True)

                temp_ARI = adjusted_rand_score(temp_obs['DeepST_refine_domain'], temp_obs['original_clusters'])
                FMS = fowlkes_mallows_score(temp_obs['DeepST_refine_domain'], temp_obs['original_clusters'])
                PS = purity_score(temp_obs['DeepST_refine_domain'], temp_obs['original_clusters'])
                print('ARI of section ID %s: %.3f' %(section_id, temp_ARI))
                print('FMS of section ID %s: %.3f' %(section_id, FMS))
                print('PS of section ID %s: %.3f' %(section_id, PS))
                if 'ARI' not in dict_out[sec1 + '_' + sec2][section_id].keys():
                    dict_out[sec1 + '_' + sec2][section_id]['ARI'] = [temp_ARI]
                else:
                    dict_out[sec1 + '_' + sec2][section_id]['ARI'].append(temp_ARI)
                
                if 'FMS' not in dict_out[sec1 + '_' + sec2][section_id].keys():
                    dict_out[sec1 + '_' + sec2][section_id]['FMS'] = [FMS]
                else:
                    dict_out[sec1 + '_' + sec2][section_id]['FMS'].append(FMS)
                
                if 'PS' not in dict_out[sec1 + '_' + sec2][section_id].keys():
                    dict_out[sec1 + '_' + sec2][section_id]['PS'] = [PS]
                else:
                    dict_out[sec1 + '_' + sec2][section_id]['PS'].append(PS)

with open('deepst_result_bc.json', 'w', encoding='utf-8') as f:
    json.dump(dict_out, f, ensure_ascii=False, indent=4)

"""all slices result BC"""
# Nan, only 2 slices for BC

dict_out = {}
"""pair-wise result mHypo"""
section_list = ['-0.04', '-0.09', '-0.14', '-0.19', '-0.24']
for index_ in range(1, len(section_list)):
    sec1 = section_list[index_-1]
    sec2 = section_list[index_]
    temp_section_list = [sec1, sec2]

    if sec1 + '_' + sec2 not in dict_out.keys():
        dict_out[sec1 + '_' + sec2] = {}

    for section_id in temp_section_list:
        temp_section_list = [sec1, sec2]

        if sec1 + '_' + sec2 not in dict_out.keys():
            dict_out[sec1 + '_' + sec2] = {}
        

        for e in range(iter_):
            data_path = "/home/yunfei/spatial_benchmarking/benchmarking_data/mHypothalamus" 
            data_name_list = temp_section_list
            save_path = "../Results" 
            n_domains = 20 
            deepen = run(save_path = save_path, 
                pca_n_comps = 200,
                pre_epochs = 800, #### According to your own hardware, choose the number of training
                epochs = 1000, #### According to your own hardware, choose the number of training
                platform = "Visium",
                )
            adata, graph_dict, domains = deepen._get_multiple_adata(data_path, data_name_list)
            adata = deepen._fit(adata, graph_dict, domains, pretrain = True)
            adata = deepen._get_cluster_data(adata, n_domains = n_domains, priori=True)


            np.save(os.path.join(save_dir_r, 'deepst_' + sec1 + '_' + sec2 + '_iter_' + str(e) + '_deepst_embedding.npy'), adata.obsm['DeepST_embed'])
                        
            # get/save ARI, purity, FMS
            for section_id in temp_section_list:
                if section_id not in dict_out[sec1 + '_' + sec2].keys():
                    dict_out[sec1 + '_' + sec2][section_id] = {}
                temp_adata = adata.obs[adata.obs['batch_name'] == section_id]
                temp_obs = temp_adata.dropna()
                save_df = temp_obs[['original_clusters', 'DeepST_refine_domain']]
                save_df.to_csv(os.path.join(save_dir_r, 'deepst_' + section_id + '_in_' + sec1 + '_' + sec2 + '_iter_' + str(e) + '_deepst_out.csv'), index=True)

                temp_ARI = adjusted_rand_score(temp_obs['DeepST_refine_domain'], temp_obs['original_clusters'])
                FMS = fowlkes_mallows_score(temp_obs['DeepST_refine_domain'], temp_obs['original_clusters'])
                PS = purity_score(temp_obs['DeepST_refine_domain'], temp_obs['original_clusters'])
                print('ARI of section ID %s: %.3f' %(section_id, temp_ARI))
                print('FMS of section ID %s: %.3f' %(section_id, FMS))
                print('PS of section ID %s: %.3f' %(section_id, PS))
                if 'ARI' not in dict_out[sec1 + '_' + sec2][section_id].keys():
                    dict_out[sec1 + '_' + sec2][section_id]['ARI'] = [temp_ARI]
                else:
                    dict_out[sec1 + '_' + sec2][section_id]['ARI'].append(temp_ARI)
                
                if 'FMS' not in dict_out[sec1 + '_' + sec2][section_id].keys():
                    dict_out[sec1 + '_' + sec2][section_id]['FMS'] = [FMS]
                else:
                    dict_out[sec1 + '_' + sec2][section_id]['FMS'].append(FMS)
                
                if 'PS' not in dict_out[sec1 + '_' + sec2][section_id].keys():
                    dict_out[sec1 + '_' + sec2][section_id]['PS'] = [PS]
                else:
                    dict_out[sec1 + '_' + sec2][section_id]['PS'].append(PS)
with open('deepst_result_mhypo_pairwise.json', 'w', encoding='utf-8') as f:
    json.dump(dict_out, f, ensure_ascii=False, indent=4)


"""all slices result mHypo"""
sets_ = [['-0.04', '-0.09', '-0.14', '-0.19', '-0.24']]
for section_list in sets_:
    temp_section_list = section_list
    sec1, sec2, sec3, sec4, sec5 = temp_section_list

    if sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 + '_' + sec5 not in dict_out.keys():
        dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 + '_' + sec5] = {}

    for section_id in temp_section_list:
        temp_section_list = [sec1, sec2]

        if sec1 + '_' + sec2 not in dict_out.keys():
            dict_out[sec1 + '_' + sec2] = {}
        

        for e in range(iter_):
            data_path = "/home/yunfei/spatial_benchmarking/benchmarking_data/mHypothalamus" 
            data_name_list = temp_section_list
            save_path = "../Results" 
            n_domains = 20 
            deepen = run(save_path = save_path, 
                pca_n_comps = 200,
                pre_epochs = 800, #### According to your own hardware, choose the number of training
                epochs = 1000, #### According to your own hardware, choose the number of training
                platform = "Visium",
                )
            adata, graph_dict, domains = deepen._get_multiple_adata(data_path, data_name_list)
            adata = deepen._fit(adata, graph_dict, domains, pretrain = True)
            adata = deepen._get_cluster_data(adata, n_domains = n_domains, priori=True)

            np.save(os.path.join(save_dir_r, 'deepst_' + sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 + '_' + sec5 + '_iter_' + str(e) + '_deepst_embedding.npy'), adata.obsm['DeepST_embed'])
                        
            # get/save ARI, purity, FMS
            for section_id in temp_section_list:
                if section_id not in dict_out[sec1 + '_' + sec2].keys():
                    dict_out[sec1 + '_' + sec2][section_id] = {}
                temp_adata = adata.obs[adata.obs['batch_name'] == section_id]
                temp_obs = temp_adata.dropna()
                save_df = temp_obs[['original_clusters', 'DeepST_refine_domain']]
                save_df.to_csv(os.path.join(save_dir_r, 'deepst_' + sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 + '_' + sec5 + '_iter_' + str(e) + '_deepst_out.csv'), index=True)

                temp_ARI = adjusted_rand_score(temp_obs['DeepST_refine_domain'], temp_obs['original_clusters'])
                FMS = fowlkes_mallows_score(temp_obs['DeepST_refine_domain'], temp_obs['original_clusters'])
                PS = purity_score(temp_obs['DeepST_refine_domain'], temp_obs['original_clusters'])
                print('ARI of section ID %s: %.3f' %(section_id, temp_ARI))
                print('FMS of section ID %s: %.3f' %(section_id, FMS))
                print('PS of section ID %s: %.3f' %(section_id, PS))
                if 'ARI' not in dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 + '_' + sec5][section_id].keys():
                    dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 + '_' + sec5][section_id]['ARI'] = [temp_ARI]
                else:
                    dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 + '_' + sec5][section_id]['ARI'].append(temp_ARI)
                
                if 'FMS' not in dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 + '_' + sec5][section_id].keys():
                    dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 + '_' + sec5][section_id]['FMS'] = [FMS]
                else:
                    dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 + '_' + sec5][section_id]['FMS'].append(FMS)
                
                if 'PS' not in dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 + '_' + sec5][section_id].keys():
                    dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 + '_' + sec5][section_id]['PS'] = [PS]
                else:
                    dict_out[sec1 + '_' + sec2 + '_' + sec3 + '_' + sec4 + '_' + sec5][section_id]['PS'].append(PS)
with open('stagate_result_mhypo_all.json', 'w', encoding='utf-8') as f:
    json.dump(dict_out, f, ensure_ascii=False, indent=4)
