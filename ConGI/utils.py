import os
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import cv2
from skimage.feature import greycomatrix, greycoprops
from tqdm import trange
import numba
from scipy.sparse import issparse
from metrics import eval_mclust_ari
from sklearn.metrics import adjusted_rand_score
from st_loading_utils import load_DLPFC, load_BC, load_mVC, load_mPFC, load_mHypothalamus, load_her2_tumor, load_mMAMP

def load_ST_file(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True, file_adj=None):
    adata_h5 = sc.read_visium(file_fold, load_images=load_images, count_file=count_file)
    adata_h5.var_names_make_unique()

    if load_images is False:
        if file_adj is None:
            file_adj = os.path.join(file_fold, "spatial/tissue_positions_list.csv")

        positions = pd.read_csv(file_adj, header=None)
        positions.columns = [
            'barcode',
            'in_tissue',
            'array_row',
            'array_col',
            'pxl_col_in_fullres',
            'pxl_row_in_fullres',
        ]
        positions.index = positions['barcode']
        adata_h5.obs = adata_h5.obs.join(positions, how="left")
        adata_h5.obsm['spatial'] = adata_h5.obs[['pxl_row_in_fullres', 'pxl_col_in_fullres']].to_numpy()
        adata_h5.obs.drop(columns=['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'], inplace=True)

    # print('adata: (' + str(adata_h5.shape[0]) + ', ' + str(adata_h5.shape[1]) + ')')
    return adata_h5


def build_her2st_data(path, name, size=112):
    cnt_path = os.path.join(path, 'ST-cnts', f'{name}.tsv')
    df_cnt = pd.read_csv(cnt_path, sep='\t', index_col=0)

    pos_path = os.path.join(path, 'ST-spotfiles', f'{name}_selection.tsv')
    df_pos = pd.read_csv(pos_path, sep='\t')

    lbl_path = os.path.join(path, 'ST-pat/lbl', f'{name}_labeled_coordinates.tsv')
    df_lbl = pd.read_csv(lbl_path, sep='\t')
    df_lbl = df_lbl.dropna(axis=0, how='any')
    df_lbl.loc[df_lbl['label'] == 'undetermined', 'label'] = np.nan
    df_lbl['x'] = (df_lbl['x']+0.5).astype(np.int64)
    df_lbl['y'] = (df_lbl['y']+0.5).astype(np.int64)

    x = df_pos['x'].values
    y = df_pos['y'].values
    ids = []
    for i in range(len(x)):
        ids.append(str(x[i])+'x'+str(y[i])) 
    df_pos['id'] = ids

    x = df_lbl['x'].values
    y = df_lbl['y'].values
    ids = []
    for i in range(len(x)):
        ids.append(str(x[i])+'x'+str(y[i])) 
    df_lbl['id'] = ids

    meta_pos = df_cnt.join(df_pos.set_index('id'))
    meta_lbl = df_cnt.join(df_lbl.set_index('id'))

    adata = anndata.AnnData(df_cnt, dtype=np.int64)
    adata.obsm['spatial'] = np.floor(meta_pos[['pixel_x','pixel_y']].values).astype(int)
    adata.obs['label'] = pd.Categorical(meta_lbl['label']).codes
    
    img_path = os.path.join(path, 'ST-imgs', name[0], name)
    full_image = cv2.imread(os.path.join(img_path, os.listdir(img_path)[0]))
    full_image = cv2.cvtColor(full_image, cv2.COLOR_BGR2RGB)
    patches = []
    for x, y in adata.obsm['spatial']:
        patches.append(full_image[y-size:y+size, x-size:x+size])
    patches = np.array(patches)
    
    return adata, patches


def adata_preprocess_hvg(adata, n_top_genes):
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_top_genes)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    return adata[:, adata.var['highly_variable']].X


def adata_preprocess_pca(i_adata, min_cells=3, pca_n_comps=300):
    sc.pp.filter_genes(i_adata, min_cells=min_cells)
    # for tumor data
    # if pca_n_comps > len(i_adata.var.index):
    #     pca_n_comps = len(i_adata.var.index)-1
    #     # print(params.cell_feat_dim)
    # if pca_n_comps > len(i_adata.obs.index):
    #     pca_n_comps = len(i_adata.obs.index)-1
    #     # print(params.cell_feat_dim)
    adata_X = sc.pp.normalize_total(i_adata, target_sum=1, exclude_highly_expressed=True, inplace=False)['X']
    adata_X = sc.pp.scale(adata_X)
    adata_X = sc.pp.pca(adata_X, n_comps=pca_n_comps)
    
    return adata_X


def extract_wash_patches(path, name, size=112, c=0.1):
    adata_h5 = load_ST_file(os.path.join(path, name))
    
    patch_name = f'{name}_{size}_patches.npy'
    if not os.path.exists(os.path.join(path, name, patch_name)):
        full_image = cv2.imread(os.path.join(path, name, f'{name}_full_image.tif'))
        full_image = cv2.cvtColor(full_image, cv2.COLOR_BGR2RGB)
        patches = []
        for x, y in adata_h5.obsm['spatial']:
            patches.append(full_image[y-size:y+size, x-size:x+size])
        patches = np.array(patches)
        np.save(os.path.join(path, name, patch_name), patches)
    patches = np.load(os.path.join(path, name, patch_name))
    
    idxs = []
    for i in trange(patches.shape[0], desc='Washing'):
        img = cv2.cvtColor(patches[i], cv2.COLOR_RGB2GRAY)
        glcm = graycomatrix(img, [16], [0], levels=256, symmetric=True, normed=True)
        if graycoprops(glcm, 'correlation')[0, 0] > c:
            idxs.append(i)
            
    wash_df = pd.DataFrame(adata_h5.obsm['spatial'])
    wash_df['wash'] = 0
    wash_df.loc[idxs, 'wash'] = 1
    
    for i in wash_df[wash_df['wash'] == 1].index:
        rep_idx = ((wash_df.loc[i, [0, 1]] - wash_df[[0, 1]]) ** 2).sum(1)[wash_df['wash'] == 0].sort_values().index[1]
        patches[i] = patches[rep_idx]
        
    np.save(os.path.join(path, name, f'{name}_{size}_patches_washed_{c}.npy'), patches)
    return patches


@numba.njit("f4(f4[:], f4[:])")
def euclid_dist(t1,t2):
    sum=0
    for i in range(t1.shape[0]):
        sum+=(t1[i]-t2[i])**2
    return np.sqrt(sum)

@numba.njit("f4[:,:](f4[:,:])", parallel=True, nogil=True)
def pairwise_distance(X):
    n=X.shape[0]
    adj=np.empty((n, n), dtype=np.float32)
    for i in numba.prange(n):
        for j in numba.prange(n):
            adj[i][j]=euclid_dist(X[i], X[j])
    return adj

def calculate_adj_matrix(x, y, x_pixel=None, y_pixel=None, image=None, beta=49, alpha=1, histology=True):
    #x,y,x_pixel, y_pixel are lists
    if histology:
        assert (x_pixel is not None) & (x_pixel is not None) & (image is not None)
        assert (len(x)==len(x_pixel)) & (len(y)==len(y_pixel))
        print("Calculateing adj matrix using histology image...")
        #beta to control the range of neighbourhood when calculate grey vale for one spot
        #alpha to control the color scale
        beta_half=round(beta/2)
        g=[]
        for i in range(len(x_pixel)):
            max_x=image.shape[0]
            max_y=image.shape[1]
            nbs=image[max(0,x_pixel[i]-beta_half):min(max_x,x_pixel[i]+beta_half+1),max(0,y_pixel[i]-beta_half):min(max_y,y_pixel[i]+beta_half+1)]
            g.append(np.mean(np.mean(nbs,axis=0),axis=0))
        c0, c1, c2=[], [], []
        for i in g:
            c0.append(i[0])
            c1.append(i[1])
            c2.append(i[2])
        c0=np.array(c0)
        c1=np.array(c1)
        c2=np.array(c2)
        print("Var of c0,c1,c2 = ", np.var(c0),np.var(c1),np.var(c2))
        c3=(c0*np.var(c0)+c1*np.var(c1)+c2*np.var(c2))/(np.var(c0)+np.var(c1)+np.var(c2))
        c4=(c3-np.mean(c3))/np.std(c3)
        z_scale=np.max([np.std(x), np.std(y)])*alpha
        z=c4*z_scale
        z=z.tolist()
        print("Var of x,y,z = ", np.var(x),np.var(y),np.var(z))
        X=np.array([x, y, z]).T.astype(np.float32)
    else:
        print("Calculateing adj matrix using xy only...")
        X=np.array([x, y]).T.astype(np.float32)
    return pairwise_distance(X)


def refine(sample_id, pred, dis, num_nbs):
    refined_pred=[]
    pred=pd.DataFrame({"pred": pred}, index=sample_id)
    dis_df=pd.DataFrame(dis, index=sample_id, columns=sample_id)
    for i in range(len(sample_id)):
        index=sample_id[i]
        dis_tmp=dis_df.loc[index, :].sort_values()
        nbs=dis_tmp[0:num_nbs+1]
        nbs_pred=pred.loc[nbs.index, "pred"]
        self_pred=pred.loc[index, "pred"]
        v_c=nbs_pred.value_counts()
        if (v_c.loc[self_pred]<num_nbs/2) and (np.max(v_c)>num_nbs/2):
            refined_pred.append(v_c.idxmax())
        else:           
            refined_pred.append(self_pred)
    return refined_pred

def count_nbr(target_cluster,cell_id, x, y, pred, radius):
    adj_2d=calculate_adj_matrix(x=x,y=y, histology=False)
    cluster_num = dict()
    df = {'cell_id': cell_id, 'x': x, "y":y, "pred":pred}
    df = pd.DataFrame(data=df)
    df.index=df['cell_id']
    target_df=df[df["pred"]==target_cluster]
    row_index=0
    num_nbr=[]
    for index, row in target_df.iterrows():
        x=row["x"]
        y=row["y"]
        tmp_nbr=df[((df["x"]-x)**2+(df["y"]-y)**2)<=(radius**2)]
        num_nbr.append(tmp_nbr.shape[0])
    return np.mean(num_nbr)

def search_radius(target_cluster,cell_id, x, y, pred, start, end, num_min=8, num_max=15,  max_run=100):
    run=0
    num_low=count_nbr(target_cluster,cell_id, x, y, pred, start)
    num_high=count_nbr(target_cluster,cell_id, x, y, pred, end)
    if num_min<=num_low<=num_max:
        print("recommended radius = ", str(start))
        return start
    elif num_min<=num_high<=num_max:
        print("recommended radius = ", str(end))
        return end
    elif num_low>num_max:
        print("Try smaller start.")
        return None
    elif num_high<num_min:
        print("Try bigger end.")
        return None
    while (num_low<num_min) and (num_high>num_min):
        run+=1
        print("Run "+str(run)+": radius ["+str(start)+", "+str(end)+"], num_nbr ["+str(num_low)+", "+str(num_high)+"]")
        if run >max_run:
            print("Exact radius not found, closest values are:\n"+"radius="+str(start)+": "+"num_nbr="+str(num_low)+"\nradius="+str(end)+": "+"num_nbr="+str(num_high))
            return None
        mid=(start+end)/2
        num_mid=count_nbr(target_cluster,cell_id, x, y, pred, mid)
        if num_min<=num_mid<=num_max:
            print("recommended radius = ", str(mid), "num_nbr="+str(num_mid))
            return mid
        if num_mid<num_min:
            start=mid
            num_low=num_mid
        elif num_mid>num_max:
            end=mid
            num_high=num_mid

def find_neighbor_clusters(target_cluster,cell_id, x, y, pred,radius, ratio=1/2):
    cluster_num = dict()
    for i in pred:
        cluster_num[i] = cluster_num.get(i, 0) + 1
    df = {'cell_id': cell_id, 'x': x, "y":y, "pred":pred}
    df = pd.DataFrame(data=df)
    df.index=df['cell_id']
    target_df=df[df["pred"]==target_cluster]
    nbr_num={}
    row_index=0
    num_nbr=[]
    for index, row in target_df.iterrows():
        x=row["x"]
        y=row["y"]
        tmp_nbr=df[((df["x"]-x)**2+(df["y"]-y)**2)<=(radius**2)]
        #tmp_nbr=df[(df["x"]<x+radius) & (df["x"]>x-radius) & (df["y"]<y+radius) & (df["y"]>y-radius)]
        num_nbr.append(tmp_nbr.shape[0])
        for p in tmp_nbr["pred"]:
            nbr_num[p]=nbr_num.get(p,0)+1
    del nbr_num[target_cluster]
    nbr_num=[(k, v)  for k, v in nbr_num.items() if v>(ratio*cluster_num[k])]
    nbr_num.sort(key=lambda x: -x[1])
    print("radius=", radius, "average number of neighbors for each spot is", np.mean(num_nbr))
    print(" Cluster",target_cluster, "has neighbors:")
    for t in nbr_num:
        print("Dmain ", t[0], ": ",t[1])
    ret=[t[0] for t in nbr_num]
    if len(ret)==0:
        print("No neighbor domain found, try bigger radius or smaller ratio.")
    else:
        return ret

def rank_genes_groups(input_adata, target_cluster,nbr_list, label_col, adj_nbr=True, log=False):
    if adj_nbr:
        nbr_list=nbr_list+[target_cluster]
        adata=input_adata[input_adata.obs[label_col].isin(nbr_list)]
    else:
        adata=input_adata.copy()
    adata.var_names_make_unique()
    adata.obs["target"]=((adata.obs[label_col]==target_cluster)*1).astype('category')
    sc.tl.rank_genes_groups(adata, groupby="target",reference="rest", n_genes=adata.shape[1],method='wilcoxon')
    pvals_adj=[i[0] for i in adata.uns['rank_genes_groups']["pvals_adj"]]
    genes=[i[1] for i in adata.uns['rank_genes_groups']["names"]]
    if issparse(adata.X):
        obs_tidy=pd.DataFrame(adata.X.A)
    else:
        obs_tidy=pd.DataFrame(adata.X)
    obs_tidy.index=adata.obs["target"].tolist()
    obs_tidy.columns=adata.var.index.tolist()
    obs_tidy=obs_tidy.loc[:,genes]
    # 1. compute mean value
    mean_obs = obs_tidy.groupby(level=0).mean()
    # 2. compute fraction of cells having value >0
    obs_bool = obs_tidy.astype(bool)
    fraction_obs = obs_bool.groupby(level=0).sum() / obs_bool.groupby(level=0).count()
    # compute fold change.
    if log: #The adata already logged
        fold_change=np.exp((mean_obs.loc[1] - mean_obs.loc[0]).values)
    else:
        fold_change = (mean_obs.loc[1] / (mean_obs.loc[0]+ 1e-9)).values
    df = {'genes': genes, 'in_group_fraction': fraction_obs.loc[1].tolist(), "out_group_fraction":fraction_obs.loc[0].tolist(),"in_out_group_ratio":(fraction_obs.loc[1]/fraction_obs.loc[0]).tolist(),"in_group_mean_exp": mean_obs.loc[1].tolist(), "out_group_mean_exp": mean_obs.loc[0].tolist(),"fold_change":fold_change.tolist(), "pvals_adj":pvals_adj}
    df = pd.DataFrame(data=df)
    return df


def get_predicted_results(dataset, name, path, z):

    if dataset in ["DLPFC", "BC", "MA"]:
        # adata = load_ST_file(os.path.join(path, name))
        # df_meta = pd.read_csv(os.path.join(path, name, 'metadata.tsv'), sep='\t')
        # label = pd.Categorical(df_meta['layer_guess']).codes
        if dataset == 'DLPFC':
            adata = load_DLPFC(path, section_id=name)
        elif dataset == 'BC':
            adata = load_BC(path, section_id=name)
        elif dataset == 'MA':
            adata = load_mMAMP(path, section_id=name)
        else:
            print("exception in dataset name")
            exit(-1)
        label = pd.Categorical(adata.obs['original_clusters']).codes
        n_clusters = label.max() + 1
        # adata = adata[label != -1]

        adj_2d = calculate_adj_matrix(x=adata.obs["array_row"].tolist(), y=adata.obs["array_col"].tolist(),
                                      histology=False)

        raw_preds = eval_mclust_ari(label[label != -1], z, n_clusters)

        if len(adata.obs)> 1000:
            num_nbs = 24
        else:
            num_nbs = 4

        refined_preds = refine(sample_id=adata.obs.index.tolist(), pred=raw_preds, dis=adj_2d, num_nbs=num_nbs)
        ari = adjusted_rand_score(label[label != -1], refined_preds)

        return ari, refined_preds

    elif dataset=="Her2st":
        adata, _ = build_her2st_data(path, name)
        label = adata.obs['label']
        n_clusters = label.max() + 1
        adata = adata[label != -1]

        adj_2d = calculate_adj_matrix(x=adata.obsm["spatial"][:, 0].tolist(), y=adata.obsm["spatial"][:, 1].tolist(),
                                      histology=False)

        raw_preds = eval_mclust_ari(label[label != -1], z, n_clusters)

        if len(adata.obs) > 1000:
            num_nbs = 24
        else:
            num_nbs = 4

        refined_preds = refine(sample_id=adata.obs.index.tolist(), pred=raw_preds, dis=adj_2d, num_nbs=num_nbs)
        ari = adjusted_rand_score(label[label != -1], refined_preds)

        return ari, refined_preds
    
    else:
        return None
        print("wrong dataset argument")
