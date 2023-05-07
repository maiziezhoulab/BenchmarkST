import numpy as np
import pandas as pd
from tqdm import tqdm
import scipy.sparse as sp
from sklearn.metrics.cluster import adjusted_rand_score
from .GAAE import GAAE_mod2, GAAE
from .utils import mclust_R, Transfer_pytorch_Data
import torch
import torch.backends.cudnn as cudnn
import torch.nn.functional as F
import scanpy as sc
from scipy.stats import ttest_ind
from scipy import stats
cudnn.deterministic = True
cudnn.benchmark = True
import time


def test(data1, data2):
    # mannwhitneyu
    stat, p = stats.mannwhitneyu(data1, data2)
    return p


def p_val_finder(test_array, ctrl_array, gene_test_num):
    test_vals = test_array[gene_test_num, :]
    ctrl_vals = ctrl_array[gene_test_num, :]

    p_val = test(test_vals, ctrl_vals)
    log_p_val = -np.log10(p_val)

    if np.mean(test_vals) > np.mean(ctrl_vals):  # overexpressed gene in test group
        cur = True
    else:  # underexpressed gene in test group
        cur = False
    return p_val, log_p_val, cur


def dif_gene_analysis(test_array, ctrl_array, gene_list, topK):
    gene_len = len(gene_list)
    dif_gene_array = np.zeros([gene_len, 3])
    dif_gene_list = list()
    for i in range(gene_len):
        dif_gene_list.append(gene_list[i])
        dif_gene_array[i, 0], dif_gene_array[i, 1], dif_gene_array[i, 2] = p_val_finder(test_array, ctrl_array, i)

    tot_df = pd.DataFrame([dif_gene_list, dif_gene_array[:, 0], dif_gene_array[:, 1], dif_gene_array[:, 2]]).transpose()
    tot_df = tot_df.set_axis(['Gene', 'p-val', '-log10(p-val)', 'cur'], axis=1, inplace=False)

    sorted_df = tot_df.sort_values(by=['cur', '-log10(p-val)'], ascending=False)

    top200_genes = sorted_df['Gene'][0:topK]

    # return sorted_df
    return top200_genes


def train_ADEPT_mod(adata, hidden_dims=None, n_epochs=750, lr=0.001,
                  gradient_clipping=5., weight_decay=0.0001, verbose=True,
                  random_seed=0, save_loss=False, save_reconstrction=False,
                  device=torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')):
    """\
    Training graph attention auto-encoder.

    Parameters
    ----------
    adata
        AnnData object of scanpy package.
    hidden_dims
        The dimension of the encoder.
    n_epochs
        Number of total epochs in training.
    lr
        Learning rate for AdamOptimizer.
    key_added
        The latent embeddings are saved in adata.obsm[key_added].
    gradient_clipping
        Gradient Clipping.
    weight_decay
        Weight decay for AdamOptimizer.
    save_loss
        If True, the training loss is saved in adata.uns['STAGATE_loss'].
    save_reconstrction
        If True, the reconstructed expression profiles are saved in adata.layers['STAGATE_ReX'].
    device
        See torch.device.

    Returns
    -------
    AnnData
    """

    # seed_everything()
    if hidden_dims is None:
        hidden_dims = [512, 30]
    seed = random_seed
    import random
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)

    adata.X = sp.csr_matrix(adata.X)

    if 'highly_variable' in adata.var.columns:
        adata_Vars = adata[:, adata.var['highly_variable']]
    else:
        adata_Vars = adata

    if verbose:
        print('Size of Input: ', adata_Vars.shape)
    if 'Spatial_Net' not in adata.uns.keys():
        raise ValueError("Spatial_Net is not existed! Run Cal_Spatial_Net first!")

    data = Transfer_pytorch_Data(adata_Vars)

    model = GAAE_mod2(hidden_dims=[adata_Vars.X.shape[1]]+hidden_dims, mlp_dims=7).to(device)
    # model = STAGATE_mod1(hidden_dims=hidden_dims, mlp_dims=7).to(device)
    data = data.to(device)
    print(model)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

    loss_list = []
    for epoch in tqdm(range(1, n_epochs + 1)):
        model.train()
        optimizer.zero_grad()
        embed, recon, pred = model(data.x, data.edge_index)
        loss = F.mse_loss(data.x, recon)  # F.nll_loss(out[data.train_mask], data.y[data.train_mask])
        loss_list.append(loss)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), gradient_clipping)
        optimizer.step()

    model.eval()
    embed, recon, pred = model(data.x, data.edge_index)

    temp = embed.to('cpu').detach().numpy()
    adata.obsm['test'] = temp
    adata = mclust_R(adata, num_cluster=7, used_obsm='test')
    out = adata.obs.dropna()
    print(out['mclust'], out['Ground Truth'])
    ARI_ini = adjusted_rand_score(out['mclust'], out['Ground Truth'])
    print("ARI = ", ARI_ini)
    lbs = torch.tensor(adata.obs['mclust'].to_numpy()) - 1
    lbs = lbs.to(device)

    for pseudo_epoch in range(750):
        model.train()
        optimizer.zero_grad()

        embed, recon, pred = model(data.x, data.edge_index)
        # m =
        # pseudo classification loss
        # print(pred.size(), lbs.size())
        loss = F.nll_loss(pred[data.train_mask], lbs[data.train_mask])
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), gradient_clipping)
        optimizer.step()
        if pseudo_epoch % 50 == 49:
            model.eval()
            pred = model(data.x, data.edge_index)[-1].argmax(dim=1)
            correct = torch.count_nonzero(pred[data.val_mask] == lbs[data.val_mask])
            acc = int(correct) / int(data.val_mask.sum())
            print(f'Validation Accuracy: {acc:.4f}')
    model.eval()
    pred = model(data.x, data.edge_index)[-1].argmax(dim=1)
    correct = torch.count_nonzero(pred[data.test_mask] == lbs[data.test_mask])
    acc = int(correct) / int(data.test_mask.sum())
    print(f'Test Accuracy: {acc:.4f}')

    model.eval()
    embed, recon, pred = model(data.x, data.edge_index)

    temp = pred.to('cpu').detach().numpy()
    adata.obsm['test'] = temp
    adata = mclust_R(adata, num_cluster=7, used_obsm='test')
    out = adata.obs.dropna()
    ARI = adjusted_rand_score(out['mclust'], out['Ground Truth'])
    # print("training ", sample[section_index][1], " done!")
    print("finalized ARI = ", ARI)
    # if save_loss:
    #     adata.uns['STAGATE_loss'] = loss
    # if save_reconstrction:
    #     ReX = out.to('cpu').detach().numpy()
    #     ReX[ReX < 0] = 0
    #     adata.layers['STAGATE_ReX'] = ReX

    return ARI_ini, ARI


def train_ADEPT_use_DE(adata, hidden_dims=None, n_epochs=1000, lr=0.001, num_cluster=7,
                     gradient_clipping=5., weight_decay=0.0001, verbose=True, dif_k=200,
                     random_seed=0, device_id='0'):
    """
    Training graph attention auto-encoder.

    Parameters
    ----------
    adata
        AnnData object of scanpy package.
    hidden_dims
        The dimension of the encoder.
    n_epochs
        Number of total epochs in training.
    lr
        Learning rate for AdamOptimizer.
    key_added
        The latent embeddings are saved in adata.obsm[key_added].
    gradient_clipping
        Gradient Clipping.
    weight_decay
        Weight decay for AdamOptimizer.
    save_loss
        If True, the training loss is saved in adata.uns['STAGATE_loss'].
    save_reconstrction
        If True, the reconstructed expression profiles are saved in adata.layers['STAGATE_ReX'].
    device
        See torch.device.

    Returns
    -------
    AnnData
    """
    device = torch.device('cuda:' + device_id if torch.cuda.is_available() else 'cpu')
    # seed_everything()
    if hidden_dims is None:
        hidden_dims = [512, 30]
    seed = random_seed
    import random
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)

    adata.X = sp.csr_matrix(adata.X)

    if 'highly_variable' in adata.var.columns:
        adata_Vars = adata[:, adata.var['highly_variable']]
    else:
        adata_Vars = adata

    if verbose:
        print('Size of Input: ', adata_Vars.shape)
    if 'Spatial_Net' not in adata.uns.keys():
        raise ValueError("Spatial_Net is not existed! Run Cal_Spatial_Net first!")

    data = Transfer_pytorch_Data(adata_Vars)

    model = GAAE(hidden_dims=[adata_Vars.X.shape[1]]+hidden_dims).to(device)
    data = data.to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

    loss_list = []
    for epoch in tqdm(range(1, n_epochs + 1)):
        model.train()
        optimizer.zero_grad()
        embed, recon = model(data.x, data.edge_index)
        loss = F.mse_loss(data.x, recon)  # F.nll_loss(out[data.train_mask], data.y[data.train_mask])
        loss_list.append(loss.to('cpu').detach())
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), gradient_clipping)
        optimizer.step()

    model.eval()
    embed, recon = model(data.x, data.edge_index)
    temp = embed.to('cpu').detach().numpy()
    adata.obsm['ade'] = temp
    adata = mclust_R(adata, num_cluster=num_cluster, used_obsm='ade', save_obs='mclust')
    out = adata.obs.dropna()
    # print(out['mclust'])
    ARI_ini = adjusted_rand_score(out['mclust'], out['Ground Truth'])
    print("ARI = ", ARI_ini)

    # get DE from labels
    print("selecting DEs after first round of clustering ...")
    dif_genes_by_cluster_num = []
    for i_ in range(num_cluster):
        # print(set(adata.obs['mclust'].to_list()))
        # print(adata.obs['mclust'] == i_+1)
        keep_bcs_ = adata.obs.loc[adata.obs['mclust'] == i_+1].index

        keep_bcs_control = adata.obs.loc[adata.obs['mclust'] != i_+1].index
        #     print(keep_bcs_)
        adata_new = adata[keep_bcs_].copy()
        test_array = adata_new.X.toarray().T

        adata_gt_new_control = adata[keep_bcs_control].copy()
        ctrl_array = adata_gt_new_control.X.toarray().T

        gene_list = adata.var.index[:, np.newaxis]
        # print(test_array.shape, ctrl_array.shape)
        # print(gene_list)
        print("label = ", i_+1)
        # if num_cluster == 7:
        #     dif_k = 200
        # if num_cluster > 10:
        #     dif_k = 70
        dif_genes = dif_gene_analysis(test_array, ctrl_array, gene_list, dif_k)

        dif_genes_new = []
        for e in dif_genes:
            dif_genes_new.append(e[0])
        dif_genes_by_cluster_num.append(dif_genes_new)
    union_ = set()
    for i in range(num_cluster):
        union_ = union_.union(set(dif_genes_by_cluster_num[i]))
    print("selecting DE num = ", len(union_))
    print("train a new AE with the potential DEs ...")
    adata_new = adata[:, list(union_)]
    nzr_ = np.count_nonzero(adata_new.X.todense())/(adata_new.X.shape[0]*adata_new.X.shape[1])
    print("DE expression matrix non-zero rate = ", str(nzr_))
    data2 = Transfer_pytorch_Data(adata_new)

    model2 = GAAE(hidden_dims=[adata_new.X.shape[1]] + hidden_dims).to(device)
    # model = STAGATE_mod1(hidden_dims=hidden_dims, mlp_dims=7).to(device)
    data2 = data2.to(device)
    # print(model2)

    optimizer = torch.optim.Adam(model2.parameters(), lr=lr, weight_decay=weight_decay)

    loss_list = []
    for epoch in tqdm(range(1, n_epochs + 1)):
        model2.train()
        optimizer.zero_grad()
        embed, recon = model2(data2.x, data2.edge_index)
        loss = F.mse_loss(data2.x, recon)  # F.nll_loss(out[data.train_mask], data.y[data.train_mask])
        loss_list.append(loss)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model2.parameters(), gradient_clipping)
        optimizer.step()

    model2.eval()
    embed, recon = model2(data2.x, data2.edge_index)

    temp = embed.to('cpu').detach().numpy()
    adata.obsm['ade_impute'] = temp
    adata = mclust_R(adata, num_cluster=num_cluster, used_obsm='ade_impute', save_obs='mclust_impute')
    out = adata.obs.dropna()
    # print(out['mclust_impute'])
    ARI = adjusted_rand_score(out['mclust_impute'], out['Ground Truth'])
    # print("training ", sample[section_index][1], " done!")
    print("finalized ARI = ", ARI)
    # if save_loss:
    #     adata.uns['STAGATE_loss'] = loss
    # if save_reconstrction:
    #     ReX = out.to('cpu').detach().numpy()
    #     ReX[ReX < 0] = 0
    #     adata.layers['STAGATE_ReX'] = ReX

    return ARI_ini, ARI, union_, adata.copy()


def train_ADEPT(adata, hidden_dims=None, n_epochs=1000, lr=0.001, num_cluster=7,
                gradient_clipping=5., weight_decay=0.0001, verbose=True,
                random_seed=0, device_id='0'):
    """
    Training graph attention auto-encoder.

    Parameters
    ----------
    adata
        AnnData object of scanpy package.
    hidden_dims
        The dimension of the encoder.
    n_epochs
        Number of total epochs in training.
    lr
        Learning rate for AdamOptimizer.
    key_added
        The latent embeddings are saved in adata.obsm[key_added].
    gradient_clipping
        Gradient Clipping.
    weight_decay
        Weight decay for AdamOptimizer.
    save_loss
        If True, the training loss is saved in adata.uns['STAGATE_loss'].
    save_reconstrction
        If True, the reconstructed expression profiles are saved in adata.layers['STAGATE_ReX'].
    device
        See torch.device.

    Returns
    -------
    AnnData
    """
    device = torch.device('cuda:' + device_id if torch.cuda.is_available() else 'cpu')
    # seed_everything()
    if hidden_dims is None:
        hidden_dims = [512, 30]
    seed = random_seed
    import random
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)

    adata.X = sp.csr_matrix(adata.X)

    if 'highly_variable' in adata.var.columns:
        adata_Vars = adata[:, adata.var['highly_variable']]
    else:
        adata_Vars = adata

    if verbose:
        print('Size of Input: ', adata_Vars.shape)
    if 'Spatial_Net' not in adata.uns.keys():
        raise ValueError("Spatial_Net is not existed! Run Cal_Spatial_Net first!")

    data = Transfer_pytorch_Data(adata_Vars)

    model = GAAE(hidden_dims=[adata_Vars.X.shape[1]]+hidden_dims).to(device)
    data = data.to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

    loss_list = []
    for epoch in tqdm(range(1, n_epochs + 1)):
        model.train()
        optimizer.zero_grad()
        embed, recon = model(data.x, data.edge_index)
        loss = F.mse_loss(data.x, recon)  # F.nll_loss(out[data.train_mask], data.y[data.train_mask])
        loss_list.append(loss.to('cpu').detach())
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), gradient_clipping)
        optimizer.step()

    model.eval()
    embed, recon = model(data.x, data.edge_index)
    temp = embed.to('cpu').detach().numpy()
    adata.obsm['ade'] = temp
    adata = mclust_R(adata, num_cluster=num_cluster, used_obsm='ade', save_obs='mclust')
    out = adata.obs.dropna()
    # print(out['mclust'])
    ARI_ini = adjusted_rand_score(out['mclust'], out['Ground Truth'])
    print("ARI = ", ARI_ini)

    return ARI_ini, adata.copy()


def DE_nzr(adata, hidden_dims=None, n_epochs=1000, lr=0.001, num_cluster=7,
                     gradient_clipping=5., weight_decay=0.0001, verbose=True, dif_k=200,
                     random_seed=0, device_id='0'):
    device = torch.device('cuda:' + device_id if torch.cuda.is_available() else 'cpu')
    # seed_everything()
    if hidden_dims is None:
        hidden_dims = [512, 30]
    seed = random_seed
    import random
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)

    adata.X = sp.csr_matrix(adata.X)

    if 'highly_variable' in adata.var.columns:
        adata_Vars = adata[:, adata.var['highly_variable']]
    else:
        adata_Vars = adata

    if verbose:
        print('Size of Input: ', adata_Vars.shape)
    if 'Spatial_Net' not in adata.uns.keys():
        raise ValueError("Spatial_Net is not existed! Run Cal_Spatial_Net first!")

    data = Transfer_pytorch_Data(adata_Vars)

    model = GAAE(hidden_dims=[adata_Vars.X.shape[1]]+hidden_dims).to(device)
    data = data.to(device)
    print(model)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

    loss_list = []
    for epoch in tqdm(range(1, n_epochs + 1)):
        model.train()
        optimizer.zero_grad()
        embed, recon = model(data.x, data.edge_index)
        loss = F.mse_loss(data.x, recon)  # F.nll_loss(out[data.train_mask], data.y[data.train_mask])
        loss_list.append(loss.to('cpu').detach())
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), gradient_clipping)
        optimizer.step()

    model.eval()
    embed, recon = model(data.x, data.edge_index)
    temp = embed.to('cpu').detach().numpy()
    adata.obsm['test'] = temp
    adata = mclust_R(adata, num_cluster=num_cluster, used_obsm='test')

    print("selecting DEs after first round of clustering ...")
    dif_genes_by_cluster_num = []
    for i_ in range(num_cluster):
        # print(set(adata.obs['mclust'].to_list()))
        # print(adata.obs['mclust'] == i_+1)
        keep_bcs_ = adata.obs.loc[adata.obs['mclust_impute'] == i_+1].index

        keep_bcs_control = adata.obs.loc[adata.obs['mclust_impute'] != i_+1].index
        #     print(keep_bcs_)
        adata_new = adata[keep_bcs_].copy()
        test_array = adata_new.X.toarray().T

        adata_gt_new_control = adata[keep_bcs_control].copy()
        ctrl_array = adata_gt_new_control.X.toarray().T

        gene_list = adata.var.index[:, np.newaxis]
        print(test_array.shape, ctrl_array.shape)
        # print(gene_list)
        print("label = ", i_+1)
        # if num_cluster == 7:
        #     dif_k = 200
        # if num_cluster > 10:
        #     dif_k = 70
        dif_genes = dif_gene_analysis(test_array, ctrl_array, gene_list, dif_k)

        dif_genes_new = []
        for e in dif_genes:
            dif_genes_new.append(e[0])
        dif_genes_by_cluster_num.append(dif_genes_new)
    union_ = set()
    for i in range(num_cluster):
        union_ = union_.union(set(dif_genes_by_cluster_num[i]))
    print("selecting DE num = ", len(union_))
    print("train a new AE with the potential DEs ...")
    adata_new = adata[:, list(union_)]
    nzr_ = np.count_nonzero(adata_new.X.todense())/(adata_new.X.shape[0]*adata_new.X.shape[1])
    print("DE expression matrix non-zero rate = ", str(nzr_))

    return nzr_


# def train_ADE_GO(adata, hidden_dims=None, n_epochs=1000, lr=0.001, num_cluster=4,
#                  gradient_clipping=5., weight_decay=0.0001, verbose=True, dif_k=200,
#                  random_seed=0, device=torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')):
#     """\
#     Training graph attention auto-encoder.

#     Parameters
#     ----------
#     adata
#         AnnData object of scanpy package.
#     hidden_dims
#         The dimension of the encoder.
#     n_epochs
#         Number of total epochs in training.
#     lr
#         Learning rate for AdamOptimizer.
#     key_added
#         The latent embeddings are saved in adata.obsm[key_added].
#     gradient_clipping
#         Gradient Clipping.
#     weight_decay
#         Weight decay for AdamOptimizer.
#     save_loss
#         If True, the training loss is saved in adata.uns['STAGATE_loss'].
#     save_reconstrction
#         If True, the reconstructed expression profiles are saved in adata.layers['STAGATE_ReX'].
#     device
#         See torch.device.

#     Returns
#     -------
#     AnnData
#     """

#     # seed_everything()
#     if hidden_dims is None:
#         hidden_dims = [512, 30]
#     seed = random_seed
#     import random
#     random.seed(seed)
#     torch.manual_seed(seed)
#     torch.cuda.manual_seed_all(seed)
#     np.random.seed(seed)

#     adata.X = sp.csr_matrix(adata.X)

#     if 'highly_variable' in adata.var.columns:
#         adata_Vars = adata[:, adata.var['highly_variable']]
#     else:
#         adata_Vars = adata

#     if verbose:
#         print('Size of Input: ', adata_Vars.shape)
#     if 'Spatial_Net' not in adata.uns.keys():
#         raise ValueError("Spatial_Net is not existed! Run Cal_Spatial_Net first!")

#     data = Transfer_pytorch_Data(adata_Vars)

#     model = GAAE(hidden_dims=[adata_Vars.X.shape[1]]+hidden_dims).to(device)
#     # model = STAGATE_mod1(hidden_dims=hidden_dims, mlp_dims=7).to(device)
#     data = data.to(device)
#     print(model)

#     optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

#     loss_list = []
#     for epoch in tqdm(range(1, n_epochs + 1)):
#         model.train()
#         optimizer.zero_grad()
#         embed, recon = model(data.x, data.edge_index)
#         loss = F.mse_loss(data.x, recon)  # F.nll_loss(out[data.train_mask], data.y[data.train_mask])
#         loss_list.append(loss.to('cpu').detach())
#         loss.backward()
#         torch.nn.utils.clip_grad_norm_(model.parameters(), gradient_clipping)
#         optimizer.step()

#     model.eval()
#     embed, recon = model(data.x, data.edge_index)
#     temp = embed.to('cpu').detach().numpy()
#     adata.obsm['test'] = temp
#     adata = mclust_R(adata, num_cluster=num_cluster, used_obsm='test')
#     # out = adata.obs.dropna()
#     # print(out['mclust'], out['Ground Truth'])
#     # ARI_ini = adjusted_rand_score(out['mclust'], out['Ground Truth'])
#     # print("ARI = ", ARI_ini)

#     # get DE from labels
#     print("selecting DEs after first round of clustering ...")
#     dif_genes_by_cluster_num = []
#     for i_ in range(num_cluster):
#         # print(set(adata.obs['mclust'].to_list()))
#         # print(adata.obs['mclust'] == i_+1)
#         keep_bcs_ = adata.obs.loc[adata.obs['mclust'] == i_+1].index

#         keep_bcs_control = adata.obs.loc[adata.obs['mclust'] != i_+1].index
#         #     print(keep_bcs_)
#         adata_new = adata[keep_bcs_].copy()
#         test_array = adata_new.X.toarray().T

#         adata_gt_new_control = adata[keep_bcs_control].copy()
#         ctrl_array = adata_gt_new_control.X.toarray().T

#         gene_list = adata.var.index[:, np.newaxis]
#         print(test_array.shape, ctrl_array.shape)
#         # print(gene_list)
#         print("label = ", i_+1)
#         # if num_cluster == 7:
#         #     dif_k = 200
#         # if num_cluster > 10:
#         #     dif_k = 70
#         dif_genes = dif_gene_analysis(test_array, ctrl_array, gene_list, dif_k)

#         dif_genes_new = []
#         for e in dif_genes:
#             dif_genes_new.append(e[0])
#         dif_genes_by_cluster_num.append(dif_genes_new)
#     union_ = set()
#     for i in range(num_cluster):
#         union_ = union_.union(set(dif_genes_by_cluster_num[i]))

#     print("train a new AE with the potential DEs ...")
#     adata_new = adata[:, list(union_)]
#     data2 = Transfer_pytorch_Data(adata_new)

#     model2 = GAAE(hidden_dims=[adata_new.X.shape[1]] + hidden_dims).to(device)
#     # model = GAAE_mod1(hidden_dims=hidden_dims, mlp_dims=7).to(device)
#     data2 = data2.to(device)
#     print(model2)

#     optimizer = torch.optim.Adam(model2.parameters(), lr=lr, weight_decay=weight_decay)

#     loss_list = []
#     for epoch in tqdm(range(1, n_epochs + 1)):
#         model2.train()
#         optimizer.zero_grad()
#         embed, recon = model2(data2.x, data2.edge_index)
#         loss = F.mse_loss(data2.x, recon)  # F.nll_loss(out[data.train_mask], data.y[data.train_mask])
#         loss_list.append(loss)
#         loss.backward()
#         torch.nn.utils.clip_grad_norm_(model2.parameters(), gradient_clipping)
#         optimizer.step()

#     model2.eval()
#     embed, recon = model2(data2.x, data2.edge_index)

#     temp = embed.to('cpu').detach().numpy()
#     adata.obsm['test'] = temp
#     adata = mclust_R(adata, num_cluster=num_cluster, used_obsm='test')
#     # out = adata.obs.dropna()

#     # get DE from labels
#     print("selecting DEs after first round of clustering ...")
#     dif_genes_by_cluster_num2 = []
#     for i_ in range(num_cluster):
#         # print(set(adata.obs['mclust'].to_list()))
#         # print(adata.obs['mclust'] == i_+1)
#         keep_bcs_ = adata.obs.loc[adata.obs['mclust'] == i_ + 1].index

#         keep_bcs_control = adata.obs.loc[adata.obs['mclust'] != i_ + 1].index
#         #     print(keep_bcs_)
#         adata_new = adata[keep_bcs_].copy()
#         test_array = adata_new.X.toarray().T

#         adata_gt_new_control = adata[keep_bcs_control].copy()
#         ctrl_array = adata_gt_new_control.X.toarray().T

#         gene_list = adata.var.index[:, np.newaxis]
#         print(test_array.shape, ctrl_array.shape)
#         # print(gene_list)
#         print("label = ", i_ + 1)
#         # if num_cluster == 7:
#         #     dif_k = 200
#         # if num_cluster > 10:
#         #     dif_k = 70
#         dif_genes = dif_gene_analysis(test_array, ctrl_array, gene_list, dif_k)

#         dif_genes_new = []
#         for e in dif_genes:
#             dif_genes_new.append(e[0])
#         # print(dif_genes_new)
#         dif_genes_by_cluster_num2.append(dif_genes_new)
#         # print(dif_genes_by_cluster_num2)

#     return dif_genes_by_cluster_num2
