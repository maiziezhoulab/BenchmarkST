import os
import numpy as np
import pandas as pd
from copy import deepcopy
import torch
from torch.utils import data
from torchvision import transforms
from torchtoolbox.transform import Cutout
import cv2
import scanpy as sc

from utils import load_ST_file, adata_preprocess_pca, adata_preprocess_hvg, extract_wash_patches, build_her2st_data
from st_loading_utils import load_DLPFC, load_mMAMP, load_BC


class Dataset(data.Dataset):
    def __init__(self, dataset, path, name, gene_preprocess='pca', n_genes=3000,
                 prob_mask=0.5, pct_mask=0.2, prob_noise=0.5, pct_noise=0.8, sigma_noise=0.5,
                 prob_swap=0.5, pct_swap=0.1, img_size=112, train=True):
        super(Dataset, self).__init__()
        
        # adata, label, image
        if dataset in ["MA", "DLPFC", "BC"]:
            if dataset == "DLPFC":
                adata = load_DLPFC(path, section_id=name)
                full_image = cv2.imread(os.path.join(path, name, f'{name}_full_image.tif'))
                full_image = cv2.cvtColor(full_image, cv2.COLOR_BGR2RGB)
                patches = []
                for x, y in adata.obsm['spatial']:
                    patches.append(full_image[y-img_size:y+img_size, x-img_size:x+img_size])
                patches = np.array(patches)
                self.label = adata.obs['original_clusters']
                self.image = patches
            elif dataset == "MA":
                adata = load_mMAMP(path, section_id=name)
                full_image = cv2.imread(os.path.join(path, name, f'{name}_full_image.tif'))
                full_image = cv2.cvtColor(full_image, cv2.COLOR_BGR2RGB)
                patches = []
                for x, y in adata.obsm['spatial']:
                    patches.append(full_image[y-img_size:y+img_size, x-img_size:x+img_size])
                patches = np.array(patches)
                self.label = adata.obs['original_clusters']
                self.image = patches
            elif dataset == "BC":
                adata = load_BC(path, section_id=name)
                full_image = cv2.imread(os.path.join(path, name, f'{name}_full_image.tif'))
                full_image = cv2.cvtColor(full_image, cv2.COLOR_BGR2RGB)
                patches = []
                for x, y in adata.obsm['spatial']:
                    patches.append(full_image[y-img_size:y+img_size, x-img_size:x+img_size])
                patches = np.array(patches)
                self.label = adata.obs['original_clusters']
                self.image = patches
            else:
                print("something unexpected ")
                exit(-1)
                

        elif dataset == "Her2st":
            adata, patches = build_her2st_data(path, name, img_size)
            self.label = adata.obs['label']
            self.image = patches

        # elif dataset in ["Mouse_brain_anterior"]:
        #     adata = sc.read_h5ad(os.path.join(path, f"{name}.h5ad"))
        #     adata.X = adata.X.A
        #     x = adata.obs['x4'].values
        #     y = adata.obs['x5'].values
        #     adata.obsm['spatial'] = np.stack((y, x), 1)
        #     # label
        #     self.label = np.zeros(adata.shape[0], dtype=int)
        #     # image
        #     full_image = cv2.imread(os.path.join(path, f"{name}_histology.tif"))
        #     full_image = cv2.cvtColor(full_image, cv2.COLOR_BGR2RGB)
        #     patches = []
        #     for x, y in adata.obsm['spatial']:
        #         patches.append(full_image[y-img_size:y+img_size, x-img_size:x+img_size])
        #     patches = np.array(patches)
        #     self.image = patches

        if dataset == "IDC":
            adata = load_ST_file(os.path.join(path, name))
            adata.X = adata.X.A
            self.label = np.zeros(adata.shape[0], dtype=int)

            # image
            full_image = cv2.imread(os.path.join(path, name, f'{name}.tif'))
            full_image = cv2.cvtColor(full_image, cv2.COLOR_BGR2RGB)
            patches = []
            for x, y in adata.obsm['spatial']:
                patches.append(full_image[y - img_size:y + img_size, x - img_size:x + img_size])
            patches = np.array(patches)
            self.image = patches

        # print(self.label.nunique())
        self.n_clusters = self.label.nunique()
        self.spatial = adata.obsm['spatial']
        self.n_pos = self.spatial.max() + 1
        
        # preprocess
        if gene_preprocess == 'pca':
            self.gene = adata_preprocess_pca(adata, pca_n_comps=n_genes).astype(np.float32)
        elif gene_preprocess == 'hvg':
            self.gene = np.array(adata_preprocess_hvg(adata, n_top_genes=n_genes)).astype(np.float32)
        
        self.train = train
        self.img_train_transform = transforms.Compose([
            Cutout(0.5),
            transforms.ToTensor(),
            transforms.RandomHorizontalFlip(p=0.5),
            transforms.RandomVerticalFlip(p=0.5),
            transforms.RandomApply([transforms.ColorJitter(0.4, 0.4, 0.4, 0.1)], p=0.8),
            transforms.RandomGrayscale(p=0.2),
            transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])
        ])
        self.img_test_transform = transforms.Compose([
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])
        ])
        
        self.gene_train_transform = GeneTransforms(n_genes, 
                                                    prob_mask=0.5, pct_mask=0.2,
                                                    prob_noise=0.5, pct_noise=0.8, sigma_noise=0.5,
                                                    prob_swap=0.5, pct_swap=0.1)
        
        self.gene = self.gene[self.label != -1]
        self.image = self.image[self.label != -1]
        self.label = self.label[self.label != -1]

    def __len__(self):
        return len(self.label)
    
    def __getitem__(self, idx):
        spatial = torch.from_numpy(self.spatial[idx])
        y = self.label[idx]
        
        if self.train:
            xg = self.gene[idx]
            xg_u = self.gene_train_transform(deepcopy(xg))
            xg_v = self.gene_train_transform(deepcopy(xg))
            
            xg = torch.from_numpy(xg)
            xg_u = torch.from_numpy(xg_u)
            xg_v = torch.from_numpy(xg_v)
            
            xi_u = self.img_train_transform(self.image[idx])
            xi_v = self.img_train_transform(self.image[idx])
            
            return xg, xg_u, xg_v, xi_u, xi_v, spatial, y, idx
        
        else:
            xg = self.gene[idx]
            xg = torch.from_numpy(xg)
            xi = self.img_test_transform(self.image[idx])
            
            return xg, xi, spatial, y, idx

    
class GeneTransforms(torch.nn.Module):
    def __init__(self, n_genes,
                prob_mask, pct_mask,
                prob_noise, pct_noise, sigma_noise,
                prob_swap, pct_swap):
        super(GeneTransforms, self).__init__()
        
        self.n_genes = n_genes
        self.prob_mask = prob_mask
        self.pct_mask = pct_mask
        self.prob_noise = prob_noise
        self.pct_noise = pct_noise
        self.sigma_noise = sigma_noise
        self.prob_swap = prob_swap
        self.pct_swap = pct_swap
        
    def build_mask(self, pct_mask):
        mask = np.concatenate([np.ones(int(self.n_genes * pct_mask), dtype=bool), 
                               np.zeros(self.n_genes - int(self.n_genes * pct_mask), dtype=bool)])
        np.random.shuffle(mask)
        return mask
        
    def forward(self, xg):
        if np.random.uniform(0, 1) < self.prob_mask:
            mask = self.build_mask(self.pct_mask)
            xg[mask] = 0
        
        if np.random.uniform(0, 1) < self.prob_noise:
            mask = self.build_mask(self.pct_noise)
            noise = np.random.normal(0, self.sigma_noise, int(self.n_genes * self.pct_noise))
            xg[mask] += noise
        
        if np.random.uniform(0, 1) < self.prob_swap:
            swap_pairs = np.random.randint(self.n_genes, size=(int(self.n_genes * self.pct_swap / 2), 2))
            xg[swap_pairs[:, 0]], xg[swap_pairs[:, 1]] = xg[swap_pairs[:, 1]], xg[swap_pairs[:, 0]]
            
        return xg
    