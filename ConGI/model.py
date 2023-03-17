import os
import numpy as np
from tqdm import tqdm
from collections import OrderedDict
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from st_loading_utils import load_DLPFC
import torch
from torch import nn
from torch.utils.tensorboard import SummaryWriter
from torchvision.models import resnet50, densenet121
from torchtoolbox.tools import mixup_data, cutmix_data, mixup_criterion

from utils import load_ST_file, calculate_adj_matrix, refine, build_her2st_data
from metrics import  eval_mclust_ari
from loss import NT_Xent
from loss import DCL
from loss import DCLW


def LinearBlock(input_dim, output_dim, p_drop):
    return nn.Sequential(
        nn.Linear(input_dim, output_dim),
        nn.BatchNorm1d(output_dim),
        nn.ELU(),
        nn.Dropout(p=p_drop),
    )


class SpaCLR(nn.Module):
    def __init__(self, gene_dims, image_dims, p_drop, n_pos, backbone='densenet', projection_dims=[64, 64]):
        super(SpaCLR, self).__init__()
        
        gene_dims.append(projection_dims[0])
        self.gene_encoder = nn.Sequential(OrderedDict([
            (f'gene_block{i+1}', LinearBlock(gene_dims[i], gene_dims[i+1], p_drop)) for i, _ in enumerate(gene_dims[:-1])
        ]))
        self.mse_loss = nn.MSELoss()
        self.gene_decoder = nn.Linear(projection_dims[0], gene_dims[0])

        if backbone == 'densenet':
            self.image_encoder = densenet121(pretrained=True)
            n_features = self.image_encoder.classifier.in_features
            self.image_encoder.classifier = nn.Identity()
        elif backbone == 'resnet':
            self.image_encoder = resnet50(pretrained=True)
            n_features = self.image_encoder.fc.in_features
            self.image_encoder.fc = nn.Identity()
        
        self.x_embedding = nn.Embedding(n_pos, n_features)
        self.y_embedding = nn.Embedding(n_pos, n_features)

        image_dims[0] = n_features
        image_dims.append(projection_dims[0])
        self.image_linear = nn.Sequential(OrderedDict([
            (f'image_block{i+1}', LinearBlock(image_dims[i], image_dims[i+1], p_drop)) for i, _ in enumerate(image_dims[:-1])
        ]))

        self.projector = nn.Sequential(
            nn.Linear(projection_dims[0], projection_dims[0]),
            nn.ReLU(),
            nn.Linear(projection_dims[0], projection_dims[1]),
        )

    def forward_image(self, xi, spatial):
        xi = self.image_encoder(xi)

        #x = self.x_embedding(spatial[:, 0])
        #y = self.y_embedding(spatial[:, 1])
        #xi = xi + x + y
        
        xi = self.image_linear(xi)
        hi = self.projector(xi)

        return xi, hi

    def forward_gene(self, xg):
        xg = self.gene_encoder(xg)
        hg = self.projector(xg)

        return xg, hg

    def forward(self, xg, xi, spatial):
        xg, hg = self.forward_gene(xg)
        xi, hi = self.forward_image(xi, spatial)

        return xg, xi, hg, hi
    
    def recon_loss(self, zg, xg):
        zg = self.gene_decoder(zg)
        return self.mse_loss(zg, xg)


class TrainerSpaCLR:
    def __init__(self, args, n_clusters, network, optimizer, log_dir, device='cuda'):
        self.n_clusters = n_clusters
        self.network = network
        self.optimizer = optimizer

        self.train_writer = SummaryWriter(log_dir+'_train')
        self.valid_writer = SummaryWriter(log_dir+'_valid')
        self.device = device

        self.args = args
        if args.dataset == "SpatialLIBD":
            adata = load_DLPFC(args.path, section_id=args.name)
            # df_meta = pd.read_csv(os.path.join(args.path, args.name, 'metadata.tsv'), sep='\t')
            # label = pd.Categorical(df_meta['layer_guess']).codes
            # adata = adata[label != -1]
            self.sample_id = adata.obs.index.tolist()
            self.adj_2d = calculate_adj_matrix(x=adata.obs["array_row"].tolist(), y=adata.obs["array_col"].tolist(), histology=False)
        elif args.dataset == "Her2st":
            adata, _ = build_her2st_data(args.path, args.name, args.img_size)
            label = adata.obs['label']
            adata = adata[label != -1]
            self.sample_id = adata.obs.index.tolist()
            self.adj_2d = calculate_adj_matrix(x=adata.obsm["spatial"][:, 0].tolist(), y=adata.obsm["spatial"][:, 1].tolist(), histology=False)
        
        self.w_g2g = args.w_g2g
        self.w_i2i = args.w_i2i
        self.w_recon = args.w_recon
    
    def eval_mclust_refined_ari(self, label, z):
        if z.shape[0] < 1000:
            num_nbs = 4
        else:
            num_nbs = 24
        ari, preds = eval_mclust_ari(label, z, self.n_clusters)
        refined_preds = refine(sample_id=self.sample_id, pred=preds, dis=self.adj_2d, num_nbs=num_nbs)
        ari = adjusted_rand_score(label, refined_preds)
        return ari

    def train(self, trainloader, epoch):
        with tqdm(total=len(trainloader)) as t:
            self.network.train()
            train_loss = 0
            train_cnt = 0

            for i, batch in enumerate(trainloader):
                t.set_description(f'Epoch {epoch} train')
                
                self.optimizer.zero_grad()
                xg, xg_u, xg_v, xi_u, xi_v, spatial, y, _ = batch
                xg = xg.to(self.device)
                xg_u = xg_u.to(self.device)
                xg_v = xg_v.to(self.device)
                xi_u = xi_u.to(self.device)
                xi_v = xi_v.to(self.device)
                spatial = spatial.to(self.device)
                
                criterion = NT_Xent(xg.shape[0])

                xg, xi_a, xi_b, lam = mixup_data(xg, xi_u)
                zg, hg = self.network.forward_gene(xg)
                zi_a, hi_a = self.network.forward_image(xi_a, spatial)
                zi_b, hi_b = self.network.forward_image(xi_b, spatial)
                g2i_loss = mixup_criterion(criterion, hg, hi_a, hi_b, lam)
                
                xg_u, xg_a, xg_b, lam = mixup_data(xg_u, xg_v)
                zg_u, hg_u = self.network.forward_gene(xg_u)
                zg_a, hg_a = self.network.forward_gene(xg_a)
                zg_b, hg_b = self.network.forward_gene(xg_b)
                g2g_loss = mixup_criterion(criterion, hg_u, hg_a, hg_b, lam) * self.w_g2g

                
                zi_c, hi_c = self.network.forward_image(xi_v, spatial)
                i2i_loss = criterion(hi_a, hi_c) * self.w_i2i
                
                recon_loss = self.network.recon_loss(zg, xg) * self.w_recon
                
                loss = g2i_loss + g2g_loss + i2i_loss + recon_loss
                loss.backward()
                self.optimizer.step()

                train_cnt += 1
                train_loss += loss.item()

                t.set_postfix(loss=f'{(train_loss/train_cnt):.3f}', 
                              g2i_loss=f'{g2i_loss.item():.3f}', 
                              g2g_loss=f'{g2g_loss.item():.3f}',
                              i2i_loss=f'{i2i_loss.item():.3f}', 
                              recon_loss=f'{recon_loss.item():.3f}')
                t.update(1)

            self.train_writer.add_scalar('loss', (train_loss/train_cnt), epoch)
            self.train_writer.flush()

    def valid(self, validloader, epoch=0):
        Xg = []
        Xi = []
        Y = []
        with torch.no_grad():
            with tqdm(total=len(validloader)) as t:
                self.network.eval()
                valid_loss = 0
                valid_cnt = 0

                for i, batch in enumerate(validloader):

                    xg, xi, spatial, y, _ = batch
                    xg = xg.to(self.device)
                    xi = xi.to(self.device)
                    spatial = spatial.to(self.device)

                    xg, xi, hg, hi = self.network(xg, xi, spatial)
                    criterion = NT_Xent(xg.shape[0])
                    loss = criterion(hg, hi)

                    valid_cnt += 1
                    valid_loss += loss.item()

                    Xg.append(xg.detach().cpu().numpy())
                    Xi.append(xi.detach().cpu().numpy())
                    Y.append(y)

                    t.set_postfix(loss=f'{(valid_loss/valid_cnt):.3f}')
                    t.update(1)
        
                Xg = np.vstack(Xg)
                Xi = np.vstack(Xi)
                Y = np.concatenate(Y, 0)

        return Xg, Xi, Y

    def fit(self, trainloader, epochs):
        self.network = self.network.to(self.device)

        for epoch in range(epochs):
            self.train(trainloader, epoch+1)

    def get_embeddings(self, validloader, save_name):

        xg, xi, _  = self.valid(validloader)
        np.save(os.path.join('preds', f'{save_name}_xg.npy'), xg)
        np.save(os.path.join('preds', f'{save_name}_xi.npy'), xi)


    def encode(self, batch):
        xg, xi, spatial, y, _ = batch
        xg = xg.to(self.device)
        xi = xi.to(self.device)
        spatial = spatial.to(self.device)
        xg, xi, hg, hi = self.network(xg, xi, spatial)
        return xg + 0.1 * xi
    
    def save_model(self, ckpt_path):
        torch.save(self.network.state_dict(), ckpt_path)

    def load_model(self, ckpt_path):
        self.network.load_state_dict(torch.load(ckpt_path))

