import numpy as np
from torch_geometric.nn.dense.linear import Linear
import torch
import torch.nn as nn
import torch.backends.cudnn as cudnn
cudnn.deterministic = True
cudnn.benchmark = True
import torch.nn.functional as F
from .gat_conv import GATConv


class GAAE(torch.nn.Module):
    def __init__(self, hidden_dims):
        super(GAAE, self).__init__()

        [in_dim, num_hidden, out_dim] = hidden_dims
        self.conv1 = GATConv(in_dim, num_hidden, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)
        self.conv2 = GATConv(num_hidden, out_dim, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)
        self.conv3 = GATConv(out_dim, num_hidden, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)
        self.conv4 = GATConv(num_hidden, in_dim, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)

    def forward(self, features, edge_index):

        h1 = F.elu(self.conv1(features, edge_index))
        h2 = self.conv2(h1, edge_index, attention=False)
        self.conv3.lin_src.data = self.conv2.lin_src.transpose(0, 1)
        self.conv3.lin_dst.data = self.conv2.lin_dst.transpose(0, 1)
        self.conv4.lin_src.data = self.conv1.lin_src.transpose(0, 1)
        self.conv4.lin_dst.data = self.conv1.lin_dst.transpose(0, 1)
        h3 = F.elu(self.conv3(h2, edge_index, attention=True,
                              tied_attention=self.conv1.attentions))
        h4 = self.conv4(h3, edge_index, attention=False)

        return h2, h4  # F.log_softmax(x, dim=-1)


class GAAE_mod1(torch.nn.Module):
    def __init__(self, in_dim, deg_dim, gae_hidden_dims, mlp_dims):
        super(GAAE_mod1, self).__init__()

        num_hidden, out_dim = gae_hidden_dims
        pred_out = mlp_dims
        self.linear_enc = Linear(in_dim, deg_dim)
        self.conv1 = GATConv(deg_dim, num_hidden, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)
        self.conv2 = GATConv(num_hidden, out_dim, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)
        self.conv3 = GATConv(out_dim, num_hidden, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)
        self.conv4 = GATConv(num_hidden, deg_dim, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)
        self.linear_dec = Linear(deg_dim, in_dim)
        self.linear_pred = Linear(out_dim, pred_out)

    def forward(self, features, edge_index):
        de = self.linear_enc(features)
        h1 = F.elu(self.conv1(de, edge_index))
        h2 = self.conv2(h1, edge_index, attention=False)

        pred = self.linear_pred(h2)

        self.conv3.lin_src.data = self.conv2.lin_src.transpose(0, 1)
        self.conv3.lin_dst.data = self.conv2.lin_dst.transpose(0, 1)
        self.conv4.lin_src.data = self.conv1.lin_src.transpose(0, 1)
        self.conv4.lin_dst.data = self.conv1.lin_dst.transpose(0, 1)
        h3 = F.elu(self.conv3(h2, edge_index, attention=True, tied_attention=self.conv1.attentions))
        h4 = self.conv4(h3, edge_index, attention=False)

        out = self.linear_dec(F.elu(h4))

        return h2, out, F.log_softmax(pred, dim=-1)


class GAAE_mod2(torch.nn.Module):
    def __init__(self, hidden_dims, mlp_dims):
        super(GAAE_mod2, self).__init__()

        [in_dim, num_hidden, out_dim] = hidden_dims
        pred_out = mlp_dims
        self.conv1 = GATConv(in_dim, num_hidden, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)
        self.conv2 = GATConv(num_hidden, out_dim, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)
        self.conv3 = GATConv(out_dim, num_hidden, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)
        self.conv4 = GATConv(num_hidden, in_dim, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)
        self.linear_pred = Linear(in_dim, pred_out)

    def forward(self, features, edge_index):
        h1 = F.elu(self.conv1(features, edge_index))
        h2 = self.conv2(h1, edge_index, attention=False)
        self.conv3.lin_src.data = self.conv2.lin_src.transpose(0, 1)
        self.conv3.lin_dst.data = self.conv2.lin_dst.transpose(0, 1)
        self.conv4.lin_src.data = self.conv1.lin_src.transpose(0, 1)
        self.conv4.lin_dst.data = self.conv1.lin_dst.transpose(0, 1)
        h3 = F.elu(self.conv3(h2, edge_index, attention=True,
                              tied_attention=self.conv1.attentions))
        h4 = self.conv4(h3, edge_index, attention=False)

        pred = self.linear_pred(h4)

        return h2, h4, F.log_softmax(pred, dim=-1)  # F.log_softmax(x, dim=-1)


class GAAE_mod3(torch.nn.Module):
    def __init__(self, hidden_dims, mlp_dims):
        super(GAAE_mod3, self).__init__()

        [in_dim, num_hidden, out_dim] = hidden_dims
        pred_out = mlp_dims
        self.conv1 = GATConv(in_dim, num_hidden, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)
        self.conv2 = GATConv(num_hidden, out_dim, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)
        self.conv3 = GATConv(out_dim, num_hidden, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)
        self.conv4 = GATConv(num_hidden, in_dim, heads=1, concat=False,
                             dropout=0, add_self_loops=False, bias=False)
        self.linear_pred = Linear(out_dim, pred_out)
        self.gene_attention_layer = nn.Parameter(torch.ones(1, in_dim))

        # for i in range(mlp_dims):
        #     self.gene_attention_layer_list.append(torch.nn.parameter.Parameter(data=torch.ones(in_dim),
        #                                                                        requires_grad=True).cuda())
        self.final_coeff = torch.zeros(in_dim).cuda()

    def forward_init(self, features, edge_index):
        h1 = F.elu(self.conv1(features, edge_index))
        h2 = self.conv2(h1, edge_index, attention=False)
        self.conv3.lin_src.data = self.conv2.lin_src.transpose(0, 1)
        self.conv3.lin_dst.data = self.conv2.lin_dst.transpose(0, 1)
        self.conv4.lin_src.data = self.conv1.lin_src.transpose(0, 1)
        self.conv4.lin_dst.data = self.conv1.lin_dst.transpose(0, 1)
        h3 = F.elu(self.conv3(h2, edge_index, attention=True,
                              tied_attention=self.conv1.attentions))
        h4 = self.conv4(h3, edge_index, attention=False)

        # pred = self.linear_pred(h4)

        return h2, h4  # , F.log_softmax(pred, dim=-1)  # F.log_softmax(x, dim=-1)

    def forward_refine(self, features, edge_index):
        att = features*self.gene_attention_layer
        h1 = F.elu(self.conv1(att, edge_index))
        h2 = self.conv2(h1, edge_index, attention=False)
        self.conv3.lin_src.data = self.conv2.lin_src.transpose(0, 1)
        self.conv3.lin_dst.data = self.conv2.lin_dst.transpose(0, 1)
        self.conv4.lin_src.data = self.conv1.lin_src.transpose(0, 1)
        self.conv4.lin_dst.data = self.conv1.lin_dst.transpose(0, 1)
        h3 = F.elu(self.conv3(h2, edge_index, attention=True,
                              tied_attention=self.conv1.attentions))
        h4 = self.conv4(h3, edge_index, attention=False)

        pred = self.linear_pred(h2)

        return h2, h4, F.log_softmax(pred, dim=-1)  # F.log_softmax(x, dim=-1)

    def predict(self, features, edge_index):
        # aggregate all the attentions
        # for e in self.gene_attention_layer_list:
        #     self.final_coeff += e
        att = features * self.gene_attention_layer
        h1 = F.elu(self.conv1(att, edge_index))
        h2 = self.conv2(h1, edge_index, attention=False)
        self.conv3.lin_src.data = self.conv2.lin_src.transpose(0, 1)
        self.conv3.lin_dst.data = self.conv2.lin_dst.transpose(0, 1)
        self.conv4.lin_src.data = self.conv1.lin_src.transpose(0, 1)
        self.conv4.lin_dst.data = self.conv1.lin_dst.transpose(0, 1)
        h3 = F.elu(self.conv3(h2, edge_index, attention=True,
                              tied_attention=self.conv1.attentions))
        h4 = self.conv4(h3, edge_index, attention=False)

        pred = self.linear_pred(h2)

        return h2, h4, F.log_softmax(pred, dim=-1)

    def get_de_attention(self):
        return self.gene_attention_layer  # .to('cpu').detach().numpy()
