
from scipy.spatial import *
from sklearn.preprocessing import *

from sklearn.metrics import *
from scipy.spatial.distance import *
import scanpy as sc
import os
import pandas as pd
import numpy as np
import anndata
import matplotlib.pyplot as plt
import random
import networkx as nx
from scipy.spatial.distance import cdist
import math
import pickle


def create_graph(adata, degree = 4):
    """
    Converts spatial coordinates into graph using networkx library.
    
    param: adata - ST Slice 
    param: degree - number of edges per vertex

    return: 1) G - networkx graph
            2) node_dict - dictionary mapping nodes to spots
    """
    D = cdist(adata.obsm['spatial'], adata.obsm['spatial'])
    # Get column indexes of the degree+1 lowest values per row
    idx = np.argsort(D, 1)[:, 0:degree+1]
    # Remove first column since it results in self loops
    idx = idx[:, 1:]

    G = nx.Graph()
    for r in range(len(idx)):
        for c in idx[r]:
            G.add_edge(r, c)

    node_dict = dict(zip(range(adata.shape[0]), adata.obs.index))
    return G, node_dict

def generate_graph_from_labels(adata, labels_dict):
    """
    Creates and returns the graph and dictionary {node: cluster_label} for specified layer
    """
    
    g, node_to_spot = create_graph(adata)
    spot_to_cluster = labels_dict

    # remove any nodes that are not mapped to a cluster
    removed_nodes = []
    for node in node_to_spot.keys():
        if (node_to_spot[node] not in spot_to_cluster.keys()):
            removed_nodes.append(node)

    for node in removed_nodes:
        del node_to_spot[node]
        g.remove_node(node)
        
    labels = dict(zip(g.nodes(), [spot_to_cluster[node_to_spot[node]] for node in g.nodes()]))
    return g, labels

def spatial_coherence_score(graph, labels):
    g, l = graph, labels
    true_entropy = spatial_entropy(g, l)
    entropies = []
    for i in range(1000):
        new_l = list(l.values())
        random.shuffle(new_l)
        labels = dict(zip(l.keys(), new_l))
        entropies.append(spatial_entropy(g, labels))
        
    return (true_entropy - np.mean(entropies))/np.std(entropies)

def spatial_entropy(g, labels):
    """
    Calculates spatial entropy of graph  
    """
    # construct contiguity matrix C which counts pairs of cluster edges
    cluster_names = np.unique(list(labels.values()))
    C = pd.DataFrame(0,index=cluster_names, columns=cluster_names)

    for e in g.edges():
        C[labels[e[0]]][labels[e[1]]] += 1

    # calculate entropy from C
    C_sum = C.values.sum()
    H = 0
    for i in range(len(cluster_names)):
        for j in range(i, len(cluster_names)):
            if (i == j):
                z = C[cluster_names[i]][cluster_names[j]]
            else:
                z = C[cluster_names[i]][cluster_names[j]] + C[cluster_names[j]][cluster_names[i]]
            if z != 0:
                H += -(z/C_sum)*math.log(z/C_sum)
    return H

def fx_1NN(i,location_in):
    location_in = np.array(location_in)
    dist_array = distance_matrix(location_in[i,:][None,:],location_in)[0,:]
    dist_array[i] = np.inf
    return np.min(dist_array)


def fx_kNN(i,location_in,k,cluster_in):

    location_in = np.array(location_in)
    cluster_in = np.array(cluster_in)


    dist_array = distance_matrix(location_in[i,:][None,:],location_in)[0,:]
    dist_array[i] = np.inf
    ind = np.argsort(dist_array)[:k]
    cluster_use = np.array(cluster_in)
    # print(cluster_use[ind])
    # print(cluster_in[i])
    if np.sum(cluster_use[ind]!=cluster_in[i])>(k/2):
        return 1
    else:
        return 0

def _compute_CHAOS(clusterlabel, location):

    clusterlabel = np.array(clusterlabel)
    location = np.array(location)
    matched_location = StandardScaler().fit_transform(location)

    clusterlabel_unique = np.unique(clusterlabel)
    dist_val = np.zeros(len(clusterlabel_unique))
    count = 0
    for k in clusterlabel_unique:
        location_cluster = matched_location[clusterlabel==k,:]
        if len(location_cluster)<=2:
            continue
        n_location_cluster = len(location_cluster)
        results = [fx_1NN(i,location_cluster) for i in range(n_location_cluster)]
        dist_val[count] = np.sum(results)
        count = count + 1

    return np.sum(dist_val)/len(clusterlabel)

def _compute_PAS(clusterlabel,location):
        
    clusterlabel = np.array(clusterlabel)
    location = np.array(location)
    matched_location = location
    k=4
    # print(k)
    results = [fx_kNN(i,matched_location,k=k,cluster_in=clusterlabel) for i in range(matched_location.shape[0])]
    # print(results)
    # print(len(clusterlabel))
    # print(np.sum(results))
    return np.sum(results)/len(clusterlabel)

def compute_ASW(clusterlabel, location):
    d = squareform(pdist(location))
    return silhouette_score(X=d,labels=clusterlabel,metric='precomputed')