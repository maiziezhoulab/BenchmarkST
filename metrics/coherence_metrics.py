import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from tqdm import trange
from scipy.spatial import *
from sklearn.preprocessing import *

from sklearn.metrics import *
from scipy.spatial.distance import *

"""
SCS
PAS
CHAOS
ASW
"""

def _get_spatial_entropy(C, C_sum):
    H = 0
    for i in range(len(C)):
        for j in range(len(C)):
            z = C[i, j]
            if z != 0:
                H += -(z / C_sum) * np.log(z / C_sum)
    return H

# def spatial_entropy(g, labels):
def spatial_entropy(k_neighbors, labels, degree=4):
    """
    Calculates spatial entropy of graph
    """
    # construct contiguity matrix C which counts pairs of cluster edges
    # nx.set_node_attributes(g, labels, "labels")
    # cluster_nums = len(np.unique(list(labels.values())))
    # C = np.zeros((cluster_nums, cluster_nums))
    # for e in g.edges():
    #     C[labels[e[0]]][labels[e[1]]] += 1


    # S = np.repeat(adata.obs[annotation_key].values[:, None], degree, axis=1)
    S = np.broadcast_to(labels[:, None], (len(labels), degree))
    N = labels[k_neighbors]
    cluster_names = np.unique(labels)
    cluster_nums = len(cluster_names)
    C = np.zeros((cluster_nums, cluster_nums))
    for i in range(cluster_nums):
        for j in range(cluster_nums):
            # C[i, j] = np.sum(np.logical_and(N == i, S == j))
            C[i, j] = np.sum(np.logical_and(S == cluster_names[i], N == cluster_names[j]))
    # cluster_names = np.unique(list(labels.values()))
    # C = pd.DataFrame(0,index=cluster_names, columns=cluster_names)
    # C = np.zeros((len(cluster_names), len(cluster_names)))

    # calculate entropy from C
    # C_sum = C.values.sum()
    C_sum = C.sum()
    # print("C_sum", C_sum)
    # H = 0
    # # for i in range(len(cluster_names)):
    # #     for j in range(i, len(cluster_names)):
    # #         if (i == j):
    # #             z = C[cluster_names[i]][cluster_names[j]]
    # #         else:
    # #             z = C[cluster_names[i]][cluster_names[j]] + C[cluster_names[j]][cluster_names[i]]
    # #         if z != 0:
    # #             H += -(z/C_sum)*math.log(z/C_sum)
    # for i in range(len(C)):
    #     for j in range(len(C)):
    #         z = C[i, j]
    #         if z != 0:
    #             H += -(z / C_sum) * np.log(z / C_sum)
    # return H
    return _get_spatial_entropy(C, C_sum)

def spatial_coherence_score(adata, annotation_key, degree=4, rep_time=1000, seed=0):
    spatial_coords = adata.obsm['spatial']
    origin_labels = adata.obs[annotation_key].values
    # Use kneighbors_graph to get the adjacency matrix
    neigh = NearestNeighbors(n_neighbors=degree, metric='euclidean').fit(spatial_coords)
    # adjacency_matrix = neigh.kneighbors_graph(n_neighbors=degree, mode='connectivity').toarray().astype(np.int32)
    k_neighbors = neigh.kneighbors(n_neighbors=degree, return_distance=False)
    true_entropy = spatial_entropy(k_neighbors, origin_labels, degree=degree)
    entropies = []
    rng = np.random.default_rng(seed)
    shuffled_labels = origin_labels.copy()
    # for _ in trange(1000):
    for _ in trange(rep_time):
        rng.shuffle(shuffled_labels)
        entropies.append(spatial_entropy(k_neighbors, shuffled_labels, degree=degree))

    return (true_entropy - np.mean(entropies)) / np.std(entropies), true_entropy, entropies


def CHAOS_score(X, pred_labels):
    """
    Calculate the CHAOS score for a given set of spatial coordinates and predicted labels.

    param: X - spatial coordinates
    param: pred_labels - predicted labels

    return: CHAOS score
    """
    # Standardize the spatial coordinates
    X = StandardScaler().fit_transform(X)

    # Get the unique cluster labels
    cluster_labels = np.unique(pred_labels)

    # Initialize the distance value and count
    dist_val = 0.
    count = 0

    # Iterate through each cluster
    for k in cluster_labels:
        # Get the spatial coordinates for the current cluster
        cluster_coords = X[pred_labels == k, :]

        # Check if there are at least 2 spatial coordinates in the cluster
        if len(cluster_coords) <= 2:
            continue
        else:
            count += len(cluster_coords)

        # Calculate the distance to the nearest neighbor for each spatial coordinate in the cluster
        nbrs = NearestNeighbors(n_neighbors=1).fit(cluster_coords)
        distances, _ = nbrs.kneighbors()

        # Sum the distances
        dist_val = dist_val + np.sum(distances)

    # Calculate the CHAOS score
    return dist_val / count


def PAS_score(X, pred_labels, k=6):
    """
    Calculate the PAS score for a given set of spatial coordinates and predicted labels.

    param: X - spatial coordinates
    param: pred_labels - predicted labels
    param: k - number of nearest neighbors to consider

    return: PAS score
    """
    # Use NearestNeighbors to find the nearest neighbors
    nbrs = NearestNeighbors(n_neighbors=k).fit(X)
    indices = nbrs.kneighbors(return_distance=False)
    # print("type(indices)", type(indices))
    # print("type(pred_labels)", type(pred_labels))
    # print("indices.shape", indices.shape)
    # print("pred_labels.shape", pred_labels.shape)
    # Calculate the PAS score
    return ((pred_labels.reshape(-1, 1) != pred_labels[indices]).sum(1) > k / 2).mean()


def ASW_score(X, pred_labels):
    d = squareform(pdist(X))
    return silhouette_score(X=d,labels=pred_labels,metric='precomputed')
