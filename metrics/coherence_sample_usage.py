# import dependencies correctly


"""SCS"""

adata = load_DLPFC(section_id="151673")
g, l = generate_graph_from_labels(adata, adata.obs['original_clusters'])
scs = spatial_coherence_score(g, l)
print(scs)

"""CHAOS/PAS/ASW"""
adata = load_DLPFC(section_id=sec)
X = adata.obsm['spatial']
pred_labels = adata.obs['original_clusters'].values
# g, l = generate_graph_from_labels(adata, adata.obs['original_clusters'])
# scs = spatial_coherence_score(g, l)
# print(clusterlabel)
# print(location)
chaos = CHAOS_score(X=X, pred_labels=pred_labels)
pas = PAS_score(X=X, pred_labels=pred_labels)
asw = ASW_score(X=X, pred_labels=pred_labels)

print(chaos, pas, asw)
