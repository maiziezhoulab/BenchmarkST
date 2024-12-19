from sklearn import metrics

ad = load_DLPFC(section_id="151673")

# calculate all metrics
ari = metrics.adjusted_rand_score(ad.obs['original_clusters'], ad.obs['pred'])
nmi = metrics.normalized_mutual_info_score(ad.obs['original_clusters'], ad.obs['pred'])
ami = metrics.adjusted_mutual_info_score(ad.obs['original_clusters'], ad.obs['pred'])
homogeneity = metrics.homogeneity_score(ad.obs['original_clusters'], ad.obs['pred'])

print(ari, nmi, ami, homogeneity)
