library(PCAtools)

vsd_PCA <- pca(assay(vsd), metadata = colData(vsd), removeVar = 0.1)

screeplot(vsd_PCA, axisLabSize = 18, titleLabSize = 22)
biplot(vsd_PCA, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)





zscore_matrix <- t(as.matrix(merged_RF[,-c("rn","Diagnose","Cohort")]))
colnames(zscore_matrix) <- merged_RF$rn
rownames(zscore_matrix) <- colnames(merged_RF[,-c("rn","Diagnose","Cohort")])  


zs_PCA <- pca(zscore_matrix, metadata = colData(vsd), removeVar = 0.1)

screeplot(zs_PCA, axisLabSize = 18, titleLabSize = 22)
biplot(zs_PCA, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
