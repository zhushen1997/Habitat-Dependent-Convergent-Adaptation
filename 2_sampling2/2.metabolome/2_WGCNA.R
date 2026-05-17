###======================================================================================================================================###
# step 2: WGCNA
###======================================================================================================================================###

library(flashClust)
library(WGCNA)
enableWGCNAThreads()

# load data
metabolome_matrix <- read.table('1_metabolome_matrix/metabolome_LOESS_PQN_missForest_matrix_noQC.tsv', row.names = 1, header = TRUE, sep = "\t") %>% t()

# select power by soft threshold
powers <- c(seq(1, 10, by=1), seq(12, 40, by=2))

sft <- pickSoftThreshold(
  metabolome_matrix, 
  powerVector = powers, 
  networkType = "signed", 
  verbose = 5)

# Scatter plot of fitted index and power value
par(mfrow = c(1, 2))

plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers,
     cex = 0.85,
     col = "red");
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], 
     sft$fitIndices[,5], 
     labels=powers, 
     cex=0.85,
     col="red")

# Select the appropriate power value
power = sft$powerEstimate

# TOM matrix
adjacency <- adjacency(metabolome_matrix, 
                       power = power, 
                       type = "signed",  
                       corFnc = "cor", 
                       corOptions = "use = 'p', method = 'spearman'")

# Calculate the similarity matrix
tom_sim <- TOMsimilarity(adjacency, TOMType = "signed")
rownames(tom_sim) <- rownames(adjacency)
colnames(tom_sim) <- colnames(adjacency)

# Calculate the dissimilarity matrix
tom_dis  <- 1 - tom_sim

# Hierarchical clustering tree
geneTree <- flashClust(as.dist(tom_dis), 
                       method = "complete")

plot(geneTree, 
     xlab = '', 
     sub = '', 
     main = 'metabolome Dendrogram', 
     labels = FALSE, 
     hang = 0.04)

# Use the dynamic shear tree mining module
minModuleSize <- 10
dynamicMods <- cutreeDynamic(dendro = geneTree, 
                             distM = tom_dis, 
                             method = "hybrid", 
                             deepSplit = 4, 
                             pamRespectsDendro = FALSE, 
                             minClusterSize = minModuleSize)

# Set the color of the module
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, 
                    dynamicColors, 
                    'Dynamic Tree Cut',
                    dendroLabels = FALSE, 
                    addGuide = TRUE, 
                    hang = 0.03, 
                    guideHang = 0.05,
                    main = 'metabolome dendrogram')


# Draw a heatmap of the correlation matrix
plot_sim <- -(1-tom_sim)
TOMplot(plot_sim, geneTree, dynamicColors, main = 'Network heatmap plot')

# Calculate the expression matrix of the original module
raw_MEList <- moduleEigengenes(metabolome_matrix, colors = dynamicColors)
raw_MEs <- as.data.frame(t(raw_MEList$eigengenes))
raw_MEs_table <- tibble::rownames_to_column(raw_MEs, "module")

# Metabolite mapping to Original Module
raw_metabolome_module <- data.frame(metabolome_name = colnames(metabolome_matrix), module = dynamicColors, stringsAsFactors = FALSE)

# Characterizing the similarity between modules
ME_cor <- cor(raw_MEs)
ME_diss <- 1-ME_cor
METree <- flashClust(as.dist(ME_diss), method = "complete")

# Draw a clustering tree
plot(METree, main = 'Clustering of module eigengenes', xlab = '', sub = '')

# Draw a clustering heatmap
plotEigengeneNetworks(raw_MEList$eigengenes, "Eigengene adjacency heatmap", marDendro = c(3,3,2,4), marHeatmap = c(3,4,2,2), plotDendrograms = TRUE, xLabelsAngle = 90)

# Merge similar modules
merge_module <- mergeCloseModules(metabolome_matrix, dynamicColors, corFnc = "cor", corOptions = list (use = 'p', method = 'spearman'), cutHeight = 0.1, verbose = 3)
mergedColors <- merge_module$colors
table(mergedColors)

# Draw the clustering tree and the original modules as well as the grouped merged modules
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged Tree Cut"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# Gene and Merge Module Correspondence Table
merge_metabolome_module <- data.frame(metabolome_name = colnames(metabolome_matrix), module = mergedColors, stringsAsFactors = FALSE)

# Calculate the expression quantity of the consolidation module
merge_MEList = moduleEigengenes(metabolome_matrix, colors = mergedColors)
merge_MEs = as.data.frame(t(merge_MEList$eigengenes))
merge_MEs_table <- tibble::rownames_to_column(merge_MEs, "module")

# raw_MEs & raw_metabolome_module are expression matrix of the original module and the mapping of module genes
write.table(raw_metabolome_module, "2_WGCNA/raw_module/raw_module_gene.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(raw_MEs_table, "2_WGCNA/raw_module/raw_module_Eigengenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# merge_MEs_table & merge_metabolome_module are merge the expression matrix of the module and the module gene mapping
write.table(merge_metabolome_module, "2_WGCNA/merge_module/merge_module_gene.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(merge_MEs_table, "2_WGCNA/merge_module/merge_module_Eigengenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
