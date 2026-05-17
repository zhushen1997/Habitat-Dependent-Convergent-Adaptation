### For Figure 2A, Figure 3A, Figure 4A, Figure 6A, Figure 6B, Figure 6C, Figure S3 

rm(list=ls())
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(ape)
library(Biostrings)
library(picante)
library(GUniFrac)
library(sna)
library(igraph)
library(ggClusterNet)
library(ggrepel)
library(network)
library(MetaNet)
library(SpiecEasi)
library(rgexf)
library(patchwork)
library(extrafont)
library(ggpolypath)
library(dplyr)
library(vegan)
library(pheatmap)
library(ggside)
library(flashClust)
library(ggtree)
library(dendextend)
library(ape)
library(ggtern)
library(aplot)
library(ggnested)
library(phyloseq)
library(ALDEx2)
library(readxl)
library(ggsignif)
library(ggbeeswarm)
library(ggh4x)
library(ggpubr)
library(Ternary)
library(extrafont)
library(monochromeR)
library(paletteer)
library(reshape2)
library(ggvenn)
library(eoffice)
library(Cairo)


# load Arial font
font_import(pattern = "arial")
loadfonts()
windowsFonts()

# Set a color for each group
color1 <- paletteer_d("nationalparkcolors::Badlands")
color2 <- paletteer_d("PrettyCols::Bright")
color3 <- paletteer_d("ggthemes::Classic_10_Medium")
color4 <- paletteer_d("MoMAColors::Althoff")

Esin_color <- color1[1]
Bger_color <- color1[4]
Pame_color <- generate_palette(paletteer_d("fishualize::Acanthurus_sohal")[1], modification = "go_lighter",  n_colours = 5, view_palette = TRUE)[2]
Pful_color <- color1[3]


###########################
#### KEGG KO abundance ####
###########################

# load KO abundance table
metagenome_function_KEGG_data <- t(read.table("function/2_downstream_analysis/1_Total_KO_TPM_filt.tsv", header = TRUE, sep = "\t", row.names = 1))

# load sample grouping table
group <- read.table('function/2_downstream_analysis/2_group.tsv', sep = '\t', header =TRUE)

# Merge KO abundance table and the sample grouping table
metagenome_function_KEGG_data_grouped <- merge(metagenome_function_KEGG_data, group, by.x = "row.names", by.y = "ID")

# Calculate the group mean
metagenome_function_KEGG_mean_abundance <- 
  metagenome_function_KEGG_data_grouped %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% 
  column_to_rownames(var = "group") %>%
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "KO")

write.table(metagenome_function_KEGG_mean_abundance, "function/2_downstream_analysis/3_metagenome_function_KEGG_mean_abundance.tsv", quote = F, row.names = F)


### (PCA) Principal Components Analysis
metagenome_function_KEGG_pca <- prcomp(log2(metagenome_function_KEGG_data + 1), scale. = TRUE)

# # Extract the coordinates (principal component scores) from the PCA results
metagenome_function_KEGG_pca_coords <- metagenome_function_KEGG_pca$x
metagenome_function_KEGG_pca_result <- merge(metagenome_function_KEGG_pca_coords, group, by.x = "row.names", by.y = "ID")

# # Extract all principal components
metagenome_function_KEGG_pca_scores <- metagenome_function_KEGG_pca$x[, 1:3] 

# Calculate the distance matrix based on PCA scores
metagenome_function_KEGG_distance <- dist(metagenome_function_KEGG_pca_scores)
metagenome_function_KEGG_distance_matrix <- as.matrix(metagenome_function_KEGG_distance)
write.table(metagenome_function_KEGG_distance_matrix, "function/2_downstream_analysis/4_metagenome_function_KEGG_sample_distance.txt", sep = "\t", quote = FALSE)

# The "pca$sdev" records the eigenvalues of the main sorting axes in the PCA sorting results (dividing each eigenvalue by the total sum of eigenvalues gives the explanatory power of each axis)
metagenome_function_KEGG_pca_eig = sum(pmax(metagenome_function_KEGG_pca$sdev[1:3]), 0)
metagenome_function_KEGG_pca_eig_percent <- round(metagenome_function_KEGG_pca$sdev[1:3]/metagenome_function_KEGG_pca_eig*100, 3) 

metagenome_function_KEGG_pca_eig_percent

# Conduct a permutation multivariate (factorial) variance analysis (PERMANOVA/adonis2)
metagenome_function_KEGG_pca_permanova_result <- adonis2(metagenome_function_KEGG_distance_matrix ~ group, data = group , permutations = 999)
metagenome_function_KEGG_pca_dune_adonis <- paste("R2", " = ", round(metagenome_function_KEGG_pca_permanova_result$R2, 3), "\nP = ", round(metagenome_function_KEGG_pca_permanova_result$`Pr(>F)`, 3))
metagenome_function_KEGG_pca_dune_adonis


### PCA plot
# Calculate the center point of each group
metagenome_function_KEGG_pca_center_points <- 
  metagenome_function_KEGG_pca_result %>%
  group_by(group) %>%
  summarise(mean_wt = mean(PC1),mean_mpg = mean(PC2))

metagenome_function_KEGG_pca_result <- 
  metagenome_function_KEGG_pca_result %>% 
  left_join(metagenome_function_KEGG_pca_center_points, by = "group")

metagenome_function_KEGG_pca_result$group <- factor(metagenome_function_KEGG_pca_result$group, level = c("Esin", "Bger", "Pame", "Pful"))

# Basic scatter plot
metagenome_function_KEGG_pca_p0 <- 
  ggplot(metagenome_function_KEGG_pca_result, 
         aes(x = PC1, y = PC2, color = group)
  ) +
  geom_point(aes(color = group, shape = group), 
             size = 2,
             show.legend = T,
  ) +
  labs(x = paste("PC 1 (", round(metagenome_function_KEGG_pca_eig_percent[1], 2), "%)", sep = ""), 
       y = paste("PC 2 (", round(metagenome_function_KEGG_pca_eig_percent[2], 2), "%)", sep = ""), 
       tag = metagenome_function_KEGG_pca_dune_adonis
  ) +
  scale_color_manual(values = c("Bger" = Bger_color, "Esin" = Esin_color, "Pame" = Pame_color, "Pful" = Pful_color)) +
  scale_shape_manual(values = c("Bger" = 19, "Esin" = 19, "Pame" = 19, "Pful" = 19)) +
  scale_x_continuous(position = "bottom") +
  scale_y_continuous(position = "left") +
  theme(
    axis.title.x.bottom = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.25, vjust = 0.5, angle = 0, lineheight = 2, color = "black"),
    axis.title.y.left = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.25, vjust = 1, angle = 90, lineheight = 2, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.25, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(2, "pt"),
    legend.frame = element_blank(),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.text.position = "right",
    legend.title = element_blank(),
    legend.title.position = "top",
    plot.tag.location = "panel",
    plot.tag = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0, angle = 0, lineheight = 1, color = "black"), 
    plot.tag.position = c(0.45, 0.60)
  )

metagenome_function_KEGG_pca_p0

# Add confidence ellipse + center point connection
metagenome_function_KEGG_pca_p1 <- 
  metagenome_function_KEGG_pca_p0 + 
  stat_ellipse(
    data = metagenome_function_KEGG_pca_result,
    aes(fill = group),
    geom = "polygon", 
    type = "t",
    level = 0.95, 
    linetype = 2, 
    linewidth = 0.5, 
    alpha=0.3, 
    show.legend = F
  ) +
  scale_fill_manual(values = c("Bger" = Bger_color, "Esin" = Esin_color, "Pame" = Pame_color, "Pful" = Pful_color)) +
  geom_segment(
    data = metagenome_function_KEGG_pca_result,
    aes(x = PC1, y = PC2, xend = mean_wt, yend = mean_mpg, color = group), 
    linetype = "dashed", 
    linewidth = 1, 
    alpha = 0.5,
    show.legend = F)

metagenome_function_KEGG_pca_p1

# Use ggside package to add marginal box plots
library(ggside)
metagenome_function_KEGG_pca_p_final <- 
  metagenome_function_KEGG_pca_p1 +
  geom_xsideboxplot(aes(y = group, fill = group),
                    orientation = "y",
                    alpha = 0.3,
                    outliers = FALSE,
                    staplewidth = 0.6,
                    linewidth = 0.25,
                    show.legend = F
  ) +
  geom_xsidepoint(aes(y = group, color = group, fill = group), 
                  position = "jitter",
                  size = 1.5,
                  alpha = 0.5,
                  show.legend = F
  ) +
  geom_ysideboxplot(aes(x = group, fill = group), 
                    orientation = "x", 
                    alpha = 0.3, 
                    outliers = FALSE, 
                    staplewidth = 0.6,
                    linewidth = 0.25,
                    show.legend = F
  ) +
  geom_ysidepoint(aes(x = group, color = group, fill = group), 
                  position = "jitter",
                  size = 1.5, 
                  alpha = 0.5,
                  show.legend = F
  ) +
  scale_xsidey_discrete() +
  scale_ysidex_discrete() +
  theme(ggside.panel.scale.x = 0.32,
        ggside.panel.scale.y = 0.32,
        legend.position = c(0.88, 0.88),
  )                                                                                                                 

print(metagenome_function_KEGG_pca_p_final)

ggsave("function/2_downstream_analysis/5_metagenome_function_KEGG_PCA_PC1~PC2.svg", plot = metagenome_function_KEGG_pca_p_final, width = 5, height = 5, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/5_metagenome_function_KEGG_PCA_PC1~PC2.tiff", plot = metagenome_function_KEGG_pca_p_final, width = 5, height = 5, units = "cm", dpi = 300)


###############################
### Hierarchical clustering ###
###############################
metagenome_function_spearman <- as.dist(1 -cor(t(metagenome_function_KEGG_data), method = "spearman"))
metagenome_function_spearman_matrix <- as.matrix(metagenome_function_spearman)

write.table(metagenome_function_spearman_matrix, "function/2_downstream_analysis/6_metagenome_function_spearman_matrix.tsv", quote = F, sep = "\t")

metagenome_function_KEGG_hc <- hclust(
  metagenome_function_spearman,
  method = "average",
  members = NULL)

#  The hclust object Convert to a phylo object
library(ape)
metagenome_function_KEGG_hc_tree <- ape::as.phylo(metagenome_function_KEGG_hc)
plot(metagenome_function_KEGG_hc_tree)

# Export Newick format tree file
write.tree(metagenome_function_KEGG_hc_tree, file = "function/2_downstream_analysis/7_metagenome_function_KEGG_spearman_average_tree.newick")

# Convert the hclust object to a dendrogram object
metagenome_function_KEGG_dend <- as.dendrogram(metagenome_function_KEGG_hc)

# Customize the dendrogram
metagenome_function_KEGG_dend <- metagenome_function_KEGG_dend %>% set("labels_col", c(rep(Esin_color, 6), rep(Bger_color, 6 ), rep(Pame_color ,6), rep(Pful_color, 6)))       # Set the color of the labels
metagenome_function_KEGG_dend <- metagenome_function_KEGG_dend %>% set("branches_k_color",c(rep(Esin_color, 6), rep(Bger_color, 6 ), rep(Pame_color ,6), rep(Pful_color, 6)))  # Set branch color
metagenome_function_KEGG_dend <- metagenome_function_KEGG_dend %>% set("branches_lwd", 4)                                                                                      # Set the branch width
metagenome_function_KEGG_dend <- metagenome_function_KEGG_dend %>% set("branches_lty", 1)                                                                                      # Set the branch line type
metagenome_function_KEGG_dend <- metagenome_function_KEGG_dend %>% set("nodes_pch", 19)                                                                                        # Set the shape of the node point
metagenome_function_KEGG_dend <- metagenome_function_KEGG_dend %>% set("nodes_col", "black")                                                                                   # Set the color of the node point
metagenome_function_KEGG_dend <- metagenome_function_KEGG_dend %>% set("nodes_cex", 0)                                                                                         # Set the size of the node point
metagenome_function_KEGG_dend <- metagenome_function_KEGG_dend %>% set("labels_cex", 1.5)                                                                                      # Set the label size
metagenome_function_KEGG_dend <- metagenome_function_KEGG_dend %>% set("leaves_pch", c(rep(17,6), rep(18,6), rep(19,6), rep(15,6)))                                            # Set the shape of the label point
metagenome_function_KEGG_dend <- metagenome_function_KEGG_dend %>% set("leaves_col", c(rep(Esin_color, 6), rep(Bger_color, 6 ), rep(Pame_color ,6), rep(Pful_color, 6)))       # Set the color of the label point
metagenome_function_KEGG_dend <- metagenome_function_KEGG_dend %>% set("leaves_cex", 2)                                                                                        # Set the size of the label point

# plot dendrogram
dev.off()
pdf("function/2_downstream_analysis/8_metagenome_function_spearman_hclust_tree.pdf")
plot(metagenome_function_KEGG_dend, main = "KO Distance Hierarchical Dendrogram")
dev.off()


### merge tree 
library(ape)
library(ggtree)
library(patchwork)
library(phytools)

groupInfo <- split(group$ID, group$group)

p_tree <-
  ggtree(groupOTU(metagenome_function_KEGG_hc_tree, groupInfo), 
         aes(color = group, fill = group, shape = group),
         size = 0.5,
         show.legend = F
  )+ 
  geom_tippoint(size = 1, 
                show.legend = F
  ) + 
  geom_text(aes(x = max(x) + 0.02, y = y, label = label),
            hjust = 0, 
            vjust = 0.5,
            size = 2
  ) + 
  scale_shape_manual(values = c("Bger" = 21, "Esin" = 22, "Pame" = 23, "Pful" = 24)) +
  scale_fill_manual(values = c("Bger" = Bger_color, "Esin" = Esin_color, "Pame" = Pame_color, "Pful" = Pful_color)) +
  scale_color_manual(values = c("Bger" = Bger_color, "Esin" = Esin_color, "Pame" = Pame_color, "Pful" = Pful_color)) +
  scale_y_reverse() +
  scale_x_continuous(
    expand = c(0, 0.02)) +
  theme(axis.text.y = element_blank(),
        axis.line.x =  element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
        axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1),
        axis.ticks.x  =  element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
        axis.ticks.length.x = unit(2, "pt"),
        legend.position = "none"
  )

p_tree

ggsave("function/2_downstream_analysis/9_metagenome_function_Hierarchical_Dendrogram.svg", plot = p_tree, width = 5, height = 5, units = "cm")
# ggsave("function/2_downstream_analysis/8_metagenome_function_Hierarchical_Dendrogram.tiff", plot = p_tree, width = 5, height = 5, units = "cm", dpi = 300)




### Hierarchical heatmap
row_annotation <- group[,c(1,2,4)]
row.names(row_annotation) <- row_annotation[, 1]
colnames(row_annotation)[2] <- "Host_group"
colnames(row_annotation)[3] <- "Host_habitat"

col_annotation <- group[,c(1,2,4)]
row.names(col_annotation) <- col_annotation[, 1]
colnames(col_annotation)[2] <- "Host_group"
colnames(col_annotation)[3] <- "Host_habitat"

annot_colors <- list(Host_group = c(Esin = Esin_color, Bger = Bger_color, Pame = Pame_color, Pful = Pful_color),
                     Host_habitat = c(Esin = "#453947FF", other = paletteer::paletteer_d("nationalparkcolors::Acadia")[1]))

min <- min(1- metagenome_function_spearman_matrix)
max <- max(1- metagenome_function_spearman_matrix)

dev.off()

pdf("function/2_downstream_analysis/10_metagenome_function_KO_spearman_Hierarchical_heatmap.pdf",
    width  = 10,
    height = 10)

pheatmap(
  as.matrix(1- metagenome_function_spearman_matrix),
  color = colorRampPalette(c( "lightyellow", "lightblue", "black"))(100),
  ### Global parameter
  main = "Spearman Similarity Heatmap",
  scale = "none",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 15,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 15,
  legend = TRUE,
  legend_breaks	= c(min, max),
  legend_labels = c(min = round(min, 1), max = round(max, 1)),
  width	= 15,
  height = 15,
  # Cluster tree and gap parameter
  cluster_rows = metagenome_function_KEGG_hc,
  cluster_cols = metagenome_function_KEGG_hc,
  cutree_rows = 2,
  cutree_cols = 2,
  treeheight_row = 40,
  treeheight_col = 40,
  # Annotation parameter
  annotation_row = row_annotation[, c(2), drop = FALSE],
  annotation_col = row_annotation[, c(2), drop = FALSE],
  annotation_colors = annot_colors,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  drop_levels = TRUE
)

dev.off()


### total_KO_heatmap
KO_dataframe <- 
  read.table("function/2_downstream_analysis/1_Total_KO_TPM_filt.tsv", header = T, row.names = 1) %>%
  t() %>%
  scale() %>%
  t()

col_annotation <- group[,c(1,2,4)]
row.names(col_annotation) <- col_annotation[, 1]
colnames(col_annotation)[2] <- "Host_group"
colnames(col_annotation)[3] <- "Host_habitat"

annot_colors <- list(Host_group = c(Esin = Esin_color, Bger = Bger_color, Pame = Pame_color, Pful = Pful_color),
                     Host_habitat = c(Esin = "#453947FF", other = paletteer::paletteer_d("nationalparkcolors::Acadia")[1]))


pdf("function/2_downstream_analysis/11_metagenome_KO_heatmap.pdf",
    width  = 10,
    height = 10)

pheatmap(KO_dataframe,
         color = colorRampPalette(c( "#648C16FF", "white","#FF7200FF"))(30),
         # color = colorRampPalette(c( "lightblue", "white", "pink"))(100),
         border_color = NA,
         cluster_rows = T,
         cluster_cols = T,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_names_row = F,
         show_rownames = F,
         annotation_col = row_annotation[, c(2), drop = FALSE],
         annotation_colors = annot_colors)

dev.off()


#########################################
### dimensionality reduction analysis ###
#########################################

### module_heatmap
KO_module_dataframe <- 
  read.table("function/1_KEGG/diff_module/ssGSEA/Total_metagenome_KEGG_ssGSEA_module.tsv", header = T, row.names = 1) %>%
  t() %>%
  scale() %>%
  t()

col_annotation <- group[,c(1,2,4)]
row.names(col_annotation) <- col_annotation[, 1]
colnames(col_annotation)[2] <- "Host_group"
colnames(col_annotation)[3] <- "Host_habitat"

annot_colors <- list(Host_group = c(Esin = Esin_color, Bger = Bger_color, Pame = Pame_color, Pful = Pful_color),
                     Host_habitat = c(Esin = "#453947FF", other = paletteer::paletteer_d("nationalparkcolors::Acadia")[1]))

pdf("function/2_downstream_analysis/12_metagenome_KO_module_heatmap.pdf",
    width  = 10,
    height = 10)

pheatmap(KO_module_dataframe,
         color = colorRampPalette(c( "#648C16FF", "white","#FF7200FF"))(30),
         # color = colorRampPalette(c( "lightblue", "white", "pink"))(100),
         border_color = NA,
         cluster_rows = T,
         cluster_cols = T,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_names_row = F,
         show_rownames = F,
         annotation_col = row_annotation[, c(2), drop = FALSE],
         annotation_colors = annot_colors)

dev.off()


### KEGG module regulation
# enriched KEGG modules heatmap
metagenome_up_module <- read.table("function/1_KEGG/diff_module/ssGSEA/1_all_up_module.tsv", header = TRUE, row.names = 1 ,sep = "\t") 

metagenome_up_module_zscore <- metagenome_up_module %>% dplyr::select(-c(1:8)) %>% t() %>% scale() %>% t() %>% as.data.frame()

metagenome_up_module[9:32] <- metagenome_up_module_zscore

metagenome_up_module <- 
  metagenome_up_module %>%
  mutate(log2trans = metagenome_up_module$Directionality.x..log10.p.adj.) %>%
  dplyr::select(log2trans, everything())

metagenome_up_module <- 
  metagenome_up_module %>%
  mutate(module = rownames(metagenome_up_module)) %>%
  dplyr::select(module, everything())

up_module_mratrix <- 
  metagenome_up_module %>%
  dplyr::select(-c(1:10)) %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "module")

up_module_mratrix$module <- factor(up_module_mratrix$module, level = rev(up_module_mratrix$module))

# wide to long
up_module_mratrix_long <-
  up_module_mratrix %>%
  pivot_longer(cols = -module,
               names_to = "sample", 
               values_to = "Enrichment_score")

sample_order <- 
  up_module_mratrix_long %>%
  distinct(sample) %>%
  mutate(
    group = str_extract(sample, "^[A-Za-z]+"),
    number = as.numeric(str_extract(sample, "\\d+"))
  ) %>%
  arrange(
    match(group, rev(c("Esin", "Bger", "Pame", "Pful"))),
    number
  ) %>%
  pull(sample)

up_module_mratrix_long$sample <- factor(up_module_mratrix_long$sample, levels = sample_order)

min <- min(up_module_mratrix_long$Enrichment_score)
max <- max(up_module_mratrix_long$Enrichment_score)

# plot main heatmap
up_module_heatmap_p <- 
  ggplot(up_module_mratrix_long, 
         aes(x = module, y = sample, fill = Enrichment_score)
  ) +
  geom_tile(color = NA,
            linewidth = 0,
            lineend = "square",
            height = 1,
            show.legend = T
  ) +
  annotate("rect", ymin = 18.5, ymax = 24.5, xmin = 0.5, xmax = 105.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 12.5, ymax = 18.5, xmin = 105.5, xmax = 150.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 6.5, ymax = 12.5, xmin = 150.5, xmax = 199.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 0.5, ymax = 6.5, xmin = 199.5, xmax = 203.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF", 
                       limit = c(min, max),
                       breaks =c(-2, 0, 2)
  ) +
  scale_x_discrete(position = "top") +
  labs(x = "Enriched metagenome modules of four species cockroaches", 
       fill = "Enrich Z-score"
  ) +
  theme_void() +
  theme(
    axis.title.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.position = "none",
    legend.key.size = unit(10, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "left",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = -90, lineheight = 1, color = "black"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

up_module_heatmap_p

# plot sample grouping annotatin
group$group <- factor(group$group, levels = c("Esin", "Bger", "Pame", "Pful"))

up_module_group_p <-
  ggplot(group, 
         aes(x = 0, y = ID, fill = group, color = group), 
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.05,
            lineend = "square",
            show.legend = F
  ) +
  labs(fill = "Host group", 
       color = "Host group"
  ) +
  scale_fill_manual(values = c("Esin" = Esin_color, "Bger" = Bger_color, "Pame" = Pame_color, "Pful" = Pful_color)) +
  theme_void()

up_module_group_p

# plot module pvalue annotatin
up_module_module_p <-
  ggplot(metagenome_up_module, 
         aes(x = module, y = 0, fill = log2trans , color = log2trans), 
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.05,
            lineend = "square",
            show.legend = T
  ) +
  labs(fill = "[Dir] x -log10(P.adj)", 
       color = "[Dir] x -log10(P.adj)"
  ) +
  scale_fill_gradient(low = paletteer::paletteer_d("MetBrewer::Benedictus")[5], 
                      high = paletteer::paletteer_d("MetBrewer::Benedictus")[3],
                      breaks = c(4, 8, 12)
  ) +
  scale_color_gradient(low = paletteer::paletteer_d("MetBrewer::Benedictus")[5], 
                       high = paletteer::paletteer_d("MetBrewer::Benedictus")[3],
                       breaks = c(4, 8, 12)
  ) +
  theme_void() + 
  theme(
    legend.position = "none",
    legend.key.size = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "left",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = -90, lineheight = 2, color = "black"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

up_module_module_p

# plot module category annotatin
up_module_module_info <- read.table("function/1_KEGG/diff_module/ssGSEA/1_all_up_module_information.tsv", header = TRUE, sep = "\t") %>% as.data.frame()

table(up_module_module_info$category_B)

up_module_module_info$category_B

up_module_module_info <- 
  up_module_module_info %>%
  mutate(Module.category = case_when(
    up_module_module_info$category_B == "Amino acid metabolism" ~ "AA. MET.",
    up_module_module_info$category_B == "Biosynthesis of other secondary metabolites" ~ "SECO. SYN.",
    up_module_module_info$category_B == "Biosynthesis of terpenoids and polyketides" ~ "SECO. SYN.",
    up_module_module_info$category_B == "Energy metabolism" ~ "ENER. MET.",
    up_module_module_info$category_B == "Lipid metabolism" ~ "LIPI. MET.",
    up_module_module_info$category_B == "Metabolism of cofactors and vitamins" ~ "VITA. MET.",
    up_module_module_info$category_B == "Xenobiotics biodegradation" ~ "XENO. DEG.",
    up_module_module_info$category_B == "Carbohydrate metabolism" ~ "CARB. MET.",
    up_module_module_info$category_B == "Glycan metabolism" ~ "GLYC. MET.",
    up_module_module_info$category_B == "Drug resistance" ~ "DRUG. RES.",
    TRUE ~ "Other"
  ))

up_module_module_info$Module.category <- 
  factor(up_module_module_info$Module.category,
         levels = c("AA. MET.", "CARB. MET.", "ENER. MET.", "GLYC. MET.", "LIPI. MET.", 
                    "VITA. MET.", "SECO. SYN.", "XENO. DEG.", "DRUG. RES.", "Other"))

up_module_module_info_p <-
  ggplot(up_module_module_info, 
         aes(x = module, y = 0, fill = Module.category , color = Module.category), 
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.05,
            lineend = "square",
            show.legend = T
  ) +
  labs(fill = "Module category", 
       color = "Module category"
  ) +
  guides(fill = guide_legend(ncol = 5, nrow = 2)) + 
  scale_fill_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:9],"grey"), drop = FALSE) +
  scale_color_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:9],"grey"), drop = FALSE) +
  theme_void() +
  theme(
    legend.position = "none",
    legend.key.size = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.spacing.x = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 2, color = "black"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "grey60")
  )

up_module_module_info_p

# merge plots 
up_module <- up_module_heatmap_p %>% 
  aplot::insert_left(up_module_group_p, width = 0.02) %>% 
  aplot::insert_top(up_module_module_p, height = 0.05) %>%
  aplot::insert_bottom(up_module_module_info_p, height = 0.05)

up_module

# save image
ggsave("function/2_downstream_analysis/13_metagenome_function_KEGG_enrich_module.svg", plot = up_module, width = 10, height = 4.5, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/13_metagenome_function_KEGG_enrich_module.tiff", plot = up_module, width = 10, height = 4.5, units = "cm", dpi = 300)


# save all legend
up_module_heatmap_p_legend <- cowplot::get_legend(up_module_heatmap_p + theme(legend.position = "right"))
up_module_module_p_legend <- cowplot::get_legend(up_module_module_p + theme(legend.position = "right"))
up_module_module_info_p_legend <- cowplot::get_legend(up_module_module_info_p + theme(legend.position = "right"))
up_module_combined_legend <- cowplot::plot_grid(up_module_module_p_legend, up_module_heatmap_p_legend, up_module_module_info_p_legend, nrow = 3)

up_module_combined_legend

ggsave("function/2_downstream_analysis/14_metagenome_function_KEGG_enrich_module_legend.svg", plot = up_module_combined_legend, width = 10, height = 10, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/14_metagenome_function_KEGG_enrich_module_legend.tiff", plot = up_module_combined_legend, width = 10, height = 10, units = "cm", dpi = 300)


# depleted KEGG modules heatmap
metagenome_down_module <- read.table("function/1_KEGG/diff_module/ssGSEA/2_all_down_module.tsv", header = TRUE, row.names = 1 ,sep = "\t") 

metagenome_down_module_zscore <- metagenome_down_module %>% dplyr::select(-c(1:8)) %>% t() %>% scale() %>% t() %>% as.data.frame()

metagenome_down_module[9:32] <- metagenome_down_module_zscore

metagenome_down_module <- 
  metagenome_down_module %>%
  mutate(log2trans = metagenome_down_module$Directionality.x..log10.p.adj.) %>%
  dplyr::select(log2trans, everything())

metagenome_down_module <- 
  metagenome_down_module %>%
  mutate(module = rownames(metagenome_down_module)) %>%
  dplyr::select(module, everything())

down_module_mratrix <- 
  metagenome_down_module %>%
  dplyr::select(-c(1:10)) %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "module")

down_module_mratrix$module <- factor(down_module_mratrix$module, level = rev(down_module_mratrix$module))

# wide to long
down_module_mratrix_long <-
  down_module_mratrix %>%
  pivot_longer(cols = -module,
               names_to = "sample", 
               values_to = "Enrichment_score")

sample_order <- 
  down_module_mratrix_long %>%
  distinct(sample) %>%
  mutate(
    group = str_extract(sample, "^[A-Za-z]+"),
    number = as.numeric(str_extract(sample, "\\d+"))
  ) %>%
  arrange(
    match(group, c("Esin", "Bger", "Pame", "Pful")),
    number
  ) %>%
  pull(sample)

down_module_mratrix_long$sample <- 
  factor(
    down_module_mratrix_long$sample,
    levels = sample_order
  )

min <- min(down_module_mratrix_long$Enrichment_score)
max <- max(down_module_mratrix_long$Enrichment_score)

# plot main heatmap
down_module_heatmap_p <- 
  ggplot(down_module_mratrix_long, 
         aes(x = module, y = sample, fill = Enrichment_score)
  ) +
  geom_tile(height = 1,
            color = NA,
            linewidth = 0,
            lineend = "square",
            show.legend = T
  ) +
  annotate("rect", ymin = 0.5, ymax = 6.5, xmin = 0.5 , xmax = 58.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 6.5, ymax = 12.5, xmin = 58.5, xmax = 126.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 12.5, ymax = 18.5, xmin = 126.5, xmax = 162.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 18.5, ymax = 24.5, xmin = 162.5, xmax = 177.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF", 
                       limit = c(min, max),
                       breaks = c(-2,0 ,2)
  ) +
  scale_x_discrete(position = "bottom") + 
  labs(x = "Depleted KEGG modules of four species cockroaches", 
       fill = "Enrich Z-score"
  ) +
  theme_void() + 
  theme(
    axis.title.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.position = "none",
    legend.key.size = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(10, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title.position = "left",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = -90, lineheight = 1, color = "black"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

down_module_heatmap_p

# plot sample grouping annotatin
group$group <- factor(group$group, levels = c("Esin", "Bger", "Pame", "Pful"))

down_module_group_p <-
  ggplot(group, 
         aes(x = 0, y = ID, fill = group, color = group), 
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.05,
            lineend = "square",
            show.legend = F
  ) +
  labs(fill = "Host group", 
       color = "Host group"
  ) +
  scale_fill_manual(values = c("Esin" = Esin_color, "Bger" = Bger_color, "Pame" = Pame_color, "Pful" = Pful_color)) +
  theme_void()

up_module_group_p 

# plot module pvalue annotatin
down_module_module_p <-
  ggplot(metagenome_down_module, 
         aes(x = module, y = 0, fill = log2trans , color = log2trans), 
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.05,
            lineend = "square",
            show.legend = T
  ) +
  labs(fill = "[Dir] x -log10(P.adj)", 
       color = "[Dir] x -log10(P.adj)"
  ) +
  scale_fill_gradient(low = paletteer::paletteer_d("MetBrewer::Benedictus")[11], 
                      high = paletteer::paletteer_d("MetBrewer::Benedictus")[9],
                      breaks = c(-4, -8, -12)
  ) +
  scale_color_gradient(low = paletteer::paletteer_d("MetBrewer::Benedictus")[11], 
                       high = paletteer::paletteer_d("MetBrewer::Benedictus")[9],
                       breaks = c(-4, -8, -12)
  ) +
  theme_void() + 
  theme(
    legend.position = "none",
    legend.key.size = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"), 
    legend.key.width = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(10, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title.position = "left",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = -90, lineheight = 2, color = "black"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white"),
  )

down_module_module_p

# plot module category annotatin
down_module_module_info <- read.table("function/1_KEGG/diff_module/ssGSEA/2_all_down_module_information.tsv", header = TRUE, sep = "\t") %>% as.data.frame()

table(down_module_module_info$category_B)

down_module_module_info <- 
  down_module_module_info %>%
  mutate(Module.category = case_when(
    down_module_module_info$category_B == "Amino acid metabolism" ~ "AA. MET.",
    down_module_module_info$category_B == "Biosynthesis of other secondary metabolites" ~ "SECO. SYN.",
    down_module_module_info$category_B == "Biosynthesis of terpenoids and polyketides" ~ "SECO. SYN.",
    down_module_module_info$category_B == "Energy metabolism" ~ "ENER. MET.",
    down_module_module_info$category_B == "Lipid metabolism" ~ "LIPI. MET.",
    down_module_module_info$category_B == "Metabolism of cofactors and vitamins" ~ "VITA. MET.",
    down_module_module_info$category_B == "Xenobiotics biodegradation" ~ "XENO. DEG.",
    down_module_module_info$category_B == "Carbohydrate metabolism" ~ "CARB. MET.",
    down_module_module_info$category_B == "Glycan metabolism" ~ "GLYC. MET.",
    down_module_module_info$category_B == "Drug resistance" ~ "DRUG. RES.",
    TRUE ~ "Other"
  ))

down_module_module_info$Module.category <- 
  factor(down_module_module_info$Module.category,
         levels = c("AA. MET.", "CARB. MET.", "ENER. MET.", "GLYC. MET.", "LIPI. MET.", 
                    "VITA. MET.", "SECO. SYN.", "XENO. DEG.", "DRUG. RES.", "Other"))

down_module_module_info_p <-
  ggplot(down_module_module_info, 
         aes(x = module, y = 0, fill = Module.category , color = Module.category), 
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.05,
            lineend = "square",
            show.legend = T
  ) +
  labs(fill = "Module category", 
       color = "Module category"
  ) +
  scale_fill_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:9],"grey"), drop = FALSE) +
  scale_color_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:9],"grey"), drop = FALSE) +
  guides(fill = guide_legend(ncol = 5, nrow = 2)) + 
  theme_void() + 
  theme(
    legend.position = "none",
    legend.key.size = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"), 
    legend.key.width = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(10, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 2, color = "black"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white"),
  )

down_module_module_info_p

# merge plots 
down_module <- down_module_heatmap_p %>% 
  aplot::insert_left(down_module_group_p, width = 0.02) %>% 
  aplot::insert_bottom(down_module_module_p, height = 0.05) %>%
  aplot::insert_top(down_module_module_info_p, height = 0.05)

down_module

# save image
ggsave("function/2_downstream_analysis/15_metagenome_function_KEGG_deplete_module.svg", plot = down_module, width = 10, height = 4.5, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/15_metagenome_function_KEGG_deplete_module.tiff", plot = down_module, width = 10, height = 4.5, units = "cm", dpi = 300)


# save all legend
down_module_heatmap_p_legend <- cowplot::get_legend(down_module_heatmap_p + theme(legend.position = "right"))
down_module_module_p_legend <- cowplot::get_legend(down_module_module_p + theme(legend.position = "right"))
down_module_module_info_p_legend <- cowplot::get_legend(down_module_module_info_p + theme(legend.position = "right"))
down_module_combined_legend <- cowplot::plot_grid(down_module_heatmap_p_legend, down_module_module_p_legend, down_module_module_info_p_legend, nrow = 3)

down_module_combined_legend

ggsave("function/2_downstream_analysis/16_metagenome_function_KEGG_deplete_module_legend.svg", plot = down_module_combined_legend, width = 10, height = 10, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/16_metagenome_function_KEGG_deplete_module_legend.tiff", plot = down_module_combined_legend, width = 10, height = 10, units = "cm", dpi = 300)




### regulation module percent
library(WeightedTreemaps)

up_module_module_percent <- 
  metagenome_up_module %>% left_join(up_module_module_info, by = "module") %>% 
  group_by(regulation, Module.category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(regulation) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

up_module_module_percent$Module.category <- factor(up_module_module_percent$Module.category, 
                                                   levels = rev(c("AA. MET.", "CARB. MET.", "ENER. MET.", "GLYC. MET.", "LIPI. MET.", 
                                                                  "VITA. MET.", "SECO. SYN.", "XENO. DEG.", "DRUG. RES.", "Other")))

up_module_module_percent$regulation <- factor(up_module_module_percent$regulation,
                                              levels = c("Esin_up", "Bger_up", "Pame_up", "Pful_up"))

up_module_module_percent_bar_plot <-
  ggplot(up_module_module_percent) + 
  geom_bar(
    aes(x = regulation, y = percentage, fill = Module.category),
    stat = "identity",
    position = "stack",
    width = 0.7
  ) +
  labs(
    x = "Regulation",
    y = "Enrich Module Percentage (%)",
    fill = "Category"
  ) +
  scale_fill_manual(values = c("grey", paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[9:1]), drop = FALSE) +
  scale_color_manual(values = c("grey", paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[9:1]), drop = FALSE) +
  theme_minimal() +
  theme(
    axis.line.x = element_blank(), 
    axis.line.y = element_blank(),  
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 1, vjust = 0.5, angle = 0, lineheight = 1, color = "black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x = unit(0, "pt"),
    axis.ticks.length.y = unit(2, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25, linetype = 1), 
    legend.position = "right",
    legend.key.size = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"), 
    legend.key.width = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(10, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 2, color = "black"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white"),
  )

print(up_module_module_percent_bar_plot)

up_module_module_total <- 
  up_module_module_percent %>%
  group_by(regulation) %>%
  summarise(sum = sum(count), .groups = 'drop')

up_module_module_total_bar_plot <- 
  ggplot(up_module_module_total) + 
  geom_bar(aes(x = regulation, y = sum, fill = regulation),
           stat = "identity",
           width = 0.7) +
  geom_label(aes(label = sum, x = regulation, y = 50), size = 3) + 
  scale_fill_manual(values = c("Esin_up" = Esin_color, "Bger_up" = Bger_color, "Pame_up" = Pame_color, "Pful_up" = Pful_color)) + 
  theme_minimal() +
  theme(
    axis.line.x = element_blank(), 
    axis.line.y = element_blank(),  
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 1, vjust = 0.5, angle = 0, lineheight = 1, color = "black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x = unit(0, "pt"),
    axis.ticks.length.y = unit(2, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25, linetype = 1), 
    legend.position = "none",
    legend.key.size = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"), 
    legend.key.width = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(10, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 2, color = "black"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white"),
  )

print(up_module_module_total_bar_plot)


down_module_module_percent <- 
  metagenome_down_module %>% left_join(down_module_module_info, by = "module") %>% 
  group_by(regulation, Module.category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(regulation) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

down_module_module_percent$Module.category <- factor(down_module_module_percent$Module.category, 
                                                     levels = rev(c("AA. MET.", "CARB. MET.", "ENER. MET.", "GLYC. MET.", "LIPI. MET.", 
                                                                    "VITA. MET.", "SECO. SYN.", "XENO. DEG.", "DRUG. RES.", "Other")))

down_module_module_percent$regulation <- factor(down_module_module_percent$regulation,
                                                levels = c("Esin_down", "Bger_down", "Pame_down", "Pful_down"))

down_module_module_percent_bar_plot <-
  ggplot(down_module_module_percent) + 
  geom_bar(
    aes(x = regulation, y = percentage, fill = Module.category),
    stat = "identity",
    position = "stack",
    width = 0.7
  ) +
  labs(
    x = "Regulation",
    y = "Deplete Module Percentage (%)",
    fill = "Category"
  ) +
  scale_fill_manual(values = c("grey", paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[9:1]), drop = FALSE) +
  scale_color_manual(values = c("grey", paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[9:1]), drop = FALSE) +
  scale_y_reverse() +
  scale_x_discrete(position = "top") + 
  theme_minimal() +
  theme(
    axis.line.x = element_blank(), 
    axis.line.y = element_blank(),  
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 1, vjust = 0.5, angle = 0, lineheight = 1, color = "black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x = unit(0, "pt"),
    axis.ticks.length.y = unit(2, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25, linetype = 1), 
    legend.position = "none",
    legend.key.size = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"), 
    legend.key.width = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(10, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 2, color = "black"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white"),
  )

down_module_module_percent_bar_plot

down_module_module_total <- 
  down_module_module_percent %>%
  group_by(regulation) %>%
  summarise(sum = sum(count), .groups = 'drop')

down_module_module_total_bar_plot <- 
  ggplot(down_module_module_total) + 
  geom_bar(aes(x = regulation, y = sum, fill = regulation),
           stat = "identity",
           width = 0.7) +
  geom_label(aes(label = sum, x = regulation, y = 35), size = 3) + 
  scale_fill_manual(values = c("Esin_down" = Esin_color, "Bger_down" = Bger_color, "Pame_down" = Pame_color, "Pful_down" = Pful_color)) + 
  scale_y_reverse() + 
  scale_x_discrete(position = "top") + 
  theme_minimal() +
  theme(
    axis.line.x = element_blank(), 
    axis.line.y = element_blank(),  
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 1, vjust = 0.5, angle = 0, lineheight = 1, color = "black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x = unit(0, "pt"),
    axis.ticks.length.y = unit(2, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25, linetype = 1), 
    legend.position = "none",
    legend.key.size = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"), 
    legend.key.width = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(10, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 2, color = "black"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white"),
  )

down_module_module_total_bar_plot


module_percent_bar_plot <- up_module_module_total_bar_plot / up_module_module_percent_bar_plot / down_module_module_percent_bar_plot / down_module_module_total_bar_plot + plot_layout(heights = c(0.2, 1, 1, 0.2))

module_percent_bar_plot

ggsave("function/2_downstream_analysis/17_module_percent_bar_plot.svg", plot = module_percent_bar_plot, width = 6, height = 12, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/17_module_percent_bar_plot.tiff", plot = module_percent_bar_plot, width = 6, height = 12, units = "cm", dpi = 300)


### highlight top 10 enriched/depleted KEGG modules

# top 10 enriched KEGG module
top10_enriched_module <- read.table("function/1_KEGG/diff_module/ssGSEA/3_each_top10_up_module.tsv", header = TRUE, sep = "\t") 
top10_enriched_module$log2trans <- ifelse(top10_enriched_module$Directionality.x..log10.p.adj. > 6, 6, top10_enriched_module$Directionality.x..log10.p.adj.)
top10_enriched_module$mark <- ifelse(top10_enriched_module$Directionality.x..log10.p.adj. > 6, "over", "not")
top10_enriched_module <- top10_enriched_module %>% dplyr::select(log2trans, mark, everything())
top10_enriched_module$regulation <- factor(top10_enriched_module$regulation, levels = c("Esin_up", "Bger_up", "Pame_up", "Pful_up"))
top10_enriched_module$module <- factor(top10_enriched_module$module, levels = top10_enriched_module$module)
top10_enriched_module_info <- read.table("function/1_KEGG/diff_module/ssGSEA/3_each_top10_up_module_information.tsv", header = TRUE, sep = "\t") 

top10_enriched_module_info <- 
  top10_enriched_module_info %>%
  mutate(Module.category = case_when(
    top10_enriched_module_info$category_B == "Amino acid metabolism" ~ "AA. MET.",
    top10_enriched_module_info$category_B == "Biosynthesis of other secondary metabolites" ~ "SECO. SYN.",
    top10_enriched_module_info$category_B == "Biosynthesis of terpenoids and polyketides" ~ "SECO. SYN.",
    top10_enriched_module_info$category_B == "Energy metabolism" ~ "ENER. MET.",
    top10_enriched_module_info$category_B == "Lipid metabolism" ~ "LIPI. MET.",
    top10_enriched_module_info$category_B == "Metabolism of cofactors and vitamins" ~ "VITA. MET.",
    top10_enriched_module_info$category_B == "Xenobiotics biodegradation" ~ "XENO. DEG.",
    top10_enriched_module_info$category_B == "Carbohydrate metabolism" ~ "CARB. MET.",
    top10_enriched_module_info$category_B == "Glycan metabolism" ~ "GLYC. MET.",
    top10_enriched_module_info$category_B == "Drug resistance" ~ "DRUG. RES.",
    TRUE ~ "Other"
  ))

top10_enriched_module_info$Module.category <- 
  factor(top10_enriched_module_info$Module.category,
         levels = c("AA. MET.", "CARB. MET.", "ENER. MET.", "GLYC. MET.", "LIPI. MET.", 
                    "VITA. MET.", "SECO. SYN.", "XENO. DEG.", "DRUG. RES.", "Other"))

top10_enriched_module_info <- 
  top10_enriched_module_info %>%
  mutate(merge.name = paste0(top10_enriched_module_info$module, ":", top10_enriched_module_info$module_name)) %>%
  separate(merge.name, into = c("first_part", "rest"), sep = ",", remove = FALSE) %>%
  mutate(merge.name = first_part) %>%
  dplyr::select(-first_part, -rest)

top10_enriched_module_info$merge.name <- gsub("M00968:\\s*Pentose bisphosphate pathway \\(nucleoside degradation\\)", "M00968:Pentose bisphosphate pathway", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("ko00571:\\s*Lipoarabinomannan \\(LAM)\\ biosynthesis", "ko00571:Lipoarabinomannan biosynthesis", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("M00631:\\s*D-Galacturonate degradation \\(bacteria\\)", "M00631:D-Galacturonate degradation", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("ko01051:\\s*Biosynthesis of ansamycins", "ko01051:Ansamycins SYN", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("biosynthesis", "SYN.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("metabolism", "MET.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("degradation", "DEG.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("oxidoreductase", "REDOX.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("oxidase", "OXID.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("cycle", "CYC.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("resistance", "RES.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("reduction", "RED.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("pathway", "PWY.", top10_enriched_module_info$merge.name)

top10_enriched_module_info$module <- factor(top10_enriched_module_info$module, levels = top10_enriched_module$module)

# top10_enriched_module bar plot
top10_enriched_module_p0 <- 
  ggplot(top10_enriched_module, 
         aes(x = log2trans, y = module , fill = regulation)
  ) +
  geom_hline(yintercept = seq_along(top10_enriched_module$module), 
             color = "grey", 
             linetype = "dotted",
             lineend = "square",
             linewidth = 0.15
  ) +
  geom_point(aes(fill = regulation, color = regulation, shape = mark),
             size = 1.5,
             show.legend = T
  )+
  scale_fill_manual(values = c("Bger_up" = Bger_color, "Esin_up" = Esin_color, "Pame_up" = Pame_color, "Pful_up" = Pful_color)) +
  scale_color_manual(values = c("Bger_up" = Bger_color, "Esin_up" = Esin_color, "Pame_up" = Pame_color, "Pful_up" = Pful_color)) +
  scale_shape_manual(values = c("over" = 22, "not" = 21)) +
  labs(x = "[Directionality] x -log10(P.adj)") +
  scale_y_discrete(position = "right",
                   labels = top10_enriched_module_info$merge.name
  )  +
  scale_x_continuous(breaks = c(0, -log10(0.05), 6), 
                     labels = c(0, "*", 6)
  ) + 
  coord_cartesian(xlim = c(0, 6), 
                  expand = TRUE
  ) +
  geom_vline(xintercept = -log10(0.05), 
             color = "black", 
             linetype = "dashed", 
             linewidth = 0.25
  ) +
  theme_void() + 
  theme(
    axis.title.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 1, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, 
                               color = c(rep(Pful_color, 4), rep(Pame_color,10), rep(Bger_color,10), rep(Esin_color, 10))), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x = unit(2, "pt"),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.25, linetype = "solid"),
    legend.key.size = unit(10, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.spacing.x = unit(5, "pt"),
    legend.spacing.y = unit(5, "pt"),
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.text.position = "right",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
  )

top10_enriched_module_p0

# top10_enriched_module module category annotation
top10_enriched_module_p1 <- 
  ggplot(top10_enriched_module_info, 
         aes(x = 0, y = module, fill = Module.category, color = Module.category)
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.25,
            lineend = "square",
            show.legend = T
  ) +  
  scale_fill_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:9],"grey"), drop = FALSE) +
  scale_color_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:9],"grey"), drop = FALSE) +
  labs(fill = "Module categary") + 
  guides(fill = guide_legend(ncol =1)) +
  theme_void() + 
  theme(
    legend.key.size = unit(10, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.direction = "vertical",
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.text.position = "right",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
  )

top10_enriched_module_p1

# plot heatmap
sample_levels <- c(paste("Esin_",c(1:6), sep = ""), paste("Bger_",c(1:6), sep = ""), paste("Pame_",c(1:6), sep = ""), paste("Pful_",c(1:6), sep = ""))

up_module_z_scores <- top10_enriched_module[, 12:35] %>% t() %>% scale() %>% t()

top10_enriched_module[, 12:35] = up_module_z_scores

# wide to long
top10_enriched_module_long <- top10_enriched_module %>%
  pivot_longer(cols = -c(1:11), names_to = "sample", values_to = "Enrich_score")

top10_enriched_module_long$sample <- factor(top10_enriched_module_long$sample, levels = sample_levels)

min <- min(top10_enriched_module_long$Enrich_score)
max <- max(top10_enriched_module_long$Enrich_score)

top10_enriched_module_p2 <- 
  ggplot(top10_enriched_module_long, 
         aes(x = sample, y = module, fill = Enrich_score)
  ) +
  geom_tile(height = 1,
            color = NA,
            linewidth = 0,
            lineend = "square",
            show.legend = T
  ) +
  annotate("rect", xmin = 0.5, xmax = 6.5, ymin = 24.5, ymax = 34.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 6.5, xmax = 12.5, ymin = 14.5, ymax = 24.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 12.5, xmax = 18.5, ymin = 4.5, ymax = 14.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 18.5, xmax = 24.5, ymin = 0.5, ymax = 4.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF",
                       limit = c(min, max),
                       breaks = c(-2, 0 ,2)
  ) +
  labs(fill = "Enrich Z-score") +
  theme_void() +
  theme(
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 1, angle = 0, lineheight = 1, color = "black"),
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0, angle = 0, lineheight = 1, color = "black"),
    legend.key.size = unit(10, "pt"), 
    legend.key.width = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.spacing.x = unit(5, "pt"),
    legend.spacing.y = unit(5, "pt"),
    legend.direction = "horizontal",
    legend.title.position = "top",
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

top10_enriched_module_p2

library(aplot)

top10_enriched_module_p <- 
  top10_enriched_module_p0 %>%
  insert_left(top10_enriched_module_p1, width = 0.2) %>%
  insert_left(top10_enriched_module_p2, width = 1)

top10_enriched_module_p

ggsave("function/2_downstream_analysis/18_top10_enriched_module.svg", plot = top10_enriched_module_p, width = 14, height = 9.5, units = "cm", dpi = 300)
#ggsave("function/2_downstream_analysis/18_top10_enriched_module.tiff", plot = top10_enriched_module_p, width = 14, height = 9.5, units = "cm", dpi = 300)


# top 10 depleted KEGG modules
top10_depleted_module <- read.table("function/1_KEGG/diff_module/ssGSEA/4_each_top10_down_module.tsv", header = TRUE, sep = "\t") 
top10_depleted_module$log2trans <- ifelse(top10_depleted_module$Directionality.x..log10.p.adj. < -6, -6, top10_depleted_module$Directionality.x..log10.p.adj.)
top10_depleted_module$mark <- ifelse(top10_depleted_module$Directionality.x..log10.p.adj. < -6, "over", "not")
top10_depleted_module <- top10_depleted_module %>% dplyr::select(log2trans, mark, everything())
top10_enriched_module$regulation <- factor(top10_enriched_module$regulation, levels = c("Esin_down", "Bger_down", "Pame_down", "Pful_down"))
top10_depleted_module$module <- factor(top10_depleted_module$module, levels = top10_depleted_module$module)
top10_depleted_module_info <- read.table("function/1_KEGG/diff_module/ssGSEA/4_each_top10_down_module_information.tsv", header = TRUE, sep = "\t") 

top10_depleted_module_info <- 
  top10_depleted_module_info %>%
  mutate(Module.category = case_when(
    top10_depleted_module_info$category_B == "Amino acid metabolism" ~ "AA. MET.",
    top10_depleted_module_info$category_B == "Biosynthesis of other secondary metabolites" ~ "SECO. SYN.",
    top10_depleted_module_info$category_B == "Biosynthesis of terpenoids and polyketides" ~ "SECO. SYN.",
    top10_depleted_module_info$category_B == "Energy metabolism" ~ "ENER. MET.",
    top10_depleted_module_info$category_B == "Lipid metabolism" ~ "LIPI. MET.",
    top10_depleted_module_info$category_B == "Metabolism of cofactors and vitamins" ~ "VITA. MET.",
    top10_depleted_module_info$category_B == "Xenobiotics biodegradation" ~ "XENO. DEG.",
    top10_depleted_module_info$category_B == "Carbohydrate metabolism" ~ "CARB. MET.",
    top10_depleted_module_info$category_B == "Glycan metabolism" ~ "GLYC. MET.",
    top10_depleted_module_info$category_B == "Drug resistance" ~ "DRUG. RES.",
    TRUE ~ "Other"
  ))

top10_depleted_module_info$Module.category <- 
  factor(top10_depleted_module_info$Module.category,
         levels = c("AA. MET.", "CARB. MET.", "ENER. MET.", "GLYC. MET.", "LIPI. MET.", 
                    "VITA. MET.", "SECO. SYN.", "XENO. DEG.", "DRUG. RES.", "Other"))

top10_depleted_module_info <- 
  top10_depleted_module_info %>%
  mutate(merge.name = paste0(top10_depleted_module_info$module, ":", top10_depleted_module_info$module_name)) %>%
  separate(merge.name, into = c("first_part", "rest"), sep = ",", remove = FALSE) %>%
  mutate(merge.name = first_part) %>%
  dplyr::select(-first_part, -rest)

top10_depleted_module_info$merge.name <- gsub("M00168:\\s*CAM \\(Crassulacean acid metabolism\\)", "M00168:Crassulacean acid metabolism", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("ko00980:\\s*Metabolism of xenobiotics by cytochrome P450", "ko00980:Xenobiotics MET by P450", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("M00415:\\s*Fatty acid elongation in endoplasmic reticulum", "M00415:Fatty acid elongation", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("M00856:\\s*Salmonella enterica pathogenicity signature", "M00856:Pathogenicity signature", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("M00375:\\s*Hydroxypropionate-hydroxybutylate", "M00375:HP-HB", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("ko00524:\\s*Neomycin", "ko00524:Neo/Kana/Genta mycin SYN", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("biosynthesis", "SYN", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("metabolism", "MET", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("degradation", "DEG", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("oxidoreductase", "REDOX", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("oxidase", "OXID", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("cycle", "CYC", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("resistance", "RES", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("reduction", "RED", top10_depleted_module_info$merge.name)

top10_depleted_module_info$module <- factor(top10_depleted_module_info$module, levels = top10_depleted_module$module)

# top10_depleted_module bar plot
top10_depleted_module_p0 <- 
  ggplot(top10_depleted_module, 
         aes(x = log2trans, y = module, fill = regulation)
  ) +
  geom_hline(yintercept = seq_along(top10_depleted_module$module), 
             color = "grey", 
             linetype = "dotted", 
             linewidth = 0.15
  ) +
  geom_point(aes(fill = regulation, color = regulation, shape = mark),
             size = 1.5,
             show.legend = T
  )+
  scale_fill_manual(values = c("Bger_down" = Bger_color, "Esin_down" = Esin_color, "Pame_down" = Pame_color, "Pful_down" = Pful_color)) +
  scale_color_manual(values = c("Bger_down" = Bger_color, "Esin_down" = Esin_color, "Pame_down" = Pame_color, "Pful_down" = Pful_color)) +
  scale_shape_manual(values = c("over" = 22, "not" = 21)) +
  labs(x = "[Directionality] x -log10(P.adj)") +
  scale_y_discrete(position = "left",
                   labels = top10_depleted_module_info$merge.name
  ) +
  scale_x_continuous(breaks = c(-6, log10(0.05), 0), 
                     labels = c(-6, "*", 0)
  ) + 
  coord_cartesian(xlim = c(-6, 0), 
                  expand = TRUE
  ) +
  geom_vline(xintercept = log10(0.05), 
             color = "black", 
             linetype = "dashed", 
             linewidth = 0.25
  ) +
  theme_void() + 
  theme(
    axis.title.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, 
                               color = c(rep(Pful_color,10), rep(Pame_color,10), rep(Bger_color,10), rep(Esin_color, 10))
                               # color = "black"
    ), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x = unit(2, "pt"),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.25, linetype = "solid"),
    legend.text.position = "right",
    legend.key.height = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.size = unit(10, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.direction = "vertical",
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
  )

top10_depleted_module_p0

# top10_enriched_module module category annotation
top10_depleted_module_p1 <- 
  ggplot(top10_depleted_module_info, 
         aes(x = 0, y = module, fill = Module.category, color = Module.category), 
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.25,
            lineend = "square",
            show.legend = T
  ) + 
  scale_fill_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:9],"grey"), drop = FALSE) +
  scale_color_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:9],"grey"), drop = FALSE) +
  labs(fill = "Module categary") + 
  guides(fill = guide_legend(ncol =1)) +
  theme_void() +
  theme(
    legend.key.size = unit(10, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.direction = "vertical",
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.text.position = "right",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
  )

top10_depleted_module_p1

# plot heatmap
sample_levels <- c(paste("Esin_",c(1:6), sep = ""), paste("Bger_",c(1:6), sep = ""), paste("Pame_",c(1:6), sep = ""), paste("Pful_",c(1:6), sep = ""))

down_module_z_scores <- top10_depleted_module[, 12:35] %>% t() %>% scale() %>% t()

top10_depleted_module[, 12:35] = down_module_z_scores

# wide to long
top10_depleted_module_long <- top10_depleted_module %>%
  pivot_longer(cols = -c(1:11), names_to = "sample", values_to = "Enrich_score")

top10_depleted_module_long$sample <- factor(top10_depleted_module_long$sample, levels = rev(sample_levels))

min <- min(top10_depleted_module_long$Enrich_score)
max <- max(top10_depleted_module_long$Enrich_score)

top10_depleted_module_p2 <- 
  ggplot(top10_depleted_module_long, 
         aes(x = sample, y = module, fill = Enrich_score)
  ) +
  geom_tile(height = 1,
            color = NA,
            linewidth = 0,
            lineend = "square",
            show.legend = T
  ) +
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF", 
                       limit = c(min, max),
                       breaks = c(-2, 0 ,2)
  ) +
  annotate("rect", xmin = 18.5, xmax = 24.5, ymin = 30.5, ymax = 40.5, fill = NA, color = "blue", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 12.5, xmax = 18.5, ymin = 20.5, ymax = 30.5, fill = NA, color = "blue", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 6.5, xmax = 12.5, ymin = 10.5, ymax = 20.5, fill = NA, color = "blue", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 0.5, xmax = 6.5, ymin = 0.5, ymax = 10.5, fill = NA, color = "blue", linewidth = 0.25, lineend = "square")+
  labs(fill = "Enrich Z-score") +
  theme_void() + 
  theme(
    legend.position = "none", 
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 1, angle = 0, lineheight = 1, color = "black"),
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0, angle = 0, lineheight = 1, color = "black"),
    legend.key.size = unit(10, "pt"), 
    legend.key.width = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.direction = "horizontal",
    legend.title.position = "top",
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

top10_depleted_module_p2

library(aplot)

top10_depleted_module_p <- 
  top10_depleted_module_p0 %>%
  insert_right(top10_depleted_module_p1, width = 0.2) %>%
  insert_right(top10_depleted_module_p2, width = 1)

top10_depleted_module_p

ggsave("function/2_downstream_analysis/19_top10_depleted_module.svg", plot = top10_depleted_module_p, width = 14, height = 9.5, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/19_top10_depleted_module.tiff", plot = top10_depleted_module_p, width =14, height = 9.5, units = "cm", dpi = 300)


# merge top10_depleted_module_p & top10_enriched_module_p
top10_delepeted_enriched_module_plot <- top10_depleted_module_p0 + top10_depleted_module_p1 + top10_depleted_module_p2 +
  top10_enriched_module_p2 +top10_enriched_module_p1 + top10_enriched_module_p0 +  
  plot_layout(ncol = 6, widths = c(3, 0.5, 4, 4, 0.5, 3)) & 
  theme(panel.spacing = unit(1, "line"))

top10_delepeted_enriched_module_plot

top10_delepeted_enriched_module_plot_0 <- top10_depleted_module_p0 + theme(axis.text.y = element_blank()) + top10_depleted_module_p1 + top10_depleted_module_p2 +
  top10_enriched_module_p2 + top10_enriched_module_p1 + top10_enriched_module_p0 + theme(axis.text.y = element_blank()) +
  plot_layout(ncol = 6, widths = c(4, 0.35, 4, 4, 0.35, 4)) & 
  theme(panel.spacing = unit(1, "line"))

top10_delepeted_enriched_module_plot_0

ggsave("function/2_downstream_analysis/20_top10_delepeted_enriched_module.svg", plot = top10_delepeted_enriched_module_plot, width = 20, height = 8, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/20_top10_delepeted_enriched_module.tiff", plot = top10_delepeted_enriched_module_plot, width = 20, height = 8, units = "cm", dpi = 300)

ggsave("function/2_downstream_analysis/21_top10_delepeted_enriched_module_0.svg", plot = top10_delepeted_enriched_module_plot_0, width = 7, height = 8, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/21_top10_delepeted_enriched_module_0.tiff", plot = top10_delepeted_enriched_module_plot_0, width = 7, height = 8, units = "cm", dpi = 300)


# save all legend 
top10_enriched_module_p0_legend <- cowplot::get_legend(top10_enriched_module_p0 + theme(legend.position = "right"))
top10_enriched_module_p1_legend <- cowplot::get_legend(top10_enriched_module_p1 + theme(legend.position = "right"))
top10_enriched_module_p2_legend <- cowplot::get_legend(top10_enriched_module_p2 + theme(legend.position = "right"))
top10_enriched_module_plot_legend <- cowplot::plot_grid(top10_enriched_module_p0_legend, top10_enriched_module_p1_legend, top10_enriched_module_p2_legend, nrow = 3)

top10_depleted_module_p0_legend <- cowplot::get_legend(top10_depleted_module_p0 + theme(legend.position = "right"))
top10_depleted_module_p1_legend <- cowplot::get_legend(top10_depleted_module_p1 + theme(legend.position = "right"))
top10_depleted_module_p2_legend <- cowplot::get_legend(top10_depleted_module_p2 + theme(legend.position = "right"))
top10_depleted_module_plot_legend <- cowplot::plot_grid(top10_depleted_module_p0_legend, top10_depleted_module_p1_legend, top10_depleted_module_p2_legend, nrow = 3)

top10_depleted_enriched_module_plot_legend <- cowplot::plot_grid(top10_depleted_module_plot_legend, top10_enriched_module_plot_legend, ncol = 2)

ggsave("function/2_downstream_analysis/22_top10_delepeted_enriched_module_legend.svg", plot = top10_depleted_enriched_module_plot_legend, width = 20, height = 20, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/22_top10_delepeted_enriched_module_legend.tiff", plot = top10_depleted_enriched_module_plot_legend, width = 20, height = 20, units = "cm", dpi = 300)




########################
### Elastic net ###
########################

### metagenome KEGG module importance ranking by Elastic Net

# 1. load data
library(glmnet)
set.seed(20241102)
lasso_data <- read.table("function/1_KEGG/diff_module/ssGSEA/5_lasso_input.tsv", header = TRUE, row.names = 1, check.names = FALSE)
x <- lasso_data %>% dplyr::select(-1) %>% as.matrix()
y <- factor(lasso_data$group)

# 2. screen alpha(0–1, step 0.05)
alpha_values <- seq(0, 1, by = 0.05)
results <- tibble(
  alpha      = alpha_values,
  lambda.min = numeric(length(alpha_values)),
  cv.error   = numeric(length(alpha_values))
)

for (i in seq_along(alpha_values)) {
  cv_fit <- cv.glmnet(x, y,
                      family  = "binomial",
                      alpha   = alpha_values[i],
                      grouped = FALSE)
  results$lambda.min[i] <- cv_fit$lambda.min
  results$cv.error[i]   <- min(cv_fit$cvm)
}

# select alpha
best_alpha <- results$alpha[which.min(results$cv.error)]

# Re-run the CV with the optimal alpha and select λ
cv_final <- cv.glmnet(x, y,
                      family  = "binomial",
                      alpha   = best_alpha)

lambda_min  <- cv_final$lambda.min
lambda_1se  <- cv_final$lambda.1se
print(tibble(lambda.min = lambda_min, lambda.1se = lambda_1se))

plot(results$alpha, results$cv.error, type = "b", xlab = "Alpha", ylab = "CV Error", main = "CV Error vs Alpha")

pdf("function/2_downstream_analysis/23_Cross-validation.pdf",
    width  = 6,
    height = 6)

plot(cv_final)

dev.off()

# cross validation plot
cv_plot <-
  ggplot(results, 
         aes(x = alpha, y = cv.error)
  ) +
  geom_line(colour = "black", 
            linewidth = 0.25,
  ) +
  geom_point(colour = "steelblue",
             shape = 19,
             size = 1) +
  geom_vline(xintercept = best_alpha, 
             linetype = 2, 
             colour = "red"
  ) +
  labs(title = "CV Error vs Alpha", 
       x = "Alpha", 
       y = "CV Error"
  ) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.title.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 00, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1, color = "black"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_line(color ="black", linewidth = 0.25, linetype = 1),
    axis.ticks.y = element_line(color ="black", linewidth = 0.25, linetype = 1),,
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25),
    panel.background = element_blank(),
    legend.position = "none"
  )

cv_plot

ggsave("function/2_downstream_analysis/24_cv_plot.svg", plot = cv_plot, width = 6, height = 6, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/24_cv_plot.tiff", plot = cv_plot, width = 6, height = 6, units = "cm", dpi = 300)


# Extract the modules with non-zero importance(lambda.min)
coef_mat <- as.matrix(coef(cv_final, s = "lambda.min"))
active.idx   <- which(coef_mat[, 1] != 0)
active.genes <- rownames(coef_mat)[active.idx] %>% setdiff("(Intercept)")

# export result
coef.df <- coef(cv_final, s = "lambda.min") %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("module") %>%
  mutate(coef = lambda.min) %>%
  filter(module != "(Intercept)") %>%
  mutate(absCoef = abs(coef)) %>%
  filter(coef != 0) %>%
  arrange(desc(coef))

write_csv(coef.df, "function/2_downstream_analysis/25_Lasso_module_ranked_by_coef.csv")


### top 20 lasso selected feature
top20_lasso_metagenome <- rbind(head(coef.df, 10), tail(coef.df, 10))

metagenome_module_info <- 
  rbind(read.table("function/1_KEGG/diff_module/ssGSEA/1_all_up_module_information.tsv", header = TRUE, sep = "\t"),
        read.table("function/1_KEGG/diff_module/ssGSEA/2_all_down_module_information.tsv", header = TRUE, sep = "\t")) %>%
  distinct(module, .keep_all = TRUE)

top20_lasso_metagenome <- 
  top20_lasso_metagenome %>% 
  left_join(metagenome_module_info, by = c("module" = "module")) %>%
  mutate(label = paste0(module, ":", module_name))

top20_lasso_metagenome$label <- gsub("ko00524:\\s*Neomycin, kanamycin and gentamicin biosynthesis", "ko00524:Neomycin/Kanamycin/Gentamicin biosynthesis", top20_lasso_metagenome$label)

top20_lasso_metagenome <- top20_lasso_metagenome %>%
  separate(label, into = c("first_part", "rest"), sep = ",", remove = FALSE) %>%
  mutate(label = first_part) %>%
  dplyr::select(-first_part, -rest)

top20_lasso_metagenome$label <- gsub("ko00625:\\s*Chloroalkane and chloroalkene degradation", "ko00625:Chloroalkane/Chloroalkene degradation", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("ko00983:\\s*Drug metabolism - other enzymes", "ko00983:Drug metabolism", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("M00613:\\s*Anoxygenic photosynthesis in green nonsulfur bacteria", "M00613:Anoxygenic photosynthesis(bacteria)", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("biosynthesis", "SYN.", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("resistance", "RES.", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("cycle", "CYC.", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("metabolism", "MET.", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("oxidase", "OXID.", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("oxidoreductase", "REDOX.", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("degradation", "DEG.", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("photosynthesis", "PS.", top20_lasso_metagenome$label)

top20_lasso_metagenome$label <- factor(top20_lasso_metagenome$label, levels = rev(top20_lasso_metagenome$label))

top20_lasso_metagenome <- top20_lasso_metagenome %>% dplyr::select(-9)

write.table(top20_lasso_metagenome, "../4.association/1_multi_omic_association/2_top10_e_d_module_information_MetaG.tsv", sep = "\t", quote = F, row.names = F)

metagenome_all_module_ssGSEA_data <- read.table("function/1_KEGG/diff_module/ssGSEA/Total_metagenome_KEGG_ssGSEA_module.tsv", sep = "\t", header = T, row.names = 1)

top20_lasso_metagenome_module_ssGSEA <- metagenome_all_module_ssGSEA_data[top20_lasso_metagenome$module, ] %>% rownames_to_column(var = "module")

write.table(top20_lasso_metagenome_module_ssGSEA, "../4.association/1_multi_omic_association/1_top10_e_d_metagenome_module.tsv", sep = "\t", quote = F, row.names = F)


### total lasso selected feature
total_lasso_metagenome <- 
  coef.df %>% 
  left_join(metagenome_module_info, by = c("module" = "module")) %>%
  mutate(label = paste0(module, ":", module_name))

total_lasso_metagenome$label <- gsub("ko00524:\\s*Neomycin, kanamycin and gentamicin biosynthesis", "ko00524:Neomycin/Kanamycin/Gentamicin biosynthesis", total_lasso_metagenome$label)

total_lasso_metagenome <- total_lasso_metagenome %>%
  separate(label, into = c("first_part", "rest"), sep = ",", remove = FALSE) %>%
  mutate(label = first_part) %>%
  dplyr::select(-first_part, -rest)

total_lasso_metagenome$label <- gsub("ko00625:\\s*Chloroalkane and chloroalkene degradation", "ko00625:Chloroalkane/Chloroalkene degradation", total_lasso_metagenome$label)
total_lasso_metagenome$label <- gsub("ko00983:\\s*Drug metabolism - other enzymes", "ko00983:Drug metabolism", total_lasso_metagenome$label)
total_lasso_metagenome$label <- gsub("M00613:\\s*Anoxygenic photosynthesis in green nonsulfur bacteria", "M00613:Anoxygenic photosynthesis(bacteria)", total_lasso_metagenome$label)
total_lasso_metagenome$label <- gsub("biosynthesis", "SYN.", total_lasso_metagenome$label)
total_lasso_metagenome$label <- gsub("resistance", "RES.", total_lasso_metagenome$label)
total_lasso_metagenome$label <- gsub("cycle", "CYC.", total_lasso_metagenome$label)
total_lasso_metagenome$label <- gsub("metabolism", "MET.", total_lasso_metagenome$label)
total_lasso_metagenome$label <- gsub("oxidase", "OXID.", total_lasso_metagenome$label)
total_lasso_metagenome$label <- gsub("oxidoreductase", "REDOX.", total_lasso_metagenome$label)
total_lasso_metagenome$label <- gsub("degradation", "DEG.", total_lasso_metagenome$label)
total_lasso_metagenome$label <- gsub("photosynthesis", "PS.", total_lasso_metagenome$label)

total_lasso_metagenome$label <- factor(total_lasso_metagenome$label, levels = rev(total_lasso_metagenome$label))

total_lasso_metagenome <- total_lasso_metagenome %>% dplyr::select(-9)

write.table(total_lasso_metagenome, "../4.association/1_multi_omic_association/17_total_e_d_module_information_MetaG.tsv", sep = "\t", quote = F, row.names = F)

total_lasso_metagenome_module_ssGSEA <- metagenome_all_module_ssGSEA_data[coef.df$module, ] %>% rownames_to_column(var = "module")

write.table(total_lasso_metagenome_module_ssGSEA, "../4.association/1_multi_omic_association/16_total_e_d_metagenome_module.tsv", sep = "\t", quote = F, row.names = F)



### lasso ranking plot
lasso_importance_rank <- 
  ggplot(top20_lasso_metagenome,
         aes(x = coef, y = module)
  )+
  geom_point(aes(x = coef, y = label, color = label), 
             size = 1,
             show.legend = F
  ) +
  geom_hline(yintercept = seq_along(top20_lasso_metagenome$label), 
             color = "grey", 
             linetype = "dotted", 
             linewidth = 0.15
  ) +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             lineend = "butt", 
             color = "black", 
             linewidth = 0.15
  ) +
  scale_color_manual(values = c(rep("lightblue", 10), rep("pink", 10))) +
  scale_y_discrete(position = "right") +
  scale_x_continuous(limits = c(-16, 16), breaks = seq(-15, 15, by = 15), labels = seq(-15, 15, by = 15)) +
  theme_void() + 
  theme(
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
    axis.text.y.right = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
    axis.title.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x = unit(2, "pt"),
    panel.border =  element_rect(color = "black", fill = NA, linewidth = 0.25, linetype = "solid"),
  )

lasso_importance_rank

lasso_importance_rank_0 <- lasso_importance_rank + theme(axis.text.y.right = element_blank())

lasso_importance_rank_0

ggsave("function/2_downstream_analysis/26_lasso_importance_rank_plot.svg", plot = lasso_importance_rank, width = 10, height = 4.5, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/26_lasso_importance_rank_plot.tiff", plot = lasso_importance_rank, width = 10, height = 4.5, units = "cm", dpi = 300)

ggsave("function/2_downstream_analysis/27_lasso_importance_rank_plot_0.svg", plot = lasso_importance_rank_0, width = 1, height = 4.5, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/26_lasso_importance_rank_plot_0.tiff", plot = lasso_importance_rank_0, width = 1, height = 4.5, units = "cm", dpi = 300)



### double volcano plot

# read file
double_volcano_data <- read.table("function/1_KEGG/diff_module/ssGSEA/6_double_volcano_log2trans.tsv", header = TRUE)

scater_color <- paletteer_d("ggthemes::Classic_Green_Orange_12")

# Quadrant division
x1 <-   log10(0.05)
x2 <-  -log10(0.05)
y1 <-   log10(0.05)
y2 <-  -log10(0.05)

double_volcano_data <- 
  double_volcano_data %>% 
  mutate(
    quadrant = case_when(
      Bger_vs_Esin < x1 & Peri_vs_Esin >= y2 ~ 1,
      Bger_vs_Esin >= x1 & Bger_vs_Esin < x2 & Peri_vs_Esin >= y2 ~ 2,
      Bger_vs_Esin >= x2 & Peri_vs_Esin >= y2 ~ 3,
      Bger_vs_Esin < x1 & Peri_vs_Esin >= y1 & Peri_vs_Esin < y2 ~ 4,
      Bger_vs_Esin >= x1 & Bger_vs_Esin < x2 & Peri_vs_Esin >= y1 & Peri_vs_Esin < y2 ~ 5,
      Bger_vs_Esin >= x2 & Peri_vs_Esin >= y1 & Peri_vs_Esin < y2 ~ 6,
      Bger_vs_Esin < x1 & Peri_vs_Esin < y1 ~ 7,
      Bger_vs_Esin >= x1 & Bger_vs_Esin < x2 & Peri_vs_Esin < y1 ~ 8,
      Bger_vs_Esin >= x2 & Peri_vs_Esin < y1 ~ 9,
      TRUE ~ NA_integer_
    )
  )

double_volcano_data$quadrant <- as.factor(double_volcano_data$quadrant)

double_volcano_data$size <- abs((double_volcano_data$Bger_vs_Esin + double_volcano_data$Peri_vs_Esin) / 2)

quadrant_counts <- double_volcano_data %>%
  group_by(quadrant) %>%
  summarise(count = n(), .groups = 'drop')

quadrant_counts$quadrant <- factor(quadrant_counts$quadrant, levels = rev(c("3", "7", "2", "8", "4", "6", "1", "9", "5")))

quadrant_counts <- 
  quadrant_counts %>% 
  mutate(type = case_when(
    quadrant == 1 ~ "Peri_enrich & Bger_deplete",
    quadrant == 2 ~ "Peri_enrich & Bger_nosig",
    quadrant == 3 ~ "Peri_enrich & Bger_enrich",
    quadrant == 4 ~ "Peri_nosig & Bger_deplete",
    quadrant == 5 ~ "Both_nosig",
    quadrant == 6 ~ "Peri_nosig & Bger_enrich",
    quadrant == 7 ~ "Peri_deplete & Bger_deplete",
    quadrant == 8 ~ "Peri_deplete & Bger_nosig",
    quadrant == 9 ~ "Peri_deplete & Bger_enrich"
  ))

quadrant_counts$percent <- quadrant_counts$count / sum(quadrant_counts$count)

quadrant_counts$percent_label <- paste(quadrant_counts$count, " (", round(quadrant_counts$percent * 100, 2), "%)", sep = "")

quadrant_counts <- 
  quadrant_counts %>% 
  mutate(Correlation = case_when(
    quadrant == 1 ~ "Negative",
    quadrant == 2 ~ "Mid",
    quadrant == 3 ~ "Positive",
    quadrant == 4 ~ "Mid",
    quadrant == 5 ~ "Weak",
    quadrant == 6 ~ "Mid",
    quadrant == 7 ~ "Positive",
    quadrant == 8 ~ "Mid",
    quadrant == 9 ~ "Negative"
  ))

quadrant_percent_pie_plot <- 
  ggplot(quadrant_counts,
         aes(x = "", y = percent, fill = quadrant)
  ) + 
  geom_bar(stat = "identity",
           position = "stack",
           width = 0.45,
           show.legend = F,
           color = NA,
           linewidth = 0
  ) +
  geom_text(aes(label = paste0(round(percent * 100, 1), "%")), 
            position = position_stack(vjust = 0.5),
            size = 1.5
  ) +
  coord_polar(theta = "y", start = -90) +
  guides(fill = guide_legend(title = "", nrow = 9, ncol = 1, reverse = TRUE))+ 
  scale_fill_manual(values = rev(c("pink", "lightblue", scater_color[3], scater_color[4], scater_color[9], scater_color[10], scater_color[1], scater_color[2], "grey80")),
                    labels = c("1" = "Peri_enrich & Bger_deplete",
                               "2" = "Peri_enrich & Bger_nosig",
                               "3" = "Peri_enrich & Bger_enrich",
                               "4" = "Peri_nosig & Bger_deplete",
                               "5" = "Both_nosig",
                               "6" = "Peri_nosig & Bger_enrich",
                               "7" = "Peri_deplete & Bger_deplete",
                               "8" = "Peri_deplete & Bger_nosig",
                               "9" = "Peri_deplete & Bger_enrich"),
                    drop = FALSE
  ) +
  theme_void()

quadrant_percent_pie_plot

ggsave("function/2_downstream_analysis/28_quadrant_percent_pie_plot.svg", plot = quadrant_percent_pie_plot, width = 4, height = 4, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/28_quadrant_percent_pie_plot.tiff", plot = quadrant_percent_pie_plot, width = 1, height = 1, units = "cm", dpi = 300)

# highlight module
highlight_module <- c(
  M00019   = "M00019", 
  M00570   = "M00570",
  M00432   = "M00432",
  ko00524  = "ko00524",
  M00535   = "M00535",
  M00029   = "M00029",
  ko00365  = "ko00365",
  M00145   = "M00145",
  M00153   = "M00153",
  M00718   = "M00718",
  M00376   = "M00376",
  M00613   = "M00613",
  M00736   = "M00736",
  ko00908  = "ko00908",
  M00033   = "M00033",
  ko00625  = "ko00625",
  M00843   = "M00843",
  M00842   = "M00842",
  ko00983  = "ko00983",
  M00845   = "M00845"
)

double_volcano_data$label <- highlight_module[double_volcano_data$gene]

quadrant_colors <- c(
  "1" = "#ffffd2", 
  "2" = "white",
  "3" = "#ffffd2",
  "4" = "white",
  "5" = "white",
  "6" = "white",
  "7" = "#ffffd2",
  "8" = "white",
  "9" = "#ffffd2"
)

# plot scatter
double_volcano_plot <- 
  ggplot(double_volcano_data, 
         aes(x = Bger_vs_Esin, y = Peri_vs_Esin , fill = quadrant, color = quadrant)
  ) +
  # xlim(c(-20, 20)) + 
  # ylim(c(-20, 20)) +
  # plot scatter
  geom_point(aes(size = size),
             alpha = 0.5,
             stroke = 0, 
             shape = 19, 
             show.legend = F
  ) +
  # plot vertical line
  geom_vline(xintercept = c(-log10(0.05), log10(0.05)), 
             linetype = 2, 
             linewidth = 0.25, 
             alpha = 1, 
             color = "grey"
  ) +
  geom_vline(xintercept = 0, 
             linetype = 2, 
             linewidth = 0.1, 
             alpha = 0.5, 
             color = "grey30"
  ) +
  # plot Horizontal line
  geom_hline(yintercept = c(-log10(0.05), log10(0.05)), 
             linetype = 2, 
             linewidth = 0.25, 
             alpha = 1, 
             color = "grey"
  ) +
  geom_hline(yintercept = 0, 
             linetype = 2, 
             linewidth = 0.1, 
             alpha = 0.5, 
             color = "grey30"
  ) +
  # highlight modules
  geom_text_repel(aes(label = label),
                  color = "black",
                  size = 1.5,
                  min.segment.length = 0,
                  # arrow = arrow(angle = 15, length = unit(0.01, "inches"), ends = "last", type = "closed"),
                  box.padding = 0.05,
                  point.padding = 0.05,
                  segment.size = 0,
                  segment.color = "black",
                  show.legend = FALSE,
                  nudge_x = 0.02,
                  nudge_y = 0.02,
                  direction = "both",
                  force = 10,
                  force_pull = 0.5,
                  max.overlaps = Inf
  ) +
  # set X / Y axis title 
  labs(x = "Bger vs Esin: [Dir] x -log10(P.adj)", 
       y = "Peri vs Esin: [Dir] x -log10(P.adj)"
  ) +
  scale_color_manual(values = c(scater_color[1], scater_color[3], "pink", scater_color[9], "grey80", scater_color[10], "lightblue", scater_color[4] ,scater_color[2]), drop = FALSE) +
  scale_fill_manual(values = c(scater_color[1], scater_color[3], "pink", scater_color[9], "grey80", scater_color[10], "lightblue", scater_color[4] ,scater_color[2]), drop = FALSE) +
  scale_size_continuous(range = c(1, 2))+
  scale_alpha_continuous(range = c(0.4 ,1))+
  theme_void() + 
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25, linetype = "solid"),
    axis.title.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
    axis.title.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1.2, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1.2, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length = unit(2, "pt"),
  )

print(double_volcano_plot)

ggsave("function/2_downstream_analysis/29_double_volcano_plot.svg", plot = double_volcano_plot, width = 4.5, height = 4.5, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/29_double_volcano_plot.tiff", plot = double_volcano_plot, width = 4.5, height = 4.5, units = "cm", dpi = 300)




### venn
# set color 
Esin_venn_color <- generate_palette(Esin_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, view_labels = FALSE)
Bger_venn_color <- generate_palette(Bger_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, view_labels = FALSE)
Pame_venn_color <- generate_palette(Pame_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, view_labels = FALSE)
Pful_venn_color <- generate_palette(Pful_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, view_labels = FALSE)

# Esin enrich
Esin_venn_dat_up  <- read.delim("function/1_KEGG/diff_module/ssGSEA/Esin_metagenome_KEGG_ssGSEA_module_each_up.tsv")
Esin_venn_up_list <- as.list(Esin_venn_dat_up)

Esin_venn_up <- ggvenn(
  data = Esin_venn_up_list,
  show_elements = F,
  label_sep = "\n",
  show_percentage = F,
  digits = 2,
  fill_color = Esin_venn_color,
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_alpha = 0,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = NA,
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
) + 
  theme_void()

Esin_venn_up

# Bger enrich
Bger_venn_dat_up  <- read.delim("function/1_KEGG/diff_module/ssGSEA/Bger_metagenome_KEGG_ssGSEA_module_each_up.tsv")
Bger_venn_up_list <- as.list(Bger_venn_dat_up )

Bger_venn_up <- ggvenn(
  data = Bger_venn_up_list,
  show_elements = F,
  label_sep = "\n",
  show_percentage = F,
  digits = 2,
  fill_color = Bger_venn_color,
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_alpha = 0,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = NA,
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
) + 
  theme_void()

Bger_venn_up

# Pame enrich
Pame_venn_dat_up  <- read.delim("function/1_KEGG/diff_module/ssGSEA/Pame_metagenome_KEGG_ssGSEA_module_each_up.tsv")
Pame_venn_up_list <- as.list(Pame_venn_dat_up)

Pame_venn_up <- ggvenn(
  data = Pame_venn_up_list,
  show_elements = F,
  label_sep = "\n",
  show_percentage = F,
  digits = 2,
  fill_color = Pame_venn_color,
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_alpha = 0,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = NA,
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
) + 
  theme_void()

Pame_venn_up

# Pful enrich
Pful_venn_dat_up  <- read.delim("function/1_KEGG/diff_module/ssGSEA/Pful_metagenome_KEGG_ssGSEA_module_each_up.tsv")
Pful_venn_up_list <- as.list(Pful_venn_dat_up )

Pful_venn_up <- ggvenn(
  data = Pful_venn_up_list,
  show_elements = F,
  label_sep = "\n",
  show_percentage = F,
  digits = 2,
  fill_color = Pful_venn_color,
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_alpha = 0,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = NA,
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
) + 
  theme_void()

Pful_venn_up

library(aplot)

venn_up <- Esin_venn_up + Bger_venn_up + Pame_venn_up + Pful_venn_up + plot_layout(ncol = 4)

ggsave("function/2_downstream_analysis/30_venn_up.svg", plot = venn_up, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/30_venn_up.tiff", plot = venn_up, width = 10, height = 5, units = "cm", dpi = 300)




# Esin deplete
Esin_venn_dat_down  <- read.delim("function/1_KEGG/diff_module/ssGSEA/Esin_metagenome_KEGG_ssGSEA_module_each_down.tsv")
Esin_venn_down_list <- as.list(Esin_venn_dat_down)

Esin_venn_down <- ggvenn(
  data = Esin_venn_down_list,
  show_elements = F,
  label_sep = "\n",
  show_percentage = F,
  digits = 2,
  fill_color = Esin_venn_color,
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_alpha = 0,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = NA,
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
) + 
  theme_void()

Esin_venn_down

# Bger deplete
Bger_venn_dat_down  <- read.delim("function/1_KEGG/diff_module/ssGSEA/Bger_metagenome_KEGG_ssGSEA_module_each_down.tsv")
Bger_venn_down_list <- as.list(Bger_venn_dat_down )

Bger_venn_down <- ggvenn(
  data = Bger_venn_down_list,
  show_elements = F,
  label_sep = "\n",
  show_percentage = F,
  digits = 2,
  fill_color = Bger_venn_color,
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_alpha = 0,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = NA,
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
) + 
  theme_void()

Bger_venn_down

# Pame deplete
Pame_venn_dat_down  <- read.delim("function/1_KEGG/diff_module/ssGSEA/Pame_metagenome_KEGG_ssGSEA_module_each_down.tsv")
Pame_venn_down_list <- as.list(Pame_venn_dat_down )

Pame_venn_down <- ggvenn(
  data = Pame_venn_down_list,
  show_elements = F,
  label_sep = "\n",
  show_percentage = F,
  digits = 2,
  fill_color = Pame_venn_color,
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_alpha = 0,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = NA,
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
) + 
  theme_void()

Pame_venn_down

# Pful deplete
Pful_venn_dat_down  <- read.delim("function/1_KEGG/diff_module/ssGSEA/Pful_metagenome_KEGG_ssGSEA_module_each_down.tsv")
Pful_venn_down_list <- as.list(Pful_venn_dat_down )

Pful_venn_down <- ggvenn(
  data = Pful_venn_down_list,
  show_elements = F,
  label_sep = "\n",
  show_percentage = F,
  digits = 2,
  fill_color = Pful_venn_color,
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_alpha = 0,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = NA,
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
) + 
  theme_void()

Pful_venn_down

library(aplot)

venn_down <- Esin_venn_down + Bger_venn_down  + Pame_venn_down + Pful_venn_down + plot_layout(ncol =4)

ggsave("function/2_downstream_analysis/31_venn_down.svg", plot = venn_down, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/31_venn_down.tiff", plot = venn_down, width = 10, height = 5, units = "cm", dpi = 300)


# plot regulation statistics bar plot
regulation <- read.table("function/1_KEGG/diff_module/ssGSEA/regulation_stat.txt",header = TRUE)
regulation$X <- factor(regulation$X, levels = unique(regulation$X))
regulation$group <- factor(regulation$group, levels = c("Esin", "Bger", "Pame", "Pful"))

regulation_stat <- 
  ggplot(regulation, 
         aes(x = X, y = number, fill = X)
  ) +
  facet_wrap(~ group, 
             scale = "free_x", 
             ncol = 4, 
             nrow = 1
  ) +
  scale_fill_manual(values = c(Esin_venn_color, Bger_venn_color, Pame_venn_color, Pful_venn_color)) + 
  geom_bar(stat = "identity", 
           position = position_stack(), 
           width = 0.95,
           show.legend = F
  ) +
  geom_hline(yintercept = 0, 
             linewidth = 0.25, 
             linetype = 2, 
             color = "purple"
  ) + 
  labs(fill = "Comparison")+
  theme_void() +
  theme(
    axis.text.y = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
    axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.y = unit(2, "pt"),
    strip.switch.pad.wrap = unit(3, "pt"),
    strip.text.x = element_blank(),
    strip.clip = "on",
  )

regulation_stat

ggsave("function/2_downstream_analysis/32_regulation_stat.svg", plot = regulation_stat, width = 10, height = 10, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/32_regulation_stat.tiff", plot = regulation_stat, width = 10, height = 10, units = "cm", dpi = 300)


# merge plot
merge_stat_plot <- venn_up / regulation_stat / venn_down + plot_layout(nrow = 3)
print(merge_stat_plot)


ggsave("function/2_downstream_analysis/33_merge_stat_plot.svg", plot = merge_stat_plot, width = 10, height = 10, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/33_merge_stat_plot.tiff", plot = merge_stat_plot, width = 10, height = 10, units = "cm", dpi = 300)







### taxonomy contribution of BCAA gene
library(glmnet)
library(caret)
library(dplyr)

### scale_I
BCAA_KO <- read.table("function/2_downstream_analysis/34_BCAA_KO_abundance.txt", header = T, sep = "\t", row.names = 1)

K00052 <- read.table("function/2_downstream_analysis/K00052_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:31, names_to = "Sample", values_to = "Abundance") %>%
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(mean)) %>%
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:29, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% c("Esin_1", "Esin_2", "Esin_3", "Esin_4", "Esin_5", "Esin_6") ~ "Esin",
    Sample %in% c("Bger_1", "Bger_2", "Bger_3", "Bger_4", "Bger_5", "Bger_6") ~ "Bger",
    Sample %in% c("Pame_1", "Pame_2", "Pame_3", "Pame_4", "Pame_5", "Pame_6") ~ "Pame",
    Sample %in% c("Pful_1", "Pful_2", "Pful_3", "Pful_4", "Pful_5", "Pful_6") ~ "Pful",
    TRUE ~ NA_character_
  )) %>%
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("Esin", "Bger", "Pame", "Pful"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K00052")

K00053 <- read.table("function/2_downstream_analysis/K00053_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:31, names_to = "Sample", values_to = "Abundance") %>%
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(mean)) %>%
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:29, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% c("Esin_1", "Esin_2", "Esin_3", "Esin_4", "Esin_5", "Esin_6") ~ "Esin",
    Sample %in% c("Bger_1", "Bger_2", "Bger_3", "Bger_4", "Bger_5", "Bger_6") ~ "Bger",
    Sample %in% c("Pame_1", "Pame_2", "Pame_3", "Pame_4", "Pame_5", "Pame_6") ~ "Pame",
    Sample %in% c("Pful_1", "Pful_2", "Pful_3", "Pful_4", "Pful_5", "Pful_6") ~ "Pful",
    TRUE ~ NA_character_
  )) %>%
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("Esin", "Bger", "Pame", "Pful"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K00053")

K00826 <- read.table("function/2_downstream_analysis/K00826_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:31, names_to = "Sample", values_to = "Abundance") %>%
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(mean)) %>%
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:29, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% c("Esin_1", "Esin_2", "Esin_3", "Esin_4", "Esin_5", "Esin_6") ~ "Esin",
    Sample %in% c("Bger_1", "Bger_2", "Bger_3", "Bger_4", "Bger_5", "Bger_6") ~ "Bger",
    Sample %in% c("Pame_1", "Pame_2", "Pame_3", "Pame_4", "Pame_5", "Pame_6") ~ "Pame",
    Sample %in% c("Pful_1", "Pful_2", "Pful_3", "Pful_4", "Pful_5", "Pful_6") ~ "Pful",
    TRUE ~ NA_character_
  )) %>%
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("Esin", "Bger", "Pame", "Pful"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K00826")

K01649 <- read.table("function/2_downstream_analysis/K01649_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:31, names_to = "Sample", values_to = "Abundance") %>%
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(mean)) %>%
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:29, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% c("Esin_1", "Esin_2", "Esin_3", "Esin_4", "Esin_5", "Esin_6") ~ "Esin",
    Sample %in% c("Bger_1", "Bger_2", "Bger_3", "Bger_4", "Bger_5", "Bger_6") ~ "Bger",
    Sample %in% c("Pame_1", "Pame_2", "Pame_3", "Pame_4", "Pame_5", "Pame_6") ~ "Pame",
    Sample %in% c("Pful_1", "Pful_2", "Pful_3", "Pful_4", "Pful_5", "Pful_6") ~ "Pful",
    TRUE ~ NA_character_
  )) %>%
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("Esin", "Bger", "Pame", "Pful"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K01649")

K01652 <- read.table("function/2_downstream_analysis/K01652_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:31, names_to = "Sample", values_to = "Abundance") %>%
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(mean)) %>%
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:29, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% c("Esin_1", "Esin_2", "Esin_3", "Esin_4", "Esin_5", "Esin_6") ~ "Esin",
    Sample %in% c("Bger_1", "Bger_2", "Bger_3", "Bger_4", "Bger_5", "Bger_6") ~ "Bger",
    Sample %in% c("Pame_1", "Pame_2", "Pame_3", "Pame_4", "Pame_5", "Pame_6") ~ "Pame",
    Sample %in% c("Pful_1", "Pful_2", "Pful_3", "Pful_4", "Pful_5", "Pful_6") ~ "Pful",
    TRUE ~ NA_character_
  )) %>%
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("Esin", "Bger", "Pame", "Pful"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K01652")

K01653 <- read.table("function/2_downstream_analysis/K01653_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:31, names_to = "Sample", values_to = "Abundance") %>%
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(mean)) %>%
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:29, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% c("Esin_1", "Esin_2", "Esin_3", "Esin_4", "Esin_5", "Esin_6") ~ "Esin",
    Sample %in% c("Bger_1", "Bger_2", "Bger_3", "Bger_4", "Bger_5", "Bger_6") ~ "Bger",
    Sample %in% c("Pame_1", "Pame_2", "Pame_3", "Pame_4", "Pame_5", "Pame_6") ~ "Pame",
    Sample %in% c("Pful_1", "Pful_2", "Pful_3", "Pful_4", "Pful_5", "Pful_6") ~ "Pful",
    TRUE ~ NA_character_
  )) %>%
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("Esin", "Bger", "Pame", "Pful"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K01653")

K01687 <- read.table("function/2_downstream_analysis/K01687_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:31, names_to = "Sample", values_to = "Abundance") %>%
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(mean)) %>%
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:29, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% c("Esin_1", "Esin_2", "Esin_3", "Esin_4", "Esin_5", "Esin_6") ~ "Esin",
    Sample %in% c("Bger_1", "Bger_2", "Bger_3", "Bger_4", "Bger_5", "Bger_6") ~ "Bger",
    Sample %in% c("Pame_1", "Pame_2", "Pame_3", "Pame_4", "Pame_5", "Pame_6") ~ "Pame",
    Sample %in% c("Pful_1", "Pful_2", "Pful_3", "Pful_4", "Pful_5", "Pful_6") ~ "Pful",
    TRUE ~ NA_character_
  )) %>%
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("Esin", "Bger", "Pame", "Pful"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K01687")

K01703 <- read.table("function/2_downstream_analysis/K01703_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:31, names_to = "Sample", values_to = "Abundance") %>%
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(mean)) %>%
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:29, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% c("Esin_1", "Esin_2", "Esin_3", "Esin_4", "Esin_5", "Esin_6") ~ "Esin",
    Sample %in% c("Bger_1", "Bger_2", "Bger_3", "Bger_4", "Bger_5", "Bger_6") ~ "Bger",
    Sample %in% c("Pame_1", "Pame_2", "Pame_3", "Pame_4", "Pame_5", "Pame_6") ~ "Pame",
    Sample %in% c("Pful_1", "Pful_2", "Pful_3", "Pful_4", "Pful_5", "Pful_6") ~ "Pful",
    TRUE ~ NA_character_
  )) %>%
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("Esin", "Bger", "Pame", "Pful"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K01703")


K01704 <- read.table("function/2_downstream_analysis/K01704_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:31, names_to = "Sample", values_to = "Abundance") %>%
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(mean)) %>%
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:29, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% c("Esin_1", "Esin_2", "Esin_3", "Esin_4", "Esin_5", "Esin_6") ~ "Esin",
    Sample %in% c("Bger_1", "Bger_2", "Bger_3", "Bger_4", "Bger_5", "Bger_6") ~ "Bger",
    Sample %in% c("Pame_1", "Pame_2", "Pame_3", "Pame_4", "Pame_5", "Pame_6") ~ "Pame",
    Sample %in% c("Pful_1", "Pful_2", "Pful_3", "Pful_4", "Pful_5", "Pful_6") ~ "Pful",
    TRUE ~ NA_character_
  )) %>%
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("Esin", "Bger", "Pame", "Pful"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K01704")


K01754 <- read.table("function/2_downstream_analysis/K01754_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:31, names_to = "Sample", values_to = "Abundance") %>%
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(mean)) %>%
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:29, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% c("Esin_1", "Esin_2", "Esin_3", "Esin_4", "Esin_5", "Esin_6") ~ "Esin",
    Sample %in% c("Bger_1", "Bger_2", "Bger_3", "Bger_4", "Bger_5", "Bger_6") ~ "Bger",
    Sample %in% c("Pame_1", "Pame_2", "Pame_3", "Pame_4", "Pame_5", "Pame_6") ~ "Pame",
    Sample %in% c("Pful_1", "Pful_2", "Pful_3", "Pful_4", "Pful_5", "Pful_6") ~ "Pful",
    TRUE ~ NA_character_
  )) %>%
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("Esin", "Bger", "Pame", "Pful"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K01754")

K01754 <- read.table("function/2_downstream_analysis/K01754_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:31, names_to = "Sample", values_to = "Abundance") %>%
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(mean)) %>%
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:29, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% c("Esin_1", "Esin_2", "Esin_3", "Esin_4", "Esin_5", "Esin_6") ~ "Esin",
    Sample %in% c("Bger_1", "Bger_2", "Bger_3", "Bger_4", "Bger_5", "Bger_6") ~ "Bger",
    Sample %in% c("Pame_1", "Pame_2", "Pame_3", "Pame_4", "Pame_5", "Pame_6") ~ "Pame",
    Sample %in% c("Pful_1", "Pful_2", "Pful_3", "Pful_4", "Pful_5", "Pful_6") ~ "Pful",
    TRUE ~ NA_character_
  )) %>%
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("Esin", "Bger", "Pame", "Pful"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K01754")

K09011 <- read.table("function/2_downstream_analysis/K09011_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:31, names_to = "Sample", values_to = "Abundance") %>%
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(mean)) %>%
  # head(.,5) %>% 
  mutate(genus = factor(genus, levels = rev(genus))) %>%
  dplyr::select(-mean) %>% 
  pivot_longer(6:29, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% c("Esin_1", "Esin_2", "Esin_3", "Esin_4", "Esin_5", "Esin_6") ~ "Esin",
    Sample %in% c("Bger_1", "Bger_2", "Bger_3", "Bger_4", "Bger_5", "Bger_6") ~ "Bger",
    Sample %in% c("Pame_1", "Pame_2", "Pame_3", "Pame_4", "Pame_5", "Pame_6") ~ "Pame",
    Sample %in% c("Pful_1", "Pful_2", "Pful_3", "Pful_4", "Pful_5", "Pful_6") ~ "Pful",
    TRUE ~ NA_character_
  )) %>%
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("Esin", "Bger", "Pame", "Pful"))) %>%
  mutate(KO = "K09011")


K17989 <- read.table("function/2_downstream_analysis/K17989_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:31, names_to = "Sample", values_to = "Abundance") %>%
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(mean)) %>%
  filter(phylum != "p_unclassify_d_Bacteria") %>% 
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:29, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% c("Esin_1", "Esin_2", "Esin_3", "Esin_4", "Esin_5", "Esin_6") ~ "Esin",
    Sample %in% c("Bger_1", "Bger_2", "Bger_3", "Bger_4", "Bger_5", "Bger_6") ~ "Bger",
    Sample %in% c("Pame_1", "Pame_2", "Pame_3", "Pame_4", "Pame_5", "Pame_6") ~ "Pame",
    Sample %in% c("Pful_1", "Pful_2", "Pful_3", "Pful_4", "Pful_5", "Pful_6") ~ "Pful",
    TRUE ~ NA_character_
  )) %>%
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("Esin", "Bger", "Pame", "Pful"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K17989")

BCAA_KO_taxonomy <-  
  rbind(K00052, K00053, K00826, K01649, K01652, K01653, K01687, K01703, K01704, K01754, K09011, K17989)

top10 <-
  BCAA_KO_taxonomy %>%
  group_by(genus) %>%
  summarise(mean = mean(mean_TPM)) %>%
  arrange(desc(mean)) %>%
  head(.,10) 

BCAA_KO_taxonomy_top10 <- 
  BCAA_KO_taxonomy %>%
  filter(genus %in% top10$genus) %>%
  mutate(genus = factor(genus, levels = rev(top10$genus)))

BCAA_KO_taxonomy_top10_bubble <- 
  ggplot(BCAA_KO_taxonomy_top10,
         aes(x = group, y = genus, size = zscore, fill = phylum)
  ) +
  geom_point(shape = 21, 
             alpha = 1,
             color = "black",
             stroke = 0.1
  ) +
  scale_size_continuous(range = c(1, 4)) +
  scale_fill_manual(values = c("p_Actinomycetota" = "#8dcfbb", "p_Bacteroidota" =  "#76c4ea", "p_Bacillota_A" =  "#ea8c8c", "p_Desulfobacterota" =  "#bab8d9", "p_Planctomycetota" = "#fbc420", "p_Synergistota" = "#b13b72"))+
  facet_wrap(~KO, 
             ncol = 2,
             scales = "free_y"
  ) +
  theme_void() + 
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25, linetype = "solid"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "italic", size = 6, hjust = 1, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length = unit(2, "pt"),
    legend.position = "none",
    strip.text = element_text(family = "Arial", face = "bold", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black")
  )

BCAA_KO_taxonomy_top10_bubble

ggsave("function/2_downstream_analysis/35_BCAA_KO_taxonomy_top10_bubble.svg", plot = BCAA_KO_taxonomy_top10_bubble, width = 10, height = 16, units = "cm", dpi = 300)
  
BCAA_KO_taxonomy_top10_bubble_legend <- cowplot::get_legend(BCAA_KO_taxonomy_top10_bubble + theme(legend.position = "right"))
BCAA_KO_taxonomy_top10_bubble_legend <-  cowplot::plot_grid(BCAA_KO_taxonomy_top10_bubble_legend, ncol = 1)
BCAA_KO_taxonomy_top10_bubble_legend

ggsave("function/2_downstream_analysis/36_BCAA_KO_taxonomy_top10_bubble_legend.svg", plot = BCAA_KO_taxonomy_top10_bubble_legend, width = 10, height = 16, units = "cm", dpi = 300) 




###
library(scatterpie)
library(PieGlyph)
library(scales)

BCAA_KO_taxonomy_top10_trans <-
  BCAA_KO_taxonomy %>%
  filter(genus %in% top10$genus) %>%
  pivot_wider(id_cols = c(phylum, class, order, family, genus, KO), names_from = group, values_from = mean_TPM) %>%
  mutate(total = Esin + Bger + Pame + Pful) %>%
  mutate(peridomestic = Bger + Pame + Pful) %>% 
  mutate(genus = factor(genus, levels = top10$genus))
  
BCAA_KO_taxonomy_top10_trans <-
  BCAA_KO_taxonomy_top10_trans %>%
  mutate(KO =  factor(KO, levels = rev(unique(BCAA_KO_taxonomy_top10_trans$KO))))

p <- ggplot(BCAA_KO_taxonomy_top10_trans) +
  geom_pie_glyph(
    aes(x = genus, radius = total, y = KO),
    #slices = c( "peridomestic", "Esin"),
    slices = c( "Pame", "Pful", "Bger", "Esin"),
    colour = "white",
    linewidth = 0.1
  ) +
  scale_radius(range = c(0.1, 0.35), name = "TPM") +
  #scale_fill_manual(values = c("Esin" = alpha(Esin_color, 0.5), "peridomestic" = alpha(Pame_color, 0.5))) +
  scale_fill_manual(values = c("Esin" = alpha(Esin_color, 0.5), "Bger" = alpha(Bger_color, 0.5), "Pame" = alpha(Pame_color, 0.5), "Pful" = alpha(Pful_color, 0.5))) +
  theme_void() + 
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25, linetype = "solid"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(family = "Arial", face = "bold", size = 6, hjust = 1, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length = unit(2, "pt"),
    legend.position = "none",
    strip.text = element_text(family = "Arial", face = "bold", size=6, hjust=0, vjust=0.5, angle=0, lineheight=1, color="black")
  )

p

ggsave("function/2_downstream_analysis/37_BCAA_KO_taxonomy_slice_bubble_plot.svg", plot = p, width = 7.5, height = 7, units = "cm", dpi = 300)

label <- BCAA_KO_taxonomy_top10_trans %>% mutate(labels = sub("^g_", "", genus)) %>% dplyr::select(genus, labels) %>% unique()

genus_annotation <- 
  ggplot(BCAA_KO_taxonomy_top10_trans
  ) + 
  geom_point(aes(x = genus, fill = phylum, color = phylum),
             y = 0.25,
             size = 2,
             shape = 21,
             show.legend = F
  ) +
  geom_label(data = label,
             aes(label = labels,x = genus, y =0),
             size = 6/2.83,
             angle = 45,
             family = "Arial",
             fontface = "italic",
             color = "black",
             fill = NA,
             label.size = 0,
             hjust = 1,
             vjust = 1
             
  ) + 
  ylim(-0.5, 0.5) + 
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("p_Bacteroidota" =  "#76c4ea", "p_Bacillota_A" =  "#ea8c8c", "p_Desulfobacterota" =  "#bab8d9", "p_Planctomycetota" = "#fbc420"))+
  scale_color_manual(values = c("p_Bacteroidota" =  "#76c4ea", "p_Bacillota_A" =  "#ea8c8c", "p_Desulfobacterota" =  "#bab8d9", "p_Planctomycetota" = "#fbc420"))+
  theme_void()

genus_annotation

ggsave("function/2_downstream_analysis/38_BCAA_KO_taxonomy_slice_bubble_plot_genus_annotation.svg", plot = genus_annotation, width = 7.5, height = 3, units = "cm", dpi = 300)




Genus_total <- 
  BCAA_KO_taxonomy_top10_trans %>%
  group_by(genus) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") 

Genus_total_p <- 
  ggplot(Genus_total) +
  geom_pie_glyph(
    aes(x = genus),
    radius = 0.25,
    y = 0.5,
    slices = c( "Pame", "Pful", "Bger", "Esin"),
    colour = "white",
    linewidth = 0.1,
    show.legend = F
  ) +
  scale_radius(range = c(0.1, 0.35), name = "TPM") +
  scale_fill_manual(values = c("Esin" = alpha(Esin_color, 0.5), "Bger" = alpha(Bger_color, 0.5), "Pame" = alpha(Pame_color, 0.5), "Pful" = alpha(Pful_color, 0.5))) +
  theme_void()

Genus_total_p


KO_total <- 
  BCAA_KO_taxonomy_top10_trans %>%
  group_by(KO) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")

KO_total_p <- 
  ggplot(KO_total) +
  geom_pie_glyph(
    
    aes(y = KO),
    radius = 0.25,
    x = 0.5,
    slices = c( "Pame", "Pful", "Bger", "Esin"),
    colour = "white",
    linewidth = 0.1,
    show.legend = F
  ) +
  scale_radius(range = c(0.1, 0.35), name = "TPM") +
  scale_fill_manual(values = c("Esin" = alpha(Esin_color, 0.5), "Bger" = alpha(Bger_color, 0.5), "Pame" = alpha(Pame_color, 0.5), "Pful" = alpha(Pful_color, 0.5))) +
  theme_void()

KO_total_p


p2 <- p %>% insert_top(Genus_total_p, height = 0.1) %>% insert_right(KO_total_p, width = 0.1)


p2





### correlation plot
total <- 
  read.table("../1.metagenome/taxonomy/1_Total_coverm_contig_count_species_absolute.tsv", header = T, sep = "\t") %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>%
  t() %>%
  as.vector()

genus_abundance <- 
  read.table("../1.metagenome/taxonomy/1_Total_coverm_contig_count_species_absolute.tsv", header = T, sep = "\t")  %>%
  separate(species, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  group_by(phylum, class, order, family, genus) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>%
  filter(genus %in% c("g_QVMG01", "g_Frigididesulfovibrio", "g_Dysgonomonas", "g_Bacteroides_H")) %>%
  # rename_with(~ str_remove(., "^Metagenome_"), .cols = everything()) %>%
  dplyr::select(-c(1:4))

genus_abundance[2:25] <- genus_abundance[2:25] / total 

order <- names(genus_abundance)[2:25] %>% sort()

genus_abundance <- genus_abundance %>% dplyr::select(1, all_of(order)) %>% column_to_rownames(var = "genus") 

BCAA_metabolite <- 
  read.table("../2.metabolome/1_metabolome_matrix/6_metabolome_LOESS_PQN_missForest_matrix_noQC.tsv", header = T, sep = "\t") %>%
  filter(CPD_ID %in% c("CPD0020", "CPD0128", "CPD0158")) %>%
  mutate(MS2_name = case_when(
    CPD_ID == "CPD0158" ~ "L-Leucine",
    CPD_ID == "CPD0128" ~ "L-Isoleucine",
    CPD_ID == "CPD0020" ~ "L-Valine"
  )) %>% 
  dplyr::select(CPD_ID, MS2_name, everything()) %>%
  dplyr::select(-CPD_ID)

BCAA_metabolite <- BCAA_metabolite %>% dplyr::select(1, all_of(order)) %>% column_to_rownames(var = "MS2_name") 

HostT_pathway <- 
  read.table("../3.transcriptome/1_OG/transcriptome_KEGG_pathway_ssGSEA_matrix.tsv", header = T, sep = "\t") %>%
  filter(gene_set %in% c("ko04150", "ko00190", "ko04214")) %>%
  mutate(pathway = case_when(
    gene_set == "ko04150" ~ "ko04150:mTOR",
    gene_set == "ko00190" ~ "ko00190:OXPHOS",
    gene_set == "ko04214" ~ "ko04214:Apoptosis"
  )) %>% 
  dplyr::select(gene_set, pathway, everything()) %>%
  dplyr::select(-gene_set)

HostT_pathway <- HostT_pathway %>% dplyr::select(1, all_of(order)) %>% column_to_rownames(var = "pathway") 

Micro_MetaB_HostT <- 
  rbind(genus_abundance, BCAA_metabolite, HostT_pathway) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>%
  mutate(group = c(rep("Bger", 6),  rep("Esin", 6), rep("Pame", 6), rep("Pful", 6))) %>%
  mutate(group = factor(group, levels = c("Esin", 'Bger', "Pame", "Pful")))

data <- Micro_MetaB_HostT

micro_raw <- data %>% dplyr::select(starts_with("g_"))
meta_raw <- data %>% dplyr::select(`L-Leucine`, `L-Isoleucine`, `L-Valine`)
host_raw <- data %>% dplyr::select(starts_with("ko"))
env_info <- data %>% dplyr::select(sample, group)

pseudo <- min(micro_raw[micro_raw > 0]) / 2
micro_log <- log10(micro_raw + pseudo)

meta_scaled <- as.data.frame(scale(meta_raw))
host_scaled <- as.data.frame(scale(host_raw))

data_norm <- cbind(env_info, micro_log, meta_scaled, host_scaled)

pca_data <- cbind(meta_scaled, host_scaled)

pca_res <- prcomp(pca_data)

pca_plot_data <- data.frame(pca_res$x, group = env_info$group)

percentage <- round(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, 2)

ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) + # 添加 95% 置信椭圆
  theme_bw() +
  labs(x = paste0("PC1 (", percentage[1], "%)"),
       y = paste0("PC2 (", percentage[2], "%)"),
       title = "PCA of Metabolome and Host Pathways")

fit_anova <- aov(`L-Leucine` ~ group, data = data)
summary(fit_anova)

TukeyHSD(fit_anova)

ggplot(data, aes(x = group, y = `L-Leucine`, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_classic() +
  labs(title = "L-Leucine Concentration by Group")

micro_for_corr <- micro_log
other_for_corr <- cbind(meta_scaled, host_scaled)

corr_matrix <- cor(micro_for_corr, other_for_corr, method = "spearman")

library(psych)
corr_test <- corr.test(micro_for_corr, other_for_corr, method = "spearman", adjust = "none")
p_matrix <- corr_test$p

add_sig <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
}
text_matrix <- add_sig(p_matrix)

micro_corr_heatmap <-
pheatmap(corr_matrix,
         display_numbers = text_matrix,
         fontsize_number = 12,
         number_color = "black",
         cluster_cols = F, 
         cluster_rows = F,
         color = colorRampPalette(c("lightblue", "white", "pink"))(100),
         main = "Spearman Correlation: Microbiome corr MetaB/HostT")

dev.off()
pdf("function/2_downstream_analysis//39_microbiome_corr.pdf",
    width  = 6,
    height = 6)
print(micro_corr_heatmap)
dev.off()





### genus_abundance_box_plot
g_QVMG01_abundance_box_plot <-
  ggplot(Micro_MetaB_HostT) + 
  geom_boxplot(aes(x = group, y = g_QVMG01, fill = group, color = group),
               alpha = 0.3,
               width = 0.8,
               staplewidth = 0.5,
               linewidth = 0.25,
               show.legend = F,
               outliers = F
  )+
  # geom_point(aes(x = group, y = g_QVMG01), size = 0.3)+
  scale_fill_manual(values = c("Esin" = ggplot2::alpha(Esin_color, 0.75), 
                               "Bger" = ggplot2::alpha(Bger_color, 0.75), 
                               "Pame" = ggplot2::alpha(Pame_color, 0.75), 
                               "Pful" = ggplot2::alpha(Pful_color, 0.75))
  ) +
  scale_color_manual(values = c("Esin" = ggplot2::alpha(Esin_color, 0.75), 
                               "Bger" = ggplot2::alpha(Bger_color, 0.75), 
                               "Pame" = ggplot2::alpha(Pame_color, 0.75), 
                               "Pful" = ggplot2::alpha(Pful_color, 0.75))
  ) +
  geom_signif(aes(x = group, y = g_QVMG01),
              comparisons = list(c("Esin", "Bger"),
                                 c("Esin", "Pame"), 
                                 c("Esin", "Pful")),
              map_signif_level = function(p) {
                ifelse(p < 0.001, "***", 
                       ifelse(p < 0.01, "**", 
                              ifelse(p < 0.05, "*", " ")))
              },
              test = "wilcox.test",
              test.args = list(exact = FALSE), 
              textsize = 4,
              vjust = 0.6,
              tip_length = 0,
              size = 0.25,
              # step_increase = 0.05,
              y_position = c(0.065, 0.070, 0.075, 0.080)
              
  ) +
  coord_cartesian(ylim = c(0.0, 0.09)) + 
  theme_void() +
  theme(
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1.2, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 1, vjust = 1, angle = 45, lineheight = 1.2, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length = unit(2, "pt"),
    panel.border = element_rect(fill = NA, linewidth = 0.5, linetype = "solid", color = "black"),
    strip.switch.pad.wrap = unit(3, "pt"),
    strip.text.x = element_blank(),
    strip.clip = "on",
  )

g_QVMG01_abundance_box_plot


g_Bacteroides_abundance_box_plot <-
  ggplot(data) + 
  geom_boxplot(aes(x = group, y = g_Bacteroides_H, fill = group, color = group),
               alpha = 0.3,
               width = 0.8,
               staplewidth = 0.5,
               linewidth = 0.25,
               show.legend = F,
               outliers = F
  )+
  # geom_point(aes(x = group, y = g_Bacteroides_H), size = 0.3)+
  scale_fill_manual(values = c("Esin" = ggplot2::alpha(Esin_color, 0.75), 
                               "Bger" = ggplot2::alpha(Bger_color, 0.75), 
                               "Pame" = ggplot2::alpha(Pame_color, 0.75), 
                               "Pful" = ggplot2::alpha(Pful_color, 0.75))
  ) +
  scale_color_manual(values = c("Esin" = ggplot2::alpha(Esin_color, 0.75), 
                                "Bger" = ggplot2::alpha(Bger_color, 0.75), 
                                "Pame" = ggplot2::alpha(Pame_color, 0.75), 
                                "Pful" = ggplot2::alpha(Pful_color, 0.75))
  ) +
  geom_signif(aes(x = group, y = g_Bacteroides_H),
              comparisons = list(c("Esin", "Bger"),
                                 c("Esin", "Pame"), 
                                 c("Esin", "Pful")),
              map_signif_level = function(p) {
                ifelse(p < 0.001, "***", 
                       ifelse(p < 0.01, "**", 
                              ifelse(p < 0.05, "*", " ")))
              },
              test = "wilcox.test",
              test.args = list(exact = FALSE), 
              textsize = 4,
              vjust = 0.6,
              tip_length = 0,
              size = 0.25,
              # step_increase = 0.05,
              y_position = c(0.26, 0.27, 0.28, 0.29)
              
  ) +
  coord_cartesian(ylim = c(0.0, 0.30)) + 
  theme_void() +
  theme(
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1.2, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 1, vjust = 1, angle = 45, lineheight = 1.2, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length = unit(2, "pt"),
    panel.border = element_rect(fill = NA, linewidth = 0.5, linetype = "solid", color = "black"),
    strip.switch.pad.wrap = unit(3, "pt"),
    strip.text.x = element_blank(),
    strip.clip = "on",
  )

g_Bacteroides_abundance_box_plot


g_Frigididesulfovibrio_abundance_box_plot <-
  ggplot(data) + 
  geom_boxplot(aes(x = group, y = g_Frigididesulfovibrio, fill = group, color = group),
               alpha = 0.3,
               width = 0.8,
               staplewidth = 0.5,
               linewidth = 0.25,
               show.legend = F,
               outliers = F
  )+
  # geom_point(aes(x = group, y = g_Frigididesulfovibrio), size = 0.3)+
  scale_fill_manual(values = c("Esin" = ggplot2::alpha(Esin_color, 0.75), 
                               "Bger" = ggplot2::alpha(Bger_color, 0.75), 
                               "Pame" = ggplot2::alpha(Pame_color, 0.75), 
                               "Pful" = ggplot2::alpha(Pful_color, 0.75))
  ) +
  scale_color_manual(values = c("Esin" = ggplot2::alpha(Esin_color, 0.75), 
                                "Bger" = ggplot2::alpha(Bger_color, 0.75), 
                                "Pame" = ggplot2::alpha(Pame_color, 0.75), 
                                "Pful" = ggplot2::alpha(Pful_color, 0.75))
  ) +
  geom_signif(aes(x = group, y = g_Frigididesulfovibrio),
              comparisons = list(c("Esin", "Bger"),
                                 c("Esin", "Pame"), 
                                 c("Esin", "Pful")),
              map_signif_level = function(p) {
                ifelse(p < 0.001, "***", 
                       ifelse(p < 0.01, "**", 
                              ifelse(p < 0.05, "*", " ")))
              },
              test = "wilcox.test",
              test.args = list(exact = FALSE), 
              textsize = 4,
              vjust = 0.6,
              tip_length = 0,
              size = 0.25,
              # step_increase = 0.05,
              y_position = c(0.115, 0.120, 0.125, 0.13)
              
  ) +
  coord_cartesian(ylim = c(0.0, 0.14)) +
  theme_void() +
  theme(
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1.2, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 1, vjust = 1, angle = 45, lineheight = 1.2, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length = unit(2, "pt"),
    panel.border = element_rect(fill = NA, linewidth = 0.5, linetype = "solid", color = "black"),
    strip.switch.pad.wrap = unit(3, "pt"),
    strip.text.x = element_blank(),
    strip.clip = "on",
  )

g_Frigididesulfovibrio_abundance_box_plot



g_Dysgonomonas_abundance_box_plot <-
  ggplot(data) + 
  geom_boxplot(aes(x = group, y = g_Dysgonomonas, fill = group, color = group),
               alpha = 0.3,
               width = 0.8,
               staplewidth = 0.5,
               linewidth = 0.25,
               show.legend = F,
               outlier.shape = NA
  )+
  # geom_point(aes(x = group, y = g_Dysgonomonas), size = 0.3)+
  scale_fill_manual(values = c("Esin" = ggplot2::alpha(Esin_color, 0.75), 
                               "Bger" = ggplot2::alpha(Bger_color, 0.75), 
                               "Pame" = ggplot2::alpha(Pame_color, 0.75), 
                               "Pful" = ggplot2::alpha(Pful_color, 0.75))
  ) +
  scale_color_manual(values = c("Esin" = ggplot2::alpha(Esin_color, 0.75), 
                                "Bger" = ggplot2::alpha(Bger_color, 0.75), 
                                "Pame" = ggplot2::alpha(Pame_color, 0.75), 
                                "Pful" = ggplot2::alpha(Pful_color, 0.75))
  ) +
  geom_signif(aes(x = group, y = g_Frigididesulfovibrio),
              comparisons = list(c("Esin", "Bger"),
                                 c("Esin", "Pame"), 
                                 c("Esin", "Pful")),
              map_signif_level = function(p) {
                ifelse(p < 0.001, "***", 
                       ifelse(p < 0.01, "**", 
                              ifelse(p < 0.05, "*", " ")))
              },
              test = "wilcox.test",
              test.args = list(exact = FALSE), 
              textsize = 4,
              vjust = 0.6,
              tip_length = 0,
              size = 0.25,
              # step_increase = 0.05,
              y_position = c(0.175, 0.180, 0.185, 0.190)
              
  ) +
  coord_cartesian(ylim = c(0.0, 0.2)) + 
  theme_void() +
  theme(
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1.2, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 1, vjust = 1, angle = 45, lineheight = 1.2, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length = unit(2, "pt"),
    panel.border = element_rect(fill = NA, linewidth = 0.5, linetype = "solid", color = "black"),
    strip.switch.pad.wrap = unit(3, "pt"),
    strip.text.x = element_blank(),
    strip.clip = "on",
  )

g_Dysgonomonas_abundance_box_plot


library(patchwork)

genus_abundance_box_plot <- 
  g_QVMG01_abundance_box_plot + 
  g_Bacteroides_abundance_box_plot + 
  g_Dysgonomonas_abundance_box_plot +
  g_Frigididesulfovibrio_abundance_box_plot +
  plot_layout(
    ncol = 2,
    nrow = 2,
    widths = c(1, 1),
    heights = c(1,1),
    guides = "collect"
  )

genus_abundance_box_plot

ggsave("function/2_downstream_analysis/40_genus_abundance_box_plot.svg", plot = genus_abundance_box_plot, width = 6, height = 8, units = "cm", dpi = 300)






