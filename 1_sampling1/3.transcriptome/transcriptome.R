rm(list=ls())
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
library(ggside)
library(vegan)
library(pheatmap)
library(ggplot2)
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
library(ggchicklet)
library(ggcorrplot)


# load Arial font
font_import(pattern = "arial")
loadfonts()
windowsFonts()

# Set a color for each group
color1 <- paletteer::paletteer_d("nationalparkcolors::Badlands")
color2 <- paletteer::paletteer_d("PrettyCols::Bright")
color3 <- paletteer::paletteer_d("ggthemes::Classic_10_Medium")
color4 <- paletteer::paletteer_d("MoMAColors::Althoff")

Esin_color <- color1[1]
Bger_color <- color1[4]
Pame_color <- generate_palette(paletteer_d("fishualize::Acanthurus_sohal")[1], modification = "go_lighter",  n_colours = 5, view_palette = TRUE)[2]
Pful_color <- color1[3]

# load transcriptome data
transcriptome_data <- read.table("2_downstream_analysis/1_OG_tpm_salmon_QN.tsv", header = TRUE, sep = "\t", row.names = 1) %>% t()
# load group file
group <- read.table('2_downstream_analysis/2_group.txt', sep = '\t', header =TRUE)
# Merge the group information into the metabolite table
transcriptome_data_grouped <- merge(group, transcriptome_data, by.x = "ID", by.y = "row.names")

# Calculate the mean value of the group
library(dplyr)
transcriptome_mean_abundance <- 
  transcriptome_data_grouped %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "group") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OG")

write.table(transcriptome_mean_abundance, "2_downstream_analysis/3_transcriptome_group_mean_abundance.tsv", row.names = FALSE, quote = FALSE)


### (PCA) Principal Components Analysis
transcriptome_pca <- prcomp(log2(transcriptome_data+1), scale. = TRUE)

# # Extract the coordinates (principal component scores) from the PCA results
transcriptome_pca_coords <- transcriptome_pca$x
transcriptome_pca_result <- merge(transcriptome_pca_coords, group, by.x = "row.names", by.y = "ID")

# # Extract all principal components
transcriptome_pca_scores <- transcriptome_pca$x[, 1:3] 

# Calculate the distance matrix based on PCA scores
transcriptome_distance <- dist(transcriptome_pca_scores)
transcriptome_distance_matrix <- as.matrix(transcriptome_distance)
write.table(transcriptome_distance_matrix, "2_downstream_analysis/4_transcriptome_sample_distance.txt", sep = "\t", quote = FALSE)

# The "pca$sdev" records the eigenvalues of the main sorting axes in the PCA sorting results (dividing each eigenvalue by the total sum of eigenvalues gives the explanatory power of each axis)
transcriptome_pca_eig = sum(pmax(transcriptome_pca$sdev[1:3]), 0)
transcriptome_pca_eig_percent <- round(transcriptome_pca$sdev[1:3]/transcriptome_pca_eig*100, 3) 
transcriptome_pca_eig_percent

# Conduct a permutation multivariate (factorial) variance analysis (PERMANOVA/adonis2)
transcriptome_pca_permanova_result <- adonis2(transcriptome_distance_matrix ~ group, data = group , permutations = 999)
transcriptome_pca_dune_adonis <- paste("R2", " = ", round(transcriptome_pca_permanova_result$R2, 3), "\nP = ", round(transcriptome_pca_permanova_result$`Pr(>F)`, 3))
transcriptome_pca_dune_adonis


### PCA plot

# Calculate the center point of each group
transcriptome_pca_center_points <- 
  transcriptome_pca_result %>%
  group_by(group) %>%
  summarise(mean_wt = mean(PC1),mean_mpg = mean(PC2))

transcriptome_pca_result <- 
  transcriptome_pca_result %>% 
  left_join(transcriptome_pca_center_points, by = "group")

transcriptome_pca_result$group <- factor(transcriptome_pca_result$group, level = c("Esin", "Bger", "Pame", "Pful"))

# Basic scatter plot
transcriptome_pca_p0 <- 
  ggplot(transcriptome_pca_result, 
         aes(x = PC1, y = PC2, color = group)
  ) +
  geom_point(aes(color = group, shape = group), 
             size = 2,
             show.legend = T,
  ) +
  labs(x = paste("PC 1 (", round(transcriptome_pca_eig_percent[1], 2), "%)", sep = ""), 
       y = paste("PC 2 (", round(transcriptome_pca_eig_percent[2], 2), "%)", sep = ""), 
       tag = transcriptome_pca_dune_adonis
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

transcriptome_pca_p0

# Add confidence ellipse + center point connection
transcriptome_pca_p1 <- 
  transcriptome_pca_p0 + 
  stat_ellipse(
    data = transcriptome_pca_result,
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
    data = transcriptome_pca_result,
    aes(x = PC1, y = PC2, xend = mean_wt, yend = mean_mpg, color = group), 
    linetype = "dashed", 
    linewidth = 1, 
    alpha = 0.5,
    show.legend = F)

transcriptome_pca_p1

# Use ggside package to add marginal box plots
library(ggside)
transcriptome_pca_p_final <- 
  transcriptome_pca_p1 +
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

print(transcriptome_pca_p_final)

ggsave("2_downstream_analysis/5_transcriptome_PCA_PC1~PC2.svg", plot = transcriptome_pca_p_final, width = 5, height = 5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/5_transcriptome_PCA_PC1~PC2.tiff", plot = transcriptome_pca_p_final, width = 5, height = 5, units = "cm", dpi = 300)




################################
#### Hierarchical clustering ###
################################

# Hierarchical clustering based on Bray-Curtis distance
transcriptome_spearman <- as.dist(1 -cor(t(transcriptome_data), method = "spearman"))
transcriptome_spearman_matrix <- as.matrix(transcriptome_spearman)
write.table(transcriptome_spearman_matrix, "2_downstream_analysis/6_transcriptome_sample_spearman.txt", sep = "\t", quote = FALSE)

transcriptome_hc <- flashClust(
  transcriptome_spearman,
  method = "average",
  members = NULL)

# Convert the hclust object to a phylo object
library(ape)
transcriptome_hc_tree <- ape::as.phylo(transcriptome_hc)
plot(transcriptome_hc_tree)

# Export Newick format tree file
write.tree(transcriptome_hc_tree, file = "2_downstream_analysis/7_transcriptome_spearman_average_tree.newick")

# Convert the hclust object to a dendrogram object
transcriptome_dend <- as.dendrogram(transcriptome_hc)

# Customize the dendrogram
transcriptome_dend <- transcriptome_dend %>% set("labels_col", c(rep(Esin_color, 6), rep(Bger_color, 6 ), rep(Pful_color ,6), rep(Pame_color, 6)))       # Set the color of the labels
transcriptome_dend <- transcriptome_dend %>% set("branches_k_color",c(rep(Esin_color, 6), rep(Bger_color, 6 ), rep(Pful_color ,6), rep(Pame_color, 6)))  # Set branch color
transcriptome_dend <- transcriptome_dend %>% set("branches_lwd", 4)                                                                                      # Set the branch width
transcriptome_dend <- transcriptome_dend %>% set("branches_lty", 1)                                                                                      # Set the branch line type
transcriptome_dend <- transcriptome_dend %>% set("nodes_pch", 19)                                                                                        # Set the shape of the node point
transcriptome_dend <- transcriptome_dend %>% set("nodes_col", "black")                                                                                   # Set the color of the node point
transcriptome_dend <- transcriptome_dend %>% set("nodes_cex", 0)                                                                                         # Set the size of the node point
transcriptome_dend <- transcriptome_dend %>% set("labels_cex", 1.5)                                                                                      # Set the label size
transcriptome_dend <- transcriptome_dend %>% set("leaves_pch", c(rep(17,6), rep(18,6), rep(15,6), rep(19,6)))                                            # Set the shape of the label point
transcriptome_dend <- transcriptome_dend %>% set("leaves_col",c(rep(Esin_color, 6), rep(Bger_color, 6 ), rep(Pful_color ,6), rep(Pame_color, 6)))        # Set the color of the label point
transcriptome_dend <- transcriptome_dend %>% set("leaves_cex", 2)                                                                                        # Set the size of the label point

# Plot dendrogram
pdf("2_downstream_analysis/8_transcriptome_Hierarchical_Dendrogram.pdf")
plot(transcriptome_dend, main = "transcriptome Hierarchical Dendrogram")
dev.off()


### merge tree 
library(ape)
library(ggtree)
library(patchwork)
library(phytools)

groupInfo <- split(group$ID, group$group)

p_tree <-
  ggtree(groupOTU(transcriptome_hc_tree, groupInfo), 
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

ggsave("2_downstream_analysis/9_transcriptome_Hierarchical_Dendrogram.svg", plot = p_tree, width = 5, height = 5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/9_transcriptome_Hierarchical_Dendrogram.tiff", plot = p_tree, width = 5, height = 5, units = "cm", dpi = 300)


### Hierarchical clustering heatmap
row_annotation <- group[,c(1,2,4)]
row.names(row_annotation) <- row_annotation[, 1]
colnames(row_annotation)[2] <- "Host_group"
colnames(row_annotation)[3] <- "Host_habitat"

col_annotation <- group[,c(1,2,4)]
row.names(col_annotation) <- col_annotation[, 1]
colnames(col_annotation)[2] <- "Host_group"
colnames(col_annotation)[3] <- "Host_habitat"

# "annot_colors" is a list that contains all the factor names and their corresponding colors.
# The names of each vector in the list correspond to the column names in the row_annotation/col_annotation matrix.
annot_colors <- list(Host_group = c(Esin = Esin_color, Bger = Bger_color, Pame = Pame_color, Pful = Pful_color),
                     Host_habitat = c(Esin = "#751C6DFF", other = paletteer::paletteer_d("nationalparkcolors::Acadia")[1]))

library(pheatmap)

min = min(1- transcriptome_spearman_matrix)
max = max(1- transcriptome_spearman_matrix)

pdf("2_downstream_analysis/10_transcriptome Hierarchical heatmap.pdf",
    width  = 10,
    height = 10)

pheatmap(
  as.matrix(1- transcriptome_spearman_matrix),
  color = colorRampPalette(c("lightyellow", "lightblue", "black"))(100),
  ### Global parameter
  main = "Ortho Genes Hierarchical Clustering Heatmap",
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
  legend_labels =c(min = round(min, 1), max = round(max, 1)),
  width	= 15,
  height = 15,
  # Cluster tree and gap parameter
  cluster_rows = transcriptome_hc,
  cluster_cols = transcriptome_hc,
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


### total_OG_heatmap
OG_dataframe <- 
  read.table("1_OG/OG_tpm_salmon_QN.tsv", header = T, row.names = 1) %>%
  t() %>%
  scale() %>%
  t()

col_annotation <- group[,c(1,2,4)]
row.names(col_annotation) <- col_annotation[, 1]
colnames(col_annotation)[2] <- "Host_group"
colnames(col_annotation)[3] <- "Host_habitat"

annot_colors <- list(Host_group = c(Esin = Esin_color, Bger = Bger_color, Pame = Pame_color, Pful = Pful_color),
                     Host_habitat = c(Esin = "#453947FF", other = paletteer::paletteer_d("nationalparkcolors::Acadia")[1]))

pdf("2_downstream_analysis/11_OG_total_heatmap.pdf",
    width  = 10,
    height = 10)

pheatmap(OG_dataframe,
         color = colorRampPalette(c( "#648C16FF", "white","#FF7200FF"))(30),
         border_color = NA,
         cluster_rows = T,
         cluster_cols = T,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "average",
         annotation_names_row = F,
         show_rownames = F,
         annotation_col = row_annotation[, c(2), drop = FALSE],
         annotation_colors = annot_colors)

dev.off()


#########################################
### dimensionality reduction analysis ###
#########################################

### pathway_heatmap
HostT_pathway_dataframe <- 
  read.table("1_OG/transcriptome_KEGG_pathway_ssGSEA_matrix.tsv", header = T, row.names = 1) %>%
  t() %>%
  scale() %>%
  t()

col_annotation <- group[,c(1,2,4)]
row.names(col_annotation) <- col_annotation[, 1]
colnames(col_annotation)[2] <- "Host_group"
colnames(col_annotation)[3] <- "Host_habitat"

annot_colors <- list(Host_group = c(Esin = Esin_color, Bger = Bger_color, Pame = Pame_color, Pful = Pful_color),
                     Host_habitat = c(Esin = "#453947FF", other = paletteer::paletteer_d("nationalparkcolors::Acadia")[1]))

pdf("2_downstream_analysis/12_OG_pathway_heatmap.pdf",
    width  = 10,
    height = 10)

pheatmap(HostT_pathway_dataframe,
         color = colorRampPalette(c( "#648C16FF", "white","#FF7200FF"))(30),
         border_color = NA,
         cluster_rows = T,
         cluster_cols = T,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "average",
         annotation_names_row = F,
         show_rownames = F,
         annotation_col = row_annotation[, c(2), drop = FALSE],
         annotation_colors = annot_colors)

dev.off()


# enriched KEGG pathway heatmap
transcriptome_up_pathway <- read.table("1_OG/1_all_up_pathway.tsv", header = TRUE, row.names = 1 ,sep = "\t") 

transcriptome_up_pathway <- 
  transcriptome_up_pathway %>%
  mutate(log2trans = transcriptome_up_pathway$Directionality.x..log10.p.adj.) %>%
  dplyr::select(log2trans, everything())

transcriptome_up_pathway <- 
  transcriptome_up_pathway %>%
  mutate(pathway = rownames(transcriptome_up_pathway)) %>% 
  dplyr::select(pathway, everything())

up_pathway_mratrix <- 
  transcriptome_up_pathway %>%
  dplyr::select(-c(1:10)) %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "pathway")

up_pathway_mratrix$pathway <- factor(up_pathway_mratrix$pathway, level = rev(up_pathway_mratrix$pathway))

# wide to long
up_pathway_mratrix_long <-
  up_pathway_mratrix %>%
  pivot_longer(cols = -pathway,
               names_to = "sample", 
               values_to = "abundance")

sample_order <- 
  up_pathway_mratrix_long %>%
  distinct(sample) %>%
  mutate(
    group = str_extract(sample, "^[A-Za-z]+"),
    number = as.numeric(str_extract(sample, "\\d+"))
  ) %>%
  arrange(
    match(group, rev(c("Esin", "Bger", "Pame", "Pful"))) ,
    number
  ) %>%
  pull(sample)

up_pathway_mratrix_long$sample <- 
  factor(
    up_pathway_mratrix_long$sample,
    levels = sample_order
  )

min <- min(up_pathway_mratrix_long$abundance)
max <- max(up_pathway_mratrix_long$abundance)

# plot main heatmap
up_pathway_heatmap_p <- 
  ggplot(up_pathway_mratrix_long, 
         aes(x = pathway, y = sample, fill = abundance)
  ) +
  geom_tile(color = NA,
            linewidth = 0,
            lineend = "square",
            height = 1,
            show.legend = T
  ) +
  annotate("rect", ymin = 18.5, ymax = 24.5, xmin = 0.5, xmax = 28.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 12.5, ymax = 18.5, xmin = 28.5, xmax = 50.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 6.5, ymax = 12.5, xmin = 50.5, xmax = 59.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 0.5, ymax = 6.5, xmin = 59.5, xmax = 65.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF", 
                       limit = c(min, max),
                       breaks = c(-2, 0, 2)
  ) +
  scale_x_discrete(position = "top") +
  labs(x = "Enriched transcriptome modules of four species cockroaches", 
       fill = "Enrich Z-score"
  ) +
  theme_void() +
  theme(
    axis.title.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.position = "none",
    legend.key.size = unit(10, "pt"),
    legend.key.width = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title.position = "left",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = -90, lineheight = 1, color = "black"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

up_pathway_heatmap_p

# plot sample grouping annotatin
group$group <- factor(group$group, levels = c("Esin", "Bger", "Pame", "Pful"))

up_pathway_group_p <-
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

up_pathway_group_p

# plot pathway pvalue annotatin
up_pathway_pathway_p <-
  ggplot(transcriptome_up_pathway, 
         aes(x = pathway, y = 0, fill = log2trans , color = log2trans), 
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
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white"),
  )

up_pathway_pathway_p

# plot pathway category annotatin
up_pathway_pathway_info <- read.table("1_OG/1_all_up_pathway_information.tsv", header = TRUE, sep = "\t") %>% as.data.frame()
up_pathway_pathway_info$pathway_name <- substring(up_pathway_pathway_info$pathway_name, 7)
up_pathway_pathway_info$Subcategory <- substring(up_pathway_pathway_info$Subcategory, 7)
up_pathway_pathway_info$Category <- substring(up_pathway_pathway_info$Category, 7)

down_pathway_pathway_info <- read.table("1_OG/2_all_down_pathway_information.tsv", header = TRUE, sep = "\t") %>% as.data.frame()
down_pathway_pathway_info$pathway_name <- substring(down_pathway_pathway_info$pathway_name, 7)
down_pathway_pathway_info$Subcategory <- substring(down_pathway_pathway_info$Subcategory, 7)
down_pathway_pathway_info$Category <- substring(down_pathway_pathway_info$Category, 7)

up_pathway_pathway_info <- 
  up_pathway_pathway_info %>%
  mutate(Pathway.category = case_when(
    up_pathway_pathway_info$Subcategory == "Amino acid metabolism" ~ "AA. MET.",
    up_pathway_pathway_info$Subcategory == "Metabolism of other amino acids" ~ "AA. MET.",
    up_pathway_pathway_info$Subcategory == "Carbohydrate metabolism" ~ "CARB. MET.",
    up_pathway_pathway_info$Subcategory == "Glycan biosynthesis and metabolism" ~ "GLYC. SYN. & MET.",
    up_pathway_pathway_info$Subcategory == "Lipid metabolism" ~ "LIPI. MET.",
    up_pathway_pathway_info$Subcategory == "Signal transduction" ~ "SIGN. TDC.",
    up_pathway_pathway_info$Subcategory == "Metabolism of cofactors and vitamins" ~ "VITA. MET.",
    up_pathway_pathway_info$Subcategory == "Transport and catabolism" ~ "TRNSP. & CATAB.",
    up_pathway_pathway_info$Subcategory == "Translation" ~ "TRANSL.",
    up_pathway_pathway_info$Subcategory == "Signaling molecules and interaction" ~ "SIGN. MOL. INT.",
    TRUE ~ "Other"
  ))

up_pathway_pathway_info$Pathway.category <- 
  factor(up_pathway_pathway_info$Pathway.category,
         levels = c("AA. MET.", "CARB. MET.", "GLYC. SYN. & MET.", "LIPI. MET.", "SIGN. TDC.", 
                    "VITA. MET.", "TRNSP. & CATAB.", "TRANSL.", "SIGN. MOL. INT.", "Other"))


pathway_info_color <- c("AA. MET." = paletteer_d("ggthemes::Classic_Green_Orange_12")[1],
                        "CARB. MET." = paletteer_d("ggthemes::Classic_Green_Orange_12")[2],
                        "GLYC. SYN. & MET." = paletteer_d("ggthemes::Classic_Green_Orange_12")[3],
                        "LIPI. MET." = paletteer_d("ggthemes::Classic_Green_Orange_12")[4],
                        "SIGN. TDC." = paletteer_d("ggthemes::Classic_Green_Orange_12")[5],
                        "VITA. MET." = paletteer_d("ggthemes::Classic_Green_Orange_12")[6],
                        "TRNSP. & CATAB." = paletteer_d("ggthemes::Classic_Green_Orange_12")[7],
                        "TRANSL." = paletteer_d("ggthemes::Classic_Green_Orange_12")[8],
                        "SIGN. MOL. INT." = paletteer_d("ggthemes::Classic_Green_Orange_12")[9],
                        "Other" = "grey"
                        )

up_pathway_pathway_info_p <-
  ggplot(up_pathway_pathway_info, 
         aes(x = pathway, y = 0, fill = Pathway.category , color = Pathway.category), 
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.05,
            lineend = "square",
            show.legend = T
  ) +
  labs(fill = "Pathway category", 
       color = "Pathway category"
  ) +
  guides(fill = guide_legend(ncol = 5, nrow = 2)) + 
  scale_fill_manual(values = c(paletteer_d("ggthemes::Classic_Green_Orange_12")[1:9], "grey"), drop = FALSE) +
  scale_color_manual(values = c(paletteer_d("ggthemes::Classic_Green_Orange_12")[1:9], "grey"), drop = FALSE) +
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
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 2, color = "black")
  )

up_pathway_pathway_info_p

# merge plots 
up_pathway <- up_pathway_heatmap_p %>% 
  aplot::insert_left(up_pathway_group_p, width = 0.02) %>% 
  aplot::insert_top(up_pathway_pathway_p, height = 0.05) %>%
  aplot::insert_bottom(up_pathway_pathway_info_p, height = 0.05)

up_pathway

# save image
ggsave("2_downstream_analysis/13_transcriptome_KEGG_pathway_enrich_pathway.svg", plot = up_pathway, width = 10, height = 4.5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/13_transcriptome_KEGG_pathway_enrich_pathway.tiff", plot = up_pathway, width = 10, height = 4.5, units = "cm", dpi = 300)

# save all legend
up_pathway_heatmap_p_legend <- cowplot::get_legend(up_pathway_heatmap_p + theme(legend.position = "right"))
up_pathway_pathway_p_legend <- cowplot::get_legend(up_pathway_pathway_p + theme(legend.position = "right"))
up_pathway_pathway_info_p_legend <- cowplot::get_legend(up_pathway_pathway_info_p + theme(legend.position = "right"))
up_pathway_combined_legend <- cowplot::plot_grid(up_pathway_heatmap_p_legend, up_pathway_pathway_p_legend, up_pathway_pathway_info_p_legend, nrow = 3)

up_pathway_combined_legend

ggsave("2_downstream_analysis/14_transcriptome_KEGG_pathway_enrich_pathway_legend.svg", plot = up_pathway_combined_legend, width = 20, height = 10, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/14_transcriptome_KEGG_pathway_enrich_pathway_legend.tiff", plot = up_pathway_combined_legend, width = 20, height = 10, units = "cm", dpi = 300)


# depleted KEGG pathways heatmap
transcriptome_down_pathway <- read.table("1_OG/2_all_down_pathway.tsv", header = TRUE, row.names = 1 ,sep = "\t") 

transcriptome_down_pathway <- 
  transcriptome_down_pathway %>%
  mutate(log2trans = transcriptome_down_pathway$Directionality.x..log10.p.adj.) %>%
  dplyr::select(log2trans, everything())

transcriptome_down_pathway <- 
  transcriptome_down_pathway %>%
  mutate(pathway = rownames(transcriptome_down_pathway)) %>%
  dplyr::select(pathway, everything())

down_pathway_mratrix <- 
  transcriptome_down_pathway %>%
  dplyr::select(-c(1:10)) %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "pathway")

down_pathway_mratrix$pathway <- factor(down_pathway_mratrix$pathway, level = rev(down_pathway_mratrix$pathway))

# wide to long
down_pathway_mratrix_long <-
  down_pathway_mratrix %>%
  pivot_longer(cols = -pathway,
               names_to = "sample", 
               values_to = "abundance")

sample_order <- 
  down_pathway_mratrix_long %>%
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

down_pathway_mratrix_long$sample <- 
  factor(
    down_pathway_mratrix_long$sample,
    levels = sample_order
  )

min <- min(down_pathway_mratrix_long$abundance)
max <- max(down_pathway_mratrix_long$abundance)

# plot main heatmap
down_pathway_heatmap_p <- 
  ggplot(down_pathway_mratrix_long, 
         aes(x = pathway, y = sample, fill = abundance)
  ) +
  geom_tile(height = 1,
            color = NA,
            linewidth = 0,
            lineend = "square",
            show.legend = T
  ) +
  annotate("rect", ymin = 0.5, ymax = 6.5, xmin = 0.5 , xmax = 31.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 6.5, ymax = 12.5, xmin = 31.5, xmax = 61.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 12.5, ymax = 18.5, xmin = 61.5, xmax = 72.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 18.5, ymax = 24.5, xmin = 72.5, xmax = 79.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF", 
                       limit = c(min, max),
                       breaks = c(-2,0 ,2)
  ) +
  scale_x_discrete(position = "bottom") + 
  labs(x = "Depleted transcriptome pathways of four species cockroaches", 
       fill = "Pathway Z-score"
  ) +
  theme_void() + 
  theme(
    axis.title.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.position = "none",
    legend.key.size = unit(10, "pt"),
    legend.key.width = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"),
    legend.key.spacing.y = unit(10, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "left",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = -90, lineheight = 1, color = "black"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

down_pathway_heatmap_p

# plot sample grouping annotatin
group$group <- factor(group$group, levels = c("Esin", "Bger", "Pame", "Pful"))

down_pathway_group_p <-
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

down_pathway_group_p

# plot pathway pvalue annotatin
down_pathway_pathway_p <-
  ggplot(transcriptome_down_pathway, 
         aes(x = pathway, y = 0, fill = log2trans , color = log2trans), 
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

down_pathway_pathway_p

# plot pathway category annotatin
down_pathway_pathway_info <- read.table("1_OG/2_all_down_pathway_information.tsv", header = TRUE, sep = "\t") %>% as.data.frame()
down_pathway_pathway_info$pathway_name <- substring(down_pathway_pathway_info$pathway_name, 7)
down_pathway_pathway_info$Subcategory <- substring(down_pathway_pathway_info$Subcategory, 7)
down_pathway_pathway_info$Category <- substring(down_pathway_pathway_info$Category, 7)

down_pathway_pathway_info <- 
  down_pathway_pathway_info %>%
  mutate(Pathway.category = case_when(
    down_pathway_pathway_info$Subcategory == "Amino acid metabolism" ~ "AA. MET.",
    down_pathway_pathway_info$Subcategory == "Metabolism of other amino acids" ~ "AA. MET.",
    down_pathway_pathway_info$Subcategory == "Carbohydrate metabolism" ~ "CARB. MET.",
    down_pathway_pathway_info$Subcategory == "Glycan biosynthesis and metabolism" ~ "GLYC. SYN. & MET.",
    down_pathway_pathway_info$Subcategory == "Lipid metabolism" ~ "LIPI. MET.",
    down_pathway_pathway_info$Subcategory == "Signal transduction" ~ "SIGN. TDC.",
    down_pathway_pathway_info$Subcategory == "Metabolism of cofactors and vitamins" ~ "VITA. MET.",
    down_pathway_pathway_info$Subcategory == "Transport and catabolism" ~ "TRNSP. & CATAB.",
    down_pathway_pathway_info$Subcategory == "Translation" ~ "TRANSL.",
    down_pathway_pathway_info$Subcategory == "Signaling molecules and interaction" ~ "SIGN. MOL. INT.",
    TRUE ~ "Other"
  ))

down_pathway_pathway_info$Pathway.category <- 
  factor(down_pathway_pathway_info$Pathway.category,
         levels = c("AA. MET.", "CARB. MET.", "GLYC. SYN. & MET.", "LIPI. MET.", "SIGN. TDC.", 
                    "VITA. MET.", "TRNSP. & CATAB.", "TRANSL.", "SIGN. MOL. INT.", "Other"))

down_pathway_pathway_info_p <-
  ggplot(down_pathway_pathway_info, 
         aes(x = pathway, y = 0, fill = Pathway.category , color = Pathway.category)
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.05,
            lineend = "square",
            show.legend = T
  ) +
  labs(fill = "Pathway category", 
       color = "Pathway category"
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

down_pathway_pathway_info_p

# merge plots 
down_pathway <- down_pathway_heatmap_p %>% 
  aplot::insert_left(down_pathway_group_p, width = 0.02) %>% 
  aplot::insert_bottom(down_pathway_pathway_p, height = 0.05) %>%
  aplot::insert_top(down_pathway_pathway_info_p, height = 0.05)

down_pathway

# save image
ggsave("2_downstream_analysis/15_transcriptome_KEGG_pathway_deplete_pathway.svg", plot = down_pathway, width = 10, height = 4.5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/15_transcriptome_KEGG_pathway_deplete_pathway.tiff", plot = down_pathway, width = 10, height = 4.5, units = "cm", dpi = 300)

# save all legend
down_pathway_heatmap_p_legend <- cowplot::get_legend(down_pathway_heatmap_p + theme(legend.position = "right"))
down_pathway_pathway_p_legend <- cowplot::get_legend(down_pathway_pathway_p + theme(legend.position = "right"))
down_pathway_pathway_info_p_legend <- cowplot::get_legend(down_pathway_pathway_info_p + theme(legend.position = "right"))
down_pathway_combined_legend <- cowplot::plot_grid(down_pathway_heatmap_p_legend, down_pathway_pathway_p_legend, down_pathway_pathway_info_p_legend, nrow = 3)

down_pathway_combined_legend

ggsave("2_downstream_analysis/16_transcriptome_KEGG_pathway_deplete_pathway_legend.svg", plot = down_pathway_combined_legend, width = 20, height = 10, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/16_transcriptome_KEGG_pathway_deplete_pathway_legend.tiff", plot = down_pathway_combined_legend, width = 20, height = 10, units = "cm", dpi = 300)


### regulation pathway percent
up_pathway_pathway_percent <- 
  transcriptome_up_pathway %>% left_join(up_pathway_pathway_info, by = "pathway") %>% 
  group_by(regulation, Pathway.category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(regulation) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

up_pathway_pathway_percent$Pathway.category <- factor(up_pathway_pathway_percent$Pathway.category,
                                                     levels = rev(c("AA. MET.", "CARB. MET.", "GLYC. SYN. & MET.", "LIPI. MET.", "SIGN. TDC.", 
                                                                "VITA. MET.", "TRNSP. & CATAB.", "TRANSL.", "SIGN. MOL. INT.", "Other")))

up_pathway_pathway_percent$regulation <- factor(up_pathway_pathway_percent$regulation,
                                              levels = c("Esin_up", "Bger_up", "Pame_up", "Pful_up"))

up_pathway_pathway_percent_bar_plot <-
  ggplot(up_pathway_pathway_percent) + 
  geom_bar(
    aes(x = regulation, y = percentage, fill = Pathway.category),
    stat = "identity",
    position = "stack",
    width = 0.7
  ) +
  labs(
    x = "Regulation",
    y = "Enrich Pathway Percentage (%)",
    fill = "Category"
  ) +
  scale_fill_manual(values = c("grey", paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[9:1]), drop =FALSE) +
  scale_color_manual(values = c("grey", paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[9:1]), drop =FALSE) +
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

up_pathway_pathway_percent_bar_plot

up_pathway_pathway_total <- 
  up_pathway_pathway_percent %>%
  group_by(regulation) %>%
  summarise(sum = sum(count), .groups = 'drop')

up_pathway_pathway_total_bar_plot <- 
  ggplot(up_pathway_pathway_total) + 
  geom_bar(aes(x = regulation, y = sum, fill = regulation),
           stat = "identity",
           width = 0.7) +
  geom_label(aes(label = sum, x = regulation, y = 15), size = 3) + 
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

up_pathway_pathway_total_bar_plot

down_pathway_pathway_percent <- 
  transcriptome_down_pathway %>% left_join(down_pathway_pathway_info, by = "pathway") %>% 
  group_by(regulation, Pathway.category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(regulation) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

down_pathway_pathway_percent$Pathway.category <- factor(down_pathway_pathway_percent$Pathway.category, 
                                                        levels = rev(c("AA. MET.", "CARB. MET.", "GLYC. SYN. & MET.", "LIPI. MET.", "SIGN. TDC.", 
                                                                       "VITA. MET.", "TRNSP. & CATAB.", "TRANSL.", "SIGN. MOL. INT.", "Other")))

down_pathway_pathway_percent$regulation <- factor(down_pathway_pathway_percent$regulation,
                                                levels = c("Esin_down", "Bger_down", "Pame_down", "Pful_down"))

down_pathway_pathway_percent_bar_plot <-
  ggplot(down_pathway_pathway_percent) + 
  geom_bar(
    aes(x = regulation, y = percentage, fill = Pathway.category),
    stat = "identity",
    position = "stack",
    width = 0.7
  ) +
  labs(
    x = "Regulation",
    y = "Deplete Pathway Percentage (%)",
    fill = "Category"
  ) +
  scale_fill_manual(values = c("grey", paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[9:1]), drop =FALSE) +
  scale_color_manual(values = c("grey", paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[9:1]), drop =FALSE) +
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

down_pathway_pathway_percent_bar_plot

down_pathway_pathway_total <- 
  down_pathway_pathway_percent %>%
  group_by(regulation) %>%
  summarise(sum = sum(count), .groups = 'drop')

down_pathway_pathway_total_bar_plot <- 
  ggplot(down_pathway_pathway_total) + 
  geom_bar(aes(x = regulation, y = sum, fill = regulation),
           stat = "identity",
           width = 0.7) +
  geom_label(aes(label = sum, x = regulation, y = 15), size = 3) + 
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

down_pathway_pathway_total_bar_plot

pathway_percent_bar_plot <- up_pathway_pathway_total_bar_plot / up_pathway_pathway_percent_bar_plot / down_pathway_pathway_percent_bar_plot / down_pathway_pathway_total_bar_plot + plot_layout(heights = c(0.2, 1, 1, 0.2))

pathway_percent_bar_plot

ggsave("2_downstream_analysis/17_pathway_percent_bar_plot.svg", plot = pathway_percent_bar_plot, width = 6, height = 12, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/17_pathway_percent_bar_plot.tiff", plot = pathway_percent_bar_plot, width = 6, height = 12, units = "cm", dpi = 300)


### highlight top 10 enriched/depleted transcriptome pathways

### top 10 enriched transcriptome pathway
top10_enriched_pathway <- read.table("1_OG/3_each_top10_up_pathway.tsv", header = TRUE, sep = "\t") 
top10_enriched_pathway$log2trans <- ifelse(top10_enriched_pathway$Directionality.x..log10.p.adj. > 6, 6, top10_enriched_pathway$Directionality.x..log10.p.adj.)
top10_enriched_pathway$mark <- ifelse(top10_enriched_pathway$Directionality.x..log10.p.adj. > 6, "over", "not")
top10_enriched_pathway <- top10_enriched_pathway %>% dplyr::select(log2trans, mark, everything())
top10_enriched_pathway$regulation <- factor(top10_enriched_pathway$regulation, levels = c("Esin_up", "Bger_up", "Pame_up", "Pful_up"))
top10_enriched_pathway$pathway <- factor(top10_enriched_pathway$pathway, levels = top10_enriched_pathway$pathway)
top10_enriched_pathway_info <- read.table("1_OG/3_each_top10_up_pathway_information.tsv.txt", header = TRUE, sep = "\t") 

top10_enriched_pathway_info$pathway_name <- substring(top10_enriched_pathway_info$pathway_name, 7)
top10_enriched_pathway_info$Subcategory <- substring(top10_enriched_pathway_info$Subcategory, 7)
top10_enriched_pathway_info$Category <- substring(top10_enriched_pathway_info$Category, 7)

top10_enriched_pathway_info <- 
  top10_enriched_pathway_info %>%
  mutate(Pathway.category = case_when(
    top10_enriched_pathway_info$Subcategory == "Amino acid metabolism" ~ "AA. MET.",
    top10_enriched_pathway_info$Subcategory == "Metabolism of other amino acids" ~ "AA. MET.",
    top10_enriched_pathway_info$Subcategory == "Carbohydrate metabolism" ~ "CARB. MET.",
    top10_enriched_pathway_info$Subcategory == "Glycan biosynthesis and metabolism" ~ "GLYC. SYN. & MET.",
    top10_enriched_pathway_info$Subcategory == "Lipid metabolism" ~ "LIPI. MET.",
    top10_enriched_pathway_info$Subcategory == "Signal transduction" ~ "SIGN. TDC.",
    top10_enriched_pathway_info$Subcategory == "Metabolism of cofactors and vitamins" ~ "VITA. MET.",
    top10_enriched_pathway_info$Subcategory == "Transport and catabolism" ~ "TRNSP. & CATAB.",
    top10_enriched_pathway_info$Subcategory == "Translation" ~ "TRANSL.",
    top10_enriched_pathway_info$Subcategory == "Signaling molecules and interaction" ~ "SIGN. MOL. INT.",
    TRUE ~ "Other"
  ))

top10_enriched_pathway_info$Pathway.category <- 
  factor(top10_enriched_pathway_info$Pathway.category,
         levels = c("AA. MET.", "CARB. MET.", "GLYC. SYN. & MET.", "LIPI. MET.", "SIGN. TDC.", 
                    "VITA. MET.", "TRNSP. & CATAB.", "TRANSL.", "SIGN. MOL. INT.", "Other"))

top10_enriched_pathway_info$pathway_name <- gsub("Apoptosis - multiple species", "Apoptosis", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Glycosylphosphatidylinositol \\(GPI)\\-anchor biosynthesis", "GPI-anchor biosynthesis", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Glycosphingolipid biosynthesis - ganglio series", "Ganglio biosynthesis", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Glycosphingolipid biosynthesis - lacto and neolacto series", "Lacto/Neolacto biosynthesis", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Hippo signaling pathway - multiple species", "Hippo signaling pathway", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Mitophagy - animal", "Mitophagy", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Glycosphingolipid biosynthesis - globo and isoglobo series", "Globo/Isoglobo biosynthesis", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Glyoxylate and dicarboxylate metabolism", "Glyoxylate/Dicarboxylate metabolism", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Pentose and glucuronate interconversions", "Pentose/Glucuronate Interconversions", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Amino sugar and nucleotide sugar metabolism", "Amino/Nucleotide sugar metabolism", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Nicotinate and nicotinamide metabolism", "Nicotinate/Nicotinamide metabolism", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Phosphonate and phosphinate metabolism", "Phosphonate/Phosphinate metabolism", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Virion - Lassa virus and SFTS virus", "Virion(Lassa/SFTS virus)", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Virion - Hepatitis viruses", "Virion(Hepatitis viruses)", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Phosphatidylinositol signaling system", "Phosphatidylinositol signaling", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Cobalamin transport and metabolism", "Cobalamin TRNSP & MET", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("Mannose type O-glycan biosynthesis", "Mannose O-glycan biosynthesis", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("biosynthesis", "SYN.", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("metabolism", "MET.", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("degradation", "DEG.", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("pathway", "PWY.", top10_enriched_pathway_info$pathway_name)
top10_enriched_pathway_info$pathway_name <- gsub("signaling", "SIG.", top10_enriched_pathway_info$pathway_name)

top10_enriched_pathway_info <- 
  top10_enriched_pathway_info %>%
  mutate(merge.name = paste0(top10_enriched_pathway_info$pathway, ":", top10_enriched_pathway_info$pathway_name))

top10_enriched_pathway_info$pathway <- factor(top10_enriched_pathway_info$pathway, levels = top10_enriched_pathway$pathway)

# top10_enriched_pathway bar plot
top10_enriched_pathway_p0 <- 
  ggplot(top10_enriched_pathway, 
         aes(x = log2trans, y = pathway , fill = regulation)
  ) +
  geom_hline(yintercept = seq_along(top10_enriched_pathway$pathway), 
             color = "grey", 
             linetype = "dotted", 
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
                   labels = top10_enriched_pathway_info$merge.name
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
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, 
                               color = c(rep(Pful_color,6), rep(Pame_color,9), rep(Bger_color,10), rep(Esin_color, 10))), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x = unit(2, "pt"),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.25, linetype = "solid"),
    legend.key.height = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.size = unit(10, "pt"),
    legend.key.spacing.y  = unit(5, "pt"),
    legend.key.spacing.x  = unit(5, "pt"),
    legend.direction = "vertical",
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.text.position = "right",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
  )

top10_enriched_pathway_p0

# top10_enriched_pathway pathway category annotation
top10_enriched_pathway_p1 <- 
  ggplot(top10_enriched_pathway_info, 
         aes(x = 0, y = pathway, fill = Pathway.category, color = Pathway.category)
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.25,
            lineend = "square",
            show.legend = T
  ) +  
  scale_fill_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:9],"grey"), drop = FALSE) +
  scale_color_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:9],"grey"), drop = FALSE) +
  labs(fill = "Pathway categary") + 
  guides(fill = guide_legend(ncol =1)) +
  theme_void() + 
  theme(
    legend.key.size = unit(10, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.direction = "vertical",
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.text.position = "right",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
  )

top10_enriched_pathway_p1

# plot heatmap
sample_levels <- c(paste("Esin_",c(1:6), sep = ""), paste("Bger_",c(1:6), sep = ""), paste("Pame_",c(1:6), sep = ""), paste("Pful_",c(1:6), sep = ""))

up_module_z_scores <- top10_enriched_pathway[, -c(1:11)] %>% t() %>% scale() %>% t()

top10_enriched_pathway[, 12:35] = up_module_z_scores

# wide to long
top10_enriched_pathway_long <- top10_enriched_pathway %>%
  pivot_longer(cols = -c(1:11), names_to = "sample", values_to = "Enrich_score")

top10_enriched_pathway_long$sample <- factor(top10_enriched_pathway_long$sample, levels = sample_levels)

min <- min(top10_enriched_pathway_long$Enrich_score)
max <- max(top10_enriched_pathway_long$Enrich_score)

top10_enriched_pathway_p2 <- 
  ggplot(top10_enriched_pathway_long, 
         aes(x = sample, y = pathway, fill = Enrich_score)
  ) +
  geom_tile(height = 1,
            color = NA,
            linewidth = 0,
            lineend = "square",
            show.legend = T
  ) +
  annotate("rect", xmin = 0.5, xmax = 6.5, ymin = 25.5, ymax = 35.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 6.5, xmax = 12.5, ymin = 15.5, ymax = 25.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 12.5, xmax = 18.5, ymin = 6.5, ymax = 15.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 18.5, xmax = 24.5, ymin = 0.5, ymax = 6.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
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
    legend.direction = "horizontal",
    legend.title.position = "top",
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

top10_enriched_pathway_p2

library(aplot)

top10_enriched_pathway_p <- 
  top10_enriched_pathway_p0 %>%
  insert_left(top10_enriched_pathway_p1, width = 0.2) %>%
  insert_left(top10_enriched_pathway_p2, width = 1)

top10_enriched_pathway_p

ggsave("2_downstream_analysis/18_top10_enriched_pathway.svg", plot = top10_enriched_pathway_p, width = 14, height = 9.5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/18_top10_enriched_pathway.tiff", plot = top10_enriched_pathway_p, width = 14, height = 9.5, units = "cm", dpi = 300)


# top 10 depleted KEGG pathways
top10_depleted_pathway <- read.table("1_OG/4_each_top10_down_pathway.tsv", header = TRUE, sep = "\t") 
top10_depleted_pathway$log2trans <- ifelse(top10_depleted_pathway$Directionality.x..log10.p.adj. < -6, -6, top10_depleted_pathway$Directionality.x..log10.p.adj.)
top10_depleted_pathway$mark <- ifelse(top10_depleted_pathway$Directionality.x..log10.p.adj. < -6, "over", "not")
top10_depleted_pathway <- top10_depleted_pathway %>% dplyr::select(log2trans, mark, everything())
top10_depleted_pathway$regulation <- factor(top10_depleted_pathway$regulation, levels = c("Esin_down", "Bger_down", "Pame_down", "Pful_down"))

top10_depleted_pathway$pathway <- factor(top10_depleted_pathway$pathway, levels = top10_depleted_pathway$pathway)
top10_depleted_pathway_info <- read.table("1_OG/4_each_top10_down_pathway_information.tsv", header = TRUE, sep = "\t") 
top10_depleted_pathway_info$pathway_name <- substring(top10_depleted_pathway_info$pathway_name, 7)
top10_depleted_pathway_info$Subcategory <- substring(top10_depleted_pathway_info$Subcategory, 7)
top10_depleted_pathway_info$Category <- substring(top10_depleted_pathway_info$Category, 7)

top10_depleted_pathway_info <- 
  top10_depleted_pathway_info %>%
  mutate(Pathway.category = case_when(
    top10_depleted_pathway_info$Subcategory == "Amino acid metabolism" ~ "AA. MET.",
    top10_depleted_pathway_info$Subcategory == "Metabolism of other amino acids" ~ "AA. MET.",
    top10_depleted_pathway_info$Subcategory == "Carbohydrate metabolism" ~ "CARB. MET.",
    top10_depleted_pathway_info$Subcategory == "Glycan biosynthesis and metabolism" ~ "GLYC. SYN. & MET.",
    top10_depleted_pathway_info$Subcategory == "Lipid metabolism" ~ "LIPI. MET.",
    top10_depleted_pathway_info$Subcategory == "Signal transduction" ~ "SIGN. TDC.",
    top10_depleted_pathway_info$Subcategory == "Metabolism of cofactors and vitamins" ~ "VITA. MET.",
    top10_depleted_pathway_info$Subcategory == "Transport and catabolism" ~ "TRNSP. & CATAB.",
    top10_depleted_pathway_info$Subcategory == "Translation" ~ "TRANSL.",
    top10_depleted_pathway_info$Subcategory == "Signaling molecules and interaction" ~ "SIGN. MOL. INT.",
    TRUE ~ "Other"
  ))

top10_depleted_pathway_info$Pathway.category <- 
  factor(top10_depleted_pathway_info$Pathway.category,
         levels = c("AA. MET.", "CARB. MET.", "GLYC. SYN. & MET.", "LIPI. MET.", "SIGN. TDC.", 
                    "VITA. MET.", "TRNSP. & CATAB.", "TRANSL.", "SIGN. MOL. INT.", "Other"))

top10_depleted_pathway_info$pathway_name <- gsub("Amino sugar and nucleotide sugar metabolism", "Amino/Nucleotide sugar metabolism", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("Glycosaminoglycan biosynthesis - heparan sulfate / heparin", "Heparan-sulfate/Heparin biosynthesis", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("Starch and sucrose metabolism", "Starch/Sucrose metabolism", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("Ribosome biogenesis in eukaryotes", "Ribosome biogenesis", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("MAPK signaling pathway - fly", "MAPK signaling pathway", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("Hippo signaling pathway - fly", "Hippo signaling pathway", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate", "Chondroitin/Dermatan sulfate biosynthesis", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("Arginine and proline metabolism", "Arginine/Proline metabolism", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("Valine, leucine and isoleucine biosynthesis", "Valine/Leucine/Isoleucine biosynthesis", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("Glycosphingolipid biosynthesis - globo and isoglobo series", "Globo/Isoglobo series biosynthesis", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("MAPK signaling pathway - fly", "MAPK signaling pathway", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("Other types of O-glycan biosynthesis", "Other O-glycan biosynthesis", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("Apoptosis - fly", "Apoptosis", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("Various types of N-glycan biosynthesis", "Various N-glycan biosynthesis", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("Mucin type O-glycan biosynthesis", "Mucin O-glycan biosynthesis", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("Cobalamin transport and metabolism", "Cobalamin TRNSP. & MET.", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("metabolism", "MET.", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("degradation", "DEG.", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("biosynthesis", "SYN.", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("pathway", "PWY.", top10_depleted_pathway_info$pathway_name)

top10_depleted_pathway_info <- 
  top10_depleted_pathway_info %>%
  mutate(merge.name = paste0(top10_depleted_pathway_info$pathway, ":", top10_depleted_pathway_info$pathway_name))

top10_depleted_pathway_info$pathway <- factor(top10_depleted_pathway_info$pathway, levels = top10_depleted_pathway$pathway)

# top10_depleted_pathway bar plot
top10_depleted_pathway_p0 <- 
  ggplot(top10_depleted_pathway, 
         aes(x = log2trans, y = pathway , fill = regulation)
  ) +
  geom_hline(yintercept = seq_along(top10_depleted_pathway$pathway), 
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
                   labels = top10_depleted_pathway_info$merge.name
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
                               color = c(rep(Pful_color,7), rep(Pame_color,10), rep(Bger_color,10), rep(Esin_color, 10))), 
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

top10_depleted_pathway_p0

# top10_depleted_pathway pathway category annotation

top10_depleted_pathway_p1 <- 
  ggplot(top10_depleted_pathway_info, 
         aes(x = 0, y = pathway, fill = Pathway.category, color = Pathway.category)
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.25,
            lineend = "square",
            show.legend = T
  ) + 
  scale_fill_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[c(1:9)],"grey"), drop = FALSE) +
  scale_color_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[c(1:9)],"grey"), drop = FALSE) +
  labs(fill = "Pathway categary") + 
  guides(fill = guide_legend(ncol =1)) +
  theme_void() + 
  theme(
    legend.key.size = unit(10, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.direction = "vertical",
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.text.position = "right",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
  )

top10_depleted_pathway_p1

# plot heatmap
sample_levels <- c(paste("Esin_",c(1:6), sep = ""), paste("Bger_",c(1:6), sep = ""), paste("Pame_",c(1:6), sep = ""), paste("Pful_",c(1:6), sep = ""))

down_module_z_scores <- top10_depleted_pathway[, 12:35] %>% t() %>% scale() %>% t()

top10_depleted_pathway[, 12:35] = down_module_z_scores

# wide to long
top10_depleted_pathway_long <- top10_depleted_pathway %>%
  pivot_longer(cols = -c(1:11), names_to = "sample", values_to = "Enrich_score")

top10_depleted_pathway_long$sample <- factor(top10_depleted_pathway_long$sample, levels = rev(sample_levels))

min <- min(top10_depleted_pathway_long$Enrich_score)
max <- max(top10_depleted_pathway_long$Enrich_score)

top10_depleted_pathway_p2 <- 
  ggplot(top10_depleted_pathway_long, 
         aes(x = sample, y = pathway, fill = Enrich_score)
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
  annotate("rect", xmin = 18.5, xmax = 24.5, ymin = 27.5, ymax = 37.5, fill = NA, color = "blue", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 12.5, xmax = 18.5, ymin = 17.5, ymax = 27.5, fill = NA, color = "blue", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 6.5, xmax = 12.5, ymin = 7.5, ymax = 17.5, fill = NA, color = "blue", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 0.5, xmax = 6.5, ymin = 0.5, ymax = 7.5, fill = NA, color = "blue", linewidth = 0.25, lineend = "square")+
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

top10_depleted_pathway_p2

library(aplot)

top10_depleted_pathway_p <- 
  top10_depleted_pathway_p0 %>%
  insert_right(top10_depleted_pathway_p1, width = 0.2) %>%
  insert_right(top10_depleted_pathway_p2, width = 1)

top10_depleted_pathway_p

ggsave("2_downstream_analysis/19_top10_depleted_pathway.svg", plot = top10_depleted_pathway_p, width = 14, height = 9.5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/19_top10_depleted_pathway.tiff", plot = top10_depleted_pathway_p, width = 14, height = 9.5, units = "cm", dpi = 300)


# merge top10_depleted_pathway_p & top10_enriched_pathway_p
top10_delepeted_enriched_pathway_plot <- top10_depleted_pathway_p0 + top10_depleted_pathway_p1 + top10_depleted_pathway_p2 +
                                         top10_enriched_pathway_p2 + top10_enriched_pathway_p1 + top10_enriched_pathway_p0 +  
                                         plot_layout(ncol = 6, widths = c(3, 0.5, 4, 4, 0.5, 3)) & 
                                         theme(panel.spacing = unit(1, "line"))

top10_delepeted_enriched_pathway_plot

top10_delepeted_enriched_pathway_plot_0 <- top10_depleted_pathway_p0 + theme(axis.text.y = element_blank()) + top10_depleted_pathway_p1 + top10_depleted_pathway_p2 +
                                           top10_enriched_pathway_p2 + top10_enriched_pathway_p1 + top10_enriched_pathway_p0 + theme(axis.text.y = element_blank()) +
                                           plot_layout(ncol = 6, widths = c(4, 0.35, 4, 4, 0.35, 4)) & 
                                           theme(panel.spacing = unit(1, "line"))

top10_delepeted_enriched_pathway_plot_0

ggsave("2_downstream_analysis/20_top10_delepeted_enriched_pathway.svg", plot = top10_delepeted_enriched_pathway_plot, width = 20, height = 8, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/20_top10_delepeted_enriched_pathway.tiff", plot = top10_delepeted_enriched_pathway_plot, width = 20, height = 8, units = "cm", dpi = 300)

ggsave("2_downstream_analysis/21_top10_delepeted_enriched_pathway_0.svg", plot = top10_delepeted_enriched_pathway_plot_0, width = 7, height = 8, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/21_top10_delepeted_enriched_pathway_0.tiff", plot = top10_delepeted_enriched_pathway_plot_0, width = 7, height = 8, units = "cm", dpi = 300)

# save all legend 
top10_enriched_pathway_p0_legend <- cowplot::get_legend(top10_enriched_pathway_p0 + theme(legend.position = "right"))
top10_enriched_pathway_p1_legend <- cowplot::get_legend(top10_enriched_pathway_p1 + theme(legend.position = "right"))
top10_enriched_pathway_p2_legend <- cowplot::get_legend(top10_enriched_pathway_p2 + theme(legend.position = "right"))
top10_enriched_pathway_plot_legend <- cowplot::plot_grid(top10_enriched_pathway_p0_legend, top10_enriched_pathway_p1_legend, top10_enriched_pathway_p2_legend, nrow = 3)

top10_depleted_pathway_p0_legend <- cowplot::get_legend(top10_depleted_pathway_p0 + theme(legend.position = "right"))
top10_depleted_pathway_p1_legend <- cowplot::get_legend(top10_depleted_pathway_p1 + theme(legend.position = "right"))
top10_depleted_pathway_p2_legend <- cowplot::get_legend(top10_depleted_pathway_p2 + theme(legend.position = "right"))
top10_depleted_pathway_plot_legend <- cowplot::plot_grid(top10_depleted_pathway_p0_legend, top10_depleted_pathway_p1_legend, top10_depleted_pathway_p2_legend, nrow = 3)

top10_depleted_enriched_pathway_plot_legend <- cowplot::plot_grid(top10_depleted_pathway_plot_legend, top10_enriched_pathway_plot_legend, ncol = 2)

top10_depleted_enriched_pathway_plot_legend 

ggsave("2_downstream_analysis/22_top10_delepeted_enriched_pathway_legend.svg", plot = top10_depleted_enriched_pathway_plot_legend, width = 20, height = 20, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/22_top10_delepeted_enriched_pathway_legend.tiff", plot = top10_depleted_enriched_pathway_plot_legend, width = 20, height = 20, units = "cm", dpi = 300)



########################
### Machine Learning ###
########################

### transcriptome KEGG pathway importance ranking by lasso

# load data
library(glmnet)
set.seed(20241102)
lasso_data <- read.table("1_OG/5_lasso_input.tsv", header = TRUE, row.names = 1, check.names = FALSE)
x <- lasso_data %>% dplyr::select(-1) %>% as.matrix()
y <- factor(lasso_data$group)

# screen alpha(0–1, step 0.05)
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
                      alpha   = 0.5)

lambda_min  <- cv_final$lambda.min
lambda_1se  <- cv_final$lambda.1se
print(tibble(lambda.min = lambda_min, lambda.1se = lambda_1se))

plot(results$alpha, results$cv.error, type = "b", xlab = "Alpha", ylab = "CV Error", main = "CV Error vs Alpha")

pdf("2_downstream_analysis/23_Cross-validation.pdf",
    width  = 6,
    height = 6)

plot(cv_final)

dev.off()

# plot
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

ggsave("2_downstream_analysis/24_cv_plot.svg", plot = cv_plot, width = 6, height = 6, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/24_cv_plot.tiff", plot = cv_plot, width = 6, height = 6, units = "cm", dpi = 300)

# Extract the modules with non-zero importance(lambda.min)
coef_mat <- as.matrix(coef(cv_final, s = "lambda.min"))
active.idx   <- which(coef_mat[, 1] != 0)
active.genes <- rownames(coef_mat)[active.idx] %>% setdiff("(Intercept)")

# export result
coef.df <- coef(cv_final, s = "lambda.min") %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pathway") %>%
  filter(pathway != "(Intercept)") %>%
  mutate(coef = lambda.min) %>%
  mutate(absCoef = abs(coef)) %>%
  filter(absCoef != 0) %>%
  arrange(desc(coef))

write_csv(coef.df, "2_downstream_analysis/25_Lasso_module_ranked_by_coef.csv")

### top 20 lasso selected feature
top20_lasso_transcriptome <- rbind(head(coef.df, 10), tail(coef.df, 10))
top20_lasso_transcriptome$pathway <- factor(top20_lasso_transcriptome$pathway, levels = rev(top20_lasso_transcriptome$pathway))

transcriptome_pathway_info <- rbind(read.table("1_OG/1_all_up_pathway_information.tsv", header = TRUE, sep = "\t"),
                                    read.table("1_OG/2_all_down_pathway_information.tsv", header = TRUE, sep = "\t")) %>%
                                    distinct(pathway, .keep_all = TRUE)

transcriptome_pathway_info$pathway_name <- substring(transcriptome_pathway_info$pathway_name, 7)
transcriptome_pathway_info$Subcategory <- substring(transcriptome_pathway_info$Subcategory, 7)
transcriptome_pathway_info$Category <- substring(transcriptome_pathway_info$Category, 7)

transcriptome_pathway_info$pathway_name <- gsub("Virion - Lassa virus and SFTS virus", "Virion(Lassa/SFTS virus)", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("Mannose type O-glycan biosynthesis", "Mannose O-glycan biosynthesis", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("Virion - Hepatitis viruses", "Virion(Hepatitis viruses)", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("Hippo signaling pathway - multiple species", "Hippo signaling pathway", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("Cobalamin transport and metabolism", "Cobalamin TRNSP. & MET.", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("Phosphatidylinositol signaling system", "Phosphatidylinositol signaling", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("signaling", "SIG.", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("pathway", "PWY.", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("biosynthesis", "SYN.", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("degradation", "DEG.", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("metabolism", "MET.", transcriptome_pathway_info$pathway_name)

top20_lasso_transcriptome <- 
  top20_lasso_transcriptome %>% 
  left_join(transcriptome_pathway_info, by = c("pathway" = "pathway")) %>%
  mutate(label = paste0(pathway, ":", pathway_name))

top20_lasso_transcriptome$label <- factor(top20_lasso_transcriptome$label, levels = rev(top20_lasso_transcriptome$label))

write.table(top20_lasso_transcriptome, "../4.association/1_multi_omic_association/6_top20_e_d_pathway_information_HostT.tsv", row.names = F, quote = F, sep = "\t")

transcriptome_all_pathway_ssGSEA_data <- read.table("1_OG/transcriptome_KEGG_pathway_ssGSEA_matrix.tsv", sep = "\t", header = T, row.names = 1)

top20_lasso_transcriptome_pathway_ssGSEA <- transcriptome_all_pathway_ssGSEA_data[top20_lasso_transcriptome$pathway, ] %>% rownames_to_column(var = "pathway")

write.table(top20_lasso_transcriptome_pathway_ssGSEA, "../4.association/1_multi_omic_association/5_top20_e_d_pathway_HostT.tsv", sep = "\t", quote = F, row.names = F)


### total lasso selected feature
total_lasso_transcriptome <- coef.df
total_lasso_transcriptome$pathway <- factor(total_lasso_transcriptome$pathway, levels = rev(total_lasso_transcriptome$pathway))

transcriptome_pathway_info <- 
  rbind(read.table("1_OG/1_all_up_pathway_information.tsv", header = TRUE, sep = "\t"),
        read.table("1_OG/2_all_down_pathway_information.tsv", header = TRUE, sep = "\t")) %>%
  distinct(pathway, .keep_all = TRUE)

transcriptome_pathway_info$pathway_name <- substring(transcriptome_pathway_info$pathway_name, 7)
transcriptome_pathway_info$Subcategory <- substring(transcriptome_pathway_info$Subcategory, 7)
transcriptome_pathway_info$Category <- substring(transcriptome_pathway_info$Category, 7)

transcriptome_pathway_info$pathway_name <- gsub("Virion - Lassa virus and SFTS virus", "Virion(Lassa/SFTS virus)", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("Mannose type O-glycan biosynthesis", "Mannose O-glycan biosynthesis", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("Virion - Hepatitis viruses", "Virion(Hepatitis viruses)", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("Hippo signaling pathway - multiple species", "Hippo signaling pathway", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("Cobalamin transport and metabolism", "Cobalamin TRNSP. & MET.", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("Phosphatidylinositol signaling system", "Phosphatidylinositol signaling", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("signaling", "SIG.", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("pathway", "PWY.", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("biosynthesis", "SYN.", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("degradation", "DEG.", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("metabolism", "MET.", transcriptome_pathway_info$pathway_name)

total_lasso_transcriptome <- 
  total_lasso_transcriptome %>% 
  left_join(transcriptome_pathway_info, by = c("pathway" = "pathway")) %>%
  mutate(label = paste0(pathway, ":", pathway_name))

total_lasso_transcriptome$label <- factor(total_lasso_transcriptome$label, levels = rev(total_lasso_transcriptome$label))

write.table(total_lasso_transcriptome, "../4.association/1_multi_omic_association/55_total_e_d_pathway_information_HostT.tsv", row.names = F, quote = F, sep = "\t")

transcriptome_all_pathway_ssGSEA_data <- read.table("1_OG/transcriptome_KEGG_pathway_ssGSEA_matrix.tsv", sep = "\t", header = T, row.names = 1)

total_lasso_transcriptome_pathway_ssGSEA <- transcriptome_all_pathway_ssGSEA_data[total_lasso_transcriptome$pathway, ] %>% rownames_to_column(var = "pathway")

write.table(total_lasso_transcriptome_pathway_ssGSEA, "../4.association/1_multi_omic_association/54_total_e_d_pathway_HostT.tsv", sep = "\t", quote = F, row.names = F)


### total lasso selected feature
lasso_importance_rank <- 
  ggplot(top20_lasso_transcriptome,
         aes(x = coef, y = pathway)
  ) +
  geom_point(aes(x = coef, y = label, color = label), 
             size = 1,
             show.legend = F
  ) +
  geom_hline(yintercept = seq_along(top20_lasso_transcriptome$label), 
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
  scale_color_manual(values = c(rep("lightblue", 10), rep("pink", 11))) +
  scale_y_discrete(position = "right") +
  scale_x_continuous(limits = c(-8, 8), breaks = c(-8, 0, 8), labels = c(-8, 0, 8)) +
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

ggsave("2_downstream_analysis/26_lasso_importance_rank_plot.svg", plot = lasso_importance_rank, width = 10, height = 4.5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/26_lasso_importance_rank_plot.tiff", plot = lasso_importance_rank, width = 10, height = 4.5, units = "cm", dpi = 300)

ggsave("2_downstream_analysis/27_lasso_importance_rank_plot_0.svg", plot = lasso_importance_rank_0, width = 1, height = 4.5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/27_lasso_importance_rank_plot_0.tiff", plot = lasso_importance_rank_0, width = 1, height = 4.5, units = "cm", dpi = 300)


### double volcano plot

# read file
double_volcano_data <- read.table("1_OG/6_double_volcano_log2trans.tsv", header = TRUE)

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
quadrant_counts$quadrant <- factor(quadrant_counts$quadrant, levels = rev(c("3", "7", "2", "8", "4", "6", "1", "9", "5")))

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
  scale_fill_manual(values = c("3" = "pink", "7" = "lightblue", "2" = scater_color[3], "8" =  scater_color[4], "4" = scater_color[9], "6" =  scater_color[10], "1" =  scater_color[1], "9" =  scater_color[2], "5" =  "grey80"),
                    labels = c("1" = "Peri_enrich & Bger_deplete",
                               "2" = "Peri_enrich & Bger_nosig",
                               "3" = "Peri_enrich & Bger_enrich",
                               "4" = "Peri_nosig & Bger_deplete",
                               "5" = "Both_nosig",
                               "6" = "Peri_nosig & Bger_enrich",
                               "7" = "Peri_deplete & Bger_deplete",
                               "8" = "Peri_deplete & Bger_nosig",
                               "9" = "Peri_deplete & Bger_enrich")
  ) +
  theme_void()

quadrant_percent_pie_plot

ggsave("2_downstream_analysis/28_quadrant_percent_pie_plot.svg", plot = quadrant_percent_pie_plot, width = 1, height = 1, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/28_quadrant_percent_pie_plot.tiff", plot = quadrant_percent_pie_plot, width = 1, height = 1, units = "cm", dpi = 300)


# highlight module
highlight_module <- c(
  ko03050 = "ko03050",
  ko03010 = "ko03010", 
  ko03018 = "ko03018",
  ko00190 = "ko00190",
  ko04150 = "ko04150",
  ko00052 = "ko00052",
  ko03060 = "ko03060",
  ko00220 = "ko00220",
  ko04122 = "ko04122",
  ko04980 = "ko04980",
  ko00515 = "ko00515",
  ko00340 = "ko00340",
  ko04082 = "ko04082",
  ko04068 = "ko04068",
  ko00562 = "ko00562",
  ko03273 = "ko03273",
  ko03272 = "ko03272",
  ko04392 = "ko04392",
  ko04081 = "ko04081",
  ko04070 = "ko04070"
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
  scale_color_manual(values = c("3" = "pink", "7" = "lightblue", "2" = scater_color[3], "8" =  scater_color[4], "4" = scater_color[9], "6" =  scater_color[10], "1" =  scater_color[1], "9" =  scater_color[2], "5" =  "grey80")) +
  scale_fill_manual(values = c("3" = "pink", "7" = "lightblue", "2" = scater_color[3], "8" =  scater_color[4], "4" = scater_color[9], "6" =  scater_color[10], "1" =  scater_color[1], "9" =  scater_color[2], "5" =  "grey80")) +
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

ggsave("2_downstream_analysis/29_double_volcano_plot.svg", plot = double_volcano_plot, width = 4.5, height = 4.5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/29_double_volcano_plot.tiff", plot = double_volcano_plot, width = 4.5, height = 4.5, units = "cm", dpi = 300)




### venn

# set color 
Esin_venn_color <- generate_palette(Esin_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, view_labels = FALSE)
Bger_venn_color <- generate_palette(Bger_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, view_labels = FALSE)
Pame_venn_color <- generate_palette(Pame_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, view_labels = FALSE)
Pful_venn_color <- generate_palette(Pful_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, view_labels = FALSE)

# Esin enrich
Esin_venn_dat_up  <- read.delim("1_OG/Esin_transcriptome_KEGG_ssGSEA_pathway_each_up.tsv")
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
  set_name_color = "white",
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
)

Esin_venn_up

# Bger enrich
Bger_venn_dat_up  <- read.delim("1_OG/Bger_transcriptome_KEGG_ssGSEA_pathway_each_up.tsv")
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
  set_name_color = "white",
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
)

Bger_venn_up

# Pame enrich
Pame_venn_dat_up  <- read.delim("1_OG/Pame_transcriptome_KEGG_ssGSEA_pathway_each_up.tsv")
Pame_venn_up_list <- as.list(Pame_venn_dat_up )

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
  set_name_color = "white",
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
)

Pame_venn_up

# Pful enrich
Pful_venn_dat_up  <- read.delim("1_OG/Pful_transcriptome_KEGG_ssGSEA_pathway_each_up.tsv")
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
  set_name_color = "white",
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
)

Pful_venn_up

library(aplot)

venn_up <-Esin_venn_up + Bger_venn_up + Pame_venn_up + Pful_venn_up + plot_layout(ncol =4)
venn_up

ggsave("2_downstream_analysis/30_venn_up.svg", plot = venn_up, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/30_venn_up.tiff", plot = venn_up, width = 10, height = 5, units = "cm", dpi = 300)


# Esin deplete
Esin_venn_dat_down  <- read.delim("1_OG/Esin_transcriptome_KEGG_ssGSEA_pathway_each_down.tsv")
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
  set_name_color = "white",
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
)

Esin_venn_down

# Bger deplete
Bger_venn_dat_down  <- read.delim("1_OG/Bger_transcriptome_KEGG_ssGSEA_pathway_each_down.tsv")
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
  set_name_color = "white",
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
)

Bger_venn_down

# Pame deplete
Pame_venn_dat_down  <- read.delim("1_OG/Pame_transcriptome_KEGG_ssGSEA_pathway_each_down.tsv")
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
  set_name_color = "white",
  set_name_size = 0, 
  text_color = "black",
  text_size = 6/2.83
)

Pame_venn_down

# Pful deplete
Pful_venn_dat_down  <- read.delim("1_OG/Pful_transcriptome_KEGG_ssGSEA_pathway_each_down.tsv")
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
  set_name_color = "white",
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
)

Pful_venn_down

library(aplot)

venn_down <- Esin_venn_down + Bger_venn_down  + Pame_venn_down + Pful_venn_down + plot_layout(ncol =4)
venn_down

ggsave("2_downstream_analysis/31_venn_down.svg", plot = venn_down, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/31_venn_down.tiff", plot = venn_down, width = 10, height = 5, units = "cm", dpi = 300)

# plot regulation statistics bar plot
regulation <- read.table("1_OG/regulation_stat.txt",header = TRUE)
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
  scale_fill_manual(values = c(
    generate_palette(Esin_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, , view_labels = FALSE), 
    generate_palette(Bger_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, , view_labels = FALSE),
    generate_palette(Pame_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, , view_labels = FALSE),
    generate_palette(Pful_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, , view_labels = FALSE))
  ) + 
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

ggsave("2_downstream_analysis/32_regulation_stat.svg", plot = regulation_stat, width = 10, height = 10, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/32_regulation_stat.tiff", plot = regulation_stat, width = 10, height = 10, units = "cm", dpi = 300)

# merge plot
merge_stat_plot <- venn_up / regulation_stat / venn_down + plot_layout(nrow = 3)
print(merge_stat_plot)

ggsave("2_downstream_analysis/33_merge_stat_plot.svg", plot = merge_stat_plot, width = 10, height = 10, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/33_merge_stat_plot.tiff", plot = merge_stat_plot, width = 10, height = 10, units = "cm", dpi = 300)













### function radar plot
library(ggradar2)

transcriptome_pathway_ssGSEA_data <- read.table("1_OG/transcriptome_KEGG_pathway_ssGSEA_matrix.tsv", header = T, row.names = 1)

transcriptome_pathway_info <- read_xlsx("1_OG/pathway_information.xlsx")

transcriptome_pathway_ssGSEA_radar_data <-
  transcriptome_pathway_ssGSEA_data %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(everything(), ~ scales::rescale(.x, to = c(0, 2)))) %>%
  mutate(group = c(rep("Other", 6), rep("Esin", 6), rep("Other", 6), rep("Other", 6))) %>%
  group_by(group) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "group") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "pathway") %>%
  left_join(., transcriptome_pathway_info, by = "pathway")

Carbohydrate_metabolism_radar <- 
  transcriptome_pathway_ssGSEA_radar_data %>%
  filter(Subcategory == "Carbohydrate metabolism") %>%
  mutate(proportion = Esin / Other) %>%
  arrange(desc(proportion)) %>%
  mutate(pathway_name = factor(pathway_name, levels = pathway_name)) %>%
  column_to_rownames(var = "pathway_name") %>%
  dplyr::select(-c(1,4:7)) %>% 
  t() %>% 
  as.data.frame() %>%
  mutate(facet1 = c("Esin", "Other"))

Carbohydrate_metabolism_radar_plot <- 
  ggradar2(plot.data = Carbohydrate_metabolism_radar %>% dplyr::select(-facet1),
           base.size = 1,
           grid.min = 0,
           grid.max = 1.1, 
           centre.y = 0.25,
           grid.line.width = 0.5,
           gridline.min.linetype = "solid",
           gridline.mid.linetype = "longdash",
           gridline.max.linetype = "solid",
           gridline.min.colour = "black", 
           gridline.mid.colour= "grey",
           gridline.max.colour = "black",
           grid.label.size = 2,
           gridline.label = c(0, 50, 100),
           gridline.label.type =  "numeric",
           axis.label.offset = 1.1,
           axis.label.size = 2,
           axis.line.colour = "grey",
           group.line.width = 0.5,
           group.point.size = 2,
           plot.legend = F,
           polygonfill = F,
           radarshape = "round"
           ) +
  coord_cartesian(clip="off") + 
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("Esin" = Esin_color, "Other" = "orange"))

Carbohydrate_metabolism_radar_plot

ggsave("2_downstream_analysis/34_Carbohydrate_metabolism_radar_plot.svg", plot = Carbohydrate_metabolism_radar_plot, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/34_Carbohydrate_metabolism_radar_plot.tiff", plot = Carbohydrate_metabolism_radar_plot, width = 10, height = 5, units = "cm", dpi = 1200)


Lipid_metabolism_radar <- 
  transcriptome_pathway_ssGSEA_radar_data %>%
  filter(Subcategory == "Lipid metabolism") %>%
  mutate(proportion = Esin / Other) %>%
  arrange(desc(proportion)) %>%
  mutate(pathway_name = factor(pathway_name, levels = pathway_name)) %>%
  column_to_rownames(var = "pathway_name") %>%
  dplyr::select(-c(1,4:7)) %>% 
  t() %>% 
  as.data.frame() %>%
  mutate(facet1 = c("Esin", "Other"))

Lipid_metabolism_radar_plot <- 
  ggradar2(plot.data = Lipid_metabolism_radar %>% dplyr::select(-facet1),
           base.size = 1,
           grid.min = 0,
           grid.max = 1.1, 
           centre.y = 0.25,
           grid.line.width = 0.5,
           gridline.min.linetype = "solid",
           gridline.mid.linetype = "longdash",
           gridline.max.linetype = "solid",
           gridline.min.colour = "black", 
           gridline.mid.colour= "grey",
           gridline.max.colour = "black",
           grid.label.size = 2,
           gridline.label = c(0, 50, 100),
           gridline.label.type =  "numeric",
           axis.label.offset = 1.1,
           axis.label.size = 2,
           axis.line.colour = "grey",
           group.line.width = 0.5,
           group.point.size = 2,
           plot.legend = F,
           polygonfill = F,
           radarshape = "round"
  ) +
  coord_cartesian(clip="off") + 
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("Esin" = Esin_color, "Other" = "orange"))

Lipid_metabolism_radar_plot

ggsave("2_downstream_analysis/35_Lipid_metabolism_radar_plot.svg", plot = Lipid_metabolism_radar_plot, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/35_Lipid_metabolism_radar_plot.tiff", plot = Lipid_metabolism_radar_plot, width = 10, height = 5, units = "cm", dpi = 1200)


Amino_acid_metabolism_radar <- 
  transcriptome_pathway_ssGSEA_radar_data %>%
  filter(Subcategory == "Amino acid metabolism") %>%
  mutate(proportion = Esin / Other) %>%
  arrange(desc(proportion)) %>%
  mutate(pathway_name = factor(pathway_name, levels = pathway_name)) %>%
  column_to_rownames(var = "pathway_name") %>%
  dplyr::select(-c(1,4:7)) %>% 
  t() %>% 
  as.data.frame() %>%
  mutate(facet1 = c("Esin", "Other"))

Amino_acid_metabolism_radar_plot <- 
  ggradar2(plot.data = Amino_acid_metabolism_radar %>% dplyr::select(-facet1),
           base.size = 1,
           grid.min = 0,
           grid.max = 1.1, 
           centre.y = 0.25,
           grid.line.width = 0.5,
           gridline.min.linetype = "solid",
           gridline.mid.linetype = "longdash",
           gridline.max.linetype = "solid",
           gridline.min.colour = "black", 
           gridline.mid.colour= "grey",
           gridline.max.colour = "black",
           grid.label.size = 2,
           gridline.label = c(0, 50, 100),
           gridline.label.type =  "numeric",
           axis.label.offset = 1.1,
           axis.label.size = 2,
           axis.line.colour = "grey",
           group.line.width = 0.5,
           group.point.size = 2,
           plot.legend = F,
           polygonfill = F,
           radarshape = "round"
  ) +
  coord_cartesian(clip="off") + 
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("Esin" = Esin_color, "Other" = "orange"))

Amino_acid_metabolism_radar_plot

ggsave("2_downstream_analysis/36_Amino_acid_metabolism_radar_plot.svg", plot = Amino_acid_metabolism_radar_plot, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/36_Amino_acid_metabolism_radar_plot.tiff", plot = Amino_acid_metabolism_radar_plot, width = 10, height = 5, units = "cm", dpi = 1200)


Signal_transduction_radar <- 
  transcriptome_pathway_ssGSEA_radar_data %>%
  filter(Subcategory %in% c("Signal transduction", "Signaling molecules and interaction")) %>%
  mutate(proportion = Esin / Other) %>%
  arrange(desc(proportion)) %>%
  mutate(pathway_name = factor(pathway_name, levels = pathway_name)) %>%
  column_to_rownames(var = "pathway_name") %>%
  dplyr::select(-c(1,4:7)) %>% 
  t() %>% 
  as.data.frame() %>%
  mutate(facet1 = c("Esin", "Other"))

Signal_transduction_radar_plot <- 
  ggradar2(plot.data = Signal_transduction_radar %>% dplyr::select(-facet1),
           base.size = 1,
           grid.min = 0,
           grid.max = 1.1, 
           centre.y = 0.25,
           grid.line.width = 0.5,
           gridline.min.linetype = "solid",
           gridline.mid.linetype = "longdash",
           gridline.max.linetype = "solid",
           gridline.min.colour = "black", 
           gridline.mid.colour= "grey",
           gridline.max.colour = "black",
           grid.label.size = 2,
           gridline.label = c(0, 50, 100),
           gridline.label.type =  "numeric",
           axis.label.offset = 1.1,
           axis.label.size = 2,
           axis.line.colour = "grey",
           group.line.width = 0.5,
           group.point.size = 2,
           plot.legend = F,
           polygonfill = F,
           radarshape = "round"
  ) +
  coord_cartesian(clip="off") + 
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("Esin" = Esin_color, "Other" = "orange"))

Signal_transduction_radar_plot

ggsave("2_downstream_analysis/37_Signal_transduction_radar_plot.svg", plot = Signal_transduction_radar_plot, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/37_Signal_transduction_radar_plot.tiff", plot = Signal_transduction_radar_plot, width = 10, height = 5, units = "cm", dpi = 1200)


Glycan_biosynthesis_and_metabolism_radar <- 
  transcriptome_pathway_ssGSEA_radar_data %>%
  filter(Subcategory == "Glycan biosynthesis and metabolism") %>%
  mutate(proportion = Esin / Other) %>%
  arrange(desc(proportion)) %>%
  mutate(pathway_name = factor(pathway_name, levels = pathway_name)) %>%
  column_to_rownames(var = "pathway_name") %>%
  dplyr::select(-c(1,4:7)) %>% 
  t() %>% 
  as.data.frame() %>%
  mutate(facet1 = c("Esin", "Other"))

Glycan_biosynthesis_and_metabolism_radar_plot <- 
  ggradar2(plot.data = Glycan_biosynthesis_and_metabolism_radar %>% dplyr::select(-facet1),
           base.size = 1,
           grid.min = 0,
           grid.max = 1.1, 
           centre.y = 0.25,
           grid.line.width = 0.5,
           gridline.min.linetype = "solid",
           gridline.mid.linetype = "longdash",
           gridline.max.linetype = "solid",
           gridline.min.colour = "black", 
           gridline.mid.colour= "grey",
           gridline.max.colour = "black",
           grid.label.size = 2,
           gridline.label = c(0, 50, 100),
           gridline.label.type =  "numeric",
           axis.label.offset = 1.1,
           axis.label.size = 2,
           axis.line.colour = "grey",
           group.line.width = 0.5,
           group.point.size = 2,
           plot.legend = F,
           polygonfill = F,
           radarshape = "round"
  ) +
  coord_cartesian(clip="off") + 
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("Esin" = Esin_color, "Other" = "orange"))

Glycan_biosynthesis_and_metabolism_radar_plot

ggsave("2_downstream_analysis/38_Glycan_biosynthesis_and_metabolism_radar_plot.svg", plot = Glycan_biosynthesis_and_metabolism_radar_plot, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/38_Glycan_biosynthesis_and_metabolism_radar_plot.tiff", plot = Glycan_biosynthesis_and_metabolism_radar_plot, width = 10, height = 5, units = "cm", dpi = 1200)


Transcription_and_Translation_radar <- 
  transcriptome_pathway_ssGSEA_radar_data %>%
  filter(Subcategory %in% c("Transcription", "Translation")) %>%
  mutate(proportion = Esin / Other) %>%
  arrange(desc(proportion)) %>%
  mutate(pathway_name = factor(pathway_name, levels = pathway_name)) %>%
  column_to_rownames(var = "pathway_name") %>%
  dplyr::select(-c(1,4:7)) %>% 
  t() %>% 
  as.data.frame() %>%
  mutate(facet1 = c("Esin", "Other"))

Transcription_and_Translation_radar_plot <- 
  ggradar2(plot.data = Transcription_and_Translation_radar %>% dplyr::select(-facet1),
           base.size = 1,
           grid.min = 0,
           grid.max = 1.1, 
           centre.y = 0.25,
           grid.line.width = 0.5,
           gridline.min.linetype = "solid",
           gridline.mid.linetype = "longdash",
           gridline.max.linetype = "solid",
           gridline.min.colour = "black", 
           gridline.mid.colour= "grey",
           gridline.max.colour = "black",
           grid.label.size = 2,
           gridline.label = c(0, 50, 100),
           gridline.label.type =  "numeric",
           axis.label.offset = 1.1,
           axis.label.size = 2,
           axis.line.colour = "grey",
           group.line.width = 0.5,
           group.point.size = 2,
           plot.legend = F,
           polygonfill = F,
           radarshape = "round"
  ) +
  coord_cartesian(clip="off") + 
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("Esin" = Esin_color, "Other" = "orange"))

Transcription_and_Translation_radar_plot

ggsave("2_downstream_analysis/39_Transcription_and_Translation_radar_plot.svg", plot = Transcription_and_Translation_radar_plot, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/39_Transcription_and_Translation_radar_plot.tiff", plot = Transcription_and_Translation_radar_plot, width = 10, height = 5, units = "cm", dpi = 1200)


Cellular_Processes_radar <- 
  transcriptome_pathway_ssGSEA_radar_data %>%
  filter(Category %in% c("Cellular Processes")) %>%
  mutate(proportion = Esin / Other) %>%
  arrange(desc(proportion)) %>%
  mutate(pathway_name = factor(pathway_name, levels = pathway_name)) %>%
  column_to_rownames(var = "pathway_name") %>%
  dplyr::select(-c(1,4:7)) %>% 
  t() %>% 
  as.data.frame() %>%
  mutate(facet1 = c("Esin", "Other"))

Cellular_Processes_radar_plot <- 
  ggradar2(plot.data = Cellular_Processes_radar %>% dplyr::select(-facet1),
           base.size = 1,
           grid.min = 0,
           grid.max = 1.1, 
           centre.y = 0.25,
           grid.line.width = 0.5,
           gridline.min.linetype = "solid",
           gridline.mid.linetype = "longdash",
           gridline.max.linetype = "solid",
           gridline.min.colour = "black", 
           gridline.mid.colour= "grey",
           gridline.max.colour = "black",
           grid.label.size = 2,
           gridline.label = c(0, 50, 100),
           gridline.label.type =  "numeric",
           axis.label.offset = 1.1,
           axis.label.size = 2,
           axis.line.colour = "grey",
           group.line.width = 0.5,
           group.point.size = 2,
           plot.legend = F,
           polygonfill = F,
           radarshape = "round"
  ) +
  coord_cartesian(clip="off") + 
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("Esin" = Esin_color, "Other" = "orange"))

Cellular_Processes_radar_plot

ggsave("2_downstream_analysis/40_Cellular_Processes_radar_plot.svg", plot = Cellular_Processes_radar_plot, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/40_Cellular_Processes_radar_plot.tiff", plot = Cellular_Processes_radar_plot, width = 10, height = 5, units = "cm", dpi = 1200)


Organismal_Systems_radar <- 
  transcriptome_pathway_ssGSEA_radar_data %>%
  filter(Category %in% c("Organismal Systems")) %>%
  mutate(proportion = Esin / Other) %>%
  arrange(desc(proportion)) %>%
  mutate(pathway_name = factor(pathway_name, levels = pathway_name)) %>%
  column_to_rownames(var = "pathway_name") %>%
  dplyr::select(-c(1,4:7)) %>% 
  t() %>% 
  as.data.frame() %>%
  mutate(facet1 = c("Esin", "Other"))

Organismal_Systems_radar_plot <- 
  ggradar2(plot.data = Organismal_Systems_radar %>% dplyr::select(-facet1),
           base.size = 1,
           grid.min = 0,
           grid.max = 1.1, 
           centre.y = 0.25,
           grid.line.width = 0.5,
           gridline.min.linetype = "solid",
           gridline.mid.linetype = "longdash",
           gridline.max.linetype = "solid",
           gridline.min.colour = "black", 
           gridline.mid.colour= "grey",
           gridline.max.colour = "black",
           grid.label.size = 2,
           gridline.label = c(0, 50, 100),
           gridline.label.type =  "numeric",
           axis.label.offset = 1.1,
           axis.label.size = 2,
           axis.line.colour = "grey",
           group.line.width = 0.5,
           group.point.size = 2,
           plot.legend = F,
           polygonfill = F,
           radarshape = "round"
  ) +
  coord_cartesian(clip="off") + 
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("Esin" = Esin_color, "Other" = "orange"))

Organismal_Systems_radar_plot

ggsave("2_downstream_analysis/41_Organismal_Systems_radar_plot.svg", plot = Organismal_Systems_radar_plot, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/41_Organismal_Systems_radar_plot.tiff", plot = Organismal_Systems_radar_plot, width = 10, height = 5, units = "cm", dpi = 1200)


Other_radar <- 
  transcriptome_pathway_ssGSEA_radar_data %>%
  filter(Subcategory %in% c("Energy metabolism", "Folding, sorting and degradation")) %>%
  mutate(proportion = Esin / Other) %>%
  arrange(desc(proportion)) %>%
  mutate(pathway_name = factor(pathway_name, levels = pathway_name)) %>%
  column_to_rownames(var = "pathway_name") %>%
  dplyr::select(-c(1,4:7)) %>%
  t() %>% 
  as.data.frame() %>%
  mutate(facet1 = c("Esin", "Other"))

Other_radar_plot <- 
  ggradar2(plot.data = Other_radar %>% dplyr::select(-facet1),
           base.size = 1,
           grid.min = 0,
           grid.max = 1.1, 
           centre.y = 0.25,
           grid.line.width = 0.5,
           gridline.min.linetype = "solid",
           gridline.mid.linetype = "longdash",
           gridline.max.linetype = "solid",
           gridline.min.colour = "black", 
           gridline.mid.colour= "grey",
           gridline.max.colour = "black",
           grid.label.size = 2,
           gridline.label = c(0, 50, 100),
           gridline.label.type =  "numeric",
           axis.label.offset = 1.1,
           axis.label.size = 2,
           axis.line.colour = "grey",
           group.line.width = 0.5,
           group.point.size = 2,
           plot.legend = F,
           polygonfill = F,
           radarshape = "round"
  ) +
  coord_cartesian(clip="off") + 
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("Esin" = Esin_color, "Other" = "orange"))

Other_radar_plot

ggsave("2_downstream_analysis/42_Other_radar_plot.svg", plot = Other_radar_plot, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/42_Other_radar_plot.tiff", plot = Other_radar_plot, width = 10, height = 5, units = "cm", dpi = 1200)



















