#############################################
### Downstream analysis and Visualization ###
#############################################

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

# load Arial font
font_import(pattern = "arial")
loadfonts()
windowsFonts()

# Set color for each group
color1 <- paletteer::paletteer_d("nationalparkcolors::Badlands")
color2 <- paletteer::paletteer_d("PrettyCols::Bright")
color3 <- paletteer::paletteer_d("ggthemes::Classic_10_Medium")
color4 <- paletteer::paletteer_d("MoMAColors::Althoff")

Esin_color <- color1[1]
Bger_color <- color1[4]
Pame_color <- generate_palette(paletteer_d("fishualize::Acanthurus_sohal")[1], modification = "go_lighter",  n_colours = 5, view_palette = TRUE)[2]
Pful_color <- color1[3]

# load metabolome data
metabolome_data <- read.table("3_down_stream_analysis/1_metabolome_LOESS_PQN_missForest_matrix_noQC.tsv", header = TRUE, sep = "\t", row.names = 1) %>% t()
# load group file
group <- read.table('3_down_stream_analysis/2_group.txt', sep = '\t', header =TRUE)
# Merge the group information into the metabolite table
metabolome_data_grouped <- merge(group, metabolome_data, by.x = "ID", by.y = "row.names")

# Calculate the mean value of the group
library(dplyr)
metabolome_mean_abundance <- metabolome_data_grouped %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "group") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "metabolites")

write.table(metabolome_mean_abundance, "3_down_stream_analysis/3_metabolome_group_mean_abundance.tsv", row.names = FALSE, quote = FALSE)


##########################################
#### Dimensionality Reduction Analysis ###
##########################################

### (PCA) Principal Components Analysis
metabolome_pca <- prcomp(metabolome_data, scale. = TRUE)

# # Extract the coordinates (principal component scores) from the PCA results
metabolome_pca_coords <- metabolome_pca$x
metabolome_pca_result <- merge(metabolome_pca_coords, group, by.x = "row.names", by.y = "ID")

# # Extract all principal components
metabolome_pca_scores <- metabolome_pca$x[, 1:3] 

# Calculate the distance matrix based on PCA scores
metabolome_distance <- dist(metabolome_pca_scores)
metabolome_distance_matrix <- as.matrix(metabolome_distance)
write.table(metabolome_distance_matrix, "3_down_stream_analysis/4_metabolome_distance_matrix.tsv", sep = "\t", quote = FALSE)

# The "pca$sdev" records the eigenvalues of the main sorting axes in the PCA sorting results (dividing each eigenvalue by the total sum of eigenvalues gives the explanatory power of each axis)
metabolome_pca_eig = sum(pmax(metabolome_pca$sdev[1:3]), 0)
metabolome_pca_eig_percent <- round(metabolome_pca$sdev[1:3]/metabolome_pca_eig*100, 3) 

# Conduct a permutation multivariate (factorial) variance analysis (PERMANOVA/adonis2)
metabolome_pca_permanova_result <- adonis2(metabolome_distance_matrix ~ group, data = group , permutations = 999)
metabolome_pca_dune_adonis <- paste("R2", " = ", round(metabolome_pca_permanova_result$R2, 3), "\nP = ", round(metabolome_pca_permanova_result$`Pr(>F)`, 3))
metabolome_pca_dune_adonis

### PCA plot
# Calculate the center point of each group
metabolome_pca_center_points <- 
  metabolome_pca_result %>%
  group_by(group) %>%
  summarise(mean_wt = mean(PC1), mean_mpg = mean(PC2))

metabolome_pca_result <- 
  metabolome_pca_result %>% 
  left_join(metabolome_pca_center_points, by = "group")

metabolome_pca_result$group <- factor(metabolome_pca_result$group, level = c("Esin", "Bger", "Pame", "Pful"))

# Basic scatter plot
metabolome_pca_p0 <- 
  ggplot(metabolome_pca_result, 
         aes(x = PC1, y = PC2, color = group)
  ) +
  geom_point(aes(color = group, shape = group), 
             size = 2,
             show.legend = T,
  ) +
  labs(x = paste("PC 1 (", round(metabolome_pca_eig_percent[1], 2), "%)", sep = ""), 
       y = paste("PC 2 (", round(metabolome_pca_eig_percent[2], 2), "%)", sep = ""), 
       tag = metabolome_pca_dune_adonis
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

metabolome_pca_p0

# Add confidence ellipse + center point connection
metabolome_pca_p1 <- 
  metabolome_pca_p0 + 
  stat_ellipse(
    data = metabolome_pca_result,
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
    data = metabolome_pca_result,
    aes(x = PC1, y = PC2, xend = mean_wt, yend = mean_mpg, color = group), 
    linetype = "dashed", 
    linewidth = 1, 
    alpha = 0.5,
    show.legend = F)

metabolome_pca_p1

# Use ggside package to add marginal box plots
library(ggside)
metabolome_pca_p_final <- 
  metabolome_pca_p1 +
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

print(metabolome_pca_p_final)

ggsave("3_down_stream_analysis/5_metabolome_PCA_PC1~PC2.svg", plot = metabolome_pca_p_final, width = 5, height = 5, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/5_metabolome_PCA_PC1~PC2.tiff", plot = metabolome_pca_p_final, width = 5, height = 5, units = "cm", dpi = 300)



################################
#### Hierarchical clustering ###
################################

# Hierarchical clustering based on spearman distance

metabolome_spearman <- as.dist(1 - cor(t(metabolome_data), method = "spearman"))
metabolome_spearman_matrix <- as.matrix(metabolome_spearman)
write.table(metabolome_spearman_matrix, "3_down_stream_analysis/6_metabolome_sample_spearman.txt", sep = "\t", quote = FALSE)

metabolome_hc <- flashClust(
  metabolome_spearman,
  method = "average",
  members = NULL)

# Convert the hclust object to a phylo object
library(ape)
metabolome_hc_tree <- ape::as.phylo(metabolome_hc)

# Export Newick format tree file
write.tree(metabolome_hc_tree, file = "3_down_stream_analysis/7_metabolome_spearman_average_tree.newick")

# Convert the hclust object to a dendrogram object
metabolome_dend <- as.dendrogram(metabolome_hc)

# Customize the dendrogram
metabolome_dend <- metabolome_dend %>% set("labels_col", c(rep(Esin_color, 6), rep(Bger_color, 6 ), rep(Pful_color ,6), rep(Pame_color, 6)))       # Set the color of the labels
metabolome_dend <- metabolome_dend %>% set("branches_k_color",c(rep(Esin_color, 6), rep(Bger_color, 6 ), rep(Pful_color ,6), rep(Pame_color, 6)))  # Set branch color
metabolome_dend <- metabolome_dend %>% set("branches_lwd", 4)                                                                                      # Set the branch width
metabolome_dend <- metabolome_dend %>% set("branches_lty", 1)                                                                                      # Set the branch line type
metabolome_dend <- metabolome_dend %>% set("nodes_pch", 19)                                                                                        # Set the shape of the node point
metabolome_dend <- metabolome_dend %>% set("nodes_col", "black")                                                                                   # Set the color of the node point
metabolome_dend <- metabolome_dend %>% set("nodes_cex", 0)                                                                                         # Set the size of the node point
metabolome_dend <- metabolome_dend %>% set("labels_cex", 1.5)                                                                                      # Set the label size
metabolome_dend <- metabolome_dend %>% set("leaves_pch", c(rep(17,6), rep(18,6), rep(15,6), rep(19,6)))                                            # Set the shape of the label point
metabolome_dend <- metabolome_dend %>% set("leaves_col",c(rep(Esin_color, 6), rep(Bger_color, 6 ), rep(Pful_color ,6), rep(Pame_color, 6)))        # Set the color of the label point
metabolome_dend <- metabolome_dend %>% set("leaves_cex", 2)                                                                                        # Set the size of the label point

# Plot dendrogram
pdf("3_down_stream_analysis/8_metabolome_Hierarchical_Dendrogram.pdf")
plot(metabolome_dend, main = "Metabolome Hierarchical Dendrogram")
dev.off()



### merge tree 
library(ape)
library(ggtree)
library(patchwork)
library(phytools)

groupInfo <- split(group$ID, group$group)

p_tree <-
  ggtree(groupOTU(metabolome_hc_tree, groupInfo), 
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

ggsave("3_down_stream_analysis/9_metabolome_Hierarchical_Dendrogram.svg", plot = p_tree, width = 5, height = 5, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/9_metabolome_Hierarchical_Dendrogram.tiff", plot = p_tree, width = 5, height = 5, units = "cm", dpi = 300)


### Hierarchical clustering heatmap
row_annotation <- group[,c(1,2,4)]
row.names(row_annotation) <- row_annotation[, 1]
colnames(row_annotation)[2] <- "Host_group"
colnames(row_annotation)[3] <- "Host_habitat"

col_annotation <- group[,c(1,2,4)]
row.names(col_annotation) <- col_annotation[, 1]
colnames(col_annotation)[2] <- "Host_group"
colnames(col_annotation)[3] <- "Host_habitat"

color1 <- paletteer::paletteer_d("nationalparkcolors::Badlands")
color2 <- paletteer::paletteer_d("PrettyCols::Bright")
color3 <- paletteer::paletteer_d("ggthemes::Classic_10_Medium")
color4 <- paletteer::paletteer_d("MoMAColors::Althoff")

Esin_color <- color1[1]
Bger_color <- color1[4]
Pame_color <- generate_palette(paletteer_d("fishualize::Acanthurus_sohal")[1], modification = "go_lighter",  n_colours = 5, view_palette = TRUE)[2]
Pful_color <- color1[3]

# "annot_colors" is a list that contains all the factor names and their corresponding colors.
# The names of each vector in the list correspond to the column names in the row_annotation/col_annotation matrix.
annot_colors <- list(Host_group = c(Esin = Esin_color, Bger = Bger_color, Pame = Pame_color, Pful = Pful_color),
                     Host_habitat = c(Esin = "#751C6DFF", other = paletteer::paletteer_d("nationalparkcolors::Acadia")[1]))

library(pheatmap)

min = min(1- metabolome_spearman_matrix)
max = max(1- metabolome_spearman_matrix)

pdf("3_down_stream_analysis/10_metabolome Hierarchical heatmap.pdf",
    width  = 10,
    height = 10)

pheatmap(
  as.matrix(1- metabolome_spearman_matrix),
  color = colorRampPalette(c("lightyellow", "lightblue", "black"))(100),
  ### Global parameter
  main = "Compounds Spearman Hierarchical heatmap",
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
  cluster_rows = metabolome_hc,
  cluster_cols = metabolome_hc,
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


### total_metabolite_heatmap
metabolite_dataframe <- 
  read.table("3_down_stream_analysis/1_metabolome_LOESS_PQN_missForest_matrix_noQC.tsv", header = T, row.names = 1) %>%
  t() %>%
  scale() %>%
  t()

col_annotation <- group[,c(1,2,4)]
row.names(col_annotation) <- col_annotation[, 1]
colnames(col_annotation)[2] <- "Host_group"
colnames(col_annotation)[3] <- "Host_habitat"

annot_colors <- list(Host_group = c(Esin = Esin_color, Bger = Bger_color, Pame = Pame_color, Pful = Pful_color),
                     Host_habitat = c(Esin = "#453947FF", other = paletteer::paletteer_d("nationalparkcolors::Acadia")[1]))

pdf("3_down_stream_analysis/11_metabolome_total_heatmap.pdf",
    width  = 10,
    height = 10)

pheatmap(metabolite_dataframe,
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

### module_heatmap
MetaB_module_dataframe <- 
  read.table("2_WGCNA/merge_module/merge_module_Eigengenes.txt", header = T, row.names = 1) %>%
  t() %>%
  scale() %>%
  t()

col_annotation <- group[,c(1,2,4)]
row.names(col_annotation) <- col_annotation[, 1]
colnames(col_annotation)[2] <- "Host_group"
colnames(col_annotation)[3] <- "Host_habitat"

annot_colors <- list(Host_group = c(Esin = Esin_color, Bger = Bger_color, Pame = Pame_color, Pful = Pful_color),
                     Host_habitat = c(Esin = "#453947FF", other = paletteer::paletteer_d("nationalparkcolors::Acadia")[1]))

pdf("3_down_stream_analysis/12_metabolome_module_heatmap.pdf",
    width  = 10,
    height = 10)

pheatmap(MetaB_module_dataframe,
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


# enriched WGCNA modules heatmap
metabolome_up_module <- read.table("2_WGCNA/merge_module/1_all_up_module.tsv", header = TRUE, row.names = 1 ,sep = "\t") 

metabolome_up_module <- 
  metabolome_up_module %>%
  mutate(log2trans = metabolome_up_module$Directionality.x..log10.p.adj.) %>%
  dplyr::select(log2trans, everything())

metabolome_up_module <- 
  metabolome_up_module %>%
  mutate(module = rownames(metabolome_up_module)) %>%
  dplyr::select(module, everything())

up_module_mratrix <- 
  metabolome_up_module %>%
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
               values_to = "abundance")

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

up_module_mratrix_long$sample <- 
  factor(
    up_module_mratrix_long$sample,
    levels = sample_order
  )

min <- min(up_module_mratrix_long$abundance)
max <- max(up_module_mratrix_long$abundance)

# plot main heatmap
up_module_heatmap_p <- 
  ggplot(up_module_mratrix_long, 
         aes(x = module, y = sample, fill = abundance)
  ) +
  geom_tile(color = NA,
            linewidth = 0,
            lineend = "square",
            height = 1,
            show.legend = T
  ) +
  annotate("rect", ymin = 18.5, ymax = 24.5, xmin = 0.5, xmax = 52.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 12.5, ymax = 18.5, xmin = 52.5, xmax = 92.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 6.5, ymax = 12.5, xmin = 92.5 , xmax = 108.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 0.5, ymax = 6.5, xmin = 108.5, xmax = 130.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF", 
                       limit = c(min, max),
                       breaks = c(-2, 0 , 2),
  ) +
  scale_x_discrete(position = "top") + 
  labs(x = "Enriched metabolome modules of four species cockroaches", 
       fill = "Enrich Z-Score"
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

up_module_heatmap_p

# plot sample grouping annotatin
group$group <- factor(group$group, levels = c("Esin", "Bger", "Pame", "Pful"))

up_module_group_p <-
  ggplot(group, 
         aes(x = 0, y = ID, fill = group, color = group)
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.05,
            lineend = "square",
            show.legend = F
  ) +
  scale_fill_manual(values = c("Esin" = Esin_color, "Bger" = Bger_color, "Pame" = Pame_color, "Pful" = Pful_color)) +
  theme_void()

up_module_group_p


# plot module pvalue annotatin
up_module_module_p <-
  ggplot(metabolome_up_module, 
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
  scale_fill_gradient2(low = paletteer::paletteer_d("MetBrewer::Benedictus")[5], 
                       high = paletteer::paletteer_d("MetBrewer::Benedictus")[3],
                       breaks = c(2, 4, 6)
  ) +
  scale_color_gradient2(low = paletteer::paletteer_d("MetBrewer::Benedictus")[5], 
                        high = paletteer::paletteer_d("MetBrewer::Benedictus")[3],
                        breaks = c(2, 4, 6)
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

up_module_module_p

# plot module category annotatin
up_module_module_info <- read.table("2_WGCNA/merge_module/1_all_up_module_information.tsv", header = TRUE, sep = "\t") %>% as.data.frame()

table(up_module_module_info$Module.category)

up_module_module_info <- 
  up_module_module_info %>%
  mutate(Module.category = case_when(
    up_module_module_info$Module.category == "Benzenoids" ~ "Benzenoids",
    up_module_module_info$Module.category == "Lipids and lipid-like molecules" ~ "Lipids",
    up_module_module_info$Module.category == "Organic acids and derivatives" ~ "Organic acids",
    up_module_module_info$Module.category == "Organoheterocyclic compounds" ~ "Organoheterocyclic",
    up_module_module_info$Module.category == "Organic oxygen compounds" ~ "Organic oxygen",
    up_module_module_info$Module.category == "Phenylpropanoids and polyketides" ~ "PPs & PKs",
    up_module_module_info$Module.category == "Miscellaneous" ~ "Miscellaneous",
  ))

up_module_module_info$Module.category <- 
  factor(up_module_module_info$Module.category,
         levels = c("Benzenoids", "Lipids", "Organic acids", "Organoheterocyclic", "Organic oxygen", 
                    "PPs & PKs", "Organic nitrogen", "Nucleosides", "Miscellaneous"))

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
  scale_fill_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:8],"grey"), drop = FALSE) +
  scale_color_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:8],"grey"), drop = FALSE) +
  theme_void() + 
  theme(
    legend.position = "none",
    legend.key.size = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"), 
    legend.key.width = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 2, color = "black"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "grey60"),
  )

up_module_module_info_p

# merge plots 
up_module <- up_module_heatmap_p %>% 
  aplot::insert_left(up_module_group_p, width = 0.02) %>% 
  aplot::insert_top(up_module_module_p, height = 0.05) %>%
  aplot::insert_bottom(up_module_module_info_p, height = 0.05)

up_module

# save image
ggsave("3_down_stream_analysis/13_metabolome_WGCNA_module_enrich_module.svg", plot = up_module, width = 10, height = 4.5, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/13_metabolome_WGCNA_module_enrich_module.tiff", plot = up_module, width = 10, height = 4.5, units = "cm", dpi = 300)

# save all legend
up_module_heatmap_p_legend <- cowplot::get_legend(up_module_heatmap_p + theme(legend.position = "right"))
up_module_module_p_legend <- cowplot::get_legend(up_module_module_p + theme(legend.position = "right"))
up_module_module_info_p_legend <- cowplot::get_legend(up_module_module_info_p + theme(legend.position = "right"))
up_module_combined_legend <- cowplot::plot_grid(up_module_heatmap_p_legend, up_module_module_p_legend, up_module_module_info_p_legend, nrow = 3)

up_module_combined_legend

ggsave("3_down_stream_analysis/14_metabolome_WGCNA_module_enrich_module_legend.svg", plot = up_module_combined_legend, width = 20, height = 10, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/14_metabolome_WGCNA_module_enrich_module_legend.tiff", plot = up_module_combined_legend, width = 20, height = 10, units = "cm", dpi = 300)


# depleted WGCNA modules heatmap
metabolome_down_module <- read.table("2_WGCNA/merge_module/2_all_down_module.tsv", header = TRUE, row.names = 1 ,sep = "\t") 

metabolome_down_module <- 
  metabolome_down_module %>%
  mutate(log2trans = metabolome_down_module$Directionality.x..log10.p.adj.) %>%
  dplyr::select(log2trans, everything())

metabolome_down_module <- 
  metabolome_down_module %>%
  mutate(module = rownames(metabolome_down_module)) %>%
  dplyr::select(module, everything())

down_module_mratrix <- 
  metabolome_down_module %>%
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
               values_to = "abundance")

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

min <- min(down_module_mratrix_long$abundance)
max <- max(down_module_mratrix_long$abundance)

# plot main heatmap
down_module_heatmap_p <- 
  ggplot(down_module_mratrix_long, 
         aes(x = module, y = sample, fill = abundance)
  ) +
  geom_tile(height = 1,
            color = NA,
            linewidth = 0,
            lineend = "square",
            show.legend = T
  ) +
  annotate("rect", ymin = 0.5, ymax = 6.5, xmin = 0.5, xmax = 45.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 6.5, ymax = 12.5, xmin = 45.5 , xmax = 86.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 12.5, ymax = 18.5, xmin = 86.5, xmax = 97.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 18.5, ymax = 24.5, xmin = 97.5, xmax = 107.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF", 
                       limit = c(min, max),
                       breaks = c(-2, 0, 2)
  ) +
  scale_x_discrete(position = "bottom") + 
  labs(x = "Depleted metabolome modules of four species cockroaches", 
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

down_module_group_p 

# plot module pvalue annotatin
down_module_module_p <-
  ggplot(metabolome_down_module, 
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
  scale_fill_gradient2(low = paletteer::paletteer_d("MetBrewer::Benedictus")[11], 
                      high = paletteer::paletteer_d("MetBrewer::Benedictus")[9],
                      breaks = c(-6, -4, -2)
  ) +
  scale_color_gradient2(low = paletteer::paletteer_d("MetBrewer::Benedictus")[11], 
                       high = paletteer::paletteer_d("MetBrewer::Benedictus")[9],
                       breaks = c(-6, -4, -2)
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

down_module_module_p

# plot module category annotatin
down_module_module_info <- read.table("2_WGCNA/merge_module/2_all_down_module_information.tsv", header = TRUE, sep = "\t") %>% as.data.frame()

table(down_module_module_info$Module.category)

down_module_module_info <- 
  down_module_module_info %>%
  mutate(Module.category = case_when(
    down_module_module_info$Module.category == "Benzenoids" ~ "Benzenoids",
    down_module_module_info$Module.category == "Lipids and lipid-like molecules" ~ "Lipids",
    down_module_module_info$Module.category == "Organic acids and derivatives" ~ "Organic acids",
    down_module_module_info$Module.category == "Organoheterocyclic compounds" ~ "Organoheterocyclic",
    down_module_module_info$Module.category == "Organic oxygen compounds" ~ "Organic oxygen",
    down_module_module_info$Module.category == "Phenylpropanoids and polyketides" ~ "PPs & PKs",
    down_module_module_info$Module.category == "Organic nitrogen compounds" ~ "Organic nitrogen",
    down_module_module_info$Module.category == "Nucleosides, nucleotides, and analogues" ~ "Nucleosides",
    down_module_module_info$Module.category == "Miscellaneous" ~ "Miscellaneous",
  ))

down_module_module_info$Module.category <- 
  factor(down_module_module_info$Module.category,
         levels = c("Benzenoids", "Lipids", "Organic acids", "Organoheterocyclic", "Organic oxygen", 
                    "PPs & PKs", "Organic nitrogen", "Nucleosides", "Miscellaneous"))

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
  scale_fill_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:8],"grey"), drop = FALSE) +
  scale_color_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:8],"grey"), drop = FALSE) +
  guides(fill = guide_legend(ncol = 5, row = 2)) + 
  theme_void() + 
  theme(
    legend.position = "none",
    legend.key.size = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"), 
    legend.key.width = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 2, color = "black"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "grey60"),
  )

down_module_module_info_p

# merge plots 
down_module <- down_module_heatmap_p %>% 
  aplot::insert_left(down_module_group_p, width = 0.02) %>% 
  aplot::insert_bottom(down_module_module_p, height = 0.05) %>%
  aplot::insert_top(down_module_module_info_p, height = 0.05)

down_module

# save image
ggsave("3_down_stream_analysis/15_metabolome_WGCNA_module_deplete_module.svg", plot = down_module, width = 10, height = 4.5, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/15_metabolome_WGCNA_module_deplete_module.tiff", plot = down_module, width = 10, height = 4.5, units = "cm", dpi = 300)


# save all legend
down_module_heatmap_p_legend <- cowplot::get_legend(down_module_heatmap_p + theme(legend.position = "right"))
down_module_module_p_legend <- cowplot::get_legend(down_module_module_p + theme(legend.position = "right"))
down_module_module_info_p_legend <- cowplot::get_legend(down_module_module_info_p + theme(legend.position = "right"))
down_module_combined_legend <- cowplot::plot_grid(down_module_heatmap_p_legend, down_module_module_p_legend, down_module_module_info_p_legend, nrow = 3)

down_module_combined_legend

ggsave("3_down_stream_analysis/16_metabolome_WGCNA_module_deplete_module_legend.svg", plot = down_module_combined_legend, width = 20, height = 10, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/16_metabolome_WGCNA_module_deplete_module_legend.tiff", plot = down_module_combined_legend, width = 20, height = 10, units = "cm", dpi = 300)


### regulation module percent
up_module_module_percent <- 
  metabolome_up_module %>% left_join(up_module_module_info, by = "module") %>% 
  group_by(regulation, Module.category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(regulation) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

up_module_module_percent$Module.category <- factor(up_module_module_percent$Module.category,
                                                   levels = rev(c("Benzenoids", "Lipids", "Organic acids", "Organoheterocyclic", "Organic oxygen", 
                                                              "PPs & PKs", "Organic nitrogen", "Nucleosides", "Miscellaneous")))

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
  scale_fill_manual(values = c("grey", paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[8:1]), drop = FALSE) +
  scale_color_manual(values = c("grey", paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[8:1]), drop = FALSE) +
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

up_module_module_percent_bar_plot

up_module_module_total <- 
  up_module_module_percent %>%
  group_by(regulation) %>%
  summarise(sum = sum(count), .groups = 'drop')

up_module_module_total_bar_plot <- 
  ggplot(up_module_module_total) + 
  geom_bar(aes(x = regulation, y = sum, fill = regulation),
           stat = "identity",
           width = 0.7) +
  geom_label(aes(label = sum, x = regulation, y = 25), size = 3) + 
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

up_module_module_total_bar_plot

down_module_module_percent <- 
  metabolome_down_module %>% left_join(down_module_module_info, by = "module") %>% 
  group_by(regulation, Module.category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(regulation) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

down_module_module_percent$Module.category <- factor(down_module_module_percent$Module.category, 
                                                     levels = rev(c("Benzenoids", "Lipids", "Organic acids", "Organoheterocyclic", "Organic oxygen", 
                                                                "PPs & PKs", "Organic nitrogen", "Nucleosides", "Miscellaneous")))

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
  scale_fill_manual(values = c("grey", paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[8:1]), drop = FALSE) +
  scale_color_manual(values = c("grey", paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[8:1]), drop = FALSE) +
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
  geom_label(aes(label = sum, x = regulation, y = 20), size = 3) + 
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

ggsave("3_down_stream_analysis/17_module_percent_bar_plot.svg", plot = module_percent_bar_plot, width = 6, height = 12, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/17_module_percent_bar_plot.tiff", plot = module_percent_bar_plot, width = 6, height = 12, units = "cm", dpi = 300)


### highlight top 10 enriched/depleted metabolome modules

### top 10 enriched metabolomic modules
top10_enriched_module <- read.table("2_WGCNA/merge_module/3_each_top10_up_module.tsv", header = TRUE, sep = "\t") 
top10_enriched_module$log2trans <- ifelse(top10_enriched_module$Directionality.x..log10.p.adj. > 6, 6, top10_enriched_module$Directionality.x..log10.p.adj.)
top10_enriched_module$mark <- ifelse(top10_enriched_module$Directionality.x..log10.p.adj. > 6, "over", "not")
top10_enriched_module <- top10_enriched_module %>% dplyr::select(log2trans, mark, everything())
top10_enriched_module$regulation <- factor(top10_enriched_module$regulation, levels = c("Esin_up", "Bger_up", "Pame_up", "Pful_up"))
top10_enriched_module$module <- factor(top10_enriched_module$module, levels = top10_enriched_module$module)
top10_enriched_module_info <- read.table("2_WGCNA/merge_module/3_each_top10_up_module_information.tsv", header = TRUE, sep = "\t") 

top10_enriched_module_info <- 
  top10_enriched_module_info %>%
  mutate(Module.category.short = case_when(
    top10_enriched_module_info$Module.category == "Benzenoids" ~ "Benzenoids",
    top10_enriched_module_info$Module.category == "Lipids and lipid-like molecules" ~ "Lipids",
    top10_enriched_module_info$Module.category == "Organic acids and derivatives" ~ "Organic acids",
    top10_enriched_module_info$Module.category == "Organoheterocyclic compounds" ~ "Organoheterocyclic",
    top10_enriched_module_info$Module.category == "Organic oxygen compounds" ~ "Organic oxygen",
    top10_enriched_module_info$Module.category == "Phenylpropanoids and polyketides" ~ "PPs & PKs",
    top10_enriched_module_info$Module.category == "Organic nitrogen compounds" ~ "Organic nitrogen",
    top10_enriched_module_info$Module.category == "Nucleosides, nucleotides, and analogues" ~ "Nucleosides",
    top10_enriched_module_info$Module.category == "Miscellaneous" ~ "Miscellaneous",
  ))

top10_enriched_module_info$merge.name <- paste0(gsub("MetaB_", "", top10_enriched_module_info$module), ":", top10_enriched_module_info$Module.category)

top10_enriched_module_info$merge.name <- gsub("Nucleosides, nucleotides, and analogues", "Nucleosides/Nucleotides/Analogues", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub(" and ", " / ", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("lipid-like", "Lipid-like", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("derivatives", "Derivatives", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("polyketides", "Polyketides", top10_enriched_module_info$merge.name)

top10_enriched_module_info$Module.category.short <- 
  factor(top10_enriched_module_info$Module.category.short,
         levels = c("Benzenoids", "Lipids", "Organic acids", "Organoheterocyclic", "Organic oxygen", 
                    "PPs & PKs", "Organic nitrogen", "Nucleosides", "Miscellaneous"))

# top10_enriched_module bar plot
top10_enriched_module_p0 <- 
  ggplot(top10_enriched_module, 
         aes(x = log2trans, y = module , fill = regulation)
  ) +
  geom_hline(yintercept = seq_along(top10_enriched_module$module), 
             color = "gray", 
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
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, 
                               color = c(rep(Pful_color,10), rep(Pame_color,10), rep(Bger_color,10), rep(Esin_color, 10))), 
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

top10_enriched_module_p0

# top10_enriched_module module category annotation
top10_enriched_module_info$module <- factor(top10_enriched_module_info$module, levels = top10_enriched_module$module)

top10_enriched_module_p1 <- 
  ggplot(top10_enriched_module_info, 
         aes(x = 0, y = module, fill = Module.category.short), 
         width = 0.15
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.25,
            lineend = "square",
            show.legend = T
  ) + 
  scale_fill_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[c(1:8)], "grey"), drop =FALSE) +
  labs(fill = "Module categary") + 
  guides(fill = guide_legend(ncol = 1)) +
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

top10_enriched_module_p1

# plot heatmap
sample_levels <- c(paste("Esin_",c(1:6), sep = ""), paste("Bger_",c(1:6), sep = ""), paste("Pame_",c(1:6), sep = ""), paste("Pful_",c(1:6), sep = ""))

up_module_z_scores <- top10_enriched_module[, 12:35] %>% t() %>% scale() %>% t()

top10_enriched_module[, 12:35] <- up_module_z_scores

# wide to long
top10_enriched_module_long <- top10_enriched_module %>%
  pivot_longer(cols = -c(1:11), names_to = "sample", values_to = "abundance")


top10_enriched_module_long$sample <- factor(top10_enriched_module_long$sample, levels = sample_levels)

min <- min(top10_enriched_module_long$abundance)
max <- max(top10_enriched_module_long$abundance)

top10_enriched_module_p2 <- 
  ggplot(top10_enriched_module_long, 
         aes(x = sample, y = module, fill = abundance)
  ) +
  geom_tile(height = 1,
            color = NA,
            linewidth = 0,
            lineend = "square",
            show.legend = T
  ) +
  annotate("rect", xmin = 0.5, xmax = 6.5, ymin = 30.5, ymax = 40.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 6.5, xmax = 12.5, ymin = 20.5, ymax = 30.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 12.5, xmax = 18.5, ymin = 10.5, ymax = 20.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 18.5, xmax = 24.5, ymin = 0.5, ymax = 10.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
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
    legend.key.spacing.y  = unit(5, "pt"),
    legend.key.spacing.x  = unit(5, "pt"),
    legend.ticks.length = unit(2.5, "pt"),
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

ggsave("3_down_stream_analysis/18_top10_enriched_module.svg", plot = top10_enriched_module_p, width = 14, height = 9.5, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/18_top10_enriched_module.tiff", plot = top10_enriched_module_p, width = 14, height = 9.5, units = "cm", dpi = 300)


# top 10 depleted metabolomic modules
top10_depleted_module <- read.table("2_WGCNA/merge_module/4_each_top10_down_module.tsv", header = TRUE, sep = "\t") 
top10_depleted_module$log2trans <- ifelse(top10_depleted_module$Directionality.x..log10.p.adj. < -6, -6, top10_depleted_module$Directionality.x..log10.p.adj.)
top10_depleted_module$mark <- ifelse(top10_depleted_module$Directionality.x..log10.p.adj. < -6, "over", "not")
top10_depleted_module <- top10_depleted_module %>% dplyr::select(log2trans, mark, everything())
top10_depleted_module$regulation <- factor(top10_depleted_module$regulation, levels = c("Esin_down", "Bger_down", "Pame_down", "Pful_down"))
top10_depleted_module$module <- factor(top10_depleted_module$module, levels = top10_depleted_module$module)

top10_depleted_module_info <- read.table("2_WGCNA/merge_module/4_each_top10_down_module_information.tsv", header = TRUE, sep = "\t") 

top10_depleted_module_info <- 
  top10_depleted_module_info %>%
  mutate(Module.category.short = case_when(
    top10_depleted_module_info$Module.category == "Benzenoids" ~ "Benzenoids",
    top10_depleted_module_info$Module.category == "Lipids and lipid-like molecules" ~ "Lipids",
    top10_depleted_module_info$Module.category == "Organic acids and derivatives" ~ "Organic acids",
    top10_depleted_module_info$Module.category == "Organoheterocyclic compounds" ~ "Organoheterocyclic",
    top10_depleted_module_info$Module.category == "Organic oxygen compounds" ~ "Organic oxygen",
    top10_depleted_module_info$Module.category == "Phenylpropanoids and polyketides" ~ "PPs & PKs",
    top10_depleted_module_info$Module.category == "Organic nitrogen compounds" ~ "Organic nitrogen",
    top10_depleted_module_info$Module.category == "Nucleosides, nucleotides, and analogues" ~ "Nucleosides",
    top10_depleted_module_info$Module.category == "Miscellaneous" ~ "Miscellaneous",
  ))

top10_depleted_module_info$merge.name <- paste0(gsub("MetaB_", "", top10_depleted_module_info$module), ":", top10_depleted_module_info$Module.category)

top10_depleted_module_info$merge.name <- gsub("Nucleosides, nucleotides, and analogues", "Nucleosides/Nucleotides/Analogues", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub(" and ", " / ", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("lipid-like", "Lipid-like", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("derivatives", "Derivatives", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("polyketides", "Polyketides", top10_depleted_module_info$merge.name)

top10_depleted_module_info$merge.name

top10_depleted_module_info$Module.category.short <- 
  factor(top10_depleted_module_info$Module.category.short,
         levels = c("Benzenoids", "Lipids", "Organic acids", "Organoheterocyclic", "Organic oxygen", 
                    "PPs & PKs", "Organic nitrogen", "Nucleosides", "Miscellaneous"))

# top10_depleted_module bar plot
top10_depleted_module_p0 <- 
  ggplot(top10_depleted_module, 
         aes(x = log2trans, y = module , fill = regulation)
  ) +
  geom_hline(yintercept = seq_along(top10_depleted_module$module), 
             color = "gray", 
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
                   label = top10_depleted_module_info$merge.name
  )  +
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
                               color = c(rep(Pful_color,10), rep(Pame_color,10), rep(Bger_color,10), rep(Esin_color, 10))), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x = unit(2, "pt"),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.25, linetype = "solid"),
    legend.key.height = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.size = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.direction = "vertical",
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.text.position = "right",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
  )

top10_depleted_module_p0

# top10_depleted_module module category annotation
top10_depleted_module_info$module <- factor(top10_depleted_module_info$module, levels = top10_depleted_module$module)

top10_depleted_module_p1 <- 
  ggplot(top10_depleted_module_info, 
         aes(x = 0, y = module, fill = Module.category.short), 
         width = 0.15
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.25,
            lineend = "square",
            show.legend = T
  ) + 
  scale_fill_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[c(1:8)], "grey"), drop = FALSE) +
  labs(fill = "Module categary") + 
  guides(fill = guide_legend(ncol = 1)) +
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

top10_depleted_module_p1

# plot heatmap
sample_levels <- c(paste("Esin_",c(1:6), sep = ""), paste("Bger_",c(1:6), sep = ""), paste("Pame_",c(1:6), sep = ""), paste("Pful_",c(1:6), sep = ""))

up_module_z_scores <- top10_depleted_module[, 12:35] %>% t() %>% scale() %>% t()

top10_depleted_module[, 12:35] <- up_module_z_scores

# wide to long
top10_depleted_module_long <- top10_depleted_module %>%
  pivot_longer(cols = -c(1:11), names_to = "sample", values_to = "abundance")

top10_depleted_module_long$sample <- factor(top10_depleted_module_long$sample, levels = rev(sample_levels))

min <- min(top10_depleted_module_long$abundance)
max <- max(top10_depleted_module_long$abundance)

top10_depleted_module_p2 <- 
  ggplot(top10_depleted_module_long, 
         aes(x = sample, y = module, fill = abundance)
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
                       breaks = c(-2, 0, 2)
  ) +
  annotate("rect", xmin = 18.5, xmax = 24.5, ymin = 30.5, ymax = 40.5, fill = NA, color = "blue", linewidth = 0.25, lineend = "square") +
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

ggsave("3_down_stream_analysis/19_top10_depleted_module.svg", plot = top10_depleted_module_p, width = 14, height = 9.5, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/19_top10_depleted_module.tiff", plot = top10_depleted_module_p, width = 14, height = 9.5, units = "cm", dpi = 300)

# merge top10_depleted_module_p & top10_enriched_module_p
top10_delepeted_enriched_module_plot <- top10_depleted_module_p0 + top10_depleted_module_p1 + top10_depleted_module_p2 + 
                                        top10_enriched_module_p2 + top10_enriched_module_p1 + top10_enriched_module_p0 +  
                                        plot_layout(ncol = 6, widths = c(3, 0.5, 4, 4, 0.5, 3)) &
                                        theme(panel.spacing = unit(0.2, "lines"))

top10_delepeted_enriched_module_plot

top10_delepeted_enriched_module_plot_0 <- top10_depleted_module_p0 + theme(axis.text.y = element_blank()) + top10_depleted_module_p1 + top10_depleted_module_p2 + 
                                          top10_enriched_module_p2 + top10_enriched_module_p1 + top10_enriched_module_p0 + theme(axis.text.y = element_blank()) +  
                                          plot_layout(ncol = 6, widths = c(3, 0.4, 3, 3, 0.4, 3)) &
                                          theme(panel.spacing = unit(0.2, "lines"))

top10_delepeted_enriched_module_plot_0

ggsave("3_down_stream_analysis/20_top10_delepeted_enriched_module.svg", plot = top10_delepeted_enriched_module_plot, width = 20, height = 8, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/20_top10_delepeted_enriched_module.tiff", plot = top10_delepeted_enriched_module_plot, width = 20, height = 8, units = "cm", dpi = 300)

ggsave("3_down_stream_analysis/21_top10_delepeted_enriched_module_0.svg", plot = top10_delepeted_enriched_module_plot_0, width = 7, height = 8, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/21_top10_delepeted_enriched_module_0.tiff", plot = top10_delepeted_enriched_module_plot_0, width = 7, height = 8, units = "cm", dpi = 300)


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

top10_depleted_enriched_module_plot_legend

ggsave("3_down_stream_analysis/22_top10_delepeted_enriched_module_legend.svg", plot = top10_depleted_enriched_module_plot_legend, width = 20, height = 20, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/22_top10_delepeted_enriched_module_legend.tiff", plot = top10_depleted_enriched_module_plot_legend, width = 20, height = 20, units = "cm", dpi = 300)



########################
### Machine Learning ###
########################

### metabolome KEGG pathway importance ranking by lasso

# load data
library(glmnet)
set.seed(20241102)
lasso_data <- read.table("2_WGCNA/merge_module/5_lasso_input.tsv", header = TRUE, row.names = 1, check.names = FALSE)
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
                      alpha   = best_alpha)

lambda_min  <- cv_final$lambda.min
lambda_1se  <- cv_final$lambda.1se
print(tibble(lambda.min = lambda_min, lambda.1se = lambda_1se))

plot(results$alpha, results$cv.error, type = "b", xlab = "Alpha", ylab = "CV Error", main = "CV Error vs Alpha")

pdf("3_down_stream_analysis/23_Cross-validation.pdf",
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

ggsave("3_down_stream_analysis/24_cv_plot.svg", plot = cv_plot, width = 6, height = 6, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/24_cv_plot.tiff", plot = cv_plot, width = 6, height = 6, units = "cm", dpi = 300)

# Extract the modules with non-zero importance(lambda.min)
coef_mat <- as.matrix(coef(cv_final, s = "lambda.min"))
active.idx   <- which(coef_mat[, 1] != 0)
active.genes <- rownames(coef_mat)[active.idx] %>% setdiff("(Intercept)")

# export result
coef.df <- coef(cv_final, s = "lambda.min") %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("module") %>%
  filter(module != "(Intercept)") %>%
  mutate(coef = lambda.min) %>%
  mutate(absCoef = abs(coef)) %>%
  filter(absCoef != 0) %>%
  arrange(desc(coef))

write_csv(coef.df, "3_down_stream_analysis/25_Lasso_module_ranked_by_coef.csv")

### top 20 lasso selected feature
top20_lasso_metabolome <- rbind(head(coef.df, 10), tail(coef.df, 10))
top20_lasso_metabolome$module <- gsub("^MetaB_", "",top20_lasso_metabolome$module)
top20_lasso_metabolome$module <- factor(top20_lasso_metabolome$module, levels = rev(top20_lasso_metabolome$module))

metabolome_module_info <- rbind(read.table("2_WGCNA/merge_module/1_all_up_module_information.tsv", header = TRUE, sep = "\t"),
                                read.table("2_WGCNA/merge_module/2_all_down_module_information.tsv", header = TRUE, sep = "\t")) %>%
                                distinct(module, .keep_all = TRUE)

metabolome_module_info$module <- gsub("MetaB_", "", metabolome_module_info$module)

top20_lasso_metabolome <- 
  top20_lasso_metabolome %>% 
  left_join(metabolome_module_info, by = c("module" = "module")) %>%
  mutate(label = paste0(module, ":", Module.category))

top20_lasso_metabolome$label <- factor(top20_lasso_metabolome$label, levels = rev(top20_lasso_metabolome$label))

write.table(top20_lasso_metabolome, "../4.association/1_multi_omic_association/4_top20_e_d_module_information_MetaB.tsv", row.names = F, quote = F, sep = "\t")

metabolome_all_module_WGCNA_data <- 
  read.table("2_WGCNA/merge_module/merge_module_Eigengenes.txt", sep = "\t", header = T, row.names = 1) %>% 
  rownames_to_column(var = "module") %>%
  mutate(module = str_remove(module, "^MetaB_")) %>%
  column_to_rownames(var = "module")

top20_lasso_metabolome_module_WGCNA <- metabolome_all_module_WGCNA_data[top20_lasso_metabolome$module, ] %>% rownames_to_column(var = "module")

write.table(top20_lasso_metabolome_module_WGCNA, "../4.association/1_multi_omic_association/3_top20_e_d_metabolome_module.tsv", sep = "\t", quote = F, row.names = F)


### total lasso selected feature
total_lasso_metabolome <- coef.df
total_lasso_metabolome$module <- gsub("^MetaB_", "",total_lasso_metabolome$module)
total_lasso_metabolome$module <- factor(total_lasso_metabolome$module, levels = rev(total_lasso_metabolome$module))

metabolome_module_info <- rbind(read.table("2_WGCNA/merge_module/1_all_up_module_information.tsv", header = TRUE, sep = "\t"),
                                read.table("2_WGCNA/merge_module/2_all_down_module_information.tsv", header = TRUE, sep = "\t")) %>%
                          distinct(module, .keep_all = TRUE)

metabolome_module_info$module <- gsub("MetaB_", "", metabolome_module_info$module)

total_lasso_metabolome <- 
  total_lasso_metabolome %>% 
  left_join(metabolome_module_info, by = c("module" = "module")) %>%
  mutate(label = paste0(module, ":", Module.category))

total_lasso_metabolome$label <- factor(total_lasso_metabolome$label, levels = rev(total_lasso_metabolome$label))

write.table(total_lasso_metabolome, "../4.association/1_multi_omic_association/53_total_e_d_module_information_MetaB.tsv", row.names = F, quote = F, sep = "\t")

metabolome_all_module_WGCNA_data <- 
  read.table("2_WGCNA/merge_module/merge_module_Eigengenes.txt", sep = "\t", header = T, row.names = 1) %>% 
  rownames_to_column(var = "module") %>%
  mutate(module = str_remove(module, "^MetaB_")) %>%
  column_to_rownames(var = "module")

total_lasso_metabolome_module_WGCNA <- metabolome_all_module_WGCNA_data[total_lasso_metabolome$module, ] %>% rownames_to_column(var = "module")

write.table(total_lasso_metabolome_module_WGCNA, "../4.association/1_multi_omic_association/52_total_e_d_metabolome_module.tsv", sep = "\t", quote = F, row.names = F)


### lasso ranking plot
lasso_importance_rank <- 
  ggplot(top20_lasso_metabolome,
         aes(x = coef, y = module)
  ) +
  geom_point(aes(x = coef, y = label, color = label), 
             size = 1,
             show.legend = F
  ) +
  geom_hline(yintercept = seq_along(top20_lasso_metabolome$label), 
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
  scale_x_continuous(limits = c(-3, 3), breaks = c(-3, 0, 3), labels = c(-3, 0, 3)) +
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

ggsave("3_down_stream_analysis/26_lasso_importance_rank_plot.svg", plot = lasso_importance_rank, width = 10, height = 4.5, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/26_lasso_importance_rank_plot.tiff", plot = lasso_importance_rank, width = 10, height = 4.5, units = "cm", dpi = 300)

ggsave("3_down_stream_analysis/27_lasso_importance_rank_plot_0.svg", plot = lasso_importance_rank_0, width = 1, height = 4.5, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/27_lasso_importance_rank_plot_0.tiff", plot = lasso_importance_rank_0, width = 1, height = 4.5, units = "cm", dpi = 300)


### double volcano plot

# read file
double_volcano_data <- read.table("2_WGCNA/merge_module/6_double_volcano_log2trans.tsv", header = TRUE)

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
  scale_fill_manual(values = rev(c("pink", "lightblue", scater_color[3], scater_color[4], scater_color[9], scater_color[10], scater_color[1], scater_color[2], "grey80")),
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

ggsave("3_down_stream_analysis/28_quadrant_percent_pie_plot.svg", plot = quadrant_percent_pie_plot, width = 1, height = 1, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/23_quadrant_percent_pie_plot.tiff", plot = quadrant_percent_pie_plot, width = 1, height = 1, units = "cm", dpi = 300)


# highlight module
highlight_module <- c(
  MetaB_M174 = "M174",
  MetaB_M040 = "M040", 
  MetaB_M074 = "M074",
  MetaB_M129 = "M129",
  MetaB_M006 = "M006",
  MetaB_M064 = "M064",
  MetaB_M007 = "M007",
  MetaB_M146 = "M146",
  MetaB_M183 = "M183",
  MetaB_M134 = "M134",
  MetaB_M224 = "M224",
  MetaB_M193 = "M193",
  MetaB_M116 = "M116",
  MetaB_M075 = "M075",
  MetaB_M085 = "M085",
  MetaB_M184 = "M184",
  MetaB_M192 = "M192",
  MetaB_M005 = "M005",
  MetaB_M188 = "M188",
  MetaB_M140 = "M140"
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

ggsave("3_down_stream_analysis/29_double_volcano_plot.svg", plot = double_volcano_plot, width = 4.5, height = 4.5, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/29_double_volcano_plot.tiff", plot = double_volcano_plot, width = 4.5, height = 4.5, units = "cm", dpi = 300)




















### venn

# set color 
Esin_venn_color <- generate_palette(Esin_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, view_labels = FALSE)
Bger_venn_color <- generate_palette(Bger_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, view_labels = FALSE)
Pame_venn_color <- generate_palette(Pame_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, view_labels = FALSE)
Pful_venn_color <- generate_palette(Pful_color, modification = "go_lighter", n_colours = 4, view_palette = TRUE, view_labels = FALSE)

# Esin enrich
Esin_venn_dat_up  <- read.delim("2_WGCNA/merge_module/Esin_metabolome_WGCNA_module_each_up.tsv")
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
)

Esin_venn_up

# Bger enrich
Bger_venn_dat_up  <- read.delim("2_WGCNA/merge_module/Bger_metabolome_WGCNA_module_each_up.tsv")
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
)

Bger_venn_up

# Pame enrich
Pame_venn_dat_up  <- read.delim("2_WGCNA/merge_module/Pame_metabolome_WGCNA_module_each_up.tsv")
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
  set_name_color = NA,
  set_name_size = 0,
  text_color = "black",
  text_size = 6/2.83
)

Pame_venn_up

# Pful enrich
Pful_venn_dat_up  <- read.delim("2_WGCNA/merge_module/Pful_metabolome_WGCNA_module_each_up.tsv")
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
)

Pful_venn_up

library(aplot)

venn_up <-Esin_venn_up + Bger_venn_up + Pame_venn_up + Pful_venn_up + plot_layout(ncol =4)
venn_up

ggsave("3_down_stream_analysis/30_venn_up.svg", plot = venn_up, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/30_venn_up.tiff", plot = venn_up, width = 10, height = 5, units = "cm", dpi = 300)


# Esin deplete
Esin_venn_dat_down  <- read.delim("2_WGCNA/merge_module/Esin_metabolome_WGCNA_module_each_down.tsv")
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
)

Esin_venn_down

# Bger deplete
Bger_venn_dat_down  <- read.delim("2_WGCNA/merge_module/Bger_metabolome_WGCNA_module_each_down.tsv")
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
)

Bger_venn_down

# Pame deplete
Pame_venn_dat_down  <- read.delim("2_WGCNA/merge_module/Pame_metabolome_WGCNA_module_each_down.tsv")
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
)

Pame_venn_down

# Pful deplete
Pful_venn_dat_down  <- read.delim("2_WGCNA/merge_module/Pful_metabolome_WGCNA_module_each_down.tsv")
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
  )

Pful_venn_down

library(aplot)

venn_down <- Esin_venn_down + Bger_venn_down  + Pame_venn_down + Pful_venn_down + plot_layout(ncol =4)
venn_down


ggsave("3_down_stream_analysis/31_venn_down.svg", plot = venn_down, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/31_venn_down.tiff", plot = venn_down, width = 10, height = 5, units = "cm", dpi = 300)


# plot regulation statistics bar plot
regulation <- read.table("2_WGCNA/merge_module/regulation_stat.txt",header = TRUE)
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

ggsave("3_down_stream_analysis/32_regulation_stat.svg", plot = regulation_stat, width = 10, height = 10, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/32_regulation_stat.tiff", plot = regulation_stat, width = 10, height = 10, units = "cm", dpi = 300)

# merge plot
merge_stat_plot <- venn_up / regulation_stat / venn_down + plot_layout(nrow = 3)
print(merge_stat_plot)

ggsave("3_down_stream_analysis/33_merge_stat_plot.svg", plot = merge_stat_plot, width = 10, height = 10, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/33_merge_stat_plot.tiff", plot = merge_stat_plot, width = 10, height = 10, units = "cm", dpi = 300)






### PLS-DA
# BiocManager::install("ropls")
library(ropls)

PLS_DA_model <- opls(metabolome_data, factor(group$group), orthoI = 0) # group is a factor vector; orthoI = 0 : PLS-DA;  orthoI = 0 : OPLS-DA;

plot(PLS_DA_model)

PLS_DA_plot_data <- PLS_DA_model@scoreMN %>% as.data.frame() %>% mutate(sample = rownames(PLS_DA_model@scoreMN), group  = PLS_DA_model@suppLs$y)

PLS_DA_plot_data$group <- factor(PLS_DA_plot_data$group, levels = c("Esin", "Bger", "Pame", "Pful"))



### pls_da plot
# Calculate the center point of each group
PLS_DA_center_points <- 
  PLS_DA_plot_data %>%
  group_by(group) %>%
  summarise(mean_wt = mean(p1), mean_mpg = mean(p2))

PLS_DA_result <- 
  PLS_DA_plot_data %>% 
  left_join(PLS_DA_center_points, by = "group")


# Basic scatter plot
PLS_DA_plot_p0 <- 
  ggplot(PLS_DA_plot_data, 
         aes(x = p1, y = p2, color = group)
  ) +
  geom_point(aes(color = group, shape = group), 
             size = 2,
             show.legend = T,
  ) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.25, color = "grey") + 
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25, color = "grey") + 
  labs(x = "Component1(33%)", 
       y = "Component2(18%)", 
       # tag = metabolome_pca_dune_adonis
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

PLS_DA_plot_p0

# Add confidence ellipse + center point connection
PLS_DA_plot_p1 <- 
  PLS_DA_plot_p0 + 
  stat_ellipse(
    data = PLS_DA_result,
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
    data = PLS_DA_result,
    aes(x = p1, y = p2, xend = mean_wt, yend = mean_mpg, color = group), 
    linetype = "dashed", 
    linewidth = 1, 
    alpha = 0.5,
    show.legend = F)

PLS_DA_plot_p1

# Use ggside package to add marginal box plots
library(ggside)
PLS_DA_plot_final <- 
  PLS_DA_plot_p1 +
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

print(PLS_DA_plot_final)

ggsave("3_down_stream_analysis/34_metabolome_pls_da.svg", plot = PLS_DA_plot_final, width = 4.5, height = 4.5, units = "cm", dpi = 300)
# ggsave("3_down_stream_analysis/34_metabolome_pls_da.tiff", plot = PLS_DA_plot_final, width = 5, height = 5, units = "cm", dpi = 300)






