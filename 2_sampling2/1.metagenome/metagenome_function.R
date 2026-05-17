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

wild_color <- "#c0c000"
lab_color  <- "#fb7d80"
fam_color  <- "#2BAA92FF"
res_color  <- "#721B3EFF"
hos_color  <- "#16317DFF"
dwelling_color <- "#751C6DFF"   ### Peridomestic enviroment
human_color <- dwelling_color   ### Peridomestic enviroment


#### KEGG KO abundance ####
# load KO abundance table
metagenome_function_KEGG_data <- 
  read.table("function/2_downstream_analysis/1_Total_KO_TPM_Filt.txt", header = TRUE, sep = "\t", row.names = 1) %>% 
  as.data.frame() %>%
  t()

## CSS normalization
library(metagenomeSeq)
metagenome_function_KEGG_data_2 <- read.table("function/1_KEGG/Total_KO_Count.txt", header = TRUE, sep = "\t", row.names = 1) %>% as.data.frame() %>% t() %>%  as.data.frame() %>% dplyr::select(where(~ any(. != 0))) 
obj <- newMRexperiment(t(metagenome_function_KEGG_data_2))
obj <- cumNorm(obj, p = cumNormStatFast(obj))
metagenome_function_KEGG_data_css_norm <- MRcounts(obj, norm = TRUE, log = TRUE, sl = 1000) %>% as.data.frame()
write.table(metagenome_function_KEGG_data_css_norm %>% rownames_to_column(var = "KO"), "function/1_KEGG/metagenome_function_KEGG_data_css_norm.txt", sep = "\t", quote = FALSE, row.names = F)


# load sample grouping table
group <- read.table('function/2_downstream_analysis/2_group.txt', sep = '\t', header =TRUE)

# Merge KO abundance table and the sample grouping table
metagenome_function_KEGG_data_grouped <- merge(metagenome_function_KEGG_data, group, by.x = "row.names", by.y = "sample")


# Calculate the group mean
metagenome_function_KEGG_mean_abundance <- 
  metagenome_function_KEGG_data_grouped %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  dplyr::select(-c(1)) %>%
  as.data.frame()

rownames(metagenome_function_KEGG_mean_abundance) <- unique(group$group) 

write.table(metagenome_function_KEGG_mean_abundance, "function/2_downstream_analysis/3_metagenome_function_KEGG_mean_abundance.tsv", quote = F, row.names = F)


### (PCA) Principal Components Analysis
metagenome_function_KEGG_pca <- prcomp(log2(metagenome_function_KEGG_data + 1) %>% as.data.frame(), # %>% dplyr::select(-c(41:50)), 
                                       scale. = TRUE)

# # Extract the coordinates (principal component scores) from the PCA results
metagenome_function_KEGG_pca_coords <- metagenome_function_KEGG_pca$x
metagenome_function_KEGG_pca_result <- merge(metagenome_function_KEGG_pca_coords, group, # group %>% filter(group != "Wild"), 
                                             by.x = "row.names", by.y = "sample")

# # Extract all principal components
metagenome_function_KEGG_pca_scores <- metagenome_function_KEGG_pca$x[, 1:3] 

# Calculate the distance matrix based on PCA scores
metagenome_function_KEGG_distance <- dist(metagenome_function_KEGG_pca_scores)
metagenome_function_KEGG_distance_matrix <- as.matrix(metagenome_function_KEGG_distance)
write.table(metagenome_function_KEGG_distance_matrix, "function/2_downstream_analysis/4_metagenome_function_KEGG_sample_distance.txt", sep = "\t", quote = FALSE)

# The "pca$sdev" records the eigenvalues of the main sorting axes in the PCA sorting results (dividing each eigenvalue by the total sum of eigenvalues gives the explanatory power of each axis)
metagenome_function_KEGG_pca_eig = sum(pmax(metagenome_function_KEGG_pca$sdev[1:3]), 0)
metagenome_function_KEGG_pca_eig_percent <- round(metagenome_function_KEGG_pca$sdev[1:3]/metagenome_function_KEGG_pca_eig*100, 3) 

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

metagenome_function_KEGG_pca_result$group <- factor(metagenome_function_KEGG_pca_result$group, level = c("Wild", "Lab", "Fam", "Res", "Hos"))



# Basic scatter plot
metagenome_function_KEGG_pca_p0 <- 
  ggplot(metagenome_function_KEGG_pca_result, 
         aes(x = PC1, y = PC2, color = group)
  ) +
  geom_point(aes(color = group, shape = group), 
             size = 1,
             show.legend = T,
  ) +
  labs(x = paste("PC 1 (", round(metagenome_function_KEGG_pca_eig_percent[1], 2), "%)", sep = ""), 
       y = paste("PC 2 (", round(metagenome_function_KEGG_pca_eig_percent[2], 2), "%)", sep = ""), 
       tag = metagenome_function_KEGG_pca_dune_adonis
       # title = "PCA: Metagenome KO abundance"
  ) +
  scale_color_manual(values = c("Wild" = wild_color, "Hos" = hos_color, "Res" = res_color, "Fam" = fam_color, "Lab" = lab_color)) +
  scale_shape_manual(values = c("Wild" = 19, "Hos" = 19, "Res" = 19, "Fam" = 19, "Lab" = 19)) +
  scale_x_continuous(position = "bottom") +
  scale_y_continuous(position = "left") +
  theme(
    plot.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.title.x.bottom = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.25, vjust = 0.5, angle = 0, lineheight = 2, color = "black"),
    axis.title.y.left = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.25, vjust = 1, angle = 90, lineheight = 2, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(4, "pt"),
    legend.frame = element_blank(),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.text.position = "right",
    legend.title = element_blank(),
    legend.title.position = "top",
    plot.tag.location = "panel",
    plot.tag = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0, angle = 0, lineheight = 1, color = "black"), 
    plot.tag.position = c(0.05, 0.65)
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
    linewidth = 0.25, 
    alpha=0.3, 
    show.legend = F
  ) +
  scale_fill_manual(values = c("Wild" = wild_color, "Hos" = hos_color, "Res" = res_color, "Fam" = fam_color, "Lab" = lab_color)) +
  geom_segment(
    data = metagenome_function_KEGG_pca_result,
    aes(x = PC1, y = PC2, xend = mean_wt, yend = mean_mpg, color = group), 
    linetype = "dashed", 
    linewidth = 0.25, 
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
                  position = position_jitter(height = 0.15, width = 0.15),
                  size = 1,
                  alpha = 0.3,
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
                  position = position_jitter(height = 0.15, width = 0.15),
                  size = 1, 
                  alpha = 0.3,
                  show.legend = F
  ) +
  scale_xsidey_discrete() +
  scale_ysidex_discrete() +
  theme(ggside.panel.scale.x = 0.32,
        ggside.panel.scale.y = 0.32,
        legend.position = c(0.88, 0.88),
  )                                                                                                                 

print(metagenome_function_KEGG_pca_p_final)

ggsave("function/2_downstream_analysis/5_metagenome_function_KEGG_PCA_PC1~PC2.svg", plot = metagenome_function_KEGG_pca_p_final, width = 6, height = 6, units = "cm", dpi = 300)
# ggsave("function/2_downstream_analysis/5_metagenome_function_KEGG_PCA_PC1~PC2.tiff", plot = metagenome_function_KEGG_pca_p_final, width = 6, height = 6, units = "cm", dpi = 300)






### sample clustering

# Computational hierarchical clustering
metagenome_function_spearman <- as.dist(1 -cor(t(metagenome_function_KEGG_data), method = "spearman"))
metagenome_function_spearman_matrix <- as.matrix(metagenome_function_spearman)

write.table(metagenome_function_spearman_matrix, "function/2_downstream_analysis/6_metagenome_function_spearman_matrix.tsv", sep = "\t", quote = FALSE, row.names = F)

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

# plot dendrogram
dev.off()
pdf("function/2_downstream_analysis/8_metagenome_function_spearman_hclust_tree.pdf")
plot(metagenome_function_KEGG_dend, main = "KO Distance Hierarchical Dendrogram")
dev.off()


### Hierarchical clustering heatmap
row_annotation <- group[,c(1,2)]
row.names(row_annotation) <- row_annotation[, 1]
colnames(row_annotation)[2] <- "Host_group"


col_annotation <- group[,c(1,2)]
row.names(col_annotation) <- col_annotation[, 1]
colnames(col_annotation)[2] <- "Host_group"


# "annot_colors" is a list that contains all the factor names and their corresponding colors.
# The names of each vector in the list correspond to the column names in the row_annotation/col_annotation matrix.
annot_colors <- list(Host_group = c(Wild = wild_color, Lab = lab_color, Fam = fam_color, Res = res_color, Hos = hos_color))

min <- min(1- metagenome_function_spearman_matrix)
max <- max(1- metagenome_function_spearman_matrix)

pdf("function/2_downstream_analysis/9_metagenome_function_KO_spearman_Hierarchical_heatmap.pdf",
    width  = 20,
    height = 20)

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











##########################################
#### Dimensionality Reduction Analysis ###
##########################################

### KEGG module regulation

# enriched KEGG modules heatmap
metagenome_up_module <- read.table("function/1_KEGG/1_all_up_module.tsv", header = TRUE, row.names = 1 ,sep = "\t") 

metagenome_up_module_zscore <- metagenome_up_module %>% dplyr::select(-c(1:8)) %>% t() %>% scale() %>% t() %>% as.data.frame()

metagenome_up_module[9:58] <- metagenome_up_module_zscore

metagenome_up_module <- 
  metagenome_up_module %>%
  mutate(log2trans = metagenome_up_module$Directionality.x..log10.p.adj.) %>%
  dplyr::select(log2trans, everything())

metagenome_up_module <- 
  metagenome_up_module %>%
  mutate(module = rownames(metagenome_up_module)) %>%
  dplyr::select(module, everything())

metagenome_up_module$module <- factor(metagenome_up_module$module, levels = metagenome_up_module$module)

up_module_mratrix <- 
  metagenome_up_module %>%
  dplyr::select(-c(1:10)) %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "module")

up_module_mratrix$module <- factor(up_module_mratrix$module, level = up_module_mratrix$module)

# wide to long
up_module_mratrix_long <-
  up_module_mratrix %>%
  pivot_longer(cols = -module,
               names_to = "sample", 
               values_to = "Enrichment_score")

sample_order <- 
  c(paste0("Metagenome_", "hos", sprintf("%02d", 10:01)),
    paste0("Metagenome_", "res", sprintf("%02d", 10:01)),
    paste0("Metagenome_", "fam", sprintf("%02d", 10:01)),
    paste0("Metagenome_", "lab", sprintf("%02d", 10:01)),
    paste0("Metagenome_", "wild", sprintf("%02d", 10:01)))


up_module_mratrix_long$sample <- factor(up_module_mratrix_long$sample, levels = sample_order)

min <- min(up_module_mratrix_long$Enrichment_score)
max <- max(up_module_mratrix_long$Enrichment_score)

up_module_mratrix_long$module

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
  annotate("rect", ymin = 40.5, ymax = 50.5, xmin = 0.5, xmax = 120.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") +
  annotate("rect", ymin = 30.5, ymax = 40.5, xmin = 120.5, xmax = 126.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 0.5, ymax = 30.5, xmin = 126.5, xmax = 136.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  scale_fill_gradient2(low = "lightblue",  # low = "#648C16FF";  low = "lightblue"
                       mid = "white", 
                       high = "#FF9898FF", # high = "#FF7200FF"; high = "#FF9898FF"
                       limit = c(min, max),
                       breaks =c(-3, 0, 3)
  ) +
  scale_x_discrete(position = "top") +
  labs(x = "Enriched metagenome modules of different position cockroaches", 
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
group$group <- factor(group$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))
group$group2 <- factor(c(rep("Dwelling", 20), rep("Lab", 10), rep("Dwelling", 10), rep("Wild", 10)), levels = c("Wild", "Lab", "Dwelling"))

group$sample <- factor(group$sample, levels = group$sample)

up_module_group_p <-
  ggplot(group, 
         aes(x = 0, y = sample, fill = group, color = group), 
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
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  theme_void()

up_module_group_p

up_module_group_p2 <-
  ggplot(group, 
         aes(x = 0, y = sample, fill = group2, color = group2)
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.05,
            lineend = "square",
            show.legend = F
  ) +
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Dwelling" = dwelling_color)) +
  theme_void()

up_module_group_p2

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
                      breaks = c(4, 12, 20)
  ) +
  scale_color_gradient(low = paletteer::paletteer_d("MetBrewer::Benedictus")[5], 
                       high = paletteer::paletteer_d("MetBrewer::Benedictus")[3],
                       breaks = c(4, 12, 20)
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
up_module_module_info <- read.table("function/1_KEGG/1_all_up_module_information.tsv", header = TRUE, sep = "\t") %>% as.data.frame()

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

up_module_module_info$module <- factor(up_module_module_info$module, levels = up_module_module_info$module) 

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
  scale_fill_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[c(1:9)],"grey"), drop = FALSE) +
  scale_color_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[c(1:9)],"grey"), drop = FALSE) +
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
  aplot::insert_left(up_module_group_p, width = 0.01) %>% 
  aplot::insert_left(up_module_group_p2, width = 0.01) %>% 
  aplot::insert_top(up_module_module_p, height = 0.05) %>%
  aplot::insert_bottom(up_module_module_info_p, height = 0.05)

up_module

# save image
ggsave("function/2_downstream_analysis/10_metagenome_function_KEGG_enrich_module.svg", plot = up_module, width = 10, height = 4.5, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/10_metagenome_function_KEGG_enrich_module.tiff", plot = up_module, width = 10, height = 4.5, units = "cm", dpi = 300)

# save all legend
up_module_heatmap_p_legend <- cowplot::get_legend(up_module_heatmap_p + theme(legend.position = "right"))
up_module_module_p_legend <- cowplot::get_legend(up_module_module_p + theme(legend.position = "right"))
up_module_module_info_p_legend <- cowplot::get_legend(up_module_module_info_p + theme(legend.position = "right"))
up_module_combined_legend <- cowplot::plot_grid(up_module_module_p_legend, up_module_heatmap_p_legend, up_module_module_info_p_legend, nrow = 3)

up_module_combined_legend

ggsave("function/2_downstream_analysis/11_metagenome_function_KEGG_enrich_module_legend.svg", plot = up_module_combined_legend, width = 10, height = 10, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/11_metagenome_function_KEGG_enrich_module_legend.tiff", plot = up_module_combined_legend, width = 10, height = 10, units = "cm", dpi = 300)


# down module heatmap
metagenome_down_module <- read.table("function/1_KEGG/2_all_down_module.tsv", header = TRUE, row.names = 1 ,sep = "\t") 

metagenome_down_module_zscore <- metagenome_down_module %>% dplyr::select(-c(1:8)) %>% t() %>% scale() %>% t() %>% as.data.frame()

metagenome_down_module[9:58] <- metagenome_down_module_zscore

metagenome_down_module <- 
  metagenome_down_module %>%
  mutate(log2trans = metagenome_down_module$Directionality.x..log10.p.adj.) %>%
  dplyr::select(log2trans, everything())

metagenome_down_module <- 
  metagenome_down_module %>%
  mutate(module = rownames(metagenome_down_module)) %>%
  dplyr::select(module, everything())

metagenome_down_module$module <- factor(metagenome_down_module$module, levels = metagenome_down_module$module)

down_module_mratrix <- 
  metagenome_down_module %>%
  dplyr::select(-c(1:10)) %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "module")

down_module_mratrix$module <- factor(down_module_mratrix$module, level = down_module_mratrix$module)

# wide to long
down_module_mratrix_long <-
  down_module_mratrix %>%
  pivot_longer(cols = -module,
               names_to = "sample", 
               values_to = "Enrichment_score")

sample_order <- 
  rev(c(paste0("Metagenome_", "hos", sprintf("%02d", 10:01)),
        paste0("Metagenome_", "res", sprintf("%02d", 10:01)),
        paste0("Metagenome_", "fam", sprintf("%02d", 10:01)),
        paste0("Metagenome_", "lab", sprintf("%02d", 10:01)),
        paste0("Metagenome_", "wild", sprintf("%02d", 10:01))))


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
  annotate("rect", ymin = 0.5, ymax = 10.5, xmin = 0.5 , xmax = 100.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 10.5, ymax = 20.5, xmin = 100.5, xmax = 113.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 20.5, ymax = 50.5, xmin = 113.5, xmax = 131.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  scale_fill_gradient2(low = "lightblue",  # low = "#648C16FF";  low = "lightblue"
                       mid = "white", 
                       high = "#FF9898FF", # high = "#FF7200FF"; high = "#FF9898FF"
                       limit = c(min, max),
                       breaks = c(-3, 0 ,3)
  ) +
  scale_x_discrete(position = "bottom") + 
  labs(x = "Depleted metagenome modules of different position cockroaches", 
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
group$group <- factor(group$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))
group$group2 <- factor(c(rep("Dwelling", 20), rep("Lab", 10), rep("Dwelling", 10), rep("Wild", 10)), levels = c("Wild", "Lab", "Dwelling"))

group$sample <- factor(group$sample, levels = group$sample)

down_module_group_p <-
  ggplot(group, 
         aes(x = 0, y = sample, fill = group, color = group), 
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
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  theme_void()

up_module_group_p 


down_module_group_p2 <-
  ggplot(group, 
         aes(x = 0, y = sample, fill = group2, color = group2), 
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
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Dwelling" = dwelling_color)) +
  theme_void()

up_module_group_p2 


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
                      breaks = c(-4, -10, -16)
  ) +
  scale_color_gradient(low = paletteer::paletteer_d("MetBrewer::Benedictus")[11], 
                       high = paletteer::paletteer_d("MetBrewer::Benedictus")[9],
                       breaks = c(-4, -10, -16)
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
down_module_module_info <- read.table("function/1_KEGG/2_all_down_module_information.tsv", header = TRUE, sep = "\t") %>% as.data.frame()

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

down_module_module_info$module <- factor(down_module_module_info$module, levels = down_module_module_info$module)

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
  aplot::insert_left(down_module_group_p, width = 0.01) %>% 
  aplot::insert_left(down_module_group_p2, width = 0.01) %>% 
  aplot::insert_bottom(down_module_module_p, height = 0.05) %>%
  aplot::insert_top(down_module_module_info_p, height = 0.05)

down_module

# save image
ggsave("function/2_downstream_analysis/12_metagenome_function_KEGG_deplete_module.svg", plot = down_module, width = 10, height = 4.5, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/12_metagenome_function_KEGG_deplete_module.tiff", plot = down_module, width = 10, height = 4.5, units = "cm", dpi = 300)


# save all legend
down_module_heatmap_p_legend <- cowplot::get_legend(down_module_heatmap_p + theme(legend.position = "right"))
down_module_module_p_legend <- cowplot::get_legend(down_module_module_p + theme(legend.position = "right"))
down_module_module_info_p_legend <- cowplot::get_legend(down_module_module_info_p + theme(legend.position = "right"))
down_module_combined_legend <- cowplot::plot_grid(down_module_heatmap_p_legend, down_module_module_p_legend, down_module_module_info_p_legend, nrow = 3)

down_module_combined_legend

ggsave("function/2_downstream_analysis/13_metagenome_function_KEGG_deplete_module_legend.svg", plot = down_module_combined_legend, width = 10, height = 10, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/13_metagenome_function_KEGG_deplete_module_legend.tiff", plot = down_module_combined_legend, width = 10, height = 10, units = "cm", dpi = 300)


### regulation module percent

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
                                              levels = c("Wild_up", "Lab_up", "Human_up"))



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
  geom_label(aes(label = sum, x = regulation, y = 60), size = 3) + 
  scale_fill_manual(values = c("Wild_up" = wild_color, "Lab_up" = lab_color, "Human_up" = human_color)) + 
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
                                                levels = c("Wild_down", "Lab_down", "Human_down"))


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
  geom_label(aes(label = sum, x = regulation, y = 50), size = 3) + 
  scale_fill_manual(values = c("Wild_down" = wild_color, "Lab_down" = lab_color, "Human_down" = human_color)) + 
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

ggsave("function/2_downstream_analysis/14_module_percent_bar_plot.svg", plot = module_percent_bar_plot, width = 6, height = 12, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/14_module_percent_bar_plot.tiff", plot = module_percent_bar_plot, width = 6, height = 12, units = "cm", dpi = 300)















### highlight top 10 enriched/depleted KEGG modules

### top 10 enriched KEGG module
top10_enriched_module <- read.table("function/1_KEGG/3_each_top10_up_module.tsv", header = TRUE, sep = "\t") 
top10_enriched_module$log2trans <- ifelse(top10_enriched_module$Directionality.x..log10.p.adj. > 6, 6, top10_enriched_module$Directionality.x..log10.p.adj.)
top10_enriched_module$mark <- ifelse(top10_enriched_module$Directionality.x..log10.p.adj. > 6, "over", "not")
top10_enriched_module <- top10_enriched_module %>% dplyr::select(log2trans, mark, everything())
top10_enriched_module$regulation <- factor(top10_enriched_module$regulation, levels = c("Wild_up", "Lab_up", "Human_up"))
top10_enriched_module$module <- factor(top10_enriched_module$module, levels = top10_enriched_module$module)
top10_enriched_module_info <- read.table("function/1_KEGG/3_each_top10_up_module_information.tsv", header = TRUE, sep = "\t") 

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

top10_enriched_module_info$merge.name <- gsub("biosynthesis", "SYN.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("metabolism", "MET.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("degradation", "DEG.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("oxidoreductase", "REDOX.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("oxidase", "OXID.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("cycle", "CYC.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("resistance", "RES.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("reduction", "RED.", top10_enriched_module_info$merge.name)
top10_enriched_module_info$merge.name <- gsub("pathway", "PWY.", top10_enriched_module_info$merge.name)


top10_enriched_module <- 
  top10_enriched_module %>% 
  left_join(top10_enriched_module_info, by = "module")


top10_enriched_module$module <- factor(top10_enriched_module$module, levels = top10_enriched_module$module)


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
  scale_fill_manual(values = c("Wild_up" = wild_color, "Lab_up" = lab_color, "Human_up" = human_color)) +
  scale_color_manual(values = c("Wild_up" = wild_color, "Lab_up" = lab_color, "Human_up" = human_color)) +
  scale_shape_manual(values = c("over" = 22, "not" = 21)) +
  labs(x = "[Directionality] x -log10(P.adj)") +
  scale_y_discrete(position = "right",
                   labels = top10_enriched_module$merge.name
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
                               color = c(rep(human_color, 10), rep(lab_color,6), rep(wild_color,10))), 
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
  ggplot(top10_enriched_module, 
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
sample_levels <- c(paste0("Metagenome_", "hos", sprintf("%02d", 10:01)),
                   paste0("Metagenome_", "res", sprintf("%02d", 10:01)),
                   paste0("Metagenome_", "fam", sprintf("%02d", 10:01)),
                   paste0("Metagenome_", "lab", sprintf("%02d", 10:01)),
                   paste0("Metagenome_", "wild", sprintf("%02d", 10:01)))



up_module_z_scores <- top10_enriched_module[, 12:61] %>% t() %>% scale() %>% t()

top10_enriched_module[, 12:61] = up_module_z_scores


# wide to long
top10_enriched_module_long <- 
  top10_enriched_module %>%
  pivot_longer(cols = -c(1:11, 62:68), names_to = "sample", values_to = "Enrich_score")

top10_enriched_module_long$sample <- factor(top10_enriched_module_long$sample, levels = rev(sample_levels))


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
  annotate("rect", xmin = 0.5, xmax = 10.5, ymin = 16.5, ymax = 26.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 10.5, xmax = 20.5, ymin = 10.5, ymax = 16.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 20.5, xmax = 50.5, ymin = 0.5, ymax = 10.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF",
                       limit = c(min, max),
                       breaks = c(-2, 0 ,2)
  ) +
  labs(fill = "Enrich Z-score") +
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

top10_enriched_module_p2


library(aplot)

top10_enriched_module_p <- 
  top10_enriched_module_p0 %>%
  insert_left(top10_enriched_module_p1, width = 0.2) %>%
  insert_left(top10_enriched_module_p2, width = 1)

top10_enriched_module_p

ggsave("function/2_downstream_analysis/15_top10_enriched_module.svg", plot = top10_enriched_module_p, width = 14, height = 9.5, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/15_top10_enriched_module.tiff", plot = top10_enriched_module_p, width = 14, height = 9.5, units = "cm", dpi = 300)





# top 10 depleted KEGG modules
top10_depleted_module <- read.table("function/1_KEGG/4_each_top10_down_module.tsv", header = TRUE, sep = "\t") 
top10_depleted_module$log2trans <- ifelse(top10_depleted_module$Directionality.x..log10.p.adj. < -6, -6, top10_depleted_module$Directionality.x..log10.p.adj.)
top10_depleted_module$mark <- ifelse(top10_depleted_module$Directionality.x..log10.p.adj. < -6, "over", "not")
top10_depleted_module <- top10_depleted_module %>% dplyr::select(log2trans, mark, everything())
top10_enriched_module$regulation <- factor(top10_enriched_module$regulation, levels = c("Wild_down", "Lab_down", "Human_down"))
top10_depleted_module$module <- factor(top10_depleted_module$module, levels = top10_depleted_module$module)
top10_depleted_module_info <- read.table("function/1_KEGG/4_each_top10_down_module_information.tsv", header = TRUE, sep = "\t") 

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


top10_depleted_module_info$merge.name <- gsub("biosynthesis", "SYN", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("metabolism", "MET", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("degradation", "DEG", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("oxidoreductase", "REDOX", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("oxidase", "OXID", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("cycle", "CYC", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("resistance", "RES", top10_depleted_module_info$merge.name)
top10_depleted_module_info$merge.name <- gsub("reduction", "RED", top10_depleted_module_info$merge.name)

top10_depleted_module_info$module <- factor(top10_depleted_module_info$module, levels = top10_depleted_module$module)

top10_depleted_module <- 
  top10_depleted_module %>%
  left_join(top10_depleted_module_info, by = "module")

top10_depleted_module$module <- factor(top10_depleted_module$module, levels = top10_depleted_module$module)

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
  scale_fill_manual(values = c("Wild_down" = wild_color, "Lab_down" = lab_color, "Human_down" = human_color)) +
  scale_color_manual(values = c("Wild_down" = wild_color, "Lab_down" = lab_color, "Human_down" = human_color)) +
  scale_shape_manual(values = c("over" = 22, "not" = 21)) +
  labs(x = "[Directionality] x -log10(P.adj)") +
  scale_y_discrete(position = "left",
                   labels = top10_depleted_module$merge.name
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
    axis.title.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 1, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, 
                               color = c(rep(human_color, 10), rep(lab_color,10), rep(wild_color,10))), 
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

top10_depleted_module_p0


# top10_depleted_module module category annotation
top10_depleted_module_p1 <- 
  ggplot(top10_depleted_module, 
         aes(x = 0, y = module, fill = Module.category, color = Module.category), 
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.25,
            lineend = "square",
            show.legend = T
  ) + 
  scale_fill_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[c(1:9)],"grey"), drop = FALSE) +
  scale_color_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[c(1:9)],"grey"), drop = FALSE) +
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
sample_levels <- c(paste0("Metagenome_", "hos", sprintf("%02d", 10:01)),
                   paste0("Metagenome_", "res", sprintf("%02d", 10:01)),
                   paste0("Metagenome_", "fam", sprintf("%02d", 10:01)),
                   paste0("Metagenome_", "lab", sprintf("%02d", 10:01)),
                   paste0("Metagenome_", "wild", sprintf("%02d", 10:01)))

down_module_z_scores <- top10_depleted_module[, 12:61] %>% t() %>% scale() %>% t()

top10_depleted_module[, 12:61] = down_module_z_scores

# wide to long
top10_depleted_module_long <- 
  top10_depleted_module %>%
  pivot_longer(cols = -c(1:11, 62:68), names_to = "sample", values_to = "Enrich_score")

top10_depleted_module_long$sample <- factor(top10_depleted_module_long$sample, levels = sample_levels)

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
  annotate("rect", xmin = 40.5, xmax = 50.5, ymin = 20.5, ymax = 30.5, fill = NA, color = "blue", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 30.5, xmax = 40.5, ymin = 10.5, ymax = 20.5, fill = NA, color = "blue", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 0.5, xmax = 30.5, ymin = 0.5, ymax = 10.5, fill = NA, color = "blue", linewidth = 0.25, lineend = "square")+
  labs(fill = "Enrich Z-score") +
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

top10_depleted_module_p2

library(aplot)

top10_depleted_module_p <- 
  top10_depleted_module_p0 %>%
  insert_right(top10_depleted_module_p1, width = 0.2) %>%
  insert_right(top10_depleted_module_p2, width = 1)

top10_depleted_module_p

ggsave("function/2_downstream_analysis/16_top10_depleted_module.svg", plot = top10_depleted_module_p, width = 14, height = 9.5, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/16_top10_depleted_module.tiff", plot = top10_depleted_module_p, width =14, height = 9.5, units = "cm", dpi = 300)



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

ggsave("function/2_downstream_analysis/17_top10_delepeted_enriched_module.svg", plot = top10_delepeted_enriched_module_plot, width = 20, height = 8, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/17_top10_delepeted_enriched_module.tiff", plot = top10_delepeted_enriched_module_plot, width = 20, height = 8, units = "cm", dpi = 300)

ggsave("function/2_downstream_analysis/17_top10_delepeted_enriched_module_0.svg", plot = top10_delepeted_enriched_module_plot_0, width = 7, height = 8, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/17_top10_delepeted_enriched_module_0.tiff", plot = top10_delepeted_enriched_module_plot_0, width = 7, height = 8, units = "cm", dpi = 300)


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

ggsave("function/2_downstream_analysis/18_top10_delepeted_enriched_module_legend.svg", plot = top10_depleted_enriched_module_plot_legend, width = 20, height = 20, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/18_top10_delepeted_enriched_module_legend.tiff", plot = top10_depleted_enriched_module_plot_legend, width = 20, height = 20, units = "cm", dpi = 300)




































### metagenome KEGG module importance ranking by lasso

# 1. load data
library(glmnet)
set.seed(20241102)
lasso_data <- read.table("function/1_KEGG/5_lasso_input.tsv", header = TRUE, row.names = 1, check.names = FALSE)
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

# 3. select alpha
best_alpha <- results$alpha[which.min(results$cv.error)]

# 4. Re-run the CV with the optimal alpha and select λ
cv_final <- cv.glmnet(x, y,
                      family  = "binomial",
                      alpha   = 0.35)

lambda_min  <- cv_final$lambda.min
lambda_1se  <- cv_final$lambda.1se
print(tibble(lambda.min = lambda_min, lambda.1se = lambda_1se))


plot(results$alpha, results$cv.error, type = "b", xlab = "Alpha", ylab = "CV Error", main = "CV Error vs Alpha")

dev.off()

pdf("function/2_downstream_analysis/19_Cross-validation.pdf",
    width  = 6,
    height = 6)

plot(cv_final)

dev.off()

# 5. plot
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
  geom_vline(xintercept = 0.35, 
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

ggsave("function/2_downstream_analysis/20_cv_plot.svg", plot = cv_plot, width = 6, height = 6, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/20_cv_plot.tiff", plot = cv_plot, width = 6, height = 6, units = "cm", dpi = 300)

# 6. Extract the modules with non-zero importance(lambda.min)
coef_mat <- as.matrix(coef(cv_final, s = "lambda.min"))
active.idx   <- which(coef_mat[, 1] != 0)
active.genes <- rownames(coef_mat)[active.idx] %>% setdiff("(Intercept)")


# 7. export result
coef.df <- coef(cv_final, s = "lambda.min") %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("module") %>%
  mutate(coef = lambda.min) %>%
  filter(module != "(Intercept)") %>%
  mutate(absCoef = abs(coef)) %>%
  filter(coef != 0) %>%
  arrange(desc(coef))

write_csv(coef.df, "function/2_downstream_analysis/21_Lasso_module_ranked_by_coef.csv")

top20_lasso_metagenome <- rbind(head(coef.df, 10), tail(coef.df, 10))

metagenome_module_info <- rbind(read.table("function/1_KEGG/1_all_up_module_information.tsv", header = TRUE, sep = "\t"),
                                read.table("function/1_KEGG/2_all_down_module_information.tsv", header = TRUE, sep = "\t")) %>%
  distinct(module, .keep_all = TRUE)

top20_lasso_metagenome <- 
  top20_lasso_metagenome %>% 
  left_join(metagenome_module_info, by = c("module" = "module")) %>%
  mutate(label = paste0(module, ":", module_name))


top20_lasso_metagenome$label <- gsub("ko00907:\\s*Pinene, camphor and geraniol degradation" , "ko00907:Pinene/Camphor/Geraniol degradation", top20_lasso_metagenome$label)

top20_lasso_metagenome$label <- gsub("M00990:\\s*Dimethylsulfoniopropionate \\(DMSP\\) degradation, demethylation pathway, DMSP => methanethiol", "M00990:DMSP degradation", top20_lasso_metagenome$label)

top20_lasso_metagenome <- top20_lasso_metagenome %>%
  separate(label, into = c("first_part", "rest"), sep = ",", remove = FALSE) %>%
  mutate(label = first_part) %>%
  dplyr::select(-first_part, -rest)

top20_lasso_metagenome$label <- gsub("biosynthesis", "SYN.", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("resistance", "RES.", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("cycle", "CYC.", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("metabolism", "MET.", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("oxidase", "OXID.", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("oxidation", "OXID.", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("oxidoreductase", "REDOX.", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("degradation", "DEG.", top20_lasso_metagenome$label)
top20_lasso_metagenome$label <- gsub("photosynthesis", "PS.", top20_lasso_metagenome$label)

top20_lasso_metagenome$label <- factor(top20_lasso_metagenome$label, levels = top20_lasso_metagenome$label)

top20_lasso_metagenome <- top20_lasso_metagenome %>% dplyr::select(-9)

write.table(top20_lasso_metagenome, "function/top10_e_d_module_information_MetaG.tsv", sep = "\t", quote = F, row.names = F)


lasso_importance_rank <- 
  ggplot(top20_lasso_metagenome)+
  geom_point(aes(x = -coef, y = label, color = label), 
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
  scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 5), labels = seq(-5, 5, by = 5)) +
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

ggsave("function/2_downstream_analysis/22_lasso_importance_rank_plot.svg", plot = lasso_importance_rank, width = 10, height = 4.5, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/22_lasso_importance_rank_plot.tiff", plot = lasso_importance_rank, width = 10, height = 4.5, units = "cm", dpi = 300)

ggsave("function/2_downstream_analysis/22_lasso_importance_rank_plot_0.svg", plot = lasso_importance_rank_0, width = 1, height = 4.5, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/22_lasso_importance_rank_plot_0.tiff", plot = lasso_importance_rank_0, width = 1, height = 4.5, units = "cm", dpi = 300)



### top 20 lasso selected feature
top20_lasso_metagenome <- rbind(head(coef.df, 10), tail(coef.df, 10))

metagenome_module_info <- rbind(read.table("function/1_KEGG/1_all_up_module_information.tsv", header = TRUE, sep = "\t"),
                                read.table("function/1_KEGG/2_all_down_module_information.tsv", header = TRUE, sep = "\t")) %>%
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

write.table(top20_lasso_metagenome, "../4.association/1_multi_omic_association/2_top20_e_d_module_information_MetaG.tsv", sep = "\t", quote = F, row.names = F)

metagenome_all_module_ssGSEA_data <- read.table("function/1_KEGG/Total_metagenome_KEGG_ssGSEA_module.tsv", sep = "\t", header = T, row.names = 1)

top20_lasso_metagenome_module_ssGSEA <- metagenome_all_module_ssGSEA_data[top20_lasso_metagenome$module, ] %>% rownames_to_column(var = "module")

write.table(top20_lasso_metagenome_module_ssGSEA, "../4.association/1_multi_omic_association/1_top20_e_d_metagenome_module.tsv", sep = "\t", quote = F, row.names = F)


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

write.table(total_lasso_metagenome, "../4.association/1_multi_omic_association/51_total_e_d_module_information_MetaG.tsv", sep = "\t", quote = F, row.names = F)

total_lasso_metagenome_module_ssGSEA <- metagenome_all_module_ssGSEA_data[coef.df$module, ] %>% rownames_to_column(var = "module")

write.table(total_lasso_metagenome_module_ssGSEA, "../4.association/1_multi_omic_association/50_total_e_d_metagenome_module.tsv", sep = "\t", quote = F, row.names = F)





### KEGG module regulation

double_volcano_data <- read.table("function/1_KEGG/6_double_volcano_log2trans.tsv", header = T, sep = "\t")

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
      Lab_vs_Wild < x1 & Human_vs_Wild >= y2 ~ 1,
      Lab_vs_Wild >= x1 & Lab_vs_Wild < x2 & Human_vs_Wild >= y2 ~ 2,
      Lab_vs_Wild >= x2 & Human_vs_Wild >= y2 ~ 3,
      Lab_vs_Wild < x1 & Human_vs_Wild >= y1 & Human_vs_Wild < y2 ~ 4,
      Lab_vs_Wild >= x1 & Lab_vs_Wild < x2 & Human_vs_Wild >= y1 & Human_vs_Wild < y2 ~ 5,
      Lab_vs_Wild >= x2 & Human_vs_Wild >= y1 & Human_vs_Wild < y2 ~ 6,
      Lab_vs_Wild < x1 & Human_vs_Wild < y1 ~ 7,
      Lab_vs_Wild >= x1 & Lab_vs_Wild < x2 & Human_vs_Wild < y1 ~ 8,
      Lab_vs_Wild >= x2 & Human_vs_Wild < y1 ~ 9,
      TRUE ~ NA_integer_
    )
  )

double_volcano_data$quadrant <- as.factor(double_volcano_data$quadrant)

double_volcano_data$size <- abs((double_volcano_data$Lab_vs_Wild + double_volcano_data$Human_vs_Wild) / 2)

quadrant_counts <- double_volcano_data %>%
  group_by(quadrant) %>%
  summarise(count = n(), .groups = 'drop')

quadrant_counts$quadrant <- factor(quadrant_counts$quadrant, levels = rev(c("3", "7", "2", "8", "4", "6", "1", "9", "5")))

quadrant_counts <- 
  quadrant_counts %>% 
  mutate(type = case_when(
    quadrant == 1 ~ "Lab_enrich & Human_deplete",
    quadrant == 2 ~ "Lab_enrich & Human_nosig",
    quadrant == 3 ~ "Lab_enrich & Human_enrich",
    quadrant == 4 ~ "Lab_nosig & Human_deplete",
    quadrant == 5 ~ "Both_nosig",
    quadrant == 6 ~ "Lab_nosig & Human_enrich",
    quadrant == 7 ~ "Lab_deplete & Human_deplete",
    quadrant == 8 ~ "Lab_deplete & Human_nosig",
    quadrant == 9 ~ "Lab_deplete & Human_enrich"
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
           color = "white",
           linewidth = 0.05
  ) +
  geom_text(aes(label = paste0(round(percent * 100, 1), "%")), 
            position = position_stack(vjust = 0.5),
            size = 1.5
  ) +
  coord_polar(theta = "y", start = -90) +
  guides(fill = guide_legend(title = "", nrow = 9, ncol = 1, reverse = TRUE))+ 
  scale_fill_manual(values = rev(c("pink", "lightblue", scater_color[3], scater_color[4], scater_color[9], scater_color[10], scater_color[1], scater_color[2], "grey80")),
                    labels = c("1" = "Lab_enrich & Human_deplete",
                               "2" = "Lab_enrich & Human_nosig",
                               "3" = "Lab_enrich & Human_enrich",
                               "4" = "Lab_nosig & Human_deplete",
                               "5" = "Both_nosig",
                               "6" = "Lab_nosig & Human_enrich",
                               "7" = "Lab_deplete & Human_deplete",
                               "8" = "Lab_deplete & Human_nosig",
                               "9" = "Lab_deplete & Human_enrich"),
                    drop = FALSE
  ) +
  theme_void()

quadrant_percent_pie_plot

ggsave("function/2_downstream_analysis/23_quadrant_percent_pie_plot.svg", plot = quadrant_percent_pie_plot, width = 1, height = 1, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/23_quadrant_percent_pie_plot.tiff", plot = quadrant_percent_pie_plot, width = 1, height = 1, units = "cm", dpi = 300)


top20_lasso_metagenome$module


# highlight module
highlight_module <- c(
  M00534 = "M00534",
  M00637 = "M00637",
  M00875 = "M00875",
  M00548 = "M00548",
  M00022 = "M00022", 
  M00915 = "M00915",
  M00638 = "M00638",
  M00595 = "M00595",
  ko00907 = "P00907",
  M00990 = "M00990",
  M00730 = "M00730", 
  M00609 = "M00609",
  M00651 = "M00651",
  M00705 = "M00705", 
  M00570 = "M00570",
  M00636 = "M00636",
  M00019 = "M00019",
  M00535 = "M00535",
  M00119 = "M00119",
  M00432 = "M00432"
)

top20_lasso_metagenome$module

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
         aes(x = Lab_vs_Wild, y = Human_vs_Wild , fill = quadrant, color = quadrant)
  ) +
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
  labs(x = "Lab vs Wild: [Dir] x -log10(P.adj)", 
       y = "Dewilling vs Wild: [Dir] x -log10(P.adj)"
  ) +
  scale_color_manual(values = c(scater_color[1], scater_color[3], "pink", scater_color[9], "grey80", scater_color[10], "lightblue", scater_color[4] ,scater_color[2]), drop = FALSE) +
  scale_fill_manual(values = c(scater_color[1], scater_color[3], "pink", scater_color[9], "grey80", scater_color[10], "lightblue", scater_color[4] ,scater_color[2]), drop = FALSE) +
  scale_size_continuous(range = c(1, 2))+
  scale_alpha_continuous(range = c(0.4 ,1))+
  # theme_void() + 
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25, linetype = "solid"),
    panel.background = element_blank(),
    axis.title.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
    axis.title.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1.2, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1.2, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length = unit(2, "pt"),
  )

print(double_volcano_plot)

ggsave("function/2_downstream_analysis/24_double_volcano_plot.svg", plot = double_volcano_plot, width = 6, height = 6, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/24_double_volcano_plot.tiff", plot = double_volcano_plot, width = 6, height = 6, units = "cm", dpi = 300)













### venn

# set color 
Wild_venn_color <- generate_palette(wild_color, modification = "go_lighter", n_colours = 3, view_palette = TRUE, view_labels = FALSE)
Lab_venn_color <- generate_palette(lab_color, modification = "go_lighter", n_colours = 3, view_palette = TRUE, view_labels = FALSE)
Human_venn_color <- generate_palette(human_color, modification = "go_lighter", n_colours = 3, view_palette = TRUE, view_labels = FALSE)


# Wild enrich
Wild_venn_dat_up  <- read.delim("function/1_KEGG/Wild_metagenome_KEGG_ssGSEA_module_each_up.tsv")
Wild_venn_up_list <- as.list(Wild_venn_dat_up)

Wild_venn_up <- ggvenn(
  data = Wild_venn_up_list,
  show_elements = F,
  label_sep = "\n",
  show_percentage = F,
  digits = 2,
  fill_color = Wild_venn_color,
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_alpha = 0,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = NA,
  set_name_size = 0,
  text_color = "black",
  text_size = 2
)

Wild_venn_up

# Lab enrich
Lab_venn_dat_up  <- read.delim("function/1_KEGG/Lab_metagenome_KEGG_ssGSEA_module_each_up.tsv")
Lab_venn_up_list <- as.list(Lab_venn_dat_up )

Lab_venn_up <- ggvenn(
  data = Lab_venn_up_list,
  show_elements = F,
  label_sep = "\n",
  show_percentage = F,
  digits = 2,
  fill_color = Lab_venn_color,
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_alpha = 0,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = NA,
  set_name_size = 0,
  text_color = "black",
  text_size = 2
)

Lab_venn_up

# Human enrich
Human_venn_dat_up  <- read.delim("function/1_KEGG/Human_metagenome_KEGG_ssGSEA_module_each_up.tsv")
Human_venn_up_list <- as.list(Human_venn_dat_up)

Human_venn_up <- ggvenn(
  data = Human_venn_up_list,
  show_elements = F,
  label_sep = "\n",
  show_percentage = F,
  digits = 2,
  fill_color = Human_venn_color,
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_alpha = 0,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = NA,
  set_name_size = 0,
  text_color = "black",
  text_size = 2
)

Human_venn_up


library(aplot)

venn_up <- Wild_venn_up + Lab_venn_up + Human_venn_up + plot_layout(ncol = 3)

venn_up

ggsave("function/2_downstream_analysis/25_venn_up.svg", plot = venn_up, width = 10, height = 5, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/25_venn_up.tiff", plot = venn_up, width = 10, height = 5, units = "cm", dpi = 300)




# Wild deplete
Wild_venn_dat_down  <- read.delim("function/1_KEGG/Wild_metagenome_KEGG_ssGSEA_module_each_down.tsv")
Wild_venn_down_list <- as.list(Wild_venn_dat_down)

Wild_venn_down <- ggvenn(
  data = Wild_venn_down_list,
  show_elements = F,
  label_sep = "\n",
  show_percentage = F,
  digits = 2,
  fill_color = Wild_venn_color,
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_alpha = 0,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = NA,
  set_name_size = 0,
  text_color = "black",
  text_size = 2
)

Wild_venn_down

# Lab deplete
Lab_venn_dat_down  <- read.delim("function/1_KEGG/Lab_metagenome_KEGG_ssGSEA_module_each_down.tsv")
Lab_venn_down_list <- as.list(Lab_venn_dat_down )

Lab_venn_down <- ggvenn(
  data = Lab_venn_down_list,
  show_elements = F,
  label_sep = "\n",
  show_percentage = F,
  digits = 2,
  fill_color = Lab_venn_color,
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_alpha = 0,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = NA,
  set_name_size = 0,
  text_color = "black",
  text_size = 2
)

Lab_venn_down

# Human deplete
Human_venn_dat_down  <- read.delim("function/1_KEGG/Human_metagenome_KEGG_ssGSEA_module_each_down.tsv")
Human_venn_down_list <- as.list(Human_venn_dat_down )

Human_venn_down <- ggvenn(
  data = Human_venn_down_list,
  show_elements = F,
  label_sep = "\n",
  show_percentage = F,
  digits = 2,
  fill_color = Human_venn_color,
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_alpha = 0,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = NA,
  set_name_size = 0,
  text_color = "black",
  text_size = 2
)

Human_venn_down


library(aplot)

venn_down <- Wild_venn_down + Lab_venn_down + Human_venn_down + plot_layout(ncol = 3)

venn_down

ggsave("function/2_downstream_analysis/26_venn_down.svg", plot = venn_down, width = 10, height = 5, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/26_venn_down.tiff", plot = venn_down, width = 10, height = 5, units = "cm", dpi = 300)


# plot regulation statistics bar plot
regulation <- read.table("function/1_KEGG/regulation_stat.txt",header = TRUE)
regulation$X <- factor(regulation$X, levels = unique(regulation$X))
regulation$group <- factor(regulation$group, levels = c("Wild", "Lab", "Human"))

regulation_stat <- 
  ggplot(regulation, 
         aes(x = X, y = number, fill = X)
  ) +
  facet_wrap(~ group, 
             scale = "free_x", 
             ncol = 3, 
             nrow = 1
  ) +
  scale_fill_manual(values = c(Wild_venn_color, Lab_venn_color, Human_venn_color)) + 
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

ggsave("function/2_downstream_analysis/27_regulation_stat.svg", plot = regulation_stat, width = 10, height = 10, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/27_regulation_stat.tiff", plot = regulation_stat, width = 10, height = 10, units = "cm", dpi = 300)


# merge plot
merge_stat_plot <- venn_up / regulation_stat / venn_down + plot_layout(nrow = 3)
print(merge_stat_plot)


ggsave("function/2_downstream_analysis/28_merge_stat_plot.svg", plot = merge_stat_plot, width = 10, height = 10, units = "cm", dpi = 300)
ggsave("function/2_downstream_analysis/28_merge_stat_plot.tiff", plot = merge_stat_plot, width = 10, height = 10, units = "cm", dpi = 300)










### taxonomy contribution of BCAA gene
library(glmnet)
library(caret)
library(dplyr)

### scale_II
# BCAA_KO <- read.table("function/2_downstream_analysis/34_BCAA_KO_abundance.txt", header = T, sep = "\t", row.names = 1)

K00052 <- read.table("function/2_downstream_analysis/K00052_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:57, names_to = "Sample", values_to = "Abundance") %>% 
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>% 
  arrange(desc(mean)) %>%
  filter(!(phylum %in% c("p_unclassify", "p_unclassify_d__Bacteria"))) %>% 
  filter(genus != "g__Blattabacterium") %>% 
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:55, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% sprintf("fam%02d", 1:10) ~ "fam",
    Sample %in% sprintf("res%02d", 1:10) ~ "res",
    Sample %in% sprintf("hos%02d", 1:10) ~ "hos",
    Sample %in% sprintf("lab%02d", 1:10) ~ "lab",
    Sample %in% sprintf("wild%02d", 1:10) ~ "wild",
    TRUE ~ NA_character_
  )) %>% 
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K00052")


K00053 <- read.table("function/2_downstream_analysis/K00053_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:57, names_to = "Sample", values_to = "Abundance") %>% 
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>% 
  arrange(desc(mean)) %>%
  filter(!(phylum %in% c("p_unclassify", "p_unclassify_d__Bacteria"))) %>% 
  filter(genus != "g__Blattabacterium") %>% 
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:55, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% sprintf("fam%02d", 1:10) ~ "fam",
    Sample %in% sprintf("res%02d", 1:10) ~ "res",
    Sample %in% sprintf("hos%02d", 1:10) ~ "hos",
    Sample %in% sprintf("lab%02d", 1:10) ~ "lab",
    Sample %in% sprintf("wild%02d", 1:10) ~ "wild",
    TRUE ~ NA_character_
  )) %>% 
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K00053")


K00826 <- read.table("function/2_downstream_analysis/K00826_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:57, names_to = "Sample", values_to = "Abundance") %>% 
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>% 
  arrange(desc(mean)) %>%
  filter(!(phylum %in% c("p_unclassify", "p_unclassify_d__Bacteria"))) %>% 
  filter(genus != "g__Blattabacterium") %>% 
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:55, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% sprintf("fam%02d", 1:10) ~ "fam",
    Sample %in% sprintf("res%02d", 1:10) ~ "res",
    Sample %in% sprintf("hos%02d", 1:10) ~ "hos",
    Sample %in% sprintf("lab%02d", 1:10) ~ "lab",
    Sample %in% sprintf("wild%02d", 1:10) ~ "wild",
    TRUE ~ NA_character_
  )) %>% 
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K00826")


K01649 <- read.table("function/2_downstream_analysis/K01649_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:57, names_to = "Sample", values_to = "Abundance") %>% 
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>% 
  arrange(desc(mean)) %>%
  filter(!(phylum %in% c("p_unclassify", "p_unclassify_d__Bacteria"))) %>% 
  filter(genus != "g__Blattabacterium") %>% 
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:55, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% sprintf("fam%02d", 1:10) ~ "fam",
    Sample %in% sprintf("res%02d", 1:10) ~ "res",
    Sample %in% sprintf("hos%02d", 1:10) ~ "hos",
    Sample %in% sprintf("lab%02d", 1:10) ~ "lab",
    Sample %in% sprintf("wild%02d", 1:10) ~ "wild",
    TRUE ~ NA_character_
  )) %>% 
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K01649")


K01652 <- read.table("function/2_downstream_analysis/K01652_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:57, names_to = "Sample", values_to = "Abundance") %>% 
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>% 
  arrange(desc(mean)) %>%
  filter(!(phylum %in% c("p_unclassify", "p_unclassify_d__Bacteria"))) %>% 
  filter(genus != "g__Blattabacterium") %>% 
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:55, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% sprintf("fam%02d", 1:10) ~ "fam",
    Sample %in% sprintf("res%02d", 1:10) ~ "res",
    Sample %in% sprintf("hos%02d", 1:10) ~ "hos",
    Sample %in% sprintf("lab%02d", 1:10) ~ "lab",
    Sample %in% sprintf("wild%02d", 1:10) ~ "wild",
    TRUE ~ NA_character_
  )) %>% 
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K01652")


K01653 <- read.table("function/2_downstream_analysis/K01653_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:57, names_to = "Sample", values_to = "Abundance") %>% 
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>% 
  arrange(desc(mean)) %>%
  filter(!(phylum %in% c("p_unclassify", "p_unclassify_d__Bacteria"))) %>% 
  filter(genus != "g__Blattabacterium") %>% 
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:55, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% sprintf("fam%02d", 1:10) ~ "fam",
    Sample %in% sprintf("res%02d", 1:10) ~ "res",
    Sample %in% sprintf("hos%02d", 1:10) ~ "hos",
    Sample %in% sprintf("lab%02d", 1:10) ~ "lab",
    Sample %in% sprintf("wild%02d", 1:10) ~ "wild",
    TRUE ~ NA_character_
  )) %>% 
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K01653")


K01687 <- read.table("function/2_downstream_analysis/K01687_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:57, names_to = "Sample", values_to = "Abundance") %>% 
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>% 
  arrange(desc(mean)) %>%
  filter(!(phylum %in% c("p_unclassify", "p_unclassify_d__Bacteria"))) %>% 
  filter(genus != "g__Blattabacterium") %>% 
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:55, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% sprintf("fam%02d", 1:10) ~ "fam",
    Sample %in% sprintf("res%02d", 1:10) ~ "res",
    Sample %in% sprintf("hos%02d", 1:10) ~ "hos",
    Sample %in% sprintf("lab%02d", 1:10) ~ "lab",
    Sample %in% sprintf("wild%02d", 1:10) ~ "wild",
    TRUE ~ NA_character_
  )) %>% 
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K01687")


K01703 <- read.table("function/2_downstream_analysis/K01703_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:57, names_to = "Sample", values_to = "Abundance") %>% 
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>% 
  arrange(desc(mean)) %>%
  filter(!(phylum %in% c("p_unclassify", "p_unclassify_d__Bacteria"))) %>% 
  filter(genus != "g__Blattabacterium") %>% 
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:55, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% sprintf("fam%02d", 1:10) ~ "fam",
    Sample %in% sprintf("res%02d", 1:10) ~ "res",
    Sample %in% sprintf("hos%02d", 1:10) ~ "hos",
    Sample %in% sprintf("lab%02d", 1:10) ~ "lab",
    Sample %in% sprintf("wild%02d", 1:10) ~ "wild",
    TRUE ~ NA_character_
  )) %>% 
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K01703")


K01704 <- read.table("function/2_downstream_analysis/K01704_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:57, names_to = "Sample", values_to = "Abundance") %>% 
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>% 
  arrange(desc(mean)) %>%
  filter(!(phylum %in% c("p_unclassify", "p_unclassify_d__Bacteria"))) %>% 
  filter(genus != "g__Blattabacterium") %>% 
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:55, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% sprintf("fam%02d", 1:10) ~ "fam",
    Sample %in% sprintf("res%02d", 1:10) ~ "res",
    Sample %in% sprintf("hos%02d", 1:10) ~ "hos",
    Sample %in% sprintf("lab%02d", 1:10) ~ "lab",
    Sample %in% sprintf("wild%02d", 1:10) ~ "wild",
    TRUE ~ NA_character_
  )) %>% 
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K01704")


K01754  <- read.table("function/2_downstream_analysis/K01754_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:57, names_to = "Sample", values_to = "Abundance") %>% 
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>% 
  arrange(desc(mean)) %>%
  filter(!(phylum %in% c("p_unclassify", "p_unclassify_d__Bacteria"))) %>% 
  filter(genus != "g__Blattabacterium") %>% 
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:55, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% sprintf("fam%02d", 1:10) ~ "fam",
    Sample %in% sprintf("res%02d", 1:10) ~ "res",
    Sample %in% sprintf("hos%02d", 1:10) ~ "hos",
    Sample %in% sprintf("lab%02d", 1:10) ~ "lab",
    Sample %in% sprintf("wild%02d", 1:10) ~ "wild",
    TRUE ~ NA_character_
  )) %>% 
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K01754")


K09011 <- read.table("function/2_downstream_analysis/K09011_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:55, names_to = "Sample", values_to = "Abundance") %>% 
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>% 
  arrange(desc(mean)) %>%
  filter(!(phylum %in% c("p_unclassify", "p_unclassify_d__Bacteria"))) %>% 
  filter(genus != "g__Blattabacterium") %>% 
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:53, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% sprintf("fam%02d", 1:10) ~ "fam",
    Sample %in% sprintf("res%02d", 1:10) ~ "res",
    Sample %in% sprintf("hos%02d", 1:10) ~ "hos",
    Sample %in% sprintf("lab%02d", 1:10) ~ "lab",
    Sample %in% sprintf("wild%02d", 1:10) ~ "wild",
    TRUE ~ NA_character_
  )) %>% 
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K09011")


K17989 <- read.table("function/2_downstream_analysis/K17989_species_absolute_contribution.tsv", header = T, sep = "\t") %>%
  separate(Species.Path, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  pivot_longer(8:25, names_to = "Sample", values_to = "Abundance") %>% 
  group_by(phylum, class, order, family, genus, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>% 
  mutate(mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>% 
  arrange(desc(mean)) %>%
  filter(!(phylum %in% c("p_unclassify", "p_unclassify_d__Bacteria"))) %>% 
  filter(genus != "g__Blattabacterium") %>% 
  # head(.,5) %>% 
  dplyr::select(-mean) %>% 
  pivot_longer(6:23, names_to = "Sample", values_to = "TPM") %>%
  mutate(group = case_when(
    Sample %in% sprintf("fam%02d", 1:10) ~ "fam",
    Sample %in% sprintf("res%02d", 1:10) ~ "res",
    Sample %in% sprintf("hos%02d", 1:10) ~ "hos",
    Sample %in% sprintf("lab%02d", 1:10) ~ "lab",
    Sample %in% sprintf("wild%02d", 1:10) ~ "wild",
    TRUE ~ NA_character_
  )) %>% 
  group_by(phylum, class, order, family, genus, group) %>% 
  summarise(mean_TPM = mean(TPM, na.rm = TRUE)) %>% 
  group_by(genus) %>%
  mutate(zscore = scales::rescale(mean_TPM, to = c(-2, 2))) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos"))) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus)))) %>%
  mutate(KO = "K17989")

BCAA_KO_taxonomy <-  
  rbind(K00052, K00053, K00826, K01649, K01652, K01653, K01687, K01703, K01704, K01754, K09011)

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
  scale_fill_manual(values = c("p__Actinomycetota" = "#8dcfbb", "p__Bacteroidota" = "#76c4ea",  "p__Bacillota"  = "#ea8c8c", "p__Pseudomonadota" = "#ebae80", "p__Desulfobacterota_I" =  "#bab8d9", "#cdcbcc") )+
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

ggsave("function/2_downstream_analysis/35_BCAA_KO_taxonomy_top10_bubble.svg", plot = BCAA_KO_taxonomy_top10_bubble, width = 10, height = 24, units = "cm", dpi = 300)

BCAA_KO_taxonomy_top10_bubble_legend <- cowplot::get_legend(BCAA_KO_taxonomy_top10_bubble + theme(legend.position = "right"))
BCAA_KO_taxonomy_top10_bubble_legend <-  cowplot::plot_grid(BCAA_KO_taxonomy_top10_bubble_legend, ncol = 1)
BCAA_KO_taxonomy_top10_bubble_legend

ggsave("function/2_downstream_analysis/36_BCAA_KO_taxonomy_top10_bubble_legend.svg", plot = BCAA_KO_taxonomy_top10_bubble_legend, width = 10, height = 24, units = "cm", dpi = 300)



### 
library(scatterpie)
library(PieGlyph)
library(scales)

BCAA_KO_taxonomy_top10_trans <-
  BCAA_KO_taxonomy %>%
  filter(genus %in% top10$genus) %>%
  pivot_wider(id_cols = c(phylum, class, order, family, genus, KO), names_from = group, values_from = mean_TPM) %>%
  mutate(total = wild + lab + fam + res + hos) %>%
  mutate(genus = factor(genus, levels = top10$genus))

BCAA_KO_taxonomy_top10_trans <-
  BCAA_KO_taxonomy_top10_trans %>%
  mutate(KO =  factor(KO, levels = rev(unique(BCAA_KO_taxonomy_top10_trans$KO))))

p <- 
  ggplot(BCAA_KO_taxonomy_top10_trans) +
  geom_pie_glyph(
    aes(x = genus, radius = total, y = KO),
    slices = c( "lab", "fam", "res", "hos", "wild"),
    colour = "white",
    linewidth = 0.1
  ) +
  scale_radius(range = c(0.1, 0.35), name = "TPM") +
  scale_fill_manual(values = c("wild" = alpha(wild_color, 0.5), "lab" = alpha(lab_color, 0.5), "fam" = alpha(fam_color, 0.5), "res" = alpha(res_color, 0.5), "hos" = alpha(hos_color, 0.5))) +
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


label <- BCAA_KO_taxonomy_top10_trans %>% mutate(labels = sub("^g__", "", genus)) %>% dplyr::select(genus, labels) %>% unique()

unique(BCAA_KO_taxonomy_top10_trans$phylum)

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
  scale_fill_manual(values = c("p__Bacteroidota" =  "#76c4ea", "p__Bacillota" =  "#ea8c8c", "p__Desulfobacterota_I" =  "#bab8d9",  "p__Pseudomonadota" = "#88d1b6"))+
  scale_color_manual(values = c("p__Bacteroidota" =  "#76c4ea", "p__Bacillota" =  "#ea8c8c", "p__Desulfobacterota_I" =  "#bab8d9", "p__Pseudomonadota"= "#88d1b6"))+
  theme_void() + 
  theme(
    legend.position = "right"
  )

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
    slices = c( "lab", "fam", "res", "hos", "wild"),
    colour = "white",
    linewidth = 0.1,
    show.legend = F
  ) +
  scale_radius(range = c(0.1, 0.35), name = "TPM") +
  scale_fill_manual(values = c("wild" = alpha(wild_color, 0.75), "lab" = alpha(lab_color, 0.75), "fam" = alpha(fam_color, 0.75), 
                               "res" = alpha(res_color, 0.75), "hos" = alpha(hos_color, 0.75))) +
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
    slices = c( "lab", "fam", "res", "hos", "wild"),
    colour = "white",
    linewidth = 0.1,
    show.legend = F
  ) +
  scale_radius(range = c(0.1, 0.35), name = "TPM") +
  scale_fill_manual(values = c("wild" = alpha(wild_color, 0.75), "lab" = alpha(lab_color, 0.75), "fam" = alpha(fam_color, 0.75), 
                               "res" = alpha(res_color, 0.75), "hos" = alpha(hos_color, 0.75))) +
  theme_void()

KO_total_p


p2 <- p %>% insert_top(Genus_total_p, height = 0.1) %>% insert_right(KO_total_p, width = 0.1)


p2




### correlation scater plot
total <- 
  read.table("../1.metagenome/taxonomy/1_Total_coverm_contigs_Count_species.tsv", header = T, sep = "\t") %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>%
  t() %>%
  as.vector()

genus_abundance <- 
  read.table("../1.metagenome/taxonomy/1_Total_coverm_contigs_Count_species.tsv", header = T, sep = "\t") %>%
  separate(species, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "drop", fill = "right") %>%
  group_by(phylum, class, order, family, genus) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>%
  filter(genus %in% c("g__QVMG01", "g__Frigididesulfovibrio", "g__Dysgonomonas", "g__Bacteroides")) %>%
  rename_with(~ str_remove(., "^Metagenome_"), .cols = everything()) %>%
  dplyr::select(-c(1:4))

genus_abundance[2:51] <- genus_abundance[2:51] / total 

order <- names(genus_abundance)[2:51] %>% sort()

genus_abundance <- genus_abundance %>% dplyr::select(1, all_of(order)) %>% column_to_rownames(var = "genus") 

BCAA_metabolite <- 
  read.table("../2.metabolome/1_metabolome_matrix/metabolome_LOESS_PQN_missForest_matrix_noQC.tsv", header = T, sep = "\t") %>%
  filter(ID %in% c("pos_1388", "neg_2285", "pos_24501")) %>%
  mutate(MS2_name = case_when(
    ID == "pos_1388" ~ "L-Leucine",
    ID == "neg_2285" ~ "L-Isoleucine",
    ID == "pos_24501" ~ "L-Valine"
  )) %>% 
  dplyr::select(ID, MS2_name, everything()) %>%
  dplyr::select(-ID)

BCAA_metabolite <- BCAA_metabolite %>% dplyr::select(1, all_of(order)) %>% column_to_rownames(var = "MS2_name") 

HostT_pathway <- 
  read.table("../3.transcriptome/1_OG/Total_transcriptome_KEGG_ssGSEA_pathway.tsv", header = T, sep = "\t") %>%
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
  mutate(group = c(rep("fam", 10),  rep("hos", 10), rep("lab", 10), rep("res", 10), rep("wild", 10))) %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos")))


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
  stat_ellipse(level = 0.95) +
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

pdf("function/2_downstream_analysis//39_microbiome_corr.pdf",
    width  = 6,
    height = 6)
print(micro_corr_heatmap)
dev.off()




### genus_abundance_box_plot
g_QVMG01_abundance_box_plot <-
  ggplot(Micro_MetaB_HostT) + 
  geom_boxplot(aes(x = group, y = g__QVMG01, fill = group, color = group),
               alpha = 0.3,
               width = 0.8,
               staplewidth = 0.5,
               linewidth = 0.25,
               show.legend = F,
               outliers = F
  )+
  # geom_point(aes(x = group, y = g__QVMG01), size = 0.3)+
  scale_fill_manual(values = c("wild" = ggplot2::alpha(wild_color, 0.75), 
                               "lab" = ggplot2::alpha(lab_color, 0.75), 
                               "fam" = ggplot2::alpha(fam_color, 0.75), 
                               "res" = ggplot2::alpha(res_color, 0.75), 
                               "hos" = ggplot2::alpha(hos_color, 0.75))
  ) +
  scale_color_manual(values = c("wild" = ggplot2::alpha(wild_color, 0.75), 
                               "lab" = ggplot2::alpha(lab_color, 0.75), 
                               "fam" = ggplot2::alpha(fam_color, 0.75), 
                               "res" = ggplot2::alpha(res_color, 0.75), 
                               "hos" = ggplot2::alpha(hos_color, 0.75))
  ) +
  geom_signif(aes(x = group, y = g__QVMG01),
              comparisons = list(c("wild", "lab"),
                                 c("wild", "fam"), 
                                 c("wild", "res"), 
                                 c("wild", "hos")),
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
              y_position = c(0.080, 0.085, 0.090, 0.095)
              
  ) +
  coord_cartesian(ylim = c(0.0, 0.12)) + 
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
  geom_boxplot(aes(x = group, y = g__Bacteroides, fill = group, color = group),
               alpha = 0.3,
               width = 0.8,
               staplewidth = 0.5,
               linewidth = 0.25,
               show.legend = F,
               outliers = F
  )+
  # geom_point(aes(x = group, y = g__Bacteroides), size = 0.3)+
  scale_fill_manual(values = c("wild" = ggplot2::alpha(wild_color, 0.75), 
                               "lab" = ggplot2::alpha(lab_color, 0.75), 
                               "fam" = ggplot2::alpha(fam_color, 0.75), 
                               "res" = ggplot2::alpha(res_color, 0.75), 
                               "hos" = ggplot2::alpha(hos_color, 0.75))
  ) +
  scale_color_manual(values = c("wild" = ggplot2::alpha(wild_color, 0.75), 
                                "lab" = ggplot2::alpha(lab_color, 0.75), 
                                "fam" = ggplot2::alpha(fam_color, 0.75), 
                                "res" = ggplot2::alpha(res_color, 0.75), 
                                "hos" = ggplot2::alpha(hos_color, 0.75))
  ) +
  geom_signif(aes(x = group, y = g__Bacteroides),
              comparisons = list(c("wild", "lab"),
                                 c("wild", "fam"), 
                                 c("wild", "res"), 
                                 c("wild", "hos")),
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
              y_position = c(0.105, 0.110, 0.115, 0.120)
              
  ) +
  coord_cartesian(ylim = c(0.0, 0.13)) + 
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
  geom_boxplot(aes(x = group, y = g__Frigididesulfovibrio, fill = group, color = group),
             alpha = 0.3,
             width = 0.8,
             staplewidth = 0.5,
             linewidth = 0.25,
             show.legend = F,
             outliers = F
  )+
  # geom_point(aes(x = group, y = g__Frigididesulfovibrio), size = 0.3)+
  scale_fill_manual(values = c("wild" = ggplot2::alpha(wild_color, 0.75), 
                               "lab" = ggplot2::alpha(lab_color, 0.75), 
                               "fam" = ggplot2::alpha(fam_color, 0.75), 
                               "res" = ggplot2::alpha(res_color, 0.75), 
                               "hos" = ggplot2::alpha(hos_color, 0.75))
  ) +
  scale_color_manual(values = c("wild" = ggplot2::alpha(wild_color, 0.75), 
                                "lab" = ggplot2::alpha(lab_color, 0.75), 
                                "fam" = ggplot2::alpha(fam_color, 0.75), 
                                "res" = ggplot2::alpha(res_color, 0.75), 
                                "hos" = ggplot2::alpha(hos_color, 0.75))
  ) +
  geom_signif(aes(x = group, y = g__Frigididesulfovibrio),
              comparisons = list(c("wild", "lab"),
                                 c("wild", "fam"), 
                                 c("wild", "res"), 
                                 c("wild", "hos")),
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
              y_position = c(0.045, 0.050, 0.055, 0.060)
              
  ) +
  coord_cartesian(ylim = c(0.0, 0.065)) +
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
  geom_boxplot(aes(x = group, y = g__Dysgonomonas, fill = group, color = group),
             alpha = 0.3,
             width = 0.8,
             staplewidth = 0.5,
             linewidth = 0.25,
             show.legend = F,
             outlier.shape = NA
  )+
  # geom_point(aes(x = group, y = g__Dysgonomonas), size = 0.3)+
  scale_fill_manual(values = c("wild" = ggplot2::alpha(wild_color, 0.75), 
                               "lab" = ggplot2::alpha(lab_color, 0.75), 
                               "fam" = ggplot2::alpha(fam_color, 0.75), 
                               "res" = ggplot2::alpha(res_color, 0.75), 
                               "hos" = ggplot2::alpha(hos_color, 0.75))
  ) +
  scale_color_manual(values = c("wild" = ggplot2::alpha(wild_color, 0.75), 
                                "lab" = ggplot2::alpha(lab_color, 0.75), 
                                "fam" = ggplot2::alpha(fam_color, 0.75), 
                                "res" = ggplot2::alpha(res_color, 0.75), 
                                "hos" = ggplot2::alpha(hos_color, 0.75))
  ) +
  geom_signif(aes(x = group, y = g__Dysgonomonas),
              comparisons = list(c("wild", "lab"),
                                 c("wild", "fam"), 
                                 c("wild", "res"), 
                                 c("wild", "hos")),
              map_signif_level = function(p) {
                ifelse(p < 0.001, "***", 
                       ifelse(p < 0.01, "**", 
                              ifelse(p < 0.05, "*", "ns")))
              },
              test = "wilcox.test",
              test.args = list(exact = FALSE), 
              textsize = 4,
              vjust = 0.6,
              tip_length = 0,
              size = 0.25,
              # step_increase = 0.05,
              y_position = c(0.055, 0.060, 0.065, 0.070)
              
  ) +
  coord_cartesian(ylim = c(0.0, 0.1)) + 
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



