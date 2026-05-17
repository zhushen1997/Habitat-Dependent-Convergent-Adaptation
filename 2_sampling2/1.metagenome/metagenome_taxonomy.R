### Part 1: taxonomy analysis ###

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

# Load taxonomy abundance matrix
taxonomy_abundance <- read.table("taxonomy/1_Total_coverm_contigs_Count_species.tsv", header = TRUE, sep = "\t", row.names = 1)

taxonomy_abundance <- taxonomy_abundance[rowSums(taxonomy_abundance) > 0, ]

# Calculate taxonomy relative abundance
relative_taxonomy_abundance <- t(sweep(taxonomy_abundance, 2, colSums(taxonomy_abundance), "/")) %>% as.data.frame()
rowSums(relative_taxonomy_abundance)

# Load the group information table
group <- read.table('taxonomy/2_group.txt', sep = '\t', header =TRUE)

# Merge the relative abundance table and the sample grouping table
relative_taxonomy_abundance_grouped <- merge(relative_taxonomy_abundance, group, by.x = "row.names", by = "sample")

# Calculate the group mean
relative_taxonomy_mean_abundance <- relative_taxonomy_abundance_grouped %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  as.data.frame()

# Export the result file
relative_taxonomy_mean_abundance <- 
  relative_taxonomy_mean_abundance %>% 
  as.data.frame()  %>% 
  column_to_rownames(var = "group") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "taxonomy")

write.table(relative_taxonomy_mean_abundance, "taxonomy/3_relative_taxonomy_mean_abundance.tsv", row.names = FALSE, quote = FALSE, sep = "\t")


#################################
#### Alpha diversity analysis ###
#################################
α_diversity <- vegan::estimateR(t(taxonomy_abundance %>% as.matrix()))
richness <- vegan::specnumber(t(taxonomy_abundance))
shannon <- vegan::diversity(t(taxonomy_abundance), "shannon")
simpson <- vegan::diversity(t(taxonomy_abundance), "simpson")
inv_simpson <- vegan::diversity(t(taxonomy_abundance), index = "invsimpson")

alpha_diversity <- data.frame(
  Sample = rownames(t(taxonomy_abundance)),
  Richness = richness,
  Shannon = shannon,
  Simpson = simpson,
  Chao1 = α_diversity[2, ],
  group = relative_taxonomy_abundance_grouped$group
)

alpha_diversity$group <- factor(alpha_diversity$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

alpha_diversity_group <- 
  alpha_diversity %>%
  dplyr::select(-1) %>%
  group_by(group) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(-group, names_to = "alpha diversity index", values_to = "value")%>% 
  group_by(`alpha diversity index`) %>%
  mutate(scale_index = scales::rescale(value, to = c(-2, 2)))

alpha_diversity_group$`alpha diversity index` <- factor(alpha_diversity_group$`alpha diversity index`, levels = rev(c("Richness", "Shannon", "Simpson", "Chao1")))

# Export the result file
write.table(alpha_diversity, "taxonomy/4_alpha_diversity.tsv", row.names = FALSE, quote = FALSE)


 

### heatmap 
min(alpha_diversity_group$scale_index)
max(alpha_diversity_group$scale_index)

alpha_diversity_heatmap <- 
  ggplot(alpha_diversity_group, 
         aes(x = group, y = `alpha diversity index`, fill = scale_index)
  ) +
  geom_tile(
    height = 1,
    color = "white",
    linewidth = 0.15,
    linetype = 1,
    show.legend = T
  ) +
  scale_fill_gradient2(low = "#648C16FF", 
                       mid = "white", 
                       high = "#FF7200FF", 
                       limit = c(-2, 2),
                       breaks = c(-2, 0, 2)
  ) +
  labs(title = "alpha diversity") + 
  theme(
    plot.title = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 8, hjust = 1, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25),
    panel.background = element_blank(),
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "plain", size = 8, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title = element_text(family = "Arial", face = "plain", size = 8, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.frame = element_blank(),
    legend.key.size = unit(10, "pt"), 
    legend.key.width = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

alpha_diversity_heatmap

ggsave("taxonomy/5_alpha_diversity_heatmap.svg", plot = alpha_diversity_heatmap, width = 7, height = 7, units = "cm", dpi = 300)
# ggsave("taxonomy/5_alpha_diversity_heatmap.tiff", plot = alpha_diversity_heatmap, width = 7, height = 7, units = "cm", dpi = 300)

### box plot
alpha_diversity_plot <- 
  ggplot(alpha_diversity,
         aes(x = group, y = Shannon)) + 
  # box plot
  geom_boxplot(aes(fill = group),
               color = "black",
               outlier.shape = NA, 
               alpha = 0.3,
               width = 0.8,
               staplewidth = 0,
               linewidth = 0.25,
               show.legend = F
  ) +
  # scatter plot
  geom_point(aes(fill = group, color  = group),
             stat = "identity",
             position = position_jitter(height = 0.15, width = 0.15) ,
             shape = 21,
             size = 2, 
             stroke = 0.1,
             alpha = 0.6,
             show.legend = F
  ) +
  # ylim(c(3.5, 7.5)) +
  labs(title = "Shannon index") +
  # group color
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  # wilcox test
  geom_signif(
    comparisons = list(c("Wild", "Lab"),
                       c("Wild", "Fam"), 
                       c("Wild", "Res"), 
                       c("Wild", "Hos")),
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
    step_increase = 0.1
  ) +
  # theme parameter
  theme(
    plot.title = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    axis.text.x = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size=8, hjust=1, vjust=0.5, angle=0, lineheight=2, color="black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill=NA, color="black", linewidth=0.25, linetype="solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
  )

alpha_diversity_plot

ggsave("taxonomy/6_alpha_diversity_plot.svg", plot = alpha_diversity_plot, width = 5.5, height = 6.5, units = "cm", dpi = 300)
# ggsave("taxonomy/6_alpha_diversity_plot.tiff", plot = alpha_diversity_plot, width = 5.5, height = 6.5, units = "cm", dpi = 300)



################################
#### Beta diversity analysis ###
################################

### (PCoA)Principal Coordinates Analysis
taxonomy_bray <- vegdist(relative_taxonomy_abundance, method = 'bray', binary = FALSE, diag = TRUE)

taxonomy_bray_matrix <- as.matrix(taxonomy_bray)

write.table(taxonomy_bray_matrix, "taxonomy/7_taxonomy_sample_bray_crutis.txt", sep = "\t", quote = FALSE) 

taxonomy_pcoa <- cmdscale(
  taxonomy_bray,
  k = 3,
  eig = TRUE,
  add = FALSE,
  x.ret = TRUE)

# Extract the coordinates of the first three dimensions of the PCoA dimensionality reduction results
taxonomy_pcoa_coords <- data.frame(taxonomy_pcoa$points)
taxonomy_pcoa_coords$Sample_ID <- rownames(taxonomy_pcoa_coords)
names(taxonomy_pcoa_coords)[1:3] <- paste0("PCoA", 1:3)
taxonomy_pcoa_result <- merge(taxonomy_pcoa_coords, group, by.x = "row.names", by.y = "sample")

# The "pcoa$eig" record shows the eigenvalues of the main sorting axes in the PCoA sorting results (dividing each eigenvalue by the total sum of eigenvalues gives the explanatory power of each axis)
taxonomy_pcoa_eig = sum(pmax(taxonomy_pcoa$eig, 0))
taxonomy_pcoa_eig_percent <- round(taxonomy_pcoa$eig/taxonomy_pcoa_eig*100, 3)

# Conduct a permutation multivariate (factorial) variance analysis (PERMANOVA/adonis2)
taxonomy_pcoa_permanova_result <- vegan::adonis2(taxonomy_bray ~ group, data = group, permutations = 999, method = "bray")
taxonomy_pcoa_dune_adonis <- paste0("R2 = ", round(taxonomy_pcoa_permanova_result$R2, 3), "\nP = ", round(taxonomy_pcoa_permanova_result$`Pr(>F)`, 3))
taxonomy_pcoa_dune_adonis

### PCoA plot

# Calculate the center point of each group
taxonomy_pcoa_center_points <- 
  taxonomy_pcoa_result %>%
  group_by(group) %>%
  summarise(
    mean_wt = mean(PCoA1),
    mean_mpg = mean(PCoA2)
  )

taxonomy_pcoa_result <- taxonomy_pcoa_result %>% left_join(taxonomy_pcoa_center_points, by = "group")

taxonomy_pcoa_result$group <- factor(taxonomy_pcoa_result$group, level = c("Wild", "Lab", "Fam", "Res", "Hos"))

# Basic scatter plot
taxonomy_pcoa_p0 <- 
  ggplot(taxonomy_pcoa_result, 
         aes(x = PCoA1, y = PCoA2, color = group, fill =group
         )
  ) +
  geom_point(aes(color = group, shape = group, fill = group), 
             size = 1,
             show.legend = T
  ) +
  labs(x = paste("PCoA 1 (", round(taxonomy_pcoa_eig_percent[1], 2), "%)", sep = ""), 
       y = paste("PCoA 2 (", round(taxonomy_pcoa_eig_percent[2], 2), "%)", sep = ""), 
       tag = taxonomy_pcoa_dune_adonis
  ) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_shape_manual(values = c("Wild" = 19, "Lab" = 19, "Fam" = 19, "Res" = 19, "Hos" = 19)) +
  scale_x_continuous(position = "bottom") +
  scale_y_continuous(position = "left") +
  # labs(title = "PCoA: bray curtis distance")+
  theme(
    plot.title = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
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
    legend.key.size = unit(5, "pt"),
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

taxonomy_pcoa_p0

# Add confidence ellipse + center point connection
taxonomy_pcoa_p1 <- 
  taxonomy_pcoa_p0 + 
  stat_ellipse(
    data = taxonomy_pcoa_result,
    aes(fill = group, color = group),
    geom = "polygon", 
    type = "t",
    level = 0.85, 
    linetype = 2, 
    linewidth = 0.5, 
    alpha = 0.3, 
    show.legend = F
  ) +
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  geom_segment(
    data = taxonomy_pcoa_result,
    aes(x = PCoA1, y = PCoA2, xend = mean_wt, yend = mean_mpg, color = group), 
    linetype = "dashed", 
    linewidth = 0.5, 
    alpha = 0.5,
    show.legend = F)

taxonomy_pcoa_p1

# Use ggside to add marginal box plots
taxonomy_pcoa_final <- 
  taxonomy_pcoa_p1 +
  geom_xsideboxplot(aes(y = group, color = group, fill = group),
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
  geom_ysideboxplot(aes(x = group, color = group, fill = group), 
                    orientation = "x", 
                    alpha = 0.3, 
                    outliers = FALSE, 
                    staplewidth = 0.8,
                    linewidth = 0.25,
                    show.legend = F
  ) +
  geom_ysidepoint(aes(x = group, color = group, fill = group), 
                  position =  position_jitter(height = 0.15, width = 0.15),
                  size = 1, 
                  alpha = 0.3,
                  show.legend = F
  ) +
  scale_xsidey_discrete() +
  scale_ysidex_discrete() +
  theme(ggside.panel.scale.x = 0.32,
        ggside.panel.scale.y = 0.32,
        legend.position = c(0.86, 0.88),
  )                            

print(taxonomy_pcoa_final)

ggsave("taxonomy/8_metagenome_taxonomy_PCoA_PCoA1~PCoA2.svg", plot = taxonomy_pcoa_final, width = 6, height = 6, units = "cm")
# ggsave("taxonomy/8_metagenome_taxonomy_PCoA_PCoA1~PCoA2.tiff", plot = taxonomy_pcoa_final, width = 6, height = 6, units = "cm", dpi = 1200)



################################
#### Hierarchical clustering ###
################################

# Hierarchical clustering based on Bray-Curtis distance
taxonomy_hc <- flashClust(
  taxonomy_bray,
  method = "average",
  members = NULL)

library(dendextend)
taxonomy_hc <- as.dendrogram(taxonomy_hc) %>% rev() %>% as.hclust()

# Convert the hclust object to a phylo object
taxonomy_hc_tree <- ape::as.phylo(taxonomy_hc)
plot(taxonomy_hc_tree)

# Export Newick format tree file
write.tree(taxonomy_hc_tree, file = "taxonomy/9_metagenome_taxonomy_bray_curtis_hclust_tree.newick")

# Convert the hclust object to a dendrogram object
taxonomy_dend <- as.dendrogram(taxonomy_hc)

# Plot dendrogram
dev.off()
pdf("taxonomy/10_metagenome_taxonomy_Bray_Curtis_Hierarchical_Dendrogram.pdf")
plot(taxonomy_dend, main = "Bray Curtis Hierarchical Dendrogram")
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

min <- min(1- taxonomy_bray_matrix)
max <- max(1- taxonomy_bray_matrix)

pdf("taxonomy/11_metagenome_taxonomy_Bray_Curtis_Similarity_heatmap.pdf",
    width  = 20,
    height = 20)

pheatmap(
  as.matrix(1-taxonomy_bray_matrix),
  color = colorRampPalette(c("lightyellow" , "lightblue", "black"))(100),
  ### Global parameter
  main = "Megenome taxonomy bray curtis",
  scale = "none",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 15,
  fontsize_row = 10,
  fontsize_col = 10,
  display_numbers	= FALSE,
  cellwidth = 15,
  cellheight = 15,
  legend = TRUE,
  legend_breaks	= c(min, max),
  legend_labels =c(min = round(min, 1), max = round(max, 1)),
  width = 15,
  height = 15,
  ### Cluster tree and gap parameter
  cluster_rows = taxonomy_hc,
  cluster_cols = taxonomy_hc,
  cutree_rows = 1,
  cutree_cols = 1,
  treeheight_row = 40,
  treeheight_col = 40,
  ### Annotation parameter
  annotation_colors = annot_colors,
  annotation_row = row_annotation[, c(2), drop = FALSE],
  annotation_col = col_annotation[, c(2), drop = FALSE],
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = TRUE,
  drop_levels = TRUE
)

dev.off()



### Hierarchical clustering based on the Bray-Curtis distance of bacterial abundance at different taxonomic levels
taxonomy_data <- read.table("taxonomy/1_Total_coverm_contigs_Count_species.tsv", header = TRUE, sep = "\t", row.names = 1)

sum <- colSums(taxonomy_data)

taxonomy_data <- taxonomy_data/sum

rownames <- rownames(taxonomy_data)
split_rownames <- strsplit(rownames, ";", fixed = TRUE)
split_df <- do.call(rbind, lapply(split_rownames, function(x) {
  c(x, rep(NA, 7 - length(x)))
}))
split_df <- data.frame(split_df)
colnames(split_df) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
taxonomy_data <- cbind(taxonomy_data, split_df)


### p_Bacteroidota

# f_Dysgonomonadaceae
Dysgonomonadaceae_abundance <- t(taxonomy_data[taxonomy_data$family %in% c("f__Dysgonomonadaceae"), ] %>% dplyr::select(-c(51:57)))
Dysgonomonadaceae_taxonomy_bray <- vegdist(Dysgonomonadaceae_abundance, method = 'bray', binary = FALSE, diag = TRUE)
Dysgonomonadaceae_taxonomy_bray_matrix <- as.matrix(Dysgonomonadaceae_taxonomy_bray)

# Computational hierarchical clustering
Dysgonomonadaceae_taxonomy_hc <- flashClust(
  Dysgonomonadaceae_taxonomy_bray,
  method = "average",
  members = NULL) %>%
  as.dendrogram() %>% 
  rev() %>% 
  as.hclust()

min <- min(1-Dysgonomonadaceae_taxonomy_bray_matrix)
max <- max(1-Dysgonomonadaceae_taxonomy_bray_matrix)

pdf("taxonomy/12_f_Dysgonomonadaceae_Bray_Curtis_similarity_heatmap.pdf",
    width  = 20,
    height = 20)

pheatmap(
  as.matrix(1-Dysgonomonadaceae_taxonomy_bray_matrix),
  color = colorRampPalette(c("lightyellow" , "lightblue", "black" ))(100),
  ### Global parameter
  main = "Dysgonomonadaceae",
  scale = "none",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 15,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 15,
  silent = FALSE,
  legend = TRUE,
  legend_breaks	= c(min, max),
  legend_labels =c(min = round(min, 1), max = round(max, 1)),
  width = 15,
  height = 15,
  ### Cluster tree and gap parameter
  cluster_rows = Dysgonomonadaceae_taxonomy_hc,
  cluster_cols = Dysgonomonadaceae_taxonomy_hc,
  cutree_rows = 1,
  cutree_cols = 1,
  treeheight_row = 40,
  treeheight_col = 40,
  ### Annotation parameter
  annotation_colors = annot_colors,
  annotation_row = row_annotation[, c(2), drop = FALSE],
  annotation_col = col_annotation[, c(2), drop = FALSE],
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  drop_levels = TRUE
)

dev.off()



# f_UBA932
UBA932_abundance <- t(taxonomy_data[taxonomy_data$family %in% c("f__UBA932"), ] %>% dplyr::select(-c(51:57)))
UBA932_taxonomy_bray <- vegdist(UBA932_abundance, method = 'bray', binary = FALSE, diag = TRUE)
UBA932_taxonomy_bray_matrix <- as.matrix(UBA932_taxonomy_bray)

# Computational hierarchical clustering
UBA932_taxonomy_hc <- flashClust(
  UBA932_taxonomy_bray,
  method = "average",
  members = NULL)

min <- min(1-UBA932_taxonomy_bray_matrix)
max <- max(1-UBA932_taxonomy_bray_matrix)

pdf("taxonomy/13_f_UBA932_Bray_Curtis_Similarity_heatmap.pdf",
    width  = 20,
    height = 20)

pheatmap(
  as.matrix(1-UBA932_taxonomy_bray_matrix),
  color = colorRampPalette(c("lightyellow" , "lightblue", "black" ))(100),
  ### global parameter
  main = "UBA932",
  scale = "none",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 15,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 15,
  silent = FALSE,
  legend = TRUE,
  legend_breaks	= c(min, max),
  legend_labels =c(min = round(min, 1), max = round(max, 1)),
  width = 15,
  height = 15,
  ### Cluster tree and gap parameter
  cluster_rows = UBA932_taxonomy_hc,
  cluster_cols = UBA932_taxonomy_hc,
  cutree_rows = 1,
  cutree_cols = 1,
  treeheight_row = 40,
  treeheight_col = 40,
  ### Annotation parameter
  annotation_colors = annot_colors,
  annotation_row = row_annotation[, c(2), drop = FALSE],
  annotation_col = col_annotation[, c(2), drop = FALSE],
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  drop_levels = TRUE
)

dev.off()



# f_Bacteroidaceae
Bacteroidaceae_abundance <- t(taxonomy_data[taxonomy_data$family %in% c("f__Bacteroidaceae"), ] %>% dplyr::select(-c(51:57)))
Bacteroidaceae_taxonomy_bray <- vegdist(Bacteroidaceae_abundance, method = 'bray', binary = FALSE, diag = TRUE)
Bacteroidaceae_taxonomy_bray_matrix <- as.matrix(Bacteroidaceae_taxonomy_bray)

# Computational hierarchical clustering
Bacteroidaceae_taxonomy_hc <- flashClust(
  Bacteroidaceae_taxonomy_bray,
  method = "average",
  members = NULL)

min <- min(1-Bacteroidaceae_taxonomy_bray_matrix)
max <- max(1-Bacteroidaceae_taxonomy_bray_matrix)

pdf("taxonomy/14_f_Bacteroidaceae_Bray_Curtis_Similarity_heatmap.pdf",
    width  = 20,
    height = 20)

pheatmap(
  as.matrix(1-Bacteroidaceae_taxonomy_bray_matrix),
  color = colorRampPalette(c("lightyellow" , "lightblue", "black" ))(100),
  ### global parameter
  main = "Bacteroidaceae",
  scale = "none",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 15,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 15,
  silent = FALSE,
  legend = TRUE,
  legend_breaks	= c(min, max),
  legend_labels =c(min = round(min, 1), max = round(max, 1)),
  width = 15,
  height = 15,
  ### Cluster tree and gap parameter
  cluster_rows = Bacteroidaceae_taxonomy_hc,
  cluster_cols = Bacteroidaceae_taxonomy_hc,
  cutree_rows = 1,
  cutree_cols = 1,
  treeheight_row = 40,
  treeheight_col = 40,
  ### Annotation parameter
  annotation_colors = annot_colors,
  annotation_row = row_annotation[, c(2), drop = FALSE],
  annotation_col = col_annotation[, c(2), drop = FALSE],
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  drop_levels = TRUE
)

dev.off()



# f_Paludibacteraceae
Paludibacteraceae_abundance <- t(taxonomy_data[taxonomy_data$family %in% c("f__Paludibacteraceae"), ] %>% dplyr::select(-c(51:57)))
Paludibacteraceae_taxonomy_bray <- vegdist(Paludibacteraceae_abundance, method = 'bray', binary = FALSE, diag = TRUE)
Paludibacteraceae_taxonomy_bray_matrix <- as.matrix(Paludibacteraceae_taxonomy_bray)

# Computational hierarchical clustering
Paludibacteraceae_taxonomy_hc <- flashClust(
  Paludibacteraceae_taxonomy_bray,
  method = "average",
  members = NULL)

min <- min(1-Paludibacteraceae_taxonomy_bray_matrix)
max <- max(1-Paludibacteraceae_taxonomy_bray_matrix)

pdf("taxonomy/15_f_Paludibacteraceae_Bray_Curtis_Similarity_heatmap.pdf",
    width  = 20,
    height = 20)

pheatmap(
  as.matrix(1-Paludibacteraceae_taxonomy_bray_matrix),
  color = colorRampPalette(c("lightyellow" , "lightblue", "black" ))(100),
  ### Computational hierarchical clustering
  main = "Paludibacteraceae",
  scale = "none",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 15,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 15,
  silent = FALSE,
  legend = TRUE,
  legend_breaks	= c(min, max),
  legend_labels =c(min = round(min, 1), max = round(max, 1)),
  width = 15,
  height = 15,
  ### Cluster tree and gap parameter
  cluster_rows = Paludibacteraceae_taxonomy_hc,
  cluster_cols = Paludibacteraceae_taxonomy_hc,
  cutree_rows = 1,
  cutree_cols = 1,
  treeheight_row = 40,
  treeheight_col = 40,
  ### Annotation parameter
  annotation_colors = annot_colors,
  annotation_row = row_annotation[, c(2), drop = FALSE],
  annotation_col = col_annotation[, c(2), drop = FALSE],
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  drop_levels = TRUE
)

dev.off()



# f_Rikenellaceae
Rikenellaceae_abundance <- t(taxonomy_data[taxonomy_data$family %in% c("f__Rikenellaceae"), ] %>% dplyr::select(-c(51:57)))
Rikenellaceae_taxonomy_bray <- vegdist(Rikenellaceae_abundance, method = 'bray', binary = FALSE, diag = TRUE)
Rikenellaceae_taxonomy_bray_matrix <- as.matrix(Rikenellaceae_taxonomy_bray)

# Computational hierarchical clustering
Rikenellaceae_taxonomy_hc <- flashClust(
  Rikenellaceae_taxonomy_bray,
  method = "average",
  members = NULL)

min <- min(1-Rikenellaceae_taxonomy_bray_matrix)
max <- max(1-Rikenellaceae_taxonomy_bray_matrix)

pdf("taxonomy/16_f_Rikenellaceae_Bray_Curtis_Similarity_heatmap.pdf",
    width  = 20,
    height = 20)

pheatmap(
  1-Rikenellaceae_taxonomy_bray_matrix,
  color = colorRampPalette(c("lightyellow" , "lightblue", "black" ))(100),
  ### global parameter
  main = "Rikenellaceae",
  scale = "none",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 15,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 15,
  silent = FALSE,
  legend = TRUE,
  legend_breaks	= c(min, max),
  legend_labels =c(min = round(min, 1), max = round(max, 1)),
  width = 15,
  height = 15,
  ### Cluster tree and gap parameter
  cluster_rows = Rikenellaceae_taxonomy_hc,
  cluster_cols = Rikenellaceae_taxonomy_hc,
  cutree_rows = 1,
  cutree_cols = 1,
  treeheight_row = 15,
  treeheight_col = 15,
  ### Annotation parameter
  annotation_colors = annot_colors,  # Annotation color mapping list file
  annotation_row = row_annotation[, c(2), drop = FALSE],
  annotation_col = col_annotation[, c(2), drop = FALSE],
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  drop_levels = TRUE
)

dev.off()



# f_Tannerellaceae
Tannerellaceae_abundance <- t(taxonomy_data[taxonomy_data$family %in% c("f__Tannerellaceae"), ] %>% dplyr::select(-c(51:57)))
Tannerellaceae_taxonomy_bray <- vegdist(Tannerellaceae_abundance, method = 'bray', binary = FALSE, diag = TRUE)
Tannerellaceae_taxonomy_bray_matrix <- as.matrix(Tannerellaceae_taxonomy_bray)

# Computational hierarchical clustering
Tannerellaceae_taxonomy_hc <- flashClust(
  Tannerellaceae_taxonomy_bray,
  method = "average",
  members = NULL)

min <- min(1-Tannerellaceae_taxonomy_bray_matrix)
max <- max(1-Tannerellaceae_taxonomy_bray_matrix)

pdf("taxonomy/17_f_Tannerellaceae_Bray_Curtis_Similarity_heatmap.pdf",
    width  = 20,
    height = 20)

pheatmap(
  as.matrix(1-Tannerellaceae_taxonomy_bray_matrix),
  color = colorRampPalette(c("lightyellow" , "lightblue", "black" ))(100),
  ### global parameter
  main = "Tannerellaceae",
  scale = "none",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 15,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 15,
  silent = FALSE,
  legend = TRUE,
  legend_breaks	= c(min, max), 
  legend_labels =c(min = round(min, 1), max = round(max, 1)),
  width = 15,
  height = 15,
  ### Cluster tree and gap parameter
  cluster_rows = Tannerellaceae_taxonomy_hc,
  cluster_cols = Tannerellaceae_taxonomy_hc,
  cutree_rows = 1,
  cutree_cols = 1,
  treeheight_row = 40,
  treeheight_col = 40,
  ### Annotation parameter
  annotation_colors = annot_colors,
  annotation_row = row_annotation[, c(2), drop = FALSE],
  annotation_col = col_annotation[, c(2), drop = FALSE],
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  drop_levels = TRUE
)

dev.off()



# p_Bacteroidota
Bacteroidota_abundance <- t(taxonomy_data[taxonomy_data$phylum %in% c("p__Bacteroidota"), ] %>% dplyr::select(-c(51:57)))
Bacteroidota_taxonomy_bray <- vegdist(Bacteroidota_abundance, method = 'bray', binary = FALSE, diag = TRUE)
Bacteroidota_taxonomy_bray_matrix <- as.matrix(Bacteroidota_taxonomy_bray)

# Computational hierarchical clustering
Bacteroidota_taxonomy_hc <- flashClust(
  Bacteroidota_taxonomy_bray,
  method = "average",
  members = NULL)

min <- min(1-Bacteroidota_taxonomy_bray_matrix)
max <- max(1-Bacteroidota_taxonomy_bray_matrix)

pdf("taxonomy/18_p_Bacteroidota_BrayCurtis_Similarity_heatmap.pdf",
    width  = 20,
    height = 20)

pheatmap(
  as.matrix(1-Bacteroidota_taxonomy_bray_matrix),
  color = colorRampPalette(c("lightyellow" , "lightblue", "black" ))(100),
  ### global parameter
  main = "Bacteroidota",
  scale = "none",
  border_color = NA, 
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 15,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 15,
  silent = FALSE,
  legend = TRUE,
  legend_breaks	= c(min, max),
  legend_labels =c(min = round(min, 1), max = round(max, 1)),
  width = 15,
  height = 15,
  ### Cluster tree and gap parameter
  cluster_rows = Bacteroidota_taxonomy_hc,
  cluster_cols = Bacteroidota_taxonomy_hc,
  cutree_rows = 1,
  cutree_cols = 1,
  treeheight_row = 40,
  treeheight_col = 40,
  ### Annotation parameter
  annotation_colors = annot_colors,  # Annotation color mapping list file
  annotation_row = row_annotation[, c(2), drop = FALSE],
  annotation_col = col_annotation[, c(2), drop = FALSE],
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  drop_levels = TRUE
)

dev.off()



# p_Bacillota
Bacillota_abundance <- t(taxonomy_data[taxonomy_data$phylum %in% c("p__Bacillota_A", "p__Bacillota_B", "p__Bacillota_C", "p__Bacillota_I", "p__Bacillota"), ] %>% dplyr::select(-c(51:57)))
Bacillota_taxonomy_bray <- vegdist(Bacillota_abundance, method = 'bray', binary = FALSE, diag = TRUE)
Bacillota_taxonomy_bray_matrix <- as.matrix(Bacillota_taxonomy_bray)

# Computational hierarchical clustering
Bacillota_taxonomy_hc <- flashClust(
  Bacillota_taxonomy_bray,
  method = "average",
  members = NULL)

min <- min(1-Bacillota_taxonomy_bray_matrix)
max <- max(1-Bacillota_taxonomy_bray_matrix)

pdf("taxonomy/19_p_Bacillota_Bray_Curtis_Similarity_heatmap.pdf",
    width  = 20,
    height = 20)

pheatmap(
  as.matrix(1-Bacillota_taxonomy_bray_matrix),
  color = colorRampPalette(c("lightyellow" , "lightblue", "black" ))(100),
  ### global parameter
  main = "Bacillota",
  scale = "none",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 15,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 15,
  silent = FALSE,
  legend = TRUE,
  legend_breaks	= c(min, max),
  legend_labels =c(min = round(min, 1), max = round(max, 1)),
  width = 15,
  height = 15,
  ### Cluster tree and gap parameter
  cluster_rows = Bacillota_taxonomy_hc,
  cluster_cols = Bacillota_taxonomy_hc,
  cutree_rows = 1,
  cutree_cols = 1,
  treeheight_row = 40,
  treeheight_col = 40,
  ### Annotation parameter
  annotation_colors = annot_colors,
  annotation_row = row_annotation[, c(2), drop = FALSE],
  annotation_col = col_annotation[, c(2), drop = FALSE],
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  drop_levels = TRUE
)

dev.off()



# p_Pseudomonadota
Pseudomonadota_abundance <- t(taxonomy_data[taxonomy_data$phylum %in% c("p__Pseudomonadota"), ] %>% dplyr::select(-c(51:57)))
Pseudomonadota_taxonomy_bray <- vegdist(Pseudomonadota_abundance, method = 'bray', binary = FALSE, diag = TRUE)
Pseudomonadota_taxonomy_bray_matrix <- as.matrix(Pseudomonadota_taxonomy_bray)

# Computational hierarchical clustering
Pseudomonadota_taxonomy_hc <- flashClust(
  Pseudomonadota_taxonomy_bray,
  method = "average",
  members = NULL)

min <- min(1-Pseudomonadota_taxonomy_bray_matrix)
max <- max(1-Pseudomonadota_taxonomy_bray_matrix)

pdf("taxonomy/20_p_Pseudomonadota_Bray_Curtis_Similarity_heatmap.pdf",
    width  = 20,
    height = 20)

pheatmap(
  as.matrix(1-Pseudomonadota_taxonomy_bray_matrix),
  color = colorRampPalette(c("lightyellow" , "lightblue", "black" ))(100),
  ### global parameter
  main = "Pseudomonadota",
  scale = "none",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 15,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 15,
  silent = FALSE,
  legend = TRUE,
  legend_breaks	= c(min, max),
  legend_labels =c(min = round(min, 1), max = round(max, 1)),
  width = 15,
  height = 15,
  ### Cluster tree and gap parameter
  cluster_rows = Pseudomonadota_taxonomy_hc,
  cluster_cols = Pseudomonadota_taxonomy_hc,
  cutree_rows = 1,
  cutree_cols = 1,
  treeheight_row = 40,
  treeheight_col = 40,
  ### Annotation parameter
  annotation_colors = annot_colors,  # Annotation color mapping list file
  annotation_row = row_annotation[, c(2), drop = FALSE],
  annotation_col = col_annotation[, c(2), drop = FALSE],
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  drop_levels = TRUE
)

dev.off()



# p_Desulfobacterota
Desulfobacterota_abundance <- t(taxonomy_data[taxonomy_data$phylum %in% c("p__Desulfobacterota", "p__Desulfobacterota_B", "p__Desulfobacterota_C", "p__Desulfobacterota_D", "p__Desulfobacterota_E", "p__Desulfobacterota_G", "p__Desulfobacterota_I"), ] %>% dplyr::select(-c(51:57)))
Desulfobacterota_taxonomy_bray <- vegdist(Desulfobacterota_abundance, method = 'bray', binary = FALSE, diag = TRUE)
Desulfobacterota_taxonomy_bray_matrix <- as.matrix(Desulfobacterota_taxonomy_bray)

# Computational hierarchical clustering
Desulfobacterota_taxonomy_hc <- flashClust(
  Desulfobacterota_taxonomy_bray,
  method = "average",
  members = NULL)%>%
  as.dendrogram() %>%
  rev() %>%
  as.hclust()

min <- min(1-Desulfobacterota_taxonomy_bray_matrix)
max <- max(1-Desulfobacterota_taxonomy_bray_matrix)

pdf("taxonomy/21_p_Desulfobacterota_Bray_Curtis_Similarity_heatmap.pdf",
    width  = 20,
    height = 20)

pheatmap(
  as.matrix(1-Desulfobacterota_taxonomy_bray_matrix),
  color = colorRampPalette(c("lightyellow" , "lightblue", "black" ))(100),
  ### global parameter
  main = "Desulfobacterota",
  scale = "none",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 15,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 15,
  silent = FALSE,
  legend = TRUE,
  legend_breaks	= c(min, max), 
  legend_labels =c(min = round(min, 1), max = round(max, 1)),
  width = 15,
  height = 15,
  ### Cluster tree and gap parameter
  cluster_rows = Desulfobacterota_taxonomy_hc,
  cluster_cols = Desulfobacterota_taxonomy_hc,
  cutree_rows = 1,
  cutree_cols = 1,
  treeheight_row = 40,
  treeheight_col = 40,
  ### Annotation parameter
  annotation_colors = annot_colors,
  annotation_row = row_annotation[, c(2), drop = FALSE],
  annotation_col = col_annotation[, c(2), drop = FALSE],
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  drop_levels = TRUE
)

dev.off()



# p_Planctomycetota
Planctomycetota_abundance <- t(taxonomy_data[taxonomy_data$phylum %in% c("p__Planctomycetota"), ] %>% dplyr::select(-c(51:57)))
Planctomycetota_taxonomy_bray <- vegdist(Planctomycetota_abundance, method = 'bray', binary = FALSE, diag = TRUE)
Planctomycetota_taxonomy_bray_matrix <- as.matrix(Planctomycetota_taxonomy_bray)

# Computational hierarchical clustering
Planctomycetota_taxonomy_hc <- flashClust(
  Planctomycetota_taxonomy_bray,
  method = "average",
  members = NULL)

min <- min(1-Planctomycetota_taxonomy_bray_matrix)
max <- max(1-Planctomycetota_taxonomy_bray_matrix)

pdf("taxonomy/22_p_Planctomycetota_Bray_Curtis_Similarity_heatmap.pdf",
    width  = 20,
    height = 20)

pheatmap(
  as.matrix(1-Planctomycetota_taxonomy_bray_matrix),
  color = colorRampPalette(c("lightyellow" , "lightblue", "black" ))(100),
  ### global parameter
  main = "Planctomycetota",
  scale = "none",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 15,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 15,
  silent = FALSE,
  legend = TRUE,
  legend_breaks	= c(min, max),
  legend_labels =c(min = round(min, 1), max = round(max, 1)),
  width = 15,
  height = 15,
  ### Cluster tree and gap parameter
  cluster_rows = Planctomycetota_taxonomy_hc,
  cluster_cols = Planctomycetota_taxonomy_hc,
  cutree_rows = 1,
  cutree_cols = 1,
  treeheight_row = 40,
  treeheight_col = 40,
  ### Annotation parameter
  annotation_colors = annot_colors,  # Annotation color mapping list file
  annotation_row = row_annotation[, c(2), drop = FALSE],
  annotation_col = col_annotation[, c(2), drop = FALSE],
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  drop_levels = TRUE
)

dev.off()



# p_Actinomycetota
Actinomycetota_abundance <- t(taxonomy_data[taxonomy_data$phylum %in% c("p__Actinomycetota"), ] %>% dplyr::select(-c(51:57)))
Actinomycetota_taxonomy_bray <- vegdist(Actinomycetota_abundance, method = 'bray', binary = FALSE, diag = TRUE)
Actinomycetota_taxonomy_bray_matrix <- as.matrix(Actinomycetota_taxonomy_bray)

# Computational hierarchical clustering
Actinomycetota_taxonomy_hc <- flashClust(
  Actinomycetota_taxonomy_bray,
  method = "average",
  members = NULL)

min <- min(1-Actinomycetota_taxonomy_bray_matrix)
max <- max(1-Actinomycetota_taxonomy_bray_matrix)

pdf("taxonomy/23_p_Actinomycetota_Bray_Curtis_Similarity_heatmap.pdf",
    width  = 20,
    height = 20)

pheatmap(
  as.matrix(1-Actinomycetota_taxonomy_bray_matrix),
  color = colorRampPalette(c("lightyellow" , "lightblue", "black" ))(100),
  ### global parameter
  main = "Actinomycetota bray curtis",
  scale = "none",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 15,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 15,
  silent = FALSE,
  legend = TRUE,
  legend_breaks	= c(min, max),
  legend_labels =c(min = round(min, 1), max = round(max, 1)),
  width = 15,
  height = 15,
  ### Cluster tree and gap parameter
  cluster_rows = Actinomycetota_taxonomy_hc,
  cluster_cols = Actinomycetota_taxonomy_hc,
  cutree_rows = 1,
  cutree_cols = 1,
  treeheight_row = 40,
  treeheight_col = 40,
  ### Annotation parameter
  annotation_colors = annot_colors,
  annotation_row = row_annotation[, c(2), drop = FALSE],
  annotation_col = col_annotation[, c(2), drop = FALSE],
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  drop_levels = TRUE
)

dev.off()



# p_Verrucomicrobiota
Verrucomicrobiota_abundance <- t(taxonomy_data[taxonomy_data$phylum %in% c("p__Verrucomicrobiota"), ] %>% dplyr::select(-c(51:57)))
Verrucomicrobiota_taxonomy_bray <- vegdist(Verrucomicrobiota_abundance, method = 'bray', binary = FALSE, diag = TRUE)
Verrucomicrobiota_taxonomy_bray_matrix <- as.matrix(Verrucomicrobiota_taxonomy_bray)

# Computational hierarchical clustering
Verrucomicrobiota_taxonomy_hc <- flashClust(
  Verrucomicrobiota_taxonomy_bray,
  method = "average",
  members = NULL) 

min <- min(1-Verrucomicrobiota_taxonomy_bray_matrix)
max <- max(1-Verrucomicrobiota_taxonomy_bray_matrix)

pdf("taxonomy/24_p_Verrucomicrobiota_Bray_Curtis_Similarity_heatmap.pdf",
    width  = 20,
    height = 20)

pheatmap(
  as.matrix(1-Verrucomicrobiota_taxonomy_bray_matrix),
  color = colorRampPalette(c("lightyellow" , "lightblue", "black" ))(100),
  ### global parameter
  main = "Verrucomicrobiota",
  scale = "none",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 15,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 15,
  silent = FALSE,
  legend = TRUE,
  legend_breaks	= c(min, max),
  legend_labels =c(min = round(min, 1), max = round(max, 1)),
  width = 15,
  height = 15,
  ### Cluster tree and gap parameter
  cluster_rows = Verrucomicrobiota_taxonomy_hc,
  cluster_cols = Verrucomicrobiota_taxonomy_hc,
  cutree_rows = 1,
  cutree_cols = 1,
  treeheight_row = 40,
  treeheight_col = 40,
  ### Annotation parameter
  annotation_colors = annot_colors,  # Annotation color mapping list file
  annotation_row = row_annotation[, c(2), drop = FALSE],
  annotation_col = col_annotation[, c(2), drop = FALSE],
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  drop_levels = TRUE
)

dev.off()


### PCoA1 ridges plot

# f_Dysgonomonadaceae
# PCoA analysis based on bray curtis distance
Dysgonomonadaceae_taxonomy_pcoa <- cmdscale(
  Dysgonomonadaceae_taxonomy_bray,
  k = 3,
  eig = TRUE,
  add = FALSE,
  x.ret = TRUE)
# Extract the coordinates of the first three dimensions of the PCoA dimensionality reduction results
Dysgonomonadaceae_taxonomy_pcoa_coords <- data.frame(Dysgonomonadaceae_taxonomy_pcoa$points)
Dysgonomonadaceae_taxonomy_pcoa_coords$Sample_ID <- rownames(Dysgonomonadaceae_taxonomy_pcoa_coords)
names(Dysgonomonadaceae_taxonomy_pcoa_coords)[1:3] <- paste0("PCoA", 1:3)
Dysgonomonadaceae_taxonomy_pcoa_result <- merge(Dysgonomonadaceae_taxonomy_pcoa_coords, group, by.x = "row.names", by.y = "sample")
Dysgonomonadaceae_taxonomy_pcoa_result$facet <- c(rep("Dysgonomonadaceae", 50))
Dysgonomonadaceae_taxonomy_pcoa_result$PCoA1_zscore <- scale(Dysgonomonadaceae_taxonomy_pcoa_result$PCoA1)
Dysgonomonadaceae_taxonomy_pcoa_result$PCoA1_scaled <- 
  -3 + ((Dysgonomonadaceae_taxonomy_pcoa_result$PCoA1 - min(Dysgonomonadaceae_taxonomy_pcoa_result$PCoA1)) / (max(Dysgonomonadaceae_taxonomy_pcoa_result$PCoA1) - min(Dysgonomonadaceae_taxonomy_pcoa_result$PCoA1))) * 6

Dysgonomonadaceae_richness <- data.frame(sample = rownames(Dysgonomonadaceae_abundance),  richness = apply(Dysgonomonadaceae_abundance, 1, function(x) sum(x > 0)))

# f_UBA932
# PCoA analysis based on bray curtis distance
UBA932_taxonomy_pcoa <- cmdscale(
  UBA932_taxonomy_bray,
  k = 3,
  eig = TRUE,
  add = FALSE,
  x.ret = TRUE)
# Extract the coordinates of the first three dimensions of the PCoA dimensionality reduction results
UBA932_taxonomy_pcoa_coords <- data.frame(UBA932_taxonomy_pcoa$points)
UBA932_taxonomy_pcoa_coords$Sample_ID <- rownames(UBA932_taxonomy_pcoa_coords)
names(UBA932_taxonomy_pcoa_coords)[1:3] <- paste0("PCoA", 1:3)
UBA932_taxonomy_pcoa_result <- merge(UBA932_taxonomy_pcoa_coords, group, by.x = "row.names", by.y = "sample")
UBA932_taxonomy_pcoa_result$facet <- c(rep("UBA932", 50))
UBA932_taxonomy_pcoa_result$PCoA1_zscore <- scale(UBA932_taxonomy_pcoa_result$PCoA1)
UBA932_taxonomy_pcoa_result$PCoA1_scaled <- 
  -3 + ((UBA932_taxonomy_pcoa_result$PCoA1 - min(UBA932_taxonomy_pcoa_result$PCoA1)) / (max(UBA932_taxonomy_pcoa_result$PCoA1) - min(UBA932_taxonomy_pcoa_result$PCoA1))) * 6

UBA932_richness <- data.frame(sample = rownames(UBA932_abundance),  richness = apply(UBA932_abundance, 1, function(x) sum(x > 0)))

# f_Bacteroidaceae
# PCoA analysis based on bray curtis distance
Bacteroidaceae_taxonomy_pcoa <- cmdscale(
  Bacteroidaceae_taxonomy_bray,
  k = 3,
  eig = TRUE,
  add = FALSE,
  x.ret = TRUE)
# Extract the coordinates of the first three dimensions of the PCoA dimensionality reduction results
Bacteroidaceae_taxonomy_pcoa_coords <- data.frame(Bacteroidaceae_taxonomy_pcoa$points)
Bacteroidaceae_taxonomy_pcoa_coords$Sample_ID <- rownames(Bacteroidaceae_taxonomy_pcoa_coords)
names(Bacteroidaceae_taxonomy_pcoa_coords)[1:3] <- paste0("PCoA", 1:3)
Bacteroidaceae_taxonomy_pcoa_result <- merge(Bacteroidaceae_taxonomy_pcoa_coords, group, by.x = "row.names", by.y = "sample")
Bacteroidaceae_taxonomy_pcoa_result$facet <- c(rep("Bacteroidaceae", 50))
Bacteroidaceae_taxonomy_pcoa_result$PCoA1_zscore <- scale(Bacteroidaceae_taxonomy_pcoa_result$PCoA1)
Bacteroidaceae_taxonomy_pcoa_result$PCoA1_scaled <- 
  -3 + ((Bacteroidaceae_taxonomy_pcoa_result$PCoA1 - min(Bacteroidaceae_taxonomy_pcoa_result$PCoA1)) / (max(Bacteroidaceae_taxonomy_pcoa_result$PCoA1) - min(Bacteroidaceae_taxonomy_pcoa_result$PCoA1))) * 6

Bacteroidaceae_richness <- data.frame(sample = rownames(Bacteroidaceae_abundance),  richness = apply(Bacteroidaceae_abundance, 1, function(x) sum(x > 0)))

# f_Paludibacteraceae
# PCoA analysis based on bray curtis distance
Paludibacteraceae_taxonomy_pcoa <- cmdscale(
  Paludibacteraceae_taxonomy_bray,
  k = 3,
  eig = TRUE,
  add = FALSE,
  x.ret = TRUE)
# Extract the coordinates of the first three dimensions of the PCoA dimensionality reduction results
Paludibacteraceae_taxonomy_pcoa_coords <- data.frame(Paludibacteraceae_taxonomy_pcoa$points)
Paludibacteraceae_taxonomy_pcoa_coords$Sample_ID <- rownames(Paludibacteraceae_taxonomy_pcoa_coords)
names(Paludibacteraceae_taxonomy_pcoa_coords)[1:3] <- paste0("PCoA", 1:3)
Paludibacteraceae_taxonomy_pcoa_result <- merge(Paludibacteraceae_taxonomy_pcoa_coords, group, by.x = "row.names", by.y = "sample")
Paludibacteraceae_taxonomy_pcoa_result$facet <- c(rep("Paludibacteraceae", 50))
Paludibacteraceae_taxonomy_pcoa_result$PCoA1_zscore <- scale(Paludibacteraceae_taxonomy_pcoa_result$PCoA1)
Paludibacteraceae_taxonomy_pcoa_result$PCoA1_scaled <- 
  -(-3 + ((Paludibacteraceae_taxonomy_pcoa_result$PCoA1 - min(Paludibacteraceae_taxonomy_pcoa_result$PCoA1)) / (max(Paludibacteraceae_taxonomy_pcoa_result$PCoA1) - min(Paludibacteraceae_taxonomy_pcoa_result$PCoA1))) * 6)

Paludibacteraceae_richness <- data.frame(sample = rownames(Paludibacteraceae_abundance),  richness = apply(Paludibacteraceae_abundance, 1, function(x) sum(x > 0)))

# f_Rikenellaceae
# PCoA analysis based on bray curtis distance
Rikenellaceae_taxonomy_pcoa <- cmdscale(
  Rikenellaceae_taxonomy_bray,
  k = 3,
  eig = TRUE,
  add = FALSE,
  x.ret = TRUE)
# Extract the coordinates of the first three dimensions of the PCoA dimensionality reduction results
Rikenellaceae_taxonomy_pcoa_coords <- data.frame(Rikenellaceae_taxonomy_pcoa$points)
Rikenellaceae_taxonomy_pcoa_coords$Sample_ID <- rownames(Rikenellaceae_taxonomy_pcoa_coords)
names(Rikenellaceae_taxonomy_pcoa_coords)[1:3] <- paste0("PCoA", 1:3)
Rikenellaceae_taxonomy_pcoa_result <- merge(Rikenellaceae_taxonomy_pcoa_coords, group, by.x = "row.names", by.y = "sample")
Rikenellaceae_taxonomy_pcoa_result$facet <- c(rep("Rikenellaceae", 50))
Rikenellaceae_taxonomy_pcoa_result$PCoA1_zscore <- scale(Rikenellaceae_taxonomy_pcoa_result$PCoA1)
Rikenellaceae_taxonomy_pcoa_result$PCoA1_scaled <- 
  -3 + ((Rikenellaceae_taxonomy_pcoa_result$PCoA1 - min(Rikenellaceae_taxonomy_pcoa_result$PCoA1)) / (max(Rikenellaceae_taxonomy_pcoa_result$PCoA1) - min(Rikenellaceae_taxonomy_pcoa_result$PCoA1))) * 6

Rikenellaceae_richness <- data.frame(sample = rownames(Rikenellaceae_abundance),  richness = apply(Rikenellaceae_abundance, 1, function(x) sum(x > 0)))

# f_Tannerellaceae
# PCoA analysis based on bray curtis distance
Tannerellaceae_taxonomy_pcoa <- cmdscale(
  Tannerellaceae_taxonomy_bray,
  k = 3,
  eig = TRUE,
  add = FALSE,
  x.ret = TRUE)
# Extract the coordinates of the first three dimensions of the PCoA dimensionality reduction results
Tannerellaceae_taxonomy_pcoa_coords <- data.frame(Tannerellaceae_taxonomy_pcoa$points)
Tannerellaceae_taxonomy_pcoa_coords$Sample_ID <- rownames(Tannerellaceae_taxonomy_pcoa_coords)
names(Tannerellaceae_taxonomy_pcoa_coords)[1:3] <- paste0("PCoA", 1:3)
Tannerellaceae_taxonomy_pcoa_result <- merge(Tannerellaceae_taxonomy_pcoa_coords, group, by.x = "row.names", by.y = "sample")
Tannerellaceae_taxonomy_pcoa_result$facet <- c(rep("Tannerellaceae", 50))
Tannerellaceae_taxonomy_pcoa_result$PCoA1_zscore <- scale(Tannerellaceae_taxonomy_pcoa_result$PCoA1)
Tannerellaceae_taxonomy_pcoa_result$PCoA1_scaled <- 
  -3 + ((Tannerellaceae_taxonomy_pcoa_result$PCoA1 - min(Tannerellaceae_taxonomy_pcoa_result$PCoA1)) / (max(Tannerellaceae_taxonomy_pcoa_result$PCoA1) - min(Tannerellaceae_taxonomy_pcoa_result$PCoA1))) * 6

Tannerellaceae_richness <- data.frame(sample = rownames(Tannerellaceae_abundance),  richness = apply(Tannerellaceae_abundance, 1, function(x) sum(x > 0)))

# p_Bacteroidota
# PCoA analysis based on bray curtis distance
Bacteroidota_taxonomy_pcoa <- cmdscale(
  Bacteroidota_taxonomy_bray,
  k = 3,
  eig = TRUE,
  add = FALSE,
  x.ret = TRUE)
# Extract the coordinates of the first three dimensions of the PCoA dimensionality reduction results
Bacteroidota_taxonomy_pcoa_coords <- data.frame(Bacteroidota_taxonomy_pcoa$points)
Bacteroidota_taxonomy_pcoa_coords$Sample_ID <- rownames(Bacteroidota_taxonomy_pcoa_coords)
names(Bacteroidota_taxonomy_pcoa_coords)[1:3] <- paste0("PCoA", 1:3)
Bacteroidota_taxonomy_pcoa_result <- merge(Bacteroidota_taxonomy_pcoa_coords, group, by.x = "row.names", by.y = "sample")
Bacteroidota_taxonomy_pcoa_result$facet <- c(rep("Bacteroidota", 50))
Bacteroidota_taxonomy_pcoa_result$PCoA1_zscore <- scale(Bacteroidota_taxonomy_pcoa_result$PCoA1)
Bacteroidota_taxonomy_pcoa_result$PCoA1_scaled <- 
  -3 + ((Bacteroidota_taxonomy_pcoa_result$PCoA1 - min(Bacteroidota_taxonomy_pcoa_result$PCoA1)) / (max(Bacteroidota_taxonomy_pcoa_result$PCoA1) - min(Bacteroidota_taxonomy_pcoa_result$PCoA1))) * 6

Bacteroidota_richness <- data.frame(sample = rownames(Bacteroidota_abundance),  richness = apply(Bacteroidota_abundance, 1, function(x) sum(x > 0)))

# p_Bacillota
# PCoA analysis based on bray curtis distance
Bacillota_taxonomy_pcoa <- cmdscale(
  Bacillota_taxonomy_bray,
  k = 3,
  eig = TRUE,
  add = FALSE,
  x.ret = TRUE)
# Extract the coordinates of the first three dimensions of the PCoA dimensionality reduction results
Bacillota_taxonomy_pcoa_coords <- data.frame(Bacillota_taxonomy_pcoa$points)
Bacillota_taxonomy_pcoa_coords$Sample_ID <- rownames(Bacillota_taxonomy_pcoa_coords)
names(Bacillota_taxonomy_pcoa_coords)[1:3] <- paste0("PCoA", 1:3)
Bacillota_taxonomy_pcoa_result <- merge(Bacillota_taxonomy_pcoa_coords, group, by.x = "row.names", by.y = "sample")
Bacillota_taxonomy_pcoa_result$facet <- c(rep("Bacillota", 50))
Bacillota_taxonomy_pcoa_result$PCoA1_zscore <- scale(Bacillota_taxonomy_pcoa_result$PCoA1)
Bacillota_taxonomy_pcoa_result$PCoA1_scaled <- 
  -3 + ((Bacillota_taxonomy_pcoa_result$PCoA1 - min(Bacillota_taxonomy_pcoa_result$PCoA1)) / (max(Bacillota_taxonomy_pcoa_result$PCoA1) - min(Bacillota_taxonomy_pcoa_result$PCoA1))) * 6

Bacillota_richness <- data.frame(sample = rownames(Bacillota_abundance),  richness = apply(Bacillota_abundance, 1, function(x) sum(x > 0)))

# p_Pseudomonadota
# PCoA analysis based on bray curtis distance
Pseudomonadota_taxonomy_pcoa <- cmdscale(
  Pseudomonadota_taxonomy_bray,
  k = 3,
  eig = TRUE,
  add = FALSE,
  x.ret = TRUE)
# Extract the coordinates of the first three dimensions of the PCoA dimensionality reduction results
Pseudomonadota_taxonomy_pcoa_coords <- data.frame(Pseudomonadota_taxonomy_pcoa$points)
Pseudomonadota_taxonomy_pcoa_coords$Sample_ID <- rownames(Pseudomonadota_taxonomy_pcoa_coords)
names(Pseudomonadota_taxonomy_pcoa_coords)[1:3] <- paste0("PCoA", 1:3)
Pseudomonadota_taxonomy_pcoa_result <- merge(Pseudomonadota_taxonomy_pcoa_coords, group, by.x = "row.names", by.y = "sample")
Pseudomonadota_taxonomy_pcoa_result$facet <- c(rep("Pseudomonadota", 50))
Pseudomonadota_taxonomy_pcoa_result$PCoA1_zscore <- scale(Pseudomonadota_taxonomy_pcoa_result$PCoA1)
Pseudomonadota_taxonomy_pcoa_result$PCoA1_scaled <- 
  -3 + ((Pseudomonadota_taxonomy_pcoa_result$PCoA1 - min(Pseudomonadota_taxonomy_pcoa_result$PCoA1)) / (max(Pseudomonadota_taxonomy_pcoa_result$PCoA1) - min(Pseudomonadota_taxonomy_pcoa_result$PCoA1))) * 6

Pseudomonadota_richness <- data.frame(sample = rownames(Pseudomonadota_abundance),  richness = apply(Pseudomonadota_abundance, 1, function(x) sum(x > 0)))

# p_Desulfobacterota
# PCoA analysis based on bray curtis distance
Desulfobacterota_taxonomy_pcoa <- cmdscale(
  Desulfobacterota_taxonomy_bray,
  k = 3,
  eig = TRUE,
  add = FALSE,
  x.ret = TRUE)
# Extract the coordinates of the first three dimensions of the PCoA dimensionality reduction results
Desulfobacterota_taxonomy_pcoa_coords <- data.frame(Desulfobacterota_taxonomy_pcoa$points)
Desulfobacterota_taxonomy_pcoa_coords$Sample_ID <- rownames(Desulfobacterota_taxonomy_pcoa_coords)
names(Desulfobacterota_taxonomy_pcoa_coords)[1:3] <- paste0("PCoA", 1:3)
Desulfobacterota_taxonomy_pcoa_result <- merge(Desulfobacterota_taxonomy_pcoa_coords, group, by.x = "row.names", by.y = "sample")
Desulfobacterota_taxonomy_pcoa_result$facet <- c(rep("Desulfobacterota", 50))
Desulfobacterota_taxonomy_pcoa_result$PCoA1_zscore <- scale(Desulfobacterota_taxonomy_pcoa_result$PCoA1)
Desulfobacterota_taxonomy_pcoa_result$PCoA1_scaled <- 
  -3 + ((Desulfobacterota_taxonomy_pcoa_result$PCoA1 - min(Desulfobacterota_taxonomy_pcoa_result$PCoA1)) / (max(Desulfobacterota_taxonomy_pcoa_result$PCoA1) - min(Desulfobacterota_taxonomy_pcoa_result$PCoA1))) * 6

Desulfobacterota_richness <- data.frame(sample = rownames(Desulfobacterota_abundance),  richness = apply(Desulfobacterota_abundance, 1, function(x) sum(x > 0)))


# p_Planctomycetota
# PCoA analysis based on bray curtis distance
Planctomycetota_taxonomy_pcoa <- cmdscale(
  Planctomycetota_taxonomy_bray,
  k = 3,
  eig = TRUE,
  add = FALSE,
  x.ret = TRUE)
# Extract the coordinates of the first three dimensions of the PCoA dimensionality reduction results
Planctomycetota_taxonomy_pcoa_coords <- data.frame(Planctomycetota_taxonomy_pcoa$points)
Planctomycetota_taxonomy_pcoa_coords$Sample_ID <- rownames(Planctomycetota_taxonomy_pcoa_coords)
names(Planctomycetota_taxonomy_pcoa_coords)[1:3] <- paste0("PCoA", 1:3)
Planctomycetota_taxonomy_pcoa_result <- merge(Planctomycetota_taxonomy_pcoa_coords, group, by.x = "row.names", by.y = "sample")
Planctomycetota_taxonomy_pcoa_result$facet <- c(rep("Planctomycetota", 50))
Planctomycetota_taxonomy_pcoa_result$PCoA1_zscore <- scale(Planctomycetota_taxonomy_pcoa_result$PCoA1)
Planctomycetota_taxonomy_pcoa_result$PCoA1_scaled <- 
  -(-3 + ((Planctomycetota_taxonomy_pcoa_result$PCoA1 - min(Planctomycetota_taxonomy_pcoa_result$PCoA1)) / (max(Planctomycetota_taxonomy_pcoa_result$PCoA1) - min(Planctomycetota_taxonomy_pcoa_result$PCoA1))) * 6)

Planctomycetota_richness <- data.frame(sample = rownames(Planctomycetota_abundance),  richness = apply(Planctomycetota_abundance, 1, function(x) sum(x > 0)))

# p_Actinomycetota
# PCoA analysis based on bray curtis distance
Actinomycetota_taxonomy_pcoa <- cmdscale(
  Actinomycetota_taxonomy_bray,
  k = 3,
  eig = TRUE,
  add = FALSE,
  x.ret = TRUE)
# Extract the coordinates of the first three dimensions of the PCoA dimensionality reduction results
Actinomycetota_taxonomy_pcoa_coords <- data.frame(Actinomycetota_taxonomy_pcoa$points)
Actinomycetota_taxonomy_pcoa_coords$Sample_ID <- rownames(Actinomycetota_taxonomy_pcoa_coords)
names(Actinomycetota_taxonomy_pcoa_coords)[1:3] <- paste0("PCoA", 1:3)
Actinomycetota_taxonomy_pcoa_result <- merge(Actinomycetota_taxonomy_pcoa_coords, group, by.x = "row.names", by.y = "sample")
Actinomycetota_taxonomy_pcoa_result$facet <- c(rep("Actinomycetota", 50))
Actinomycetota_taxonomy_pcoa_result$PCoA1_zscore <- scale(Actinomycetota_taxonomy_pcoa_result$PCoA1)
Actinomycetota_taxonomy_pcoa_result$PCoA1_scaled <- 
  -(-3 + ((Actinomycetota_taxonomy_pcoa_result$PCoA1 - min(Actinomycetota_taxonomy_pcoa_result$PCoA1)) / (max(Actinomycetota_taxonomy_pcoa_result$PCoA1) - min(Actinomycetota_taxonomy_pcoa_result$PCoA1))) * 6)

Actinomycetota_richness <- data.frame(sample = rownames(Actinomycetota_abundance),  richness = apply(Actinomycetota_abundance, 1, function(x) sum(x > 0)))

# p_Verrucomicrobiota
# PCoA analysis based on bray curtis distance
Verrucomicrobiota_taxonomy_pcoa <- cmdscale(
  Verrucomicrobiota_taxonomy_bray,
  k = 3,
  eig = TRUE,
  add = FALSE,
  x.ret = TRUE)
# Extract the coordinates of the first three dimensions of the PCoA dimensionality reduction results
Verrucomicrobiota_taxonomy_pcoa_coords <- data.frame(Verrucomicrobiota_taxonomy_pcoa$points)
Verrucomicrobiota_taxonomy_pcoa_coords$Sample_ID <- rownames(Verrucomicrobiota_taxonomy_pcoa_coords)
names(Verrucomicrobiota_taxonomy_pcoa_coords)[1:3] <- paste0("PCoA", 1:3)
Verrucomicrobiota_taxonomy_pcoa_result <- merge(Verrucomicrobiota_taxonomy_pcoa_coords, group, by.x = "row.names", by.y = "sample")
Verrucomicrobiota_taxonomy_pcoa_result$facet <- c(rep("Verrucomicrobiota", 50))
Verrucomicrobiota_taxonomy_pcoa_result$PCoA1_zscore <- scale(Verrucomicrobiota_taxonomy_pcoa_result$PCoA1)
Verrucomicrobiota_taxonomy_pcoa_result$PCoA1_scaled <- 
  -3 + ((Verrucomicrobiota_taxonomy_pcoa_result$PCoA1 - min(Verrucomicrobiota_taxonomy_pcoa_result$PCoA1)) / (max(Verrucomicrobiota_taxonomy_pcoa_result$PCoA1) - min(Verrucomicrobiota_taxonomy_pcoa_result$PCoA1))) * 6

Verrucomicrobiota_richness <- data.frame(sample = rownames(Verrucomicrobiota_abundance),  richness = apply(Verrucomicrobiota_abundance, 1, function(x) sum(x > 0)))

total_taxonomy_pcoa_result <-
  rbind(Desulfobacterota_taxonomy_pcoa_result, 
        Bacillota_taxonomy_pcoa_result,
        Pseudomonadota_taxonomy_pcoa_result,  
        Dysgonomonadaceae_taxonomy_pcoa_result, 
        UBA932_taxonomy_pcoa_result, 
        Bacteroidaceae_taxonomy_pcoa_result, 
        Paludibacteraceae_taxonomy_pcoa_result, 
        Rikenellaceae_taxonomy_pcoa_result, 
        Tannerellaceae_taxonomy_pcoa_result, 
        Bacteroidota_taxonomy_pcoa_result,  
        Planctomycetota_taxonomy_pcoa_result, 
        Actinomycetota_taxonomy_pcoa_result, 
        Verrucomicrobiota_taxonomy_pcoa_result)

total_taxonomy_pcoa_result$group <- factor(total_taxonomy_pcoa_result$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

total_taxonomy_pcoa_result$facet <- factor(total_taxonomy_pcoa_result$facet, levels = c("Desulfobacterota", 
                                                                                        "Bacillota",
                                                                                        "Pseudomonadota",
                                                                                        "Actinomycetota", 
                                                                                        "Planctomycetota", 
                                                                                        "Dysgonomonadaceae", 
                                                                                        "UBA932", 
                                                                                        "Bacteroidaceae", 
                                                                                        "Paludibacteraceae", 
                                                                                        "Rikenellaceae", 
                                                                                        "Tannerellaceae", 
                                                                                        "Bacteroidota",  
                                                                                        "Verrucomicrobiota"))


library(ggridges)

ridges <- 
  ggplot(total_taxonomy_pcoa_result,
         aes(x = PCoA1_scaled, y = 1, fill = group, color = group)
  ) +
  geom_density_ridges(panel_scaling = TRUE,
                      scale = 1,
                      alpha = 0.5,
                      linewidth = 0.25,
                      show.legend = F,
                      jittered_points = TRUE, 
                      point_y_offset = -0.5,
                      position = position_points_jitter(width = 0.1, height = 0), 
                      point_shape = '|', 
                      point_size = 2, 
                      point_alpha = 1, 
                      alpha = 0.8
  ) +
  xlim(c(-5,5)) +
  coord_cartesian(xlim = c(-3, 3)) + 
  facet_wrap(~ facet, 
             scales = "free_y",
             nrow = 13, 
             ncol = 1, 
             shrink = TRUE, 
             dir = "v"
  ) +
  labs(x = "Scaled PCoA1") +
  scale_color_manual(values = c("Wild" = wild_color, "Hos" = hos_color, "Res" = res_color, "Fam" = fam_color, "Lab" = lab_color)) +
  scale_fill_manual(values = c("Wild" = wild_color, "Hos" = hos_color, "Res" = res_color, "Fam" = fam_color, "Lab" = lab_color)) +
  theme_ridges() +
  theme(
    axis.title.x.bottom = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 2, color = "black"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 2, color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.5, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x.bottom = unit(2, "pt"),
    axis.ticks.length.y.left = unit(0, "pt"),
    axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.y = element_blank(),
    panel.border =  element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(5, "pt"),
    strip.background = element_blank(),
    strip.clip = "on",
    strip.placement = "inside",
    strip.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    strip.switch.pad.wrap = unit(0, "pt")
  )

print(ridges)

ggsave("taxonomy/25_ridges_PCoA1.svg", plot = ridges, width = 3, height = 10, units = "cm", dpi = 300)
# ggsave("taxonomy/6_downstream_analysis/25_ridges_PCoA1.tiff", plot = ridges, width = 3, height = 10, units = "cm", dpi = 300) 





### alluvial plot
library(ggalluvial)
alluvial_plot <- read.table("taxonomy/26_alluvial_plot_data.txt", header = TRUE)

alluvial_plot$taxonomy <- factor(alluvial_plot$taxonomy, levels = c("Desulfobacterota", 
                                                                "Bacillota",
                                                                "Pseudomonadota",
                                                                "Actinomycetota", 
                                                                "Planctomycetota", 
                                                                "Dysgonomonadaceae", 
                                                                "UBA932", 
                                                                "Bacteroidaceae", 
                                                                "Paludibacteraceae", 
                                                                "Rikenellaceae", 
                                                                "Tannerellaceae", 
                                                                "Bacteroidota",  
                                                                "Verrucomicrobiota"))

alluvial_plot_long <- alluvial_plot %>% pivot_longer(cols = -taxonomy, names_to = "group", values_to = "abundance")

alluvial_plot_long$group <- factor(alluvial_plot_long$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

alluvial_plot <- 
  ggplot(alluvial_plot_long, 
         aes(x = group, y = abundance, alluvium = taxonomy, stratum = taxonomy)
  ) +
  geom_alluvium(aes(fill = "grey95"), 
                alpha = 0.2, 
                width = 1/3,
                show.legend = F
  ) +
  geom_stratum(aes(fill = group, color = group),
               alpha = 0.5, 
               width = 1/3,
               linetype = 1,
               linewidth = 0.25,
               show.legend = F
               
  ) +
  labs(x = "Relative\nAbundance") +
  scale_color_manual(values = c("Wild" = wild_color, "Hos" = hos_color, "Res" = res_color, "Fam" = fam_color, "Lab" = lab_color)) +
  scale_fill_manual(values = c("Wild" = wild_color, "Hos" = hos_color, "Res" = res_color, "Fam" = fam_color, "Lab" = lab_color)) +
  facet_wrap(~ taxonomy, 
             scales = "free_y",
             nrow = 13, 
             ncol = 1, 
             shrink = TRUE, 
             dir = "v"
  ) +
  theme(
    axis.title.x.bottom = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x.bottom = unit(2, "pt"),
    axis.ticks.length.y.left = unit(0, "pt"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(5, "pt"),
    strip.background = element_blank(),
    strip.clip = "on",
    strip.placement = "inside",
    strip.text = element_blank(),
    strip.switch.pad.wrap = unit(0, "pt")
  )

print(alluvial_plot)

ggsave("taxonomy/27_alluvial_plot.svg", plot = alluvial_plot, width = 3, height = 10, units = "cm", dpi = 300)
# ggsave("taxonomy/27_alluvial_plot.tiff", plot = alluvial_plot, width = 3, height = 10, units = "cm", dpi = 300) 


### richness plot 
total_richness_result <- 
  rbind(Desulfobacterota_richness, 
        Bacillota_richness,
        Pseudomonadota_richness,  
        Dysgonomonadaceae_richness, 
        UBA932_richness, 
        Bacteroidaceae_richness, 
        Paludibacteraceae_richness, 
        Rikenellaceae_richness, 
        Tannerellaceae_richness, 
        Bacteroidota_richness,  
        Planctomycetota_richness, 
        Actinomycetota_richness, 
        Verrucomicrobiota_richness)

total_richness_result$facet <- c(rep("Desulfobacterota", 50), 
                                 rep("Bacillota", 50),
                                 rep("Pseudomonadota", 50),
                                 rep("Actinomycetota", 50), 
                                 rep("Planctomycetota", 50), 
                                 rep("Dysgonomonadaceae", 50), 
                                 rep("UBA932", 50), 
                                 rep("Bacteroidaceae", 50), 
                                 rep("Paludibacteraceae", 50), 
                                 rep("Rikenellaceae", 50), 
                                 rep("Tannerellaceae", 50), 
                                 rep("Bacteroidota", 50),  
                                 rep("Verrucomicrobiota", 50))

total_richness_result$facet <- factor(total_richness_result$facet, levels = c("Desulfobacterota", 
                                                                              "Bacillota",
                                                                              "Pseudomonadota",
                                                                              "Actinomycetota", 
                                                                              "Planctomycetota", 
                                                                              "Dysgonomonadaceae", 
                                                                              "UBA932", 
                                                                              "Bacteroidaceae", 
                                                                              "Paludibacteraceae", 
                                                                              "Rikenellaceae", 
                                                                              "Tannerellaceae", 
                                                                              "Bacteroidota",  
                                                                              "Verrucomicrobiota"))

total_richness_result$group <- c(rep(c(rep("Wild", 10), rep("Lab", 10), rep("Fam", 10), rep("Res", 10), rep("Hos", 10)), 13))

total_richness_result$group <- factor(total_richness_result$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

total_richness_result <- total_richness_result %>%
  group_by(group, facet) %>%
  summarise(richness = mean(richness, na.rm = TRUE))

richness_plot <- 
  ggplot(total_richness_result, 
         aes(x = group, y = richness, alluvium = facet, stratum = facet)
  ) +
  geom_alluvium(aes(fill = "grey95"), 
                alpha = 0.2, 
                width = 1/3,
                show.legend = F
  ) +
  geom_stratum(aes(fill = group, color = group),
               alpha = 0.5, 
               width = 1/3,
               linetype = 1,
               linewidth = 0.25,
               show.legend = F
               
  ) +
  facet_wrap(~ facet, 
             scales = "free_y",
             nrow = 13, 
             ncol = 1, 
             shrink = TRUE, 
             dir = "v"
  ) +
  labs(x = "Species\nRichness") +
  scale_color_manual(values = c("Wild" = wild_color, "Hos" = hos_color, "Res" = res_color, "Fam" = fam_color, "Lab" = lab_color)) +
  scale_fill_manual(values = c("Wild" = wild_color, "Hos" = hos_color, "Res" = res_color, "Fam" = fam_color, "Lab" = lab_color)) +
  theme(
    axis.title.x.bottom = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x.bottom = unit(2, "pt"),
    axis.ticks.length.y.left = unit(0, "pt"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(5, "pt"),
    strip.background = element_blank(),
    strip.clip = "on",
    strip.placement = "inside",
    strip.text = element_blank(),
    strip.switch.pad.wrap = unit(0, "pt")
  )

richness_plot

ridges_alluvial_richness <- 
  ridges + alluvial_plot + richness_plot
plot_layout(ncol = 3, widths = c(6, 2, 2)) 

ridges_alluvial_richness

ggsave("taxonomy/28_ridges_alluvial_richness.svg", plot = ridges_alluvial_richness, width = 8, height = 15, units = "cm", dpi = 300)
# ggsave("taxonomy/28_ridges_alluvial_richness.tiff", plot = ridges_alluvial_richness, width = 8, height = 15, units = "cm", dpi = 300)



#############################
#### Taxonomy composition ###
#############################

### composition stacked bar plot

taxonomy_data_species <- read.table("taxonomy/1_Total_coverm_contigs_Count_species.tsv", header = TRUE, sep = "\t", row.names = 1) 
group <- read.table('taxonomy/2_group.txt', sep = '\t', header = TRUE)

# Calculate relative abundance
colsum <- colSums(taxonomy_data_species)
taxonomy_data_species_relative <- sweep(taxonomy_data_species, 2, colsum, "/") 
rowmean <- rowMeans(taxonomy_data_species_relative)
taxonomy_data_species_relative$average <- rowmean

# add species name
taxonomy_data_species_relative$taxonomy <- row.names(taxonomy_data_species_relative)

# Extract taxonomy information
taxonomy <- rownames(taxonomy_data_species_relative) %>% strsplit(";")
taxonomy_info <- data.frame(do.call(rbind, taxonomy), stringsAsFactors = FALSE)
rownames(taxonomy_info) <- rownames(taxonomy_data_species_relative)
colnames(taxonomy_info) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

taxonomy_data_species_relative <- cbind(taxonomy_data_species_relative, taxonomy_info)

# Convert wide to long
abundance_long <- 
  taxonomy_data_species_relative %>% 
  pivot_longer(-c("taxonomy", "kingdom", "phylum", "class", "order", "family", "genus", "species"), names_to = "sample", values_to = "Abundance")

# Add taxonomy information and group information
abundance_result <- merge(abundance_long, group, by.x = "sample", by.y = "sample")
abundance_result$group <- factor(abundance_result$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

### group catergory
group_catergory_plot <- 
  ggplot(group,
         aes(x = sample, y = 1, fill = group, color = group)
  ) +
  geom_tile(height = 1,
            color = "white",
            linewidth = 0.05,
            lineend = "square",
            show.legend = F
  ) +
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  theme_void()

group_catergory_plot


### species level
# Calculate the abundance for each species
total_abundance_species <- 
  abundance_result %>%
  group_by(species) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total))

# high light top 11 species
top11_species <- total_abundance_species$species[1:16]
top11_species <- top11_species[!grepl("^s_unclassify", top11_species, ignore.case = TRUE)] 
abundance_result_species <- 
  abundance_result %>%
  mutate(species = ifelse(species %in% top11_species, species, "others"))

# Re-calculate the abundance for top 11 species
abundance_result_species <- 
  abundance_result_species %>%
  group_by(sample, species) %>%
  summarise(Abundance = sum(Abundance), .groups = "keep")

# Reorder the list, placing "others" at the end.
abundance_result_species<- 
  abundance_result_species %>%
  ungroup() %>%
  mutate(species = factor(species, levels = c(top11_species, "others"))) %>%
  arrange(sample, desc(species))

# add sample information
abundance_result_species <- merge(abundance_result_species, group, by.x = "sample", by.y = "sample")
abundance_result_species$group <- factor(abundance_result_species$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

# stacked bar plot species
stacked_bar_color <- c("grey40", paletteer_d("ggthemes::Green_Orange_Teal")[2:12])
abundance_result_species$species <- sub("s__", "", abundance_result_species$species)
# Make sure that the "family" column is of factor type and that the level order is correct.
unique_species <- unique(abundance_result_species$species)
abundance_result_species$species <- factor(abundance_result_species$species, levels = unique_species)


### Stacked_bar_species 
Stacked_bar_species <- 
  ggplot(abundance_result_species, 
         aes(x = sample, y = Abundance, fill = species)
  ) +
  geom_bar(stat = "identity",
           position = "fill",
           width = 0.75, 
           orientation = "x",
           show.legend = T
  ) +
  facet_wrap(~ group, 
             ncol = 5, 
             nrow = 1, 
             scales = "free_x", 
             strip.position = "bottom"
  ) +
  guides(fill = guide_legend(ncol = 1, position = "bottom", reverse = TRUE)) +
  scale_fill_manual(values = stacked_bar_color) + 
  scale_y_continuous(
    expand = c(0.005, 0.005),
    sec.axis = dup_axis(name = NULL, labels = scales::percent_format(scale = 100))
  ) +
  coord_cartesian() +
  labs(y = "Species Abundance") +
  theme(
    plot.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 00, lineheight = 2, color = "black"),
    axis.title.x.bottom = element_blank(), 
    axis.title.y.left = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.title.y.right = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y.right = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.text.y.left = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y.right = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y.left = element_blank(),
    axis.ticks.length.x.bottom = unit(0, "pt"),
    axis.ticks.length.y.left = unit(2, "pt"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(3, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(family = "Arial", face = "plain", size = 8, hjust = 1, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(), 
    legend.frame = element_blank(),
    legend.key.size = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "italic", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black")
  )  

Stacked_bar_species

ggsave("taxonomy/29_taxonomy_Stacked_bar_species.svg", plot = Stacked_bar_species, width = 5.5, height = 7, units = "cm", dpi = 300)
# ggsave("taxonomy/29_taxonomy_Stacked_bar_species.tiff", plot = Stacked_bar_species, width = 5.5, height = 7, units = "cm", dpi = 300) 

Stacked_bar_species_legend <- cowplot::get_legend(Stacked_bar_species + theme(legend.position = "right"))   



### genus level
# Calculate the abundance for each genus
total_abundance_genus <- 
  abundance_result %>%
  group_by(genus) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total))

# high light top 11 genus
top11_genus <- total_abundance_genus$genus[1:15]
top11_genus <- top11_genus[!grepl("^g_unclassify", top11_genus, ignore.case = TRUE)] 
abundance_result_genus <- 
  abundance_result %>%
  mutate(genus = ifelse(genus %in% top11_genus, genus, "others"))

# Re-calculate the abundance for top 11 genus
abundance_result_genus <- 
  abundance_result_genus %>%
  group_by(sample, genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "keep")

# Reorder the list, placing "others" at the end.
abundance_result_genus<- 
  abundance_result_genus %>%
  ungroup() %>%
  mutate(genus = factor(genus, levels = c(top11_genus, "others"))) %>%
  arrange(sample, desc(genus))

# add sample information
abundance_result_genus <- merge(abundance_result_genus, group, by.x = "sample", by.y = "sample")
abundance_result_genus$group <- factor(abundance_result_genus$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

# stacked bar plot genus
stacked_bar_color <- c("grey40", paletteer_d("ggthemes::Green_Orange_Teal")[2:12])
abundance_result_genus$genus <- sub("g__", "", abundance_result_genus$genus)
# Make sure that the "family" column is of factor type and that the level order is correct.
unique_genus <- unique(abundance_result_genus$genus)
abundance_result_genus$genus <- factor(abundance_result_genus$genus, levels = unique_genus)

### Stacked_bar_genus 
Stacked_bar_genus <- 
  ggplot(abundance_result_genus, 
         aes(x = sample, y = Abundance, fill = genus)
  ) +
  geom_bar(stat = "identity",
           position = "fill",
           width = 0.75, 
           orientation = "x",
           show.legend = T
  ) +
  facet_wrap(~ group, 
             ncol = 5, 
             nrow = 1, 
             scales = "free_x", 
             strip.position = "bottom"
  ) +
  guides(fill = guide_legend(ncol = 1, position = "bottom", reverse = TRUE)) +
  scale_fill_manual(values = stacked_bar_color) + 
  scale_y_continuous(
    expand = c(0.005, 0.005),
    sec.axis = dup_axis(name = NULL, labels = scales::percent_format(scale = 100))
  ) +
  coord_cartesian() +
  labs(y = "Genus Abundance") +
  theme(
    plot.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 00, lineheight = 2, color = "black"),
    axis.title.x.bottom = element_blank(), 
    axis.title.y.left = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.title.y.right = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y.right = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.text.y.left = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y.right = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y.left = element_blank(),
    axis.ticks.length.x.bottom = unit(0, "pt"),
    axis.ticks.length.y.left = unit(2, "pt"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(3, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(family = "Arial", face = "plain", size = 8, hjust = 1, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(), 
    legend.frame = element_blank(),
    legend.key.size = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "italic", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black")
  ) 

Stacked_bar_genus

ggsave("taxonomy/30_taxonomy_Stacked_bar_genus.svg", plot = Stacked_bar_genus, width = 5.5, height = 7, units = "cm", dpi = 300)
# ggsave("taxonomy/30_taxonomy_Stacked_bar_genus.tiff", plot = Stacked_bar_genus, width = 5.5, height = 7, units = "cm", dpi = 300) 

Stacked_bar_genus_legend <- cowplot::get_legend(Stacked_bar_genus + theme(legend.position = "right"))   



### family level

# Calculate the abundance for each family
total_abundance_family <- 
  abundance_result %>%
  group_by(family) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total))

# high light top 11 family
top11_family <- total_abundance_family$family[1:13]
top11_family <- top11_family[!grepl("^f_unclassify", top11_family, ignore.case = TRUE)] 
abundance_result_family <- 
  abundance_result %>%
  mutate(family = ifelse(family %in% top11_family, family, "others"))

# Re-calculate the abundance for top 11 class
abundance_result_family <- 
  abundance_result_family %>%
  group_by(sample, family) %>%
  summarise(Abundance = sum(Abundance), .groups = "keep")

# Reorder the list, placing "others" at the end.
abundance_result_family <- 
  abundance_result_family %>%
  ungroup() %>%
  mutate(family = factor(family, levels = c(top11_family, "others"))) %>%
  arrange(sample, desc(family))

# add sample information
abundance_result_family <- merge(abundance_result_family, group, by.x = "sample", by.y = "sample")
abundance_result_family$group <- factor(abundance_result_family$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

# stacked bar plot family
stacked_bar_color <- c("grey40", paletteer_d("ggthemes::Green_Orange_Teal")[2:12])
abundance_result_family$family <- sub("f__", "", abundance_result_family$family)
# Make sure that the "family" column is of factor type and that the level order is correct.
unique_family <- unique(abundance_result_family$family)
abundance_result_family$family <- factor(abundance_result_family$family, levels = unique_family)

### Stacked_bar_family 
Stacked_bar_family <- 
  ggplot(abundance_result_family, 
         aes(x = sample, y = Abundance, fill = family)
  ) +
  geom_bar(stat = "identity",
           position = "fill",
           width = 0.75, 
           orientation = "x",
           show.legend = T
  ) +
  facet_wrap(~ group, 
             ncol = 5, 
             nrow = 1, 
             scales = "free_x", 
             strip.position = "bottom"
  ) +
  guides(fill = guide_legend(ncol = 1, position = "bottom", reverse = TRUE)) +
  scale_fill_manual(values = stacked_bar_color) + 
  scale_y_continuous(
    expand = c(0.005, 0.005),
    sec.axis = dup_axis(name = NULL, labels = scales::percent_format(scale = 100))
  ) +
  coord_cartesian() +
  labs(y = "Family Abundance") +
  theme(
    plot.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 00, lineheight = 2, color = "black"),
    axis.title.x.bottom = element_blank(), 
    axis.title.y.left = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.title.y.right = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y.right = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.text.y.left = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y.right = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y.left = element_blank(),
    axis.ticks.length.x.bottom = unit(0, "pt"),
    axis.ticks.length.y.left = unit(2, "pt"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(3, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(family = "Arial", face = "plain", size = 8, hjust = 1, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(), 
    legend.frame = element_blank(),
    legend.key.size = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "italic", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black")
  ) 

Stacked_bar_family

ggsave("taxonomy/31_taxonomy_Stacked_bar_family.svg", plot = Stacked_bar_family, width = 5.5, height = 7, units = "cm", dpi = 300)
# ggsave("taxonomy/31_taxonomy_Stacked_bar_family.tiff", plot = Stacked_bar_family, width = 5.5, height = 7, units = "cm", dpi = 300) 

Stacked_bar_family_legend <- cowplot::get_legend(Stacked_bar_family + theme(legend.position = "right"))



### order level
# Calculate the abundance for each order
total_abundance_order <- 
  abundance_result %>%
  group_by(order) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total))

# high light top 11 order
top11_order <- total_abundance_order$order[1:13]
top11_order <- top11_order[!grepl("^o_unclassify", top11_order, ignore.case = TRUE)] 
abundance_result_order <- 
  abundance_result %>%
  mutate(order = ifelse(order %in% top11_order, order, "others"))

# Re-calculate the abundance for top 10 order
abundance_result_order <- 
  abundance_result_order %>%
  group_by(sample, order) %>%
  summarise(Abundance = sum(Abundance), .groups = "keep")

# Reorder the list, placing "others" at the end.
abundance_result_order<- 
  abundance_result_order %>%
  ungroup() %>%
  mutate(order = factor(order, levels = c(top11_order, "others"))) %>%
  arrange(sample, desc(order))

# add sample information
abundance_result_order <- merge(abundance_result_order, group, by.x = "sample", by.y = "sample")
abundance_result_order$group <- factor(abundance_result_order$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

# stacked bar plot order
stacked_bar_color <- c("grey40", paletteer_d("ggthemes::Green_Orange_Teal")[2:12])
abundance_result_order$order <- sub("o__", "", abundance_result_order$order)
# Make sure that the "family" column is of factor type and that the level order is correct.
unique_order <- unique(abundance_result_order$order)
abundance_result_order$order <- factor(abundance_result_order$order, levels = unique_order)

### Stacked_bar_order 
Stacked_bar_order <- 
  ggplot(abundance_result_order, 
         aes(x = sample, y = Abundance, fill = order)
  ) +
  geom_bar(stat = "identity",
           position = "fill",
           width = 0.75, 
           orientation = "x",
           show.legend = T
  ) +
  facet_wrap(~ group, 
             ncol = 5, 
             nrow = 1, 
             scales = "free_x", 
             strip.position = "bottom"
  ) +
  guides(fill = guide_legend(ncol = 1, position = "bottom", reverse = TRUE)) +
  scale_fill_manual(values = stacked_bar_color) + 
  scale_y_continuous(
    expand = c(0.005, 0.005),
    sec.axis = dup_axis(name = NULL, labels = scales::percent_format(scale = 100))
  ) +
  coord_cartesian() +
  labs(y = "Order Abundance") +
  theme(
    plot.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 00, lineheight = 2, color = "black"),
    axis.title.x.bottom = element_blank(), 
    axis.title.y.left = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.title.y.right = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y.right = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.text.y.left = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y.right = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y.left = element_blank(),
    axis.ticks.length.x.bottom = unit(0, "pt"),
    axis.ticks.length.y.left = unit(2, "pt"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(3, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(family = "Arial", face = "plain", size = 8, hjust = 1, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(), 
    legend.frame = element_blank(),
    legend.key.size = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "italic", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black")
  ) 

Stacked_bar_order

ggsave("taxonomy/32_taxonomy_Stacked_bar_order.svg", plot = Stacked_bar_order, width = 5.5, height = 7, units = "cm", dpi = 300)
# ggsave("taxonomy/32_taxonomy_Stacked_bar_order.tiff", plot = Stacked_bar_order, width = 5.5, height = 7, units = "cm", dpi = 300) 

Stacked_bar_order_legend <- cowplot::get_legend(Stacked_bar_order + theme(legend.position = "right")) 



# Calculate the abundance for each class
total_abundance_class <- 
  abundance_result %>%
  group_by(class) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total))

# high light top 12 class
top11_classes <- total_abundance_class$class[1:13]
top11_classes <- top11_classes[!grepl("^c_unclassify", top11_classes, ignore.case = TRUE)] 
abundance_result_class <- 
  abundance_result %>%
  mutate(class = ifelse(class %in% top11_classes, class, "others"))

# Re-calculate the abundance for top 11 class
abundance_result_class <- 
  abundance_result_class %>%
  group_by(sample, class) %>%
  summarise(Abundance = sum(Abundance), .groups = "keep")

# Reorder the list, placing "others" at the end.
abundance_result_class <- 
  abundance_result_class %>%
  ungroup() %>%
  mutate(class = factor(class, levels = c(top11_classes, "others"))) %>%
  arrange(sample, desc(class))

# add sample information
abundance_result_class <- merge(abundance_result_class, group, by.x = "sample", by.y = "sample")
abundance_result_class$group <- factor(abundance_result_class$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

# stacked bar plot class
abundance_result_class$class <- sub("c__", "", abundance_result_class$class)
# Make sure that the "class" column is of factor type and that the level order is correct.
unique_classes <- unique(abundance_result_class$class)
abundance_result_class$class <- factor(abundance_result_class$class, levels = unique_classes)

stacked_bar_color <- c(others = "grey40", 
                       Synergistia = paletteer_d("ggthemes::Green_Orange_Teal")[2],
                       Bacilli_A = paletteer_d("ggthemes::Green_Orange_Teal")[3],
                       Desulfarculia = paletteer_d("ggthemes::Green_Orange_Teal")[4],
                       Actinomycetes = paletteer_d("ggthemes::Green_Orange_Teal")[5],
                       Planctomycetia = paletteer_d("ggthemes::Green_Orange_Teal")[6],
                       Gammaproteobacteria = paletteer_d("ggthemes::Green_Orange_Teal")[7],
                       `SZUA-567` = paletteer_d("ggthemes::Green_Orange_Teal")[8],
                       Bacilli = paletteer_d("ggthemes::Green_Orange_Teal")[9],
                       Desulfovibrionia = paletteer_d("ggthemes::Green_Orange_Teal")[10],
                       Clostridia = paletteer_d("ggthemes::Green_Orange_Teal")[11],
                       Bacteroidia = paletteer_d("ggthemes::Green_Orange_Teal")[12],
                       Alphaproteobacteria = paletteer_d("lisa::AndyWarhol")[1],
                       Desulfobulbia = paletteer_d("lisa::AndyWarhol")[2])

unique_classes
paletteer_d("ggthemes::Green_Orange_Teal")
paletteer_d("lisa::AndyWarhol")

### Stacked_bar_class 
Stacked_bar_class <- 
  ggplot(abundance_result_class, 
         aes(x = sample, y = Abundance, fill = class)
  ) +
  geom_bar(stat = "identity",
           position = "fill",
           width = 0.75, 
           orientation = "x",
           show.legend = T
  ) +
  facet_wrap(~ group, 
             ncol = 5, 
             nrow = 1, 
             scales = "free_x", 
             strip.position = "bottom"
  ) +
  guides(fill = guide_legend(ncol = 1, position = "bottom", reverse = TRUE)) +
  scale_fill_manual(values = stacked_bar_color, drop = F) + 
  scale_y_continuous(
    expand = c(0.005, 0.005),
    sec.axis = dup_axis(name = NULL, labels = scales::percent_format(scale = 100))
  ) +
  coord_cartesian() +
  labs(y = "Class Abundance") +
  theme(
    plot.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 00, lineheight = 2, color = "black"),
    axis.title.x.bottom = element_blank(), 
    axis.title.y.left = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.title.y.right = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y.right = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.text.y.left = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y.right = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y.left = element_blank(),
    axis.ticks.length.x.bottom = unit(0, "pt"),
    axis.ticks.length.y.left = unit(2, "pt"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(3, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(family = "Arial", face = "plain", size = 8, hjust = 1, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(), 
    legend.frame = element_blank(),
    legend.key.size = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "italic", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black")
  ) 

Stacked_bar_class

ggsave("taxonomy/33_taxonomy_Stacked_bar_class.svg", plot = Stacked_bar_class, width = 6.5, height = 5, units = "cm", dpi = 300)
# ggsave("taxonomy/33_taxonomy_Stacked_bar_class.tiff", plot = Stacked_bar_class, width = 6.5, height = 5, units = "cm", dpi = 300) 

Stacked_bar_class_legend <- cowplot::get_legend(Stacked_bar_class + theme(legend.position = "right"))


# Calculate the abundance for each phylum
total_abundance_phylum <- 
  abundance_result %>%
  group_by(phylum) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total))

# high light top 11 phylum
top11_phylums <- total_abundance_phylum$phylum[1:13]
top11_phylums <- top11_phylums[!grepl("^p_unclassify", top11_phylums, ignore.case = TRUE)]
abundance_result_phylum <- 
  abundance_result %>%
  mutate(phylum = ifelse(phylum %in% top11_phylums, phylum, "others"))

# Re-calculate the abundance for top 10 phylum
abundance_result_phylum <- 
  abundance_result_phylum %>%
  group_by(sample, phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "keep")

# Reorder the list, placing "others" at the end.
abundance_result_phylum <- 
  abundance_result_phylum %>%
  ungroup() %>%
  mutate(phylum = factor(phylum, levels = c(top11_phylums, "others"))) %>%
  arrange(sample, desc(phylum))

# add sample information
abundance_result_phylum <- merge(abundance_result_phylum, group, by.x = "sample", by.y = "sample")
abundance_result_phylum$group <- factor(abundance_result_phylum$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

### stacked bar plot phylum
stacked_bar_color <- c("grey40", paletteer_d("ggthemes::Green_Orange_Teal"))
abundance_result_phylum$phylum <- sub("p__", "", abundance_result_phylum$phylum)
# Make sure that the "phylum" column is of factor type and that the level order is correct.
unique_phylums <- unique(abundance_result_phylum$phylum)
abundance_result_phylum$phylum <- factor(abundance_result_phylum$phylum, levels = unique_phylums)

Stacked_bar_phylum <- 
  ggplot(abundance_result_phylum, 
         aes(x = sample, y = Abundance, fill = phylum)
  ) +
  geom_bar(stat = "identity",
           position = "fill",
           width = 0.75, 
           orientation = "x",
           show.legend = T
  ) +
  facet_wrap(~ group, 
             ncol = 5, 
             nrow = 1, 
             scales = "free_x", 
             strip.position = "bottom"
  ) +
  guides(fill = guide_legend(ncol = 1, position = "bottom", reverse = TRUE)) +
  scale_fill_manual(values = stacked_bar_color) + 
  scale_y_continuous(
    expand = c(0.005, 0.005),
    sec.axis = dup_axis(name = NULL, labels = scales::percent_format(scale = 100))
  ) +
  coord_cartesian() +
  labs(y = "Phylum Abundance") +
  theme(
    plot.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 00, lineheight = 2, color = "black"),
    axis.title.x.bottom = element_blank(), 
    axis.title.y.left = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.title.y.right = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y.right = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.text.y.left = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y.right = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y.left = element_blank(),
    axis.ticks.length.x.bottom = unit(0, "pt"),
    axis.ticks.length.y.left = unit(2, "pt"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(3, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(family = "Arial", face = "plain", size = 8, hjust = 1, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(), 
    legend.frame = element_blank(),
    legend.key.size = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "italic", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black")
  ) 

Stacked_bar_phylum

ggsave("taxonomy/34_taxonomy_Stacked_bar_phylum.svg", plot = Stacked_bar_phylum, width = 5.5, height = 7, units = "cm", dpi = 300)
# ggsave("taxonomy/34_taxonomy_Stacked_bar_phylum.tiff", plot = Stacked_bar_phylum, width = 5.5, height = 7, units = "cm", dpi = 300) 

Stacked_bar_phylum_legend <- cowplot::get_legend(Stacked_bar_phylum + theme(legend.position = "right"))

legend <- cowplot::plot_grid(Stacked_bar_phylum_legend, 
                             Stacked_bar_class_legend, 
                             Stacked_bar_order_legend, 
                             Stacked_bar_family_legend, 
                             Stacked_bar_genus_legend, 
                             Stacked_bar_species_legend, 
                             ncol = 6)

legend

ggsave("taxonomy/35_legend.svg", plot = legend, width = 24, height = 8, units = "cm", dpi = 300)
# ggsave("taxonomy/35_legend.tiff", plot = legend, width = 24, height = 8, units = "cm", dpi = 300)




#######################################
#### Taxonomy differential analysis ###
#######################################

# LDA score plot
LDA_score <- read_excel("taxonomy/36_lefse_res.xlsx", sheet = 2)
LDA_score$enrichment <- factor(LDA_score$enrichment, levels = c( "Wild", "Lab", "Fam", "Res", "Hos"))
LDA_score$Pvalue <- as.numeric(as.character(LDA_score$`p value`))
LDA_score$log10Pvalue <- -log10(LDA_score$Pvalue)
# group by LDA_score$enrichment & reserve top5
LDA_score <- LDA_score %>%
  group_by(enrichment) %>%
  arrange(desc(`LDA score`)) %>%
  slice_head(n = 5) %>%
  ungroup()

LDA_score$species <- gsub("^s__", "", LDA_score$species)
LDA_score$species <- factor(LDA_score$species, levels = LDA_score$species[order(LDA_score$enrichment, LDA_score$`LDA score`)])
LDA_score <- as.data.frame(LDA_score)
rownames(LDA_score) <- LDA_score$taxonomy

# LDA score
LDA_p0 <- 
  ggplot(LDA_score, 
         aes(x = `LDA score`, y = species, alpha = log10Pvalue, fill = enrichment)
  ) +
  geom_bar(stat = "identity",
           width = 0.85, 
           orientation = "y",
           show.legend = T
  ) +
  geom_vline(xintercept = 2.5, 
             color = "black", 
             linetype = "dashed", 
             linewidth = 0.25
  ) +
  labs(x = "LDA Score", 
       y = NULL, 
       fill = "Enrichment", 
       alpha = "-log10(P)"
  ) +
  guides(fill = guide_legend(ncol = 1, nrow = 4, position = "bottom", reverse = TRUE),
         alpha = guide_legend(ncol = 1, nrow = 4, position = "bottom", reverse = TRUE)
  )+
  scale_fill_manual(values = c("Bger" = Bger_color, "Esin" = Esin_color, "Pame" = Pame_color, "Pful" = Pful_color)) +
  scale_x_continuous(expand = c(0.01, 0.01), position = "bottom") +
  scale_y_discrete(position = "right", breaks = LDA_score$species, labels = LDA_score$species) + 
  theme(
    axis.title.x.bottom = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 2, color = "black"),
    axis.title.y.left = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 2, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 2, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "italic", size = 6, hjust = 0, vjust = 0.5, 
                               angle = 0, lineheight = 2, 
                               color = c(rep(Pful_color, 5), rep(Pame_color, 5), rep(Bger_color, 5), rep(Esin_color, 5))),
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x.bottom = unit(2, "pt"),
    axis.ticks.length.y.left = unit(2, "pt"),
    axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.y = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.frame = element_blank(),
    legend.key.size = unit(10, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.key.width = unit(10, "pt"),
    legend.key.spacing.y  = unit(5, "pt"),
    legend.key.spacing.x  = unit(5, "pt"),
    legend.direction = "horizontal",
    legend.position = "none",
    legend.position.inside = c(0.1, 0.8),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.text.position = "right",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title.position = "top"
  )

LDA_p0

# taxonomy classify plot
LDA_p1_color <- as.vector(paletteer::paletteer_d("ggthemes::Classic_10_Medium"))

LDA_score$class <- gsub("^c__", "", LDA_score$class)

LDA_score$class <- factor(LDA_score$class, levels = c("Bacteroidia", "Clostridia", "Desulfovibrionia", "Bacilli", "SZUA_567",
                                                      "Gammaproteobacteria", "Planctomycetia", "Actinomycetes", "Synergistia"))

LDA_score$phylum <- gsub("^p_", "", LDA_score$phylum)

LDA_score$phylum <- factor(LDA_score$phylum, levels = c("Bacteroidota", "Desulfobacterota", "Planctomycetota", "Synergistota", 
                                                        "Bacillota_A", "Bacillota", "Pseudomonadota", "Actinomycetota"))

LDA_p1 <- 
  ggplot(LDA_score, aes(x = 0, y = species, fill = class)) +
  geom_tile(color = "white",
            linewidth = 0,
            lineend = "square",
            height = 0.85,
            show.legend = T
  ) +
  guides(fill = guide_legend(ncol = 3, nrow = 4, position = "bottom", reverse = FALSE))+
  scale_fill_manual(values = c("Bacteroidia" = paletteer_d("ggthemes::Green_Orange_Teal")[12],
                               "Clostridia" = paletteer_d("ggthemes::Green_Orange_Teal")[11],
                               "Desulfovibrionia" = paletteer_d("ggthemes::Green_Orange_Teal")[10],
                               "Bacilli" = paletteer_d("ggthemes::Green_Orange_Teal")[9],
                               "SZUA_567" = paletteer_d("ggthemes::Green_Orange_Teal")[8],
                               "Gammaproteobacteria" = paletteer_d("ggthemes::Green_Orange_Teal")[7],
                               "Planctomycetia" = paletteer_d("ggthemes::Green_Orange_Teal")[6],
                               "Actinomycetes" = paletteer_d("ggthemes::Green_Orange_Teal")[5],
                               "Synergistia" = paletteer_d("ggthemes::Green_Orange_Teal")[2])
  ) +
  scale_x_discrete(expand = c(0.01, 0.01)) +
  theme_minimal() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.frame = element_blank(),
        legend.key.size = unit(10, "pt"), 
        legend.key.height = unit(10, "pt"), 
        legend.key.width = unit(10, "pt"),
        legend.key.spacing.y = unit(5, "pt"),
        legend.key.spacing.x = unit(5, "pt"),
        legend.text = element_text(family = "Arial", face = "italic", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.title.position = "top",
        legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 2, color = "black")
  )

LDA_p1

# relative abundance heatmap
relative_abundance_species_highlight <- as.data.frame(read_excel("taxonomy/36_lefse_res.xlsx", sheet = 3))
rownames(relative_abundance_species_highlight) <- relative_abundance_species_highlight$taxonomy
relative_abundance_species_highlight <- relative_abundance_species_highlight %>% dplyr:: select(-1)

# z-score
relative_abundance_species_highlight_zscore <- 
  relative_abundance_species_highlight %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame()

colnames(relative_abundance_species_highlight_zscore) <- colnames(relative_abundance_species_highlight)

# wide to long
relative_abundance_species_highlight_long <- 
  relative_abundance_species_highlight_zscore %>%
  rownames_to_column(var = "taxonomy") %>%
  pivot_longer(cols = -taxonomy, names_to = "sample", values_to = "abundance")

# Split the "taxonomy" column
relative_abundance_species_highlight_long <- relative_abundance_species_highlight_long %>%
  dplyr::mutate(sample_name = relative_abundance_species_highlight_long$sample) %>%
  tidyr::separate(sample_name, into = c("group", "subgroup"), sep = "_") %>%
  dplyr::mutate(taxonomy_sep = relative_abundance_species_highlight_long$taxonomy) %>%
  tidyr::separate(taxonomy_sep, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = "\\|", fill = "right") %>%
  dplyr::select(-subgroup)

relative_abundance_species_highlight_long$species <- gsub("^s_", "", relative_abundance_species_highlight_long$species)

# set factor levels
relative_abundance_species_highlight_long$group <- factor(relative_abundance_species_highlight_long$group, levels = c("Esin", "Bger", "Pame", "Pful"))
relative_abundance_species_highlight_long$species <- factor(relative_abundance_species_highlight_long$species, levels = LDA_score$species)
relative_abundance_species_highlight_long$sample <- factor(relative_abundance_species_highlight_long$sample, c(paste0("Esin_", 1:6), paste0("Bger_", 1:6), paste0("Pame_", 1:6), paste0("Pful_", 1:6)))

# calculate the species with the highest expression in each sample
high_expression_species <- relative_abundance_species_highlight_long %>%
  group_by(group) %>%
  top_n(5, abundance) %>%
  ungroup()

# set color gradient limits
min <- min(as.matrix(relative_abundance_species_highlight_zscore))
max <- max(as.matrix(relative_abundance_species_highlight_zscore))

# plot heatmap
paletteer_d("MoMAColors::Avedon")

LDA_p2 <- 
  ggplot(relative_abundance_species_highlight_long, 
         aes(x = sample, y = species, fill = abundance)
  ) +
  geom_tile(color = "white",
            linewidth = 0.25,
            lineend = "square",
            height = 1,
            show.legend = T
  ) +
  #annotate("rect", xmin = 0.5, xmax = 6.5, ymin = 15.5 , ymax = 20.5, fill = NA, color = "red", linetype = 1, linewidth = 0.25, lineend = "square") +
  #annotate("rect", xmin = 6.5, xmax = 12.5, ymin = 10.5 , ymax = 15.5, fill = NA, color = "red", linetype = 1, linewidth = 0.25, lineend = "square") + 
  #annotate("rect", xmin = 12.5, xmax = 18.5, ymin = 5.5 , ymax = 10.5, fill = NA, color = "red", linetype = 1, linewidth = 0.25, lineend = "square") + 
  #annotate("rect", xmin = 18.5, xmax = 24.5, ymin = 0.5 , ymax = 5.5, fill = NA, color = "red", linetype = 1, linewidth = 0.25, lineend = "square") + 
  scale_fill_gradient2(low = "#648C16FF", 
                       mid = "white", 
                       high = "#FF7200FF", 
                       limit = c(min, max),
                       breaks = c(0, 2, 4)
  ) +
  labs(x = "Sample", 
       y = "Species", 
       fill = "Abundance Z-score"
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25),
    legend.position = "none",
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.frame = element_blank(),
    legend.key.size = unit(10, "pt"), 
    legend.key.width = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

LDA_p2

LDA_plot <- LDA_p0 %>%
  insert_left(LDA_p1, width = 0.2) %>%
  insert_left(LDA_p2, width = 2)

LDA_plot

ggsave("taxonomy/37_LDA_score_plot.svg", plot = LDA_plot, width = 10, height = 6, units = "cm", dpi = 300)
# ggsave("taxonomy37_LDA_score_plot.tiff", plot = LDA_plot, width = 10, height = 6, units = "cm", dpi = 300)

### LDA legend
LDA_p0_legend <- cowplot::get_legend(LDA_p0 + theme(legend.position = "right"))
LDA_p1_legend <- cowplot::get_legend(LDA_p1+ theme(legend.position = "bottom"))
LDA_p2_legend <- cowplot::get_legend(LDA_p2+ theme(legend.position = "right"))

LDA_plot_legend <- cowplot::plot_grid(Stacked_bar_class_legend, LDA_p2_legend, LDA_p0_legend, ncol = 3)

LDA_plot_legend

ggsave("taxonomy/38_LDA_score_plot_legend.svg", plot = LDA_plot_legend, width = 30, height = 10, units = "cm", dpi = 300)
# ggsave("taxonomy/38_LDA_score_plot_legend.tiff", plot = LDA_plot_legend, width = 30, height = 10, units = "cm", dpi = 300)

# Differential taxonomy boxplot
species <- as.data.frame(read_excel("taxonomy/6_downstream_analysis/36_lefse_res.xlsx", sheet = 3))
rownames(species) <- species$taxonomy
species_2 <- merge(species, LDA_score, by = "row.names")
rownames(species_2) <- species_2$taxonomy.x

species_2 <- species_2[match(row.names(LDA_score), species_2$Row.names), ]
species_2 <- as.data.frame(t(species_2 %>% dplyr::select(c(3:26))))
species_2$group <- factor(c(rep("Pame",6), rep("Pful",6), rep("Esin",6), rep("Bger",6)),  levels = c("Esin", "Bger", "Pame", "Pful"))
species_2$sample <- row.names(species_2)

species_long <- species_2 %>%
  pivot_longer(
    cols = c(1:20),
    names_to = "species",
    values_to = "abundance")

species_long <- species_long %>%
  separate(species, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = "\\|", fill = "right")

species_long$species <- gsub("^s_", "", species_long$species)

facet_levels <- c(LDA_score$species[16:20], LDA_score$species[11:15], LDA_score$species[6:10], LDA_score$species[1:5])

species_long$species <- factor(species_long$species, levels = facet_levels)

# The Y-axis uses scientific notation. For all values except zero, it also uses scientific notation and retains two significant figures.
custom_labels <- function(x) {
  ifelse(x == 0, "0", sprintf("%.1e", x))
}

# 1.3 x [maximum value for each species] as the upper limit for the Y-axis.
species_max <- species_long %>%
  group_by(species) %>%
  summarise(max_abundance = max(abundance, na.rm = TRUE) * 1.3) %>%
  ungroup()

get_ylimits <- function(species_name) {
  max_val <- species_max$max_abundance[species_max$species == species_name]
  c(0, max_val)
}

# background color
background_data <- data.frame(
  xmin = c(-Inf, 1.5, 2.5, 3.5),
  xmax = c(1.5, 2.5, 3.5, Inf),
  group = c("Esin", "Bger", "Pame", "Pful")
)

# box plot
top5_species_box_plot <- 
  ggplot(species_long, 
         aes(x = group, y = abundance)
  ) +
  # plot 1% Y-axis hline
  geom_hline(yintercept = 0.01,
             color = "grey30", 
             linetype = "dashed",
             alpha = 0.5,
             linewidth = 0.25
  ) + 
  # scatter plot
  geom_point(aes(fill = group, color = group),
             shape = 21,
             stat = "identity",
             position = "identity",
             size = 2, 
             alpha = 0.6,
             show.legend = T
  ) +
  # box plot
  geom_boxplot(aes(fill = group, color = group),
               outlier.shape = NA, 
               alpha = 0.3,
               width = 0.8,
               staplewidth = 0.8,
               linewidth = 0.1,
               show.legend = T
  ) +
  # facet
  facet_wrap(~ species, 
             scales = "free_y", 
             nrow = 5, 
             ncol = 4, 
             shrink = TRUE, 
             dir = "v"
  ) +
  # group color
  scale_fill_manual(values = c("Bger" = Bger_color, "Esin" = Esin_color, "Pame" = Pame_color, "Pful" = Pful_color)) +
  scale_color_manual(values = c("Bger" = Bger_color, "Esin" = Esin_color, "Pame" = Pame_color, "Pful" = Pful_color)) +
  # 1.05 x [maximum value for each species] as the upper limit for the Y-axis.
  scale_y_continuous(labels = custom_labels,
                     limits = function(x) c(0, max(x) * 1.05)
  ) +
  # wilcox.test
  geom_signif(
    comparisons = list(c("Esin", "Bger"),
                       c("Bger", "Pame"), 
                       c("Pame", "Pful"), 
                       c("Esin", "Pame"),
                       c("Bger", "Pful"), 
                       c("Esin", "Pful")),
    map_signif_level = function(p) {
      ifelse(p < 0.001, "***", 
             ifelse(p < 0.01, "**", 
                    ifelse(p < 0.05, "*", " ")))
    },
    test = "wilcox.test",
    test.args = list(exact = FALSE), 
    textsize = 3,
    vjust = 0.6,
    tip_length = 0,
    size = 0.25,
    step_increase = 0.125
  ) +
  # theme parameter
  theme(
    axis.title.x.bottom = element_blank(), 
    axis.title.y.left = element_blank(),
    axis.text.x.bottom = element_blank(),
    axis.text.y = element_text(family="Arial", face="plain", size=6, hjust=1, vjust=0.5, angle=0, lineheight=2, color="black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x.bottom = unit(0, "pt"),
    axis.ticks.length.y.left = unit(2, "pt"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill=NA, color="black", linewidth=0.25, linetype="solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.frame = element_blank(),
    legend.key.size = unit(20, "pt"),
    legend.key.height = unit(20, "pt"),
    legend.key.width = unit(20, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.position = "top",                                                                                                                 
    legend.title = element_blank(),
    legend.title.position = "top",
    legend.text = element_text(family="Arial", face="plain", size=6, hjust=0, vjust=0.5, angle=0, lineheight=1, color="black"),
    legend.text.position = "right", 
    strip.background = element_blank(),
    strip.clip = "off",
    strip.text = element_text(family="Arial", face="italic", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black"),
    # strip.switch.pad.wrap = unit(8, "pt"),
  )

print(top5_species_box_plot)

ggsave("taxonomy/39_top5_species_box_plot.svg", plot = top5_species_box_plot, width = 14, height = 16, units = "cm", dpi = 300)
# ggsave("taxonomy/39_top5_species_box_plot.tiff", plot = top5_species_box_plot, width = 14, height = 16, units = "cm", dpi = 300)


### Ternary plot

# load genus-level taxonomy abundance table
genus <- read.table("taxonomy/40_ternary_plot_data.txt", sep = "\t", header = TRUE, row.names = 1) %>%  dplyr::select(-25) %>% t() %>% as.data.frame()
genus$group <- c(rep('Peri', 12), rep("Esin", 6), rep("Bger", 6))
regulation <- read.table("taxonomy/40_ternary_plot_data.txt", sep = "\t", header = TRUE, row.names = 1) %>% as.data.frame() %>%  dplyr::select(25) 

genus_group_mean <- 
  genus %>% 
  group_by(group) %>% 
  summarise(across(everything(), mean, na.rm = TRUE)) %>% 
  as.data.frame()

rownames(genus_group_mean) <- 
  genus_group_mean$group 

genus_group_mean <- 
  genus_group_mean %>%  
  dplyr::select(-1) %>% 
  t() %>% 
  as.data.frame()

genus_group_mean$regulation <- regulation$regulation

ternary_data <- genus_group_mean[, c("Peri", "Esin", "Bger")]
ternary_data$taxonomy <- rownames(ternary_data)

ternary_data <- ternary_data %>%
  separate(taxonomy, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
           sep = ";", 
           fill = "right")

ternary_data$regulation <- regulation$regulation

ternary_data$total_abundance <- rowSums(ternary_data[, 1:3])
ternary_data$Peri_trans <- ternary_data$Peri / rowSums(ternary_data[, 1:3])
ternary_data$Bger_trans <- ternary_data$Bger / rowSums(ternary_data[, 1:3])
ternary_data$Esin_trans <- ternary_data$Esin / rowSums(ternary_data[, 1:3])

ternary_data$regulation <- factor(ternary_data$regulation, levels = c("Esin", "Peri", "Bger", "no-sig"))

# set color
Peri_rgb <- (col2rgb(Pame_color) /255 + col2rgb(Pful_color) /255 ) / 2
Peri_color <-  monochromeR::generate_palette(rgb(Peri_rgb[1,], Peri_rgb[2,], Peri_rgb[3,]), modification = "go_lighter",  n_colours = 10, view_palette = TRUE)

# Calculate the minimum and maximum values
min_value <- min(ternary_data$total_abundance)
max_value <- max(ternary_data$total_abundance)

hightlighr_genus <- 
  c(g_Frigididesulfovibrio = "Frigididesulfovibrio",
    g_WRHT01 = "WRHT01",
    g_FLUQ01 = "FLUQ01",
    g_JAJBTS01 = "JAJBTS01",
    g_Bacteroides_H = "Bacteroides_H",
    g_Dysgonomonas = "Dysgonomonas",
    g_Tannerella = "Tannerella"
  )

ternary_data$label <- hightlighr_genus[ternary_data$Genus]

### ggtern ###
library(ggtern)

TernaryPlot <- 
  ggtern(data = ternary_data, 
         aes(x = Bger_trans, y = Esin_trans, z = Peri_trans)
  ) +
  geom_point(aes(color = regulation, fill = regulation, size = total_abundance), 
             shape = 19,
             alpha = 1,
             show.legend = T
  ) +
  geom_text(aes(label = label),
            color = "black",
            size = 2,
            fontface ="italic",
            vjust = 0.5,
            hjust = 0.5
  ) +
  scale_size_continuous(range = c(1, 3), 
                        name = "Abundance", 
                        breaks = seq(min(ternary_data$total_abundance), max(ternary_data$total_abundance), length.out = 4),
                        labels = c("0.1%", "1%", "10%", "33%")
  ) + 
  labs(# title = "Genus Abundance Ternary Plot",
    x = "Bger",
    y = "Esin",
    z = "Peri",
  ) + 
  scale_color_manual(values = c("Esin" = Esin_color[1], "Peri" = Peri_color[1], "Bger" = Bger_color[1], "no-sig" = "NA")) +
  scale_fill_manual(values = c("Esin" = Esin_color[1], "Peri" = Peri_color[1], "Bger" = Bger_color[1], "no-sig" = "NA")) +
  guides(color = guide_legend(override.aes = list(size = 2), nrow = 4),
         size = guide_legend(nrow = 5)
  ) +
  theme(
    # tern axis parameters
    tern.axis.clockwise = TRUE,                                                                                                                  # 轴线是否按照顺时针
    tern.axis.line.ontop = FALSE,                                                                                                                # 轴线是否总是在最上层显示
    tern.axis.line.T = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),                                    # top 轴线条设置
    tern.axis.line.R = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),                                    # right 轴线条设置
    tern.axis.line.L = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),                                    # left 轴线条设置
    tern.axis.text.show = TRUE,                                                                                                                  # 是否显示轴上文字
    tern.axis.text.T = element_text(family="Arial", face="bold", size=6, hjust=0.5, vjust=0, angle=0, lineheight=2, color="black"),              # top 轴文字设置
    tern.axis.text.R = element_text(family="Arial", face="bold", size=6, hjust=0.5, vjust=0, angle=0, lineheight=2, color="black"),              # right 轴文字设置
    tern.axis.text.L = element_text(family="Arial", face="bold", size=6, hjust=0.5, vjust=1, angle=0, lineheight=2, color="black"),              # left 轴文字设置
    tern.axis.title.show = TRUE,                                                                                                                 # 是否显示轴标题
    tern.axis.title.T = element_text(family="Arial", face="bold", size=6, hjust=0.5, vjust=0, angle=0, lineheight=2, color=Esin_color[1]),       # top 轴标题
    tern.axis.title.R = element_text(family="Arial", face="bold", size=6, hjust=0.5, vjust=0, angle=0, lineheight=2, color=Peri_color[1]),       # right 轴标题
    tern.axis.title.L = element_text(family="Arial", face="bold", size=6, hjust=0.5, vjust=0, angle=0, lineheight=2, color=Bger_color[1]),       # left 轴标题
    # tern axis ticks parameters
    tern.axis.ticks.outside = TRUE,                                                                                                              # 轴刻度是否朝外 
    tern.axis.ticks.primary.show = TRUE,                                                                                                         # 是否显示主要刻度 
    # tern.axis.ticks.length.major = unit(2, "pt"),                                                                                              # 主刻度长度 
    tern.axis.ticks.major.T = element_line(linewidth = 0.25, linetype = "solid", color = "black"),
    tern.axis.ticks.major.R = element_line(linewidth = 0.25, linetype = "solid", color = "black"),
    tern.axis.ticks.major.L = element_line(linewidth = 0.25, linetype = "solid", color = "black"),
    tern.axis.ticks.secondary.show = FALSE,                                                                                                      # 是否显示次级刻度 
    # tern.axis.ticks.length.minor = unit(0, "pt"),                                                                                              # 次级刻度长度
    tern.axis.ticks.minor.T = element_line(linewidth = 0.1, linetype = "solid", color = "black"),
    tern.axis.ticks.minor.R = element_line(linewidth = 0.1, linetype = "solid", color = "black"),
    tern.axis.ticks.minor.L = element_line(linewidth = 0.1, linetype = "solid", color = "black"),
    # tern grid parameters
    tern.plot.background = element_blank(),                                                                                                      # 图背景
    tern.plot.latex = FALSE,                                                                                                                     # 是否将文字转成letex公式
    panel.background = element_blank(),                                                                                                          # 面板背景
    tern.panel.grid.ontop = FALSE,                                                                                                               # 面板网格线是否在最上层显示
    tern.panel.grid.major.show = TRUE,                                                                                                           # 是否显示主刻度网格                                  
    tern.panel.grid.major.T = element_line(linewidth = 0.25, linetype = "dashed", color = Esin_color[1]),
    tern.panel.grid.major.R = element_line(linewidth = 0.25, linetype = "dashed", color = Peri_color[1]),
    tern.panel.grid.major.L = element_line(linewidth = 0.25, linetype = "dashed", color = Bger_color[1]),
    tern.panel.grid.minor.show = FALSE,                                                                                                          # 是否显示次级刻度网格
    tern.panel.grid.minor.T = element_line(linewidth = 0.25, linetype = "dashed", color = Esin_color[1]),
    tern.panel.grid.minor.R = element_line(linewidth = 0.25, linetype = "dashed", color = Peri_color[1]),
    tern.panel.grid.minor.L = element_line(linewidth = 0.25, linetype = "dashed", color = Bger_color[1]),
    # tern axis arrow parameters 
    tern.axis.arrow.show = TRUE,                                                                                                                 # 是否显示轴上箭头
    tern.axis.arrow.sep = 0.1,                                                                                                                   # 轴箭头与轴的距离
    tern.axis.arrow.start = 0.35,                                                                                                                # 箭头起始位置
    tern.axis.arrow.finish = 0.65,                                                                                                               # 箭头终止位置
    tern.axis.arrow.T = element_line(linewidth = 0.25, linetype = "solid", color = Esin_color[1]),                                               # top 轴箭头设置
    tern.axis.arrow.R = element_line(linewidth = 0.25, linetype = "solid", color = Peri_color[1]),                                               # right 轴箭头设置
    tern.axis.arrow.L = element_line(linewidth = 0.25, linetype = "solid", color = Bger_color[1]),                                               # left 轴箭头设置
    tern.axis.arrow.text.T = element_text(family="Arial", face="bold", size=6, hjust=0.5, vjust=0, angle=0, lineheight=2, color=Esin_color[1]),  # top 轴箭头标签
    tern.axis.arrow.text.R = element_text(family="Arial", face="bold", size=6, hjust=0.5, vjust=0, angle=0, lineheight=2, color=Peri_color[1]),  # right 轴箭头标签
    tern.axis.arrow.text.L = element_text(family="Arial", face="bold", size=6, hjust=0.5, vjust=1, angle=0, lineheight=2, color=Bger_color[1]),  # left 轴箭头标签
    # legend parameters
    legend.box.margin = margin(2, 2, 2, 2, "pt"),                                                                                                # 图例距离整个图例绘制区四周的边际
    legend.key = element_blank(),                                                                                                                # 图例键的背景设置
    legend.background = element_blank(),                                                                                                         # 图例背景设置
    legend.box.background = element_blank(),
    legend.frame = element_blank(),                                                                                                              # 用于调整整个图例区域的边框样式, 适合对图例的整体外观进行定义
    legend.key.height = unit(10, "pt"),                                                                                                          # 图例键的高度
    legend.key.width = unit(10, "pt"),                                                                                                           # 图例键的宽度
    legend.key.spacing.y = unit(5, "pt"),                                                                                                        # 图例键Y方向上的间距
    legend.key.spacing.x = unit(5, "pt"),                                                                                                        # 图例键X方向上的间距
    legend.direction = "vertical",                                                                                                               # 图例的排列方式, 水平或者垂直
    legend.box = "vertical",                                                                                                                     # 图例格子的排列方式, 水平或者垂直                                                                                                              # 图例的位置, "top" "bottom" "left" "right" "inside" "none"                                                                                                                   
    legend.title = element_text(family="Arial", face="bold", size=6, hjust=0, vjust=0.5, angle=0, lineheight=1, color="black"),                  # 图例标题字体设置                                                                                                            # 图例标题位置, "top" "bottom" "left" "right" "none" 
    legend.text = element_text(family="Arial", face="bold", size=6, hjust=0, vjust=0.5, angle=0, lineheight=1, color="black"),                   # 图例文字字体设置
  )

TernaryPlot

# Native ggsave export will lose some elements
ggtern::ggsave("taxonomy/41_ternaryplot_genus.svg", plot = TernaryPlot, width = 8, height = 10, units = "cm", dpi = 300, device = svg)
# ggtern::ggsave("taxonomy/41_ternaryplot_genus.tiff", plot = TernaryPlot, width = 8, height = 10, units = "cm", dpi = 300)





