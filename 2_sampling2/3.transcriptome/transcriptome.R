#### 
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
library(monochromeR)
library(paletteer)
library(reshape2)
library(ggvenn)
library(ggchicklet)
library(ggcorrplot)
library(mediation)
library(doParallel)
library(future.apply)
library(Hmisc)
library(aplot)
library(patchwork)
plan(multisession)
library(preprocessCore)

# load Arial font
font_import(pattern = "arial")
loadfonts()
windowsFonts()


wild_color <- "#c0c000"
lab_color  <- "#fb7d80"
fam_color  <- "#2BAA92FF"
res_color  <- "#721B3EFF"
hos_color  <- "#16317DFF"
dwelling_color <- "#751C6DFF" 
human_color <- dwelling_color




# load transcriptome data
transcriptome_data <- read.table("2_downstream_analysis/1_transcript_tpm.tsv", header = TRUE, sep = "\t", row.names = 1)

keep <- rowSums(transcriptome_data < 2) <= 25
transcriptome_filt_data <- transcriptome_data[keep, ] %>% {. + 1} %>%  log2() %>% t()

write.table(transcriptome_filt_data %>% t() %>% as.data.frame() %>% rownames_to_column(var = "gene") , "2_downstream_analysis/1_transcript_tpm_filt.tsv", quote = F, sep = "\t", row.names = F)

transcriptome_data_long <- 
  transcriptome_filt_data %>% 
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>% 
  pivot_longer(-c(sample), names_to = "gene", values_to = "TPM") 

ggplot(transcriptome_data_long, aes(x = sample, y = TPM)) + 
  geom_boxplot(outliers = F)


# load group file
group <- read.table("2_downstream_analysis/2_group.txt", sep = '\t', header =TRUE)
# Merge the group information into the metabolite table
transcriptome_data_grouped <- merge(group, transcriptome_filt_data, by.x = "sample", by.y = "row.names")


# Calculate the mean value of the group
library(dplyr)
transcriptome_mean_abundance <- 
  transcriptome_data_grouped %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "group") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene")

write.table(transcriptome_mean_abundance, "2_downstream_analysis/3_transcriptome_group_mean_abundance.tsv", row.names = FALSE, quote = FALSE)



### (PCoA)Principal Coordinates Analysis
transcriptome_bray <- vegdist(transcriptome_filt_data, method = 'bray', binary = FALSE, diag = TRUE)
transcriptome_bray_matrix <- as.matrix(transcriptome_bray)

write.table(transcriptome_bray_matrix, "2_downstream_analysis/4_transcriptome_sample_bray_crutis.txt", sep = "\t", quote = FALSE) 



# PCoA analysis based on spearman distance
transcriptome_pcoa <- cmdscale(
  transcriptome_bray,
  k = 3,
  eig = TRUE,
  add = FALSE,
  x.ret = TRUE
)

# Extract the coordinates of the first three dimensions of the PCoA dimensionality reduction results
transcriptome_pcoa_coords <- data.frame(transcriptome_pcoa$points)
transcriptome_pcoa_coords$Sample_ID <- rownames(transcriptome_pcoa_coords)
names(transcriptome_pcoa_coords)[1:3] <- paste0("PCoA", 1:3)
transcriptome_pcoa_result <- merge(transcriptome_pcoa_coords, group, by.x = "row.names", by.y = "sample")

# The "pcoa$eig" record shows the eigenvalues of the main sorting axes in the PCoA sorting results (dividing each eigenvalue by the total sum of eigenvalues gives the explanatory power of each axis)
transcriptome_pcoa_eig = sum(pmax(transcriptome_pcoa$eig, 0))
transcriptome_pcoa_eig_percent <- round(transcriptome_pcoa$eig/transcriptome_pcoa_eig*100, 3)

# Conduct a permutation multivariate (factorial) variance analysis (PERMANOVA/adonis2)
transcriptome_pcoa_permanova_result <- adonis2(transcriptome_bray  ~ group, data = group, permutations = 999, method = "bray")
transcriptome_pcoa_dune_adonis <- paste0("R2 = ", round(transcriptome_pcoa_permanova_result$R2, 3), "  P-value = ", round(transcriptome_pcoa_permanova_result$`Pr(>F)`, 3))
transcriptome_pcoa_dune_adonis


### PCoA

# Calculate the center point of each group
transcriptome_pcoa_center_points <- transcriptome_pcoa_result %>%
  group_by(group) %>%
  summarise(
    mean_wt = mean(PCoA1),
    mean_mpg = mean(PCoA2)
  )

transcriptome_pcoa_result <- transcriptome_pcoa_result %>% left_join(transcriptome_pcoa_center_points, by = "group")

transcriptome_pcoa_result$group <- factor(transcriptome_pcoa_result$group, level = c("Wild", "Lab", "Fam", "Res", "Hos"))

# Basic scatter plot
transcriptome_pcoa_p0 <- 
  ggplot(transcriptome_pcoa_result, 
         aes(x = PCoA1, y = PCoA2, fill = group, color = group)
  ) +
  geom_point(aes(color = group, shape = group, fill = group), 
             size = 1,
             show.legend = T
  ) +
  labs(x = paste("PCoA 1 (", transcriptome_pcoa_eig_percent[1], "%)", sep = ""), 
       y = paste("PCoA 2 (", transcriptome_pcoa_eig_percent[2], "%)", sep = ""), 
       tag = transcriptome_pcoa_dune_adonis
  ) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_shape_manual(values = c("Wild" = 19, "Lab" = 19, "Fam" = 19, "Res" = 19, "Hos" = 19)) +
  scale_x_continuous(position = "bottom") +
  scale_y_continuous(position = "left") +
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
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
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

transcriptome_pcoa_p0


# Add confidence ellipse + center point connection
transcriptome_pcoa_p1 <- 
  transcriptome_pcoa_p0 + 
  stat_ellipse(
    data = transcriptome_pcoa_result,
    aes(fill = group, color = group),
    geom = "polygon", 
    type = "t",
    level = 0.95, 
    linetype = 2, 
    linewidth = 0.25, 
    alpha = 0.3, 
    show.legend = F
  ) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  geom_segment(
    data = transcriptome_pcoa_result,
    aes(x = PCoA1, y = PCoA2, xend = mean_wt, yend = mean_mpg, color = group), 
    linetype = "dashed", 
    linewidth = 0.25, 
    alpha = 0.5)

transcriptome_pcoa_p1


# Use ggside to add marginal box plots
transcriptome_pcoa_final <- 
  transcriptome_pcoa_p1 +
  geom_xsideboxplot(aes(y = group, color = group, fill = group),
                    orientation = "y",
                    alpha = 0.3,
                    outliers = FALSE,
                    staplewidth = 0.6,
                    linewidth = 0.25,
                    show.legend = F
  ) +
  geom_xsidepoint(aes(y = group, color = group, fill = group), 
                  position = "jitter",
                  size = 1,
                  alpha = 0.5,
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
                  position = "jitter",
                  size = 1, 
                  alpha = 0.5,
                  show.legend = F
  ) +
  scale_xsidey_discrete() +
  scale_ysidex_discrete() +
  theme(ggside.panel.scale.x = 0.32,
        ggside.panel.scale.y = 0.32,
        legend.position = c(0.86, 0.88),
  )

print(transcriptome_pcoa_final)

ggsave("2_downstream_analysis/5_transcriptome_spearman_PCoA1~PCoA2.svg", plot = transcriptome_pcoa_final, width = 6, height = 6, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/5_transcriptome_spearman_PCoA1~PCoA2.tiff", plot = transcriptome_pcoa_final, width = 6, height = 6, units = "cm", dpi = 300)



################################
#### Hierarchical clustering ###
################################

# Hierarchical clustering based on Bray-Curtis distance
transcriptome_spearman <- as.dist(1 - cor(t(transcriptome_filt_data), method = "spearman"))
transcriptome_spearman_matrix <- transcriptome_spearman %>% as.matrix()

write.table(transcriptome_spearman_matrix, "2_downstream_analysis/6_transcriptome_sample_spearman.txt")

transcriptome_hc <- flashClust(
  transcriptome_spearman,
  method = "complete",
  members = NULL)

# Convert the hclust object to a phylo object
library(ape)
transcriptome_hc_tree <- ape::as.phylo(transcriptome_hc)
plot(transcriptome_hc_tree)

# Export Newick format tree file
write.tree(transcriptome_hc_tree, file = "2_downstream_analysis/7_transcriptome_spearman_average_tree.newick")

# Convert the hclust object to a dendrogram object
transcriptome_dend <- as.dendrogram(transcriptome_hc)

# Plot dendrogram
dev.off()
pdf("2_downstream_analysis/8_transcriptome_Hierarchical_Dendrogram.pdf")
plot(transcriptome_dend, main = "transcriptome Hierarchical Dendrogram")
dev.off()





### merge tree 
library(ape)
library(ggtree)
library(patchwork)
library(phytools)

groupInfo <- split(group$sample, group$group)

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
  scale_shape_manual(values = c("Wild" = 21, "Lab" = 22, "Fam" = 23, "Res" = 24, "Hos" = 24)) +
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color, "Human_up" = human_color)) + 
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color, "Human_up" = human_color)) + 
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
row_annotation <- group[,c(1,2)]
row.names(row_annotation) <- row_annotation[, 1]
colnames(row_annotation)[2] <- "Host_group"

col_annotation <- group[,c(1,2)]
row.names(col_annotation) <- col_annotation[, 1]
colnames(col_annotation)[2] <- "Host_group"

# "annot_colors" is a list that contains all the factor names and their corresponding colors.
# The names of each vector in the list correspond to the column names in the row_annotation/col_annotation matrix.
annot_colors <- list(Host_group = c(Wild = wild_color, Lab = lab_color, Fam = fam_color, Res = res_color, Hos = hos_color))

library(pheatmap)

min = min(1- transcriptome_spearman_matrix)
max = max(1- transcriptome_spearman_matrix)

dev.off()

pdf("2_downstream_analysis/10_transcriptome Hierarchical heatmap.pdf",
    width  = 20,
    height = 20)

pheatmap(
  as.matrix(1- transcriptome_spearman_matrix),
  color = colorRampPalette(c("lightyellow", "lightblue", "black"))(100),
  ### Global parameter
  main = "transcriptome Hierarchical heatmap",
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
  cluster_rows = transcriptome_hc,
  cluster_cols = transcriptome_hc,
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





### total_mRNA_heatmap
OG_dataframe <- 
  read.table("2_downstream_analysis/1_transcript_tpm_filt.tsv", header = T, row.names = 1) %>%
  t() %>%
  scale() %>%
  t()

row_annotation <- group[,c(1,2)]
row.names(row_annotation) <- row_annotation[, 1]
colnames(row_annotation)[2] <- "Host_group"

col_annotation <- group[,c(1,2)]
row.names(col_annotation) <- col_annotation[, 1]
colnames(col_annotation)[2] <- "Host_group"

# "annot_colors" is a list that contains all the factor names and their corresponding colors.
# The names of each vector in the list correspond to the column names in the row_annotation/col_annotation matrix.
annot_colors <- list(Host_group = c(Wild = wild_color, Lab = lab_color, Fam = fam_color, Res = res_color, Hos = hos_color))

pdf("2_downstream_analysis/11_mRNA_total_heatmap.pdf",
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
  read.table("1_OG/Total_transcriptome_KEGG_ssGSEA_pathway.tsv", header = T, row.names = 1) %>%
  t() %>%
  scale() %>%
  t()

row_annotation <- group[,c(1,2)]
row.names(row_annotation) <- row_annotation[, 1]
colnames(row_annotation)[2] <- "Host_group"

col_annotation <- group[,c(1,2)]
row.names(col_annotation) <- col_annotation[, 1]
colnames(col_annotation)[2] <- "Host_group"

# "annot_colors" is a list that contains all the factor names and their corresponding colors.
# The names of each vector in the list correspond to the column names in the row_annotation/col_annotation matrix.
annot_colors <- list(Host_group = c(Wild = wild_color, Lab = lab_color, Fam = fam_color, Res = res_color, Hos = hos_color))

pdf("2_downstream_analysis/12_mRNA_pathway_heatmap.pdf",
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







### KEGG pathway regulation

# enriched KEGG pathway heatmap
transcriptome_up_pathway <- read.table("1_OG/1_all_up_pathway.tsv", header = TRUE, row.names = 1, sep = "\t") 

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
    match(group, rev(c("wild", "lab", "fam", "res", "hos"))) ,
    number
  ) %>%
  pull(sample)


up_pathway_mratrix_long$sample <- 
  factor(
    up_pathway_mratrix_long$sample,
    levels = sample_order
  )

up_pathway_mratrix_long$pathway <- factor(up_pathway_mratrix_long$pathway, levels = rev(transcriptome_up_pathway$pathway))

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
  annotate("rect", ymin = 40.5, ymax = 50.5, xmin = 0.5, xmax = 31.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 30.5, ymax = 40.5, xmin = 31.5, xmax = 40.5, fill = NA, color = "red", linewidth = 0.1, lineend = "square") + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF", 
                       limit = c(min, max),
                       breaks = c(-2, 0, 2)
  ) +
  scale_x_discrete(position = "top") +
  labs(x = "Enriched transcriptome modules of Pame from different habitat", 
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
group$group <- factor(group$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))
group$group2 <- factor(c(rep("Dwelling", 20), rep("Lab", 10), rep("Dwelling", 10), rep("Wild", 10)), levels = c("Wild", "Lab", "Dwelling"))

group$sample <- factor(group$sample, levels = sample_order)

up_pathway_group_p <-
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

up_pathway_group_p


up_pathway_group_p2 <-
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

up_pathway_group_p2




# plot pathway pvalue annotatin
transcriptome_up_pathway$pathway <- factor(transcriptome_up_pathway$pathway, levels = rev(transcriptome_up_pathway$pathway))

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

up_pathway_pathway_info$pathway <- factor(up_pathway_pathway_info$pathway, levels = rev(transcriptome_up_pathway$pathway))

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
  aplot::insert_left(up_pathway_group_p, width = 0.01) %>% 
  aplot::insert_left(up_pathway_group_p2, width = 0.01) %>% 
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
    match(group, c("wild", "lab", "fam", "res", "hos")),
    number
  ) %>%
  pull(sample)

down_pathway_mratrix_long$sample <- 
  factor(
    down_pathway_mratrix_long$sample,
    levels = sample_order
  )

down_pathway_mratrix_long$pathway <- factor(down_pathway_mratrix_long$pathway, levels = rev(transcriptome_down_pathway$pathway))

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
  annotate("rect", ymin = 0.5, ymax = 10.5, xmin = 0.5 , xmax = 36.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 10.5, ymax = 20.5, xmin = 36.5, xmax = 40.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  annotate("rect", ymin = 20.5, ymax = 50.5, xmin = 40.5, xmax = 46.5, fill = NA, color = "blue", linewidth = 0.1, lineend = "square") + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF", 
                       limit = c(min, max),
                       breaks = c(-2,0 ,2)
  ) +
  scale_x_discrete(position = "bottom") + 
  labs(x = "Depleted transcriptome pathways of Pame from different habitat", 
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
group$group <- factor(group$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))
group$group2 <- factor(c(rep("Dwelling", 20), rep("Lab", 10), rep("Dwelling", 10), rep("Wild", 10)), levels = c("Wild", "Lab", "Dwelling"))

group$sample <- factor(group$sample, levels = sample_order)

down_pathway_group_p <-
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

down_pathway_group_p


down_pathway_group_p2 <-
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

down_pathway_group_p2

# plot pathway pvalue annotatin
transcriptome_down_pathway$pathway <- factor(transcriptome_down_pathway$pathway, levels = rev(transcriptome_down_pathway$pathway))

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
                      breaks = c(-2, -6, -10)
  ) +
  scale_color_gradient(low = paletteer::paletteer_d("MetBrewer::Benedictus")[11], 
                       high = paletteer::paletteer_d("MetBrewer::Benedictus")[9],
                       breaks = c(-2, -6, -10)
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

down_pathway_pathway_info$pathway <- factor(down_pathway_pathway_info$pathway, levels = rev(transcriptome_down_pathway$pathway))

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
  aplot::insert_left(down_pathway_group_p, width = 0.01) %>% 
  aplot::insert_left(down_pathway_group_p2, width = 0.01) %>% 
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
  transcriptome_up_pathway %>% 
  left_join(up_pathway_pathway_info, by = "pathway") %>% 
  mutate(Pathway.category = replace_na(Pathway.category, "Other")) %>% 
  group_by(regulation, Pathway.category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(regulation) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup() %>%
  add_row(regulation = "Human_up",
              Pathway.category = "AA. MET.",
              count = 0,
              percentage = 0)


up_pathway_pathway_percent$Pathway.category <- factor(up_pathway_pathway_percent$Pathway.category,
                                                      levels = rev(c("AA. MET.", "CARB. MET.", "GLYC. SYN. & MET.", "LIPI. MET.", "SIGN. TDC.", 
                                                                     "VITA. MET.", "TRNSP. & CATAB.", "TRANSL.", "SIGN. MOL. INT.", "Other")))


up_pathway_pathway_percent$regulation <- factor(up_pathway_pathway_percent$regulation,
                                                levels = c("Wild_up", "Lab_up", "Human_up"))



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
  scale_fill_manual(values = c("Wild_up" = wild_color, "Lab_up" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color, "Human_up" = human_color)) + 
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
  transcriptome_down_pathway %>% 
  left_join(down_pathway_pathway_info, by = "pathway") %>%
  mutate(Pathway.category = replace_na(Pathway.category, "Other")) %>% 
  group_by(regulation, Pathway.category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(regulation) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()


down_pathway_pathway_percent$Pathway.category <- factor(down_pathway_pathway_percent$Pathway.category, 
                                                        levels = rev(c("AA. MET.", "CARB. MET.", "GLYC. SYN. & MET.", "LIPI. MET.", "SIGN. TDC.", 
                                                                       "VITA. MET.", "TRNSP. & CATAB.", "TRANSL.", "SIGN. MOL. INT.", "Other")))


down_pathway_pathway_percent$regulation <- factor(down_pathway_pathway_percent$regulation,
                                                  levels = c("Wild_down", "Lab_down", "Human_down"))


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
  scale_fill_manual(values = c("Wild_down" = wild_color, "Lab_down" = lab_color, "Fam_down" = fam_color, "Res_down" = res_color, "Hos_down" = hos_color, "Human_down" = human_color)) + 
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
top10_enriched_pathway$regulation <- factor(top10_enriched_pathway$regulation, levels = c("Wild_up", "Lab_up", "Human_up"))
top10_enriched_pathway$pathway <- factor(top10_enriched_pathway$pathway, levels = top10_enriched_pathway$pathway)
top10_enriched_pathway_info <- read.table("1_OG/3_each_top10_up_pathway_information.tsv", header = TRUE, sep = "\t") 

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
top10_enriched_pathway_info$pathway_name <- gsub("Glycosaminoglycan biosynthesis - heparan sulfate / heparin", "Heparan sulfate/Heparin biosynthesis", top10_enriched_pathway_info$pathway_name)
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
top10_enriched_pathway_info$pathway_name <- gsub("Virion - Ebolavirus, Lyssavirus and Morbillivirus", "Virion(Ebola/Lyssa/Morbilli virus)", top10_enriched_pathway_info$pathway_name)
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
  scale_fill_manual(values = c("Wild_up" = wild_color, "Lab_up" = lab_color, "Fam_up" = lab_color, "Res" = res_color, "Hos" = hos_color, "Human_up" = human_color)) + 
  scale_color_manual(values = c("Wild_up" = wild_color, "Lab_up" = lab_color, "Fam_up" = lab_color, "Res" = res_color, "Hos" = hos_color, "Human_up" = human_color)) + 
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
                               color = c(rep(lab_color,10), rep(wild_color,10))), 
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
sample_levels <- c(
  sprintf("wild%02d", 1:10),
  sprintf("lab%02d",  1:10),
  sprintf("fam%02d",  1:10),
  sprintf("res%02d",  1:10),
  sprintf("hos%02d",  1:10)
)

up_module_z_scores <- top10_enriched_pathway[, -c(1:11)] %>% t() %>% scale() %>% t()

top10_enriched_pathway[, 12:61] = up_module_z_scores

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
  annotate("rect", xmin = 0.5, xmax = 10.5, ymin = 9.5, ymax = 19.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 10.5, xmax = 20.5, ymin = 0.5, ymax = 9.5, fill = NA, color = "red", linewidth = 0.25, lineend = "square")+
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
top10_depleted_pathway$regulation <- factor(top10_depleted_pathway$regulation, levels = c("Wild_down", "Lab_down", "Human_down"))

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
top10_depleted_pathway_info$pathway_name <- gsub("Virion - Ebolavirus, Lyssavirus and Morbillivirus", "Virion(Ebola/Lyssa/Morbilli virus)", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("Biosynthesis of unsaturated fatty acids", "Unsaturated fatty acids biosynthesis", top10_depleted_pathway_info$pathway_name)
top10_depleted_pathway_info$pathway_name <- gsub("Glycosylphosphatidylinositol \\(GPI\\)-anchor biosynthesis", "GPI-anchor biosynthesis", top10_depleted_pathway_info$pathway_name)
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
  scale_fill_manual(values = c("Wild_down" = wild_color, "Lab_down" = lab_color, "Fam_down" = lab_color, "Res_donw" = res_color, "Hos_down" = hos_color, "Human_down" = human_color)) + 
  scale_color_manual(values = c("Wild_down" = wild_color, "Lab_down" = lab_color, "Fam_down" = lab_color, "Res_down" = res_color, "Hos_down" = hos_color, "Human_down" = human_color)) + 
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
                               color = c(rep(human_color, 6), rep(lab_color, 4), rep(wild_color, 10))), 
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
sample_levels <- c(
  sprintf("wild%02d", 1:10),
  sprintf("lab%02d",  1:10),
  sprintf("fam%02d",  1:10),
  sprintf("res%02d",  1:10),
  sprintf("hos%02d",  1:10)
)

down_module_z_scores <- top10_depleted_pathway[, 12:61] %>% t() %>% scale() %>% t()

top10_depleted_pathway[, 12:61] = down_module_z_scores

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
  annotate("rect", xmin = 40.5, xmax = 50.5, ymin = 10.5, ymax = 20.5, fill = NA, color = "blue", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 30.5, xmax = 40.5, ymin = 6.5, ymax = 10.5, fill = NA, color = "blue", linewidth = 0.25, lineend = "square")+
  annotate("rect", xmin = 0.5, xmax = 30.5, ymin = 0.5, ymax = 6.5, fill = NA, color = "blue", linewidth = 0.25, lineend = "square")+
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








### transcriptome KEGG pathway importance ranking by lasso

# 1. load data
library(glmnet)
set.seed(20241102)
lasso_data <- read.table("1_OG/5_lasso_input.tsv", header = TRUE, row.names = 1, check.names = FALSE)
x <- lasso_data %>% dplyr::select(-1) %>% as.matrix()
y <- factor(lasso_data$group)

library(caret)
set.seed(123)
foldid <- createFolds(y, k = 10, list = FALSE, returnTrain = FALSE)


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
                      nfolds = 10,
                      grouped = FALSE)
  results$lambda.min[i] <- cv_fit$lambda.min
  results$cv.error[i]   <- min(cv_fit$cvm)
}


# 3. select alpha
best_alpha <- results$alpha[which.min(results$cv.error)]


# 4. Re-run the CV with the optimal alpha and select λ

cv_final <- cv.glmnet(x, y,
                      family  = "binomial",
                      alpha   = 0.5)

lambda_min  <- cv_final$lambda.min
lambda_1se  <- cv_final$lambda.1se
print(tibble(lambda.min = lambda_min, lambda.1se = lambda_1se))


plot(results$alpha, results$cv.error, type = "b", xlab = "Alpha", ylab = "CV Error", main = "CV Error vs Alpha")

dev.off()

pdf("2_downstream_analysis/23_Cross-validation.pdf",
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
  geom_vline(xintercept = 0.5, 
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


# 6. Extract the modules with non-zero importance(lambda.min)
coef_mat <- as.matrix(coef(cv_final, s = "lambda.min"))
active.idx   <- which(coef_mat[, 1] != 0)
active.genes <- rownames(coef_mat)[active.idx] %>% setdiff("(Intercept)")


# 7. export result
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
transcriptome_pathway_info$pathway_name <- gsub("Glycosaminoglycan biosynthesis - heparan sulfate / heparin ", "Heparan sulfate/Heparin biosynthesis", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("signaling", "SIG.", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("pathway", "PWY.", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("biosynthesis", "SYN.", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("degradation", "DEG.", transcriptome_pathway_info$pathway_name)
transcriptome_pathway_info$pathway_name <- gsub("metabolism", "MET.", transcriptome_pathway_info$pathway_name)

top20_lasso_transcriptome <- 
  top20_lasso_transcriptome %>% 
  left_join(transcriptome_pathway_info, by = c("pathway" = "pathway")) %>%
  mutate(label = paste0(pathway, ":", pathway_name))

top20_lasso_transcriptome$label <- factor(top20_lasso_transcriptome$label, levels = top20_lasso_transcriptome$label)

top20_lasso_transcriptome$coef <- -top20_lasso_transcriptome$coef

write.table(top20_lasso_transcriptome, "../4.association/1_multi_omic_association/6_top10_e_d_pathway_information_HostT.tsv", row.names = F, quote = F, sep = "\t")

transcriptome_all_pathway_ssGSEA_data <- read.table("1_OG/Total_transcriptome_KEGG_ssGSEA_pathway.tsv", sep = "\t", header = T, row.names = 1)

top20_lasso_transcriptome_pathway_ssGSEA <- transcriptome_all_pathway_ssGSEA_data[top20_lasso_transcriptome$pathway, ] %>% rownames_to_column(var = "pathway")

write.table(top20_lasso_transcriptome_pathway_ssGSEA, "../4.association/1_multi_omic_association/5_top10_e_d_pathway_HostT.tsv", sep = "\t", quote = F, row.names = F)



### total lasso selected feature

total_lasso_transcriptome <- coef.df
total_lasso_transcriptome$pathway <- factor(total_lasso_transcriptome$pathway, levels = rev(total_lasso_transcriptome$pathway))

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

total_lasso_transcriptome <- 
  total_lasso_transcriptome %>% 
  left_join(transcriptome_pathway_info, by = c("pathway" = "pathway")) %>%
  mutate(label = paste0(pathway, ":", pathway_name))

total_lasso_transcriptome$label <- factor(total_lasso_transcriptome$label, levels = rev(total_lasso_transcriptome$label))

write.table(total_lasso_transcriptome, "../4.association/1_multi_omic_association/21_total_e_d_pathway_information_HostT.tsv", row.names = F, quote = F, sep = "\t")

transcriptome_all_pathway_ssGSEA_data <- read.table("1_OG/Total_transcriptome_KEGG_ssGSEA_pathway.tsv", sep = "\t", header = T, row.names = 1)

total_lasso_transcriptome_pathway_ssGSEA <- transcriptome_all_pathway_ssGSEA_data[total_lasso_transcriptome$pathway, ] %>% rownames_to_column(var = "pathway")

write.table(total_lasso_transcriptome_pathway_ssGSEA, "../4.association/1_multi_omic_association/20_total_e_d_pathway_HostT.tsv", sep = "\t", quote = F, row.names = F)


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
  scale_color_manual(values = c(rep("lightblue", 10), rep("pink", 10))) +
  scale_y_discrete(position = "right") +
  scale_x_continuous(limits = c(-22, 22), breaks = c(-22, 0, 22), labels = c(-22, 0, 22)) +
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
      Human_vs_Wild < x1 & Lab_vs_Wild >= y2 ~ 1,
      Human_vs_Wild >= x1 & Human_vs_Wild < x2 & Lab_vs_Wild >= y2 ~ 2,
      Human_vs_Wild >= x2 & Lab_vs_Wild >= y2 ~ 3,
      Human_vs_Wild < x1 & Lab_vs_Wild >= y1 & Lab_vs_Wild < y2 ~ 4,
      Human_vs_Wild >= x1 & Human_vs_Wild < x2 & Lab_vs_Wild >= y1 & Lab_vs_Wild < y2 ~ 5,
      Human_vs_Wild >= x2 & Lab_vs_Wild >= y1 & Lab_vs_Wild < y2 ~ 6,
      Human_vs_Wild < x1 & Lab_vs_Wild < y1 ~ 7,
      Human_vs_Wild >= x1 & Human_vs_Wild < x2 & Lab_vs_Wild < y1 ~ 8,
      Human_vs_Wild >= x2 & Lab_vs_Wild < y1 ~ 9,
      TRUE ~ NA_integer_
    )
  )

double_volcano_data$quadrant <- as.factor(double_volcano_data$quadrant)

double_volcano_data$quadrant <- as.factor(double_volcano_data$quadrant)

double_volcano_data$size <- abs((double_volcano_data$Human_vs_Wild + double_volcano_data$Lab_vs_Wild) / 2)

quadrant_counts <- double_volcano_data %>%
  group_by(quadrant) %>%
  summarise(count = n(), .groups = 'drop')

quadrant_counts$quadrant <- factor(quadrant_counts$quadrant, levels = rev(c("3", "7", "2", "8", "4", "6", "1", "9", "5")))

quadrant_counts <- 
  quadrant_counts %>% 
  mutate(type = case_when(
    quadrant == 1 ~ "Human_enrich & Lab_deplete",
    quadrant == 2 ~ "Human_enrich & Lab_nosig",
    quadrant == 3 ~ "Human_enrich & Lab_enrich",
    quadrant == 4 ~ "Human_nosig & Lab_deplete",
    quadrant == 5 ~ "Both_nosig",
    quadrant == 6 ~ "Human_nosig & Lab_enrich",
    quadrant == 7 ~ "Human_deplete & Lab_deplete",
    quadrant == 8 ~ "Human_deplete & Lab_nosig",
    quadrant == 9 ~ "Human_deplete & Lab_enrich"
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
  scale_fill_manual(values = c("3" = "pink", "7" = "lightblue", "2" = scater_color[3], "8" =  scater_color[4], "4" = scater_color[9], "6" =  scater_color[10], "1" =  scater_color[1], "9" =  scater_color[2], "5" =  "grey80")) +
  theme_void()

quadrant_percent_pie_plot

ggsave("2_downstream_analysis/28_quadrant_percent_pie_plot.svg", plot = quadrant_percent_pie_plot, width = 1, height = 1, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/28_quadrant_percent_pie_plot.tiff", plot = quadrant_percent_pie_plot, width = 1, height = 1, units = "cm", dpi = 300)

top20_lasso_transcriptome$pathway

# highlight module
highlight_module <- c(
  ko04214 = "ko04214",
  ko04215 = "ko04215",
  ko01210 = "ko01210",
  ko00350 = "ko00350",
  ko04137 = "ko04137",
  ko00785 = "ko00785",
  ko00534 = "ko00534",
  ko00910 = "ko00910",
  ko03050 = "ko03050",
  ko03266 = "ko03266",
  ko04140 = "ko04140",
  ko04711 = "ko04711",
  ko00760 = "ko00760",
  ko04122 = "ko04122",
  ko04148 = "ko04148",
  ko00510 = "ko00510",
  ko04150 = "ko04150",
  ko00310 = "ko00310",
  ko04070 = "ko04070",
  ko04136 = "ko04136"
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
         aes(x = Lab_vs_Wild, y = Human_vs_Wild , fill = quadrant, color = quadrant)
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
  labs(x = "Lab vs Wild: [Dir] x -log10(P.adj)", 
       y = "Dwelling vs Wild: [Dir] x -log10(P.adj)"
  ) +
  scale_color_manual(values = c("3" = "pink", "7" = "lightblue", "2" = scater_color[3], "8" =  scater_color[4], "4" = scater_color[9], "6" =  scater_color[10], "1" =  scater_color[1], "9" =  scater_color[2], "5" =  "grey80")) +
  scale_fill_manual(values = c("3" = "pink", "7" = "lightblue", "2" = scater_color[3], "8" =  scater_color[4], "4" = scater_color[9], "6" =  scater_color[10], "1" =  scater_color[1], "9" =  scater_color[2], "5" =  "grey80")) +
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

ggsave("2_downstream_analysis/29_double_volcano_plot.svg", plot = double_volcano_plot, width = 6, height = 6, units = "cm", dpi = 300)
ggsave("2_downstream_analysis/29_double_volcano_plot.tiff", plot = double_volcano_plot, width = 6, height = 6, units = "cm", dpi = 300)














# set color 
Wild_venn_color <- generate_palette(wild_color, modification = "go_lighter", n_colours = 3, view_palette = TRUE, view_labels = FALSE)
Lab_venn_color <- generate_palette(lab_color, modification = "go_lighter", n_colours = 3, view_palette = TRUE, view_labels = FALSE)
Human_venn_color <- generate_palette(human_color, modification = "go_lighter", n_colours = 3, view_palette = TRUE, view_labels = FALSE)


# Wild enrich
Wild_venn_dat_up  <- read.delim("1_OG/Wild_transcriptome_KEGG_ssGSEA_pathway_each_up.tsv")
Wild_venn_up_list <- as.list(Wild_venn_dat_up)


Wild_venn_up <- 
  ggvenn(
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
Lab_venn_dat_up  <- read.delim("1_OG/Lab_transcriptome_KEGG_ssGSEA_pathway_each_up.tsv")
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
Human_venn_dat_up  <- read.delim("1_OG/Human_transcriptome_KEGG_ssGSEA_pathway_each_up.tsv")
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

ggsave("2_downstream_analysis/30_venn_up.svg", plot = venn_up, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/30_venn_up.tiff", plot = venn_up, width = 10, height = 5, units = "cm", dpi = 300)




# Wild deplete
Wild_venn_dat_down  <- read.delim("1_OG/Wild_transcriptome_KEGG_ssGSEA_pathway_each_down.tsv")
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
Lab_venn_dat_down  <- read.delim("1_OG/Lab_transcriptome_KEGG_ssGSEA_pathway_each_down.tsv")
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
Human_venn_dat_down  <- read.delim("1_OG/Human_transcriptome_KEGG_ssGSEA_pathway_each_down.tsv")
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

ggsave("2_downstream_analysis/31_venn_down.svg", plot = venn_down, width = 10, height = 5, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/31_venn_down.tiff", plot = venn_down, width = 10, height = 5, units = "cm", dpi = 300)




# plot regulation statistics bar plot
regulation <- read.table("1_OG/regulation_stat.txt",header = TRUE)
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

ggsave("2_downstream_analysis/32_regulation_stat.svg", plot = regulation_stat, width = 10, height = 10, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/32_regulation_stat.tiff", plot = regulation_stat, width = 10, height = 10, units = "cm", dpi = 300)


# merge plot
merge_stat_plot <- venn_up / regulation_stat / venn_down + plot_layout(nrow = 3)
print(merge_stat_plot)


ggsave("2_downstream_analysis/33_merge_stat_plot.svg", plot = merge_stat_plot, width = 10, height = 10, units = "cm", dpi = 300)
# ggsave("2_downstream_analysis/33_merge_stat_plot.tiff", plot = merge_stat_plot, width = 10, height = 10, units = "cm", dpi = 300)



































