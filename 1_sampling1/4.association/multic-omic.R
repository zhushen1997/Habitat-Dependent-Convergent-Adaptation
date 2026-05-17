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
library(MetaNet)
library(igraph)
library(dplyr)
library(pcutils)
library(ggfun)

# load Arial font
font_import(pattern = "arial")
loadfonts()
windowsFonts()

# Set a color for each group
wild_color <- "#c0c000"
lab_color  <- "#fb7d80"
fam_color  <- "#2BAA92FF"
res_color  <- "#721B3EFF"
hos_color  <- "#16317DFF"
dwelling_color <- "#751C6DFF"
human_color <- dwelling_color


### ================================================================================================================================================================================== ###

# top10 enrich+deplete importance module of multi-omic
metagenome <- read.table("1_multi_omic_association/1_top10_e_d_metagenome_module.tsv", row.names = 1, header = TRUE, sep = "\t")
colnames(metagenome) <- gsub("^Metagenome_", "", colnames(metagenome))
metabolome <- read.table("1_multi_omic_association/3_top10_e_d_metabolome_module.tsv", row.names = 1, header = TRUE, sep = "\t")
transcriptome <- read.table("1_multi_omic_association/5_top10_e_d_pathway_HostT.tsv", row.names = 1, header = TRUE, sep = "\t")

metagenome_info <- read.table("1_multi_omic_association/2_top10_e_d_module_information_MetaG.tsv", header = TRUE, sep = "\t")
metabolome_info <- read.table("1_multi_omic_association/4_top10_e_d_module_information_MetaB.tsv",, header = TRUE, sep = "\t")
transcriptome_info <- read.table("1_multi_omic_association/6_top10_e_d_pathway_information_HostT.tsv", header = TRUE, sep = "\t")

common_samples <- Reduce(intersect, list(
  colnames(metagenome),
  colnames(metabolome),
  colnames(transcriptome)
))

metagenome <- t(metagenome[, common_samples]) %>% as.data.frame()
metabolome <- t(metabolome[, common_samples]) %>% as.data.frame()
transcriptome <- t(transcriptome[, common_samples]) %>% as.data.frame()

MetaG_MetaB_combined_matrix <- as.matrix(cbind(metabolome, metagenome))
MetaB_HostT_combined_matrix <- as.matrix(cbind(metabolome, transcriptome))

# Use rcorr to calculate the correlation matrix and the p-value matrix
MetaG_MetaB_Correlation <- rcorr(MetaG_MetaB_combined_matrix, type = "spearman")
MetaB_HostT_Correlation <- rcorr(MetaB_HostT_combined_matrix, type = "spearman")

# Extract the correlation matrix and the p-value matrix
MetaG_MetaB_cor_matrix <- MetaG_MetaB_Correlation$r[1:21, 22:41]
MetaG_MetaB_pvalue_matrix <- MetaG_MetaB_Correlation$P[1:21, 22:41]

MetaB_HostT_cor_matrix <- MetaB_HostT_Correlation$r[1:21, 22:41]
MetaB_HostT_pvalue_matrix <- MetaB_HostT_Correlation$P[1:21, 22:41]

# BH(Benjamini-Hochberg),BY(Benjamini-Yekutieli),Bonferroni
MetaG_MetaB_pvalue_matrix <- 
  apply(MetaG_MetaB_pvalue_matrix, 2, function(x) {
    p.adjust(x, method = "BH")
  })

MetaB_HostT_pvalue_matrix <-
  apply(MetaB_HostT_pvalue_matrix, 2, function(x) {
    p.adjust(x, method = "BH")
  })

# wide to long
MetaG_MetaB_cor_matrix_long <- MetaG_MetaB_cor_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "MetaG", values_to = "Correlation")

MetaG_MetaB_pvalue_matrix_long <- MetaG_MetaB_pvalue_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "MetaG", values_to = "P_value")

MetaB_HostT_cor_matrix_long <- MetaB_HostT_cor_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "HostT", values_to = "Correlation")

MetaB_HostT_pvalue_matrix_long <- MetaB_HostT_pvalue_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "HostT", values_to = "P_value")

# merge
MetaG_MetaB_data <- MetaG_MetaB_cor_matrix_long %>%
  left_join(MetaG_MetaB_pvalue_matrix_long, by = c("MetaB", "MetaG"))

MetaB_HostT_data <- MetaB_HostT_cor_matrix_long %>%
  left_join(MetaB_HostT_pvalue_matrix_long, by = c("MetaB", "HostT"))

MetaG_MetaB_data$MetaG <- factor(MetaG_MetaB_data$MetaG, levels = rev(colnames(metagenome)))
MetaG_MetaB_data$MetaB <- factor(MetaG_MetaB_data$MetaB, levels = colnames(metabolome))

MetaB_HostT_data$HostT <- factor(MetaB_HostT_data$HostT, levels = rev(colnames(transcriptome)))
MetaB_HostT_data$MetaB <- factor(MetaB_HostT_data$MetaB, levels = colnames(metabolome))

metagenome_info$label <- gsub("Neomycin/Kanamycin/Gentamicin", "Neo/Kana/Genta mycin", metagenome_info$label)
metagenome_info$label <- paste0("MetaG-", metagenome_info$label)
metagenome_info$label <- factor(metagenome_info$label, levels = rev(metagenome_info$label))

metabolome_info$merge.name <- paste0("MetaB-", metabolome_info$label)
metabolome_info$merge.name <- gsub(" and ", " / ", metabolome_info$merge.name)
metabolome_info$merge.name <- gsub("derivatives", "Derivatives", metabolome_info$merge.name)
metabolome_info$merge.name <- gsub("lipid-like", "Lipid-like", metabolome_info$merge.name)
metabolome_info$merge.name <- gsub(" and ", " / ", metabolome_info$merge.name)
metabolome_info$merge.name <- gsub(" and ", " / ", metabolome_info$merge.name)

MetaG_MetaB_data <- MetaG_MetaB_data %>%
  mutate("-log10(P_value)" = -log10(P_value)) %>%
  mutate(Significance = case_when(
    (Correlation >= 0.45 | Correlation <= -0.45) & P_value <= 0.05 ~ "*",
    TRUE ~ NA
  )) %>%
  mutate(Fill = case_when(
    (Correlation >= 0.45 | Correlation <= -0.45) & P_value < 0.05 ~ Correlation,
    TRUE ~ NA
  )) %>%
  mutate(logtrans = case_when(
    Fill != 0 ~ -log10(P_value),
    TRUE ~ NA
  ))

MetaB_HostT_data <- MetaB_HostT_data %>%
  mutate("-log10(P_value)" = -log10(P_value)) %>%
  mutate(Significance = case_when(
    (Correlation >= 0.45 | Correlation <= -0.45) & P_value <= 0.05~ "*",
    TRUE ~ NA
  )) %>%
  mutate(Fill = case_when(
    (Correlation >= 0.45 | Correlation <= -0.45) & P_value < 0.05 ~ Correlation,
    TRUE ~ NA
  )) %>%
  mutate(logtrans = case_when(
    Fill != 0 ~ -log10(P_value),
    TRUE ~ NA
  ))

### MetaG_MetaB_heatmap
library(ggforce)

MetaG_MetaB_heatmap_p0 <- 
  ggplot(MetaG_MetaB_data, 
         aes(x = MetaB, y = MetaG, fill = Fill)
  ) +
  geom_shape(radius = unit(1.5, "pt"),
             expand = unit(3, "pt"),
             color = "white",
             linewidth = 0,
             show.legend = T
  ) +
  scale_fill_gradientn(colors = c(c("#ADD8E6", "#BDDFEB", "#CDE7F0", "#DEEFF5", "#EEF7FA"), "white", rev(c("#D9636C", "#E08289", "#E8A1A6", "#EFC0C4", "#F7DFE1"))),
                       name = "Correlation", 
                       na.value = NA,
                       breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  scale_x_discrete(position = "top",
                   labels = rev(paste0("MetaB-", metabolome_info$module))
  ) + 
  scale_y_discrete(position = "right",
                   labels = rev(metagenome_info$label)
  ) + 
  theme_void() + 
  theme(axis.text.x = element_text(family = "Arial", size = 6, hjust = 0, vjust = 0, angle = 45, lineheight = 1, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(5, "pt"),
        legend.key.spacing.y = unit(2, "pt"),
        legend.key.spacing.x = unit(2, "pt"),
        legend.text = element_text(family = "Arial", face = "plain", size = 4, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.title = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.position = "none",
        legend.ticks.length = unit(1, "pt"),
        legend.ticks = element_line(linewidth = 0.25, linetype = "solid", color = "white"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linetype = 1, linewidth = 0.25)
  )

MetaG_MetaB_heatmap_p0

# MetaG legendry
metagenome_info <- 
  metagenome_info %>% 
  mutate(type = case_when(
    metagenome_info$category_C == "Branched-chain amino acid metabolism" ~ "BCAA MET.",
    metagenome_info$category_C == "Arginine and proline metabolism" ~ "Arg/Pro MET.",
    metagenome_info$category_C == "Cofactor and vitamin metabolism" ~ "Vitamin/Cofactor MET.",
    metagenome_info$category_C == "Xenobiotics biodegradation and metabolism" ~ "Xenobiotics DEG. & MET.",
    metagenome_info$category_C == "ATP synthesis" ~ "ATP SYN.",
    TRUE ~ "Other"
  ))

metagenome_info$type <- factor(metagenome_info$type, levels = c("BCAA MET.", "Arg/Pro MET.", "Vitamin/Cofactor MET.", "Xenobiotics DEG. & MET.", "ATP SYN.", "Other"))
metagenome_info$module <- factor(metagenome_info$module, levels = rev(metagenome_info$module))

metag_color <- c("BCAA MET." = paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1],
                 "Arg/Pro MET." = paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[2],
                 "Vitamin/Cofactor MET." = paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[3],
                 "Xenobiotics DEG. & MET." = paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[4],
                 "ATP SYN." = paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[5],
                 "Other" = "grey")

MetaG_legendry <-
  ggplot(metagenome_info) + 
  geom_tile(aes(x = 0, y = module, fill = type), 
            color = "white", 
            height = 1,
            linewidth = 0.5,
            show.legend = T
  ) + 
  scale_fill_manual(values = metag_color, drop = FALSE) +
  scale_color_manual(values = metag_color, drop = FALSE) +
  theme_void() +
  theme(legend.position = "none",
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(5, "pt"),
        legend.key.spacing.y = unit(2, "pt"),
        legend.key.spacing.x = unit(2, "pt"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.text = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.title = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linetype = 1, linewidth = 0.25),
        panel.background = element_blank()
  )

MetaG_legendry

MetaG_MetaB_heatmap <- MetaG_MetaB_heatmap_p0  %>% insert_left(MetaG_legendry, width = 0.06)
MetaG_MetaB_heatmap0 <- MetaG_MetaB_heatmap_p0 + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
MetaG_MetaB_heatmap1 <- MetaG_MetaB_heatmap_p0 + theme(axis.text.x = element_blank())
MetaG_MetaB_heatmap2 <- MetaG_MetaB_heatmap_p0 + theme(axis.text.y = element_blank())

MetaG_MetaB_heatmap
MetaG_MetaB_heatmap0
MetaG_MetaB_heatmap1
MetaG_MetaB_heatmap2

ggsave("1_multi_omic_association/7_MetaG_MetaB_heatmap0.svg", plot = MetaG_MetaB_heatmap0, width = 5.5, height = 5.5, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/7_MetaG_MetaB_heatmap0.tiff", plot = MetaG_MetaB_heatmap0, width = 5.5, height = 5.5, units = "cm", dpi = 300)

ggsave("1_multi_omic_association/8_MetaG_MetaB_heatmap1.svg", plot = MetaG_MetaB_heatmap1, width = 10, height = 5.5, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/8_MetaG_MetaB_heatmap1.tiff", plot = MetaG_MetaB_heatmap1, width = 10, height = 5.5, units = "cm", dpi = 300)

ggsave("1_multi_omic_association/9_MetaG_MetaB_heatmap2.svg", plot = MetaG_MetaB_heatmap2, width = 5.5, height = 10, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/9_MetaG_MetaB_heatmap2.tiff", plot = MetaG_MetaB_heatmap2, width = 5.5, height = 10, units = "cm", dpi = 300)

ggsave("1_multi_omic_association/10_MetaG_legendry.svg", plot = MetaG_legendry, width = 0.30, height = 5.5, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/10_MetaG_legendry.tiff", plot = MetaG_legendry, width = 0.30, height = 5.5, units = "cm", dpi = 300)


### MetaB_HostT_heatmap
MetaB_HostT_heatmap_p0 <- 
  ggplot(MetaB_HostT_data, 
         aes(x = MetaB, y = HostT, fill = Fill)
  ) +
  geom_shape(radius = unit(1.5, "pt"),
             expand = unit(3, "pt"),
             color = "white",
             linewidth = 0,
             show.legend = T
  ) +
  scale_fill_gradientn(colors = c(c("#ADD8E6", "#BDDFEB", "#CDE7F0", "#DEEFF5", "#EEF7FA"), "white", rev(c("#D9636C", "#E08289", "#E8A1A6", "#EFC0C4", "#F7DFE1"))),
                       name = "Correlation", 
                       na.value = NA,
                       breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  scale_x_discrete(position = "top",
                   labels = rev(paste0("MetaB-", metabolome_info$module))
  ) + 
  scale_y_discrete(position = "right",
                   labels = paste0("HostT-", rev(transcriptome_info$label))
  ) + 
  theme_void() + 
  theme(axis.text.x = element_text(family = "Arial", size = 6, hjust = 0, vjust = 0, angle = 45, lineheight = 1, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(5, "pt"),
        legend.key.spacing.y = unit(2, "pt"),
        legend.key.spacing.x = unit(2, "pt"),
        legend.text = element_text(family = "Arial", face = "plain", size = 4, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.title = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.position = "none",
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.ticks.length = unit(1, "pt"),
        legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white"),
        panel.border = element_rect(fill = NA, color = "black", linetype = 1, linewidth = 0.25),
  )

MetaB_HostT_heatmap_p0


### legendry

# HostT_legendry
transcriptome_info <- transcriptome_info %>% 
  mutate(type = case_when(
    transcriptome_info$Subcategory == "Signal transduction" ~ "Signal transduction",
    transcriptome_info$Subcategory == "Folding, sorting and degradation" ~ "Folding/Sorting/Degradation",
    transcriptome_info$Subcategory == "Amino acid metabolism" ~ "Amino acid metabolism",
    transcriptome_info$Subcategory == "Carbohydrate metabolism" ~ "Carbohydrate metabolism",
    transcriptome_info$Subcategory == "Signaling molecules and interaction" ~ "Signaling molecules and interaction",
    TRUE ~ "Other"
  ))

transcriptome_info$pathway <- factor(transcriptome_info$pathway, levels = rev(transcriptome_info$pathway))
transcriptome_info$type <- factor(transcriptome_info$type, levels = c("Signal transduction", "Folding/Sorting/Degradation", "Amino acid metabolism", "Carbohydrate metabolism", "Signaling molecules and interaction", "Other"))

hostt_color <- c("Signal transduction" = paletteer::paletteer_d("wesanderson::AsteroidCity1")[1],
                "Folding/Sorting/Degradation" = paletteer::paletteer_d("wesanderson::AsteroidCity1")[2],
                "Amino acid metabolism" = paletteer::paletteer_d("wesanderson::AsteroidCity1")[3],
                "Carbohydrate metabolism" = paletteer::paletteer_d("wesanderson::AsteroidCity1")[4],
                "Signaling molecules and interaction" = paletteer::paletteer_d("wesanderson::AsteroidCity1")[5],
                "Other" = "grey")

HostT_legendry <-
  ggplot(transcriptome_info) + 
  geom_tile(aes(x = 0, y = pathway, fill = type), 
            color = "white", 
            height = 1,
            linewidth = 0.5,
            show.legend = T
  ) + 
  scale_fill_manual(values = hostt_color, drop = FALSE) +
  scale_color_manual(values = hostt_color, drop = FALSE) +
  theme_void() +
  theme(legend.position = "none",
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(5, "pt"),
        legend.key.spacing.y = unit(2, "pt"),
        legend.key.spacing.x = unit(2, "pt"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.text = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.title = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linetype = 1, linewidth = 0.25),
        panel.background = element_blank()
  )

HostT_legendry

# MetaB_legendry 

metabolome_info$MetaB <- paste0("MetaB_", metabolome_info$module)
metabolome_info$MetaB <- factor(metabolome_info$MetaB, levels = metabolome_info$MetaB)
metabolome_info$module <- factor(metabolome_info$module, levels = metabolome_info$module)

metabolome_info$Module.category <- factor(metabolome_info$Module.category, 
                                          levels = c("Phenylpropanoids and polyketides", "Lipids and lipid-like molecules", "Organic acids and derivatives", 
                                                     "Organoheterocyclic compounds", "Organic oxygen compounds", "Miscellaneous"))

metab_color <- c("Phenylpropanoids and polyketides" = paletteer::paletteer_d("LaCroixColoR::Apricot")[1],
              "Lipids and lipid-like molecules" = paletteer::paletteer_d("LaCroixColoR::Apricot")[2],
              "Organic acids and derivatives" = paletteer::paletteer_d("LaCroixColoR::Apricot")[3],
              "Organoheterocyclic compounds" = paletteer::paletteer_d("LaCroixColoR::Apricot")[4],
              "Organic oxygen compounds" = paletteer::paletteer_d("LaCroixColoR::Apricot")[5],
              "Miscellaneous" = "grey")

MetaB_legendry <-
  ggplot(metabolome_info) + 
  geom_tile(aes(x = rev(module), y = 0, fill = Module.category), 
            color = "white", 
            height = 1,
            linewidth = 0.5,
            show.legend = T
  ) + 
  scale_fill_manual(values = metab_color, drop = FALSE) +
  scale_color_manual(values = metab_color, drop = FALSE) +
  theme_void() +
  theme(legend.position = "none",
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(5, "pt"),
        legend.key.spacing.y = unit(2, "pt"),
        legend.key.spacing.x = unit(2, "pt"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.text = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.title = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linetype = 1, linewidth = 0.25),
        panel.background = element_blank()
  )

MetaB_legendry

MetaB_HostT_heatmap <- MetaB_HostT_heatmap_p0  %>% insert_left(HostT_legendry, width = 0.06) %>% insert_bottom(MetaB_legendry, height = 0.06)
MetaB_HostT_heatmap0 <-  MetaB_HostT_heatmap_p0 + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
MetaB_HostT_heatmap1 <- MetaB_HostT_heatmap_p0 + theme(axis.text.x = element_blank())

MetaB_HostT_heatmap
MetaB_HostT_heatmap0
MetaB_HostT_heatmap1

ggsave("1_multi_omic_association/11_MetaB_HostT_heatmap0.svg", plot = MetaB_HostT_heatmap0, width = 5.5, height = 5.5, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/11_MetaB_HostT_heatmap0.tiff", plot = MetaB_HostT_heatmap0, width = 5.5, height = 5.5, units = "cm", dpi = 300)

ggsave("1_multi_omic_association/12_MetaB_HostT_heatmap1.svg", plot = MetaB_HostT_heatmap1, width = 10, height = 5.5, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/12_MetaB_HostT_heatmap1.tiff", plot = MetaB_HostT_heatmap1, width = 10, height = 5.5, units = "cm", dpi = 300)

ggsave("1_multi_omic_association/13_HostT_legendry.svg", plot = HostT_legendry, width = 0.30, height = 5.5, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/13_HostT_legendry.tiff", plot = HostT_legendry, width = 0.30, height = 5.5, units = "cm", dpi = 300)

ggsave("1_multi_omic_association/14_MetaB_legendry.svg", plot = MetaB_legendry, width = 5.5, height = 0.30, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/14_MetaB_legendry.tiff", plot = MetaB_legendry, width = 5.5, height = 0.30, units = "cm", dpi = 300)

### save all legend
MetaG_MetaB_heatmap_p0_legend <- cowplot::get_legend(MetaG_MetaB_heatmap_p0 + theme(legend.position = "right"))
MetaG_legendry_legend <- cowplot::get_legend(MetaG_legendry + theme(legend.position = "right"))
MetaB_HostT_heatmap_p0_legend <- cowplot::get_legend(MetaB_HostT_heatmap_p0 + theme(legend.position = "right"))
MetaB_legendry_legend <- cowplot::get_legend(MetaB_legendry + theme(legend.position = "right"))
HostT_legendry_legend <- cowplot::get_legend(HostT_legendry + theme(legend.position = "right"))

MetaG_MetaB_HostT_legend <- cowplot::plot_grid(MetaG_MetaB_heatmap_p0_legend,
                                               MetaB_HostT_heatmap_p0_legend,
                                               MetaG_legendry_legend, 
                                               MetaB_legendry_legend, 
                                               HostT_legendry_legend, ncol = 5)

MetaG_MetaB_HostT_legend

ggsave("1_multi_omic_association/15_MetaG_MetaB_HostT_legend.svg", plot = MetaG_MetaB_HostT_legend, width = 20, height = 10, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/15_MetaG_MetaB_HostT_legend.tiff", plot = MetaG_MetaB_HostT_legend, width = 20, height = 10, units = "cm", dpi = 300)








### ================================================================================================================================================================================== ###

### Total feature modules of multi-omic
lasso_metagenome <- read.table("1_multi_omic_association/16_total_e_d_metagenome_module.tsv", row.names = 1, header = TRUE, sep = "\t")
colnames(lasso_metagenome) <- gsub("^Metagenome_", "", colnames(lasso_metagenome))
lasso_metabolome <- read.table("1_multi_omic_association/18_total_e_d_metabolome_module.tsv", row.names = 1, header = TRUE, sep = "\t")
lasso_transcriptome <- read.table("1_multi_omic_association/20_total_e_d_pathway_HostT.tsv", row.names = 1, header = TRUE, sep = "\t")

lasso_metagenome_info <- read.table("1_multi_omic_association/17_total_e_d_module_information_MetaG.tsv", header = TRUE, sep = "\t")
lasso_metabolome_info <- read.table("1_multi_omic_association/19_total_e_d_module_information_MetaB.tsv", header = TRUE, sep = "\t")
lasso_transcriptome_info <- read.table("1_multi_omic_association/21_total_e_d_pathway_information_HostT.tsv", header = TRUE, sep = "\t")

lasso_common_samples <- Reduce(intersect, list(
  colnames(lasso_metagenome),
  colnames(lasso_metabolome),
  colnames(lasso_transcriptome)
))

lasso_metagenome <- t(lasso_metagenome[, lasso_common_samples]) %>% as.data.frame()
lasso_metabolome <- t(lasso_metabolome[, lasso_common_samples]) %>% as.data.frame()
lasso_transcriptome <- t(lasso_transcriptome[, lasso_common_samples]) %>% as.data.frame()

lasso_MetaG_MetaB_combined_matrix <- as.matrix(cbind(lasso_metabolome, lasso_metagenome))
lasso_MetaB_HostT_combined_matrix <- as.matrix(cbind(lasso_metabolome, lasso_transcriptome))

# Use rcorr to calculate the correlation matrix and the p-value matrix
lasso_MetaG_MetaB_Correlation <- rcorr(lasso_MetaG_MetaB_combined_matrix, type = "spearman")
lasso_MetaB_HostT_Correlation <- rcorr(lasso_MetaB_HostT_combined_matrix, type = "spearman")
# lasso_MetaG_MetaB_Correlation <- rcorr(lasso_MetaG_MetaB_combined_matrix, type = "pearson")
# lasso_MetaB_HostT_Correlation <- rcorr(lasso_MetaB_HostT_combined_matrix, type = "pearson")

# Extract the correlation matrix and the p-value matrix
lasso_MetaG_MetaB_cor_matrix <- lasso_MetaG_MetaB_Correlation$r[1:57, 58:110]
lasso_MetaG_MetaB_pvalue_matrix <- lasso_MetaG_MetaB_Correlation$P[1:57, 58:110]

lasso_MetaB_HostT_cor_matrix <- lasso_MetaB_HostT_Correlation$r[1:57, 58:85]
lasso_MetaB_HostT_pvalue_matrix <- lasso_MetaB_HostT_Correlation$P[1:57, 58:85]

# BH(Benjamini-Hochberg),BY(Benjamini-Yekutieli),Bonferroni
lasso_MetaG_MetaB_pvalue_matrix <- 
  apply(lasso_MetaG_MetaB_pvalue_matrix, 2, function(x) {
    p.adjust(x, method = "BH")
  })

lasso_MetaB_HostT_pvalue_matrix <-
  apply(lasso_MetaB_HostT_pvalue_matrix, 2, function(x) {
    p.adjust(x, method = "BH")
  })

# wide to long
lasso_MetaG_MetaB_cor_matrix_long <- 
  lasso_MetaG_MetaB_cor_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "MetaG", values_to = "Correlation")

lasso_MetaG_MetaB_pvalue_matrix_long <- 
  lasso_MetaG_MetaB_pvalue_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "MetaG", values_to = "P_value")

lasso_MetaB_HostT_cor_matrix_long <- 
  lasso_MetaB_HostT_cor_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "HostT", values_to = "Correlation")

lasso_MetaB_HostT_pvalue_matrix_long <- 
  lasso_MetaB_HostT_pvalue_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "HostT", values_to = "P_value")

# merge
lasso_MetaG_MetaB_data <- 
  lasso_MetaG_MetaB_cor_matrix_long %>%
  left_join(lasso_MetaG_MetaB_pvalue_matrix_long, by = c("MetaB", "MetaG"))

lasso_MetaB_HostT_data <- 
  lasso_MetaB_HostT_cor_matrix_long %>%
  left_join(lasso_MetaB_HostT_pvalue_matrix_long, by = c("MetaB", "HostT"))

lasso_MetaG_MetaB_data$MetaG <- factor(lasso_MetaG_MetaB_data$MetaG, levels = colnames(lasso_metagenome))
lasso_MetaG_MetaB_data$MetaB <- factor(lasso_MetaG_MetaB_data$MetaB, levels = rev(colnames(lasso_metabolome)))

lasso_MetaB_HostT_data$HostT <- factor(lasso_MetaB_HostT_data$HostT, levels = colnames(lasso_transcriptome))
lasso_MetaB_HostT_data$MetaB <- factor(lasso_MetaB_HostT_data$MetaB, levels = rev(colnames(lasso_metabolome)))

lasso_metagenome_info$label <- gsub("Neomycin/Kanamycin/Gentamicin", "Neo/Kana/Genta mycin", lasso_metagenome_info$label)
lasso_metagenome_info$label <- paste0("MetaG-", lasso_metagenome_info$label)
lasso_metagenome_info$label <- factor(lasso_metagenome_info$label, levels = rev(lasso_metagenome_info$label))

lasso_metabolome_info$merge.name <- paste0("MetaB-", lasso_metabolome_info$label)
lasso_metabolome_info$merge.name <- gsub(" and ", " / ", lasso_metabolome_info$merge.name)
lasso_metabolome_info$merge.name <- gsub("derivatives", "Derivatives", lasso_metabolome_info$merge.name)
lasso_metabolome_info$merge.name <- gsub("lipid-like", "Lipid-like", lasso_metabolome_info$merge.name)
lasso_metabolome_info$merge.name <- gsub(" and ", " / ", lasso_metabolome_info$merge.name)
lasso_metabolome_info$merge.name <- gsub(" and ", " / ", lasso_metabolome_info$merge.name)

lasso_MetaG_MetaB_data <- 
  lasso_MetaG_MetaB_data %>%
  mutate("-log10(P_value)" = -log10(P_value)) %>%
  mutate(Significance = case_when(
    (Correlation >= 0.450 | Correlation <= -0.450) & P_value <= 0.05 ~ "*",
    TRUE ~ NA
  )) %>%
  mutate(Fill = case_when(
    (Correlation >= 0.450 | Correlation <= -0.450) & P_value < 0.05 ~ Correlation,
    TRUE ~ NA
  )) %>%
  mutate(logtrans = case_when(
    Fill != 0 ~ -log10(P_value),
    TRUE ~ NA
  ))

lasso_MetaB_HostT_data <- 
  lasso_MetaB_HostT_data %>%
  mutate("-log10(P_value)" = -log10(P_value)) %>%
  mutate(Significance = case_when(
    (Correlation >= 0.450 | Correlation <= -0.450) & P_value <= 0.05 ~ "*",
    TRUE ~ NA
  )) %>%
  mutate(Fill = case_when(
    (Correlation >= 0.450 | Correlation <= -0.450) & P_value < 0.05 ~ Correlation,
    TRUE ~ NA
  )) %>%
  mutate(logtrans = case_when(
    Fill != 0 ~ -log10(P_value),
    TRUE ~ NA
  ))

### MetaG_MetaB_heatmap
library(ggforce)

lasso_MetaG_MetaB_heatmap_p0 <- 
  ggplot(lasso_MetaG_MetaB_data, 
         aes(x = MetaB, y = MetaG, fill = Fill)
  ) +
  geom_shape(radius = unit(1.5, "pt"),
             expand = unit(3, "pt"),
             color = "white",
             linewidth = 0,
             show.legend = T
  ) +
  scale_fill_gradientn(colors = c(c("#ADD8E6", "#BDDFEB", "#CDE7F0", "#DEEFF5", "#EEF7FA"), "white", rev(c("#D9636C", "#E08289", "#E8A1A6", "#EFC0C4", "#F7DFE1"))),
                       name = "Correlation", 
                       na.value = NA,
                       breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  scale_x_discrete(position = "top",
                   labels = rev(paste0("MetaB-", lasso_metabolome_info$label))
  ) + 
  scale_y_discrete(position = "right",
                   labels = lasso_metagenome_info$label
  ) + 
  theme_void() + 
  theme(axis.text.x = element_text(family = "Arial", size = 6, hjust = 0, vjust = 0, angle = 45, lineheight = 1, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(5, "pt"),
        legend.key.spacing.y = unit(2, "pt"),
        legend.key.spacing.x = unit(2, "pt"),
        legend.text = element_text(family = "Arial", face = "plain", size = 4, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.title = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.position = "none",
        legend.ticks.length = unit(1, "pt"),
        legend.ticks = element_line(linewidth = 0.25, linetype = "solid", color = "white"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linetype = 1, linewidth = 0.25)
  )

lasso_MetaG_MetaB_heatmap_p0

# MetaG legendry
lasso_metagenome_info <- 
  lasso_metagenome_info %>% 
  mutate(type = case_when(
    lasso_metagenome_info$category_C == "Branched-chain amino acid metabolism" ~ "BCAA MET.",
    lasso_metagenome_info$category_C == "Arginine and proline metabolism" ~ "Arg/Pro MET.",
    lasso_metagenome_info$category_C == "Cofactor and vitamin metabolism" ~ "Vitamin/Cofactor MET.",
    lasso_metagenome_info$category_C == "Xenobiotics biodegradation and metabolism" ~ "Xenobiotics DEG. & MET.",
    lasso_metagenome_info$category_C == "ATP synthesis" ~ "ATP SYN.",
    TRUE ~ "Other"
  ))

lasso_metagenome_info$type <- factor(lasso_metagenome_info$type, levels = c("BCAA MET.", "Arg/Pro MET.", "Vitamin/Cofactor MET.", "Xenobiotics DEG. & MET.", "ATP SYN.", "Other"))
lasso_metagenome_info$module <- factor(lasso_metagenome_info$module, levels = lasso_metagenome_info$module)

lasso_MetaG_legendry <-
  ggplot(lasso_metagenome_info) + 
  geom_tile(aes(x = 0, y = module, fill = type), 
            color = "white", 
            height = 1,
            linewidth = 0.5,
            show.legend = T
  ) + 
  scale_fill_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:5], "grey")) +
  scale_color_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:5], "grey")) +
  theme_void() +
  theme(legend.position = "none",
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(5, "pt"),
        legend.key.spacing.y = unit(2, "pt"),
        legend.key.spacing.x = unit(2, "pt"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.text = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.title = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linetype = 1, linewidth = 0.25),
        panel.background = element_blank()
  )

lasso_MetaG_legendry

lasso_MetaG_MetaB_heatmap <- lasso_MetaG_MetaB_heatmap_p0  %>% insert_left(lasso_MetaG_legendry, width = 0.06)
lasso_MetaG_MetaB_heatmap0 <- lasso_MetaG_MetaB_heatmap_p0 + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
lasso_MetaG_MetaB_heatmap1 <- lasso_MetaG_MetaB_heatmap_p0 + theme(axis.text.x = element_blank())
lasso_MetaG_MetaB_heatmap2 <- lasso_MetaG_MetaB_heatmap_p0 + theme(axis.text.y = element_blank())

lasso_MetaG_MetaB_heatmap
lasso_MetaG_MetaB_heatmap0
lasso_MetaG_MetaB_heatmap1
lasso_MetaG_MetaB_heatmap2

ggsave("1_multi_omic_association/22_total_MetaG_MetaB_heatmap0.svg", plot = lasso_MetaG_MetaB_heatmap0, width = 10, height = 12, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/22_total_MetaG_MetaB_heatmap0.tiff", plot = lasso_MetaG_MetaB_heatmap0, width = 10, height = 12, units = "cm", dpi = 300)

ggsave("1_multi_omic_association/23_total_MetaG_MetaB_heatmap1.svg", plot = lasso_MetaG_MetaB_heatmap1, width = 15, height = 12, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/23_total_MetaG_MetaB_heatmap1.tiff", plot = lasso_MetaG_MetaB_heatmap1, width = 15, height = 12, units = "cm", dpi = 300)

ggsave("1_multi_omic_association/24_total_MetaG_MetaB_heatmap2.svg", plot = lasso_MetaG_MetaB_heatmap2, width = 10, height = 15, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/24_total_MetaG_MetaB_heatmap2.tiff", plot = lasso_MetaG_MetaB_heatmap2, width = 10, height = 15, units = "cm", dpi = 300)

ggsave("1_multi_omic_association/25_total_MetaG_legendry.svg", plot = lasso_MetaG_legendry, width = 0.30, height = 12, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/25_total_MetaG_legendry.tiff", plot = lasso_MetaG_legendry, width = 0.30, height = 12, units = "cm", dpi = 300)

### MetaB_HostT_heatmap
lasso_MetaB_HostT_heatmap_p0 <- 
  ggplot(lasso_MetaB_HostT_data, 
         aes(x = MetaB, y = HostT, fill = Fill)
  ) +
  geom_shape(radius = unit(1.5, "pt"),
             expand = unit(3, "pt"),
             color = "white",
             linewidth = 0,
             show.legend = T
  ) +
  scale_fill_gradientn(colors = c(c("#ADD8E6", "#BDDFEB", "#CDE7F0", "#DEEFF5", "#EEF7FA"), "white", rev(c("#D9636C", "#E08289", "#E8A1A6", "#EFC0C4", "#F7DFE1"))),
                       name = "Correlation", 
                       na.value = NA,
                       breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  scale_x_discrete(position = "top",
                   labels = rev(paste0("MetaB-", lasso_metabolome_info$module))
  ) + 
  scale_y_discrete(position = "right",
                   labels = paste0("HostT-", lasso_transcriptome_info$label)
  ) + 
  theme_void() + 
  theme(axis.text.x = element_text(family = "Arial", size = 6, hjust = 0, vjust = 0, angle = 45, lineheight = 1, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(5, "pt"),
        legend.key.spacing.y = unit(2, "pt"),
        legend.key.spacing.x = unit(2, "pt"),
        legend.text = element_text(family = "Arial", face = "plain", size = 4, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.title = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.position = "none",
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.ticks.length = unit(1, "pt"),
        legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white"),
        panel.border = element_rect(fill = NA, color = "black", linetype = 1, linewidth = 0.25),
  )

lasso_MetaB_HostT_heatmap_p0

### legendry
# HostT_legendry
lasso_transcriptome_info <- 
  lasso_transcriptome_info %>% 
  mutate(type = case_when(
    lasso_transcriptome_info$Subcategory == "Signal transduction" ~ "Signal transduction",
    lasso_transcriptome_info$Subcategory == "Folding, sorting and degradation" ~ "Folding/Sorting/Degradation",
    lasso_transcriptome_info$Subcategory == "Amino acid metabolism" ~ "Amino acid metabolism",
    lasso_transcriptome_info$Subcategory == "Carbohydrate metabolism" ~ "Carbohydrate metabolism",
    lasso_transcriptome_info$Subcategory == "Signaling molecules and interaction" ~ "Signaling molecules and interaction",
    TRUE ~ "Other"
  ))

lasso_transcriptome_info$pathway <- factor(lasso_transcriptome_info$pathway, levels = lasso_transcriptome_info$pathway)
lasso_transcriptome_info$type <- factor(lasso_transcriptome_info$type, levels = c("Signal transduction", "Folding/Sorting/Degradation", "Amino acid metabolism", "Carbohydrate metabolism", "Signaling molecules and interaction", "Other"))

lasso_HostT_legendry <-
  ggplot(lasso_transcriptome_info) + 
  geom_tile(aes(x = 0, y = pathway, fill = type), 
            color = "white", 
            height = 1,
            linewidth = 0.5,
            show.legend = T
  ) + 
  scale_fill_manual(values = c(paletteer::paletteer_d("wesanderson::AsteroidCity1"), "grey")) +
  scale_color_manual(values = c(paletteer::paletteer_d("wesanderson::AsteroidCity1"), "grey")) +
  theme_void() +
  theme(legend.position = "none",
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(5, "pt"),
        legend.key.spacing.y = unit(2, "pt"),
        legend.key.spacing.x = unit(2, "pt"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.text = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.title = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linetype = 1, linewidth = 0.25),
        panel.background = element_blank()
  )

lasso_HostT_legendry

# MetaB_legendry 
lasso_metabolome_info$MetaB <- paste0("MetaB_", lasso_metabolome_info$module)
lasso_metabolome_info$MetaB <- factor(lasso_metabolome_info$MetaB, levels = lasso_metabolome_info$MetaB)
lasso_metabolome_info$module <- factor(lasso_metabolome_info$module, levels = lasso_metabolome_info$module)
lasso_metabolome_info$Module.category <- factor(lasso_metabolome_info$Module.category, 
                                                levels = c("Benzenoids", "Lipids and lipid-like molecules", "Organic acids and derivatives", 
                                                           "Organoheterocyclic compounds", "Organic oxygen compounds", "Nucleosides, nucleotides, and analogues", "Miscellaneous"))

lasso_MetaB_legendry <-
  ggplot(lasso_metabolome_info) + 
  geom_tile(aes(x = rev(module), y = 0, fill = Module.category), 
            color = "white", 
            height = 1,
            linewidth = 0.5,
            show.legend = T
  ) + 
  scale_fill_manual(values = c(paletteer::paletteer_d("LaCroixColoR::Apricot")[1:6], "grey")) +
  # scale_color_manual(values = c(paletteer::paletteer_d("LaCroixColoR::Apricot")[1:6], "grey")) +
  theme_void() +
  theme(legend.position = "none",
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(5, "pt"),
        legend.key.spacing.y = unit(2, "pt"),
        legend.key.spacing.x = unit(2, "pt"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.text = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        legend.title = element_text(family = "Arial", face = "plain", size = 4, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linetype = 1, linewidth = 0.25),
        panel.background = element_blank()
  )

lasso_MetaB_legendry

lasso_MetaB_HostT_heatmap <- lasso_MetaB_HostT_heatmap_p0  %>% insert_left(lasso_HostT_legendry, width = 0.06) %>% insert_bottom(lasso_MetaB_legendry, height = 0.06)
lasso_MetaB_HostT_heatmap0 <-  lasso_MetaB_HostT_heatmap_p0 + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
lasso_MetaB_HostT_heatmap1 <- lasso_MetaB_HostT_heatmap_p0 + theme(axis.text.x = element_blank())

lasso_MetaB_HostT_heatmap
lasso_MetaB_HostT_heatmap0
lasso_MetaB_HostT_heatmap1

ggsave("1_multi_omic_association/26_tatol_MetaB_HostT_heatmap0.svg", plot = lasso_MetaB_HostT_heatmap0, width = 10, height = 12, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/26_tatol_MetaB_HostT_heatmap0.tiff", plot = lasso_MetaB_HostT_heatmap0, width = 10, height = 12, units = "cm", dpi = 300)

ggsave("1_multi_omic_association/27_tatol_MetaB_HostT_heatmap1.svg", plot = lasso_MetaB_HostT_heatmap1, width = 15, height = 12, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/27_tatol_MetaB_HostT_heatmap1.tiff", plot = lasso_MetaB_HostT_heatmap1, width = 15, height = 12, units = "cm", dpi = 300)

ggsave("1_multi_omic_association/28_tatol_HostT_legendry.svg", plot = lasso_HostT_legendry, width = 0.30, height = 12, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/28_tatol_HostT_legendry.tiff", plot = lasso_HostT_legendry, width = 0.30, height = 12, units = "cm", dpi = 300)

ggsave("1_multi_omic_association/29_tatol_MetaB_legendry.svg", plot = lasso_MetaB_legendry, width = 10, height = 0.30, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/29_tatol_MetaB_legendry.tiff", plot = lasso_MetaB_legendry, width = 10, height = 0.30, units = "cm", dpi = 300)


### save all legend
lasso_MetaG_MetaB_heatmap_p0_legend <- cowplot::get_legend(lasso_MetaG_MetaB_heatmap_p0 + theme(legend.position = "right"))
lasso_MetaG_legendry_legend <- cowplot::get_legend(lasso_MetaG_legendry + theme(legend.position = "right"))
lasso_MetaB_HostT_heatmap_p0_legend <- cowplot::get_legend(lasso_MetaB_HostT_heatmap_p0 + theme(legend.position = "right"))
lasso_MetaB_legendry_legend <- cowplot::get_legend(lasso_MetaB_legendry + theme(legend.position = "right"))
lasso_HostT_legendry_legend <- cowplot::get_legend(lasso_HostT_legendry + theme(legend.position = "right"))
lasso_MetaG_MetaB_HostT_legend <- cowplot::plot_grid(lasso_MetaG_MetaB_heatmap_p0_legend,
                                                     lasso_MetaB_HostT_heatmap_p0_legend,
                                                     lasso_MetaG_legendry_legend, 
                                                     lasso_MetaB_legendry_legend, 
                                                     lasso_HostT_legendry_legend, ncol = 5)
lasso_MetaG_MetaB_HostT_legend

ggsave("1_multi_omic_association/30_total_MetaG_MetaB_HostT_legend.svg", plot = lasso_MetaG_MetaB_HostT_legend, width = 20, height = 10, units = "cm", dpi = 300)
# ggsave("1_multi_omic_association/30_total_MetaG_MetaB_HostT_legend.tiff", plot = lasso_MetaG_MetaB_HostT_legend, width = 20, height = 10, units = "cm", dpi = 300)



### ================================================================================================================================================================================== ###


### correlation analysis
multi_omic <- 
  metagenome %>%
  rename_with(~ paste0("MetaG_", .x)) %>%
  rownames_to_column("sample_id") %>% 
  left_join(metabolome %>% rename_with(~ paste0("MetaB_", .x)) %>% rownames_to_column("sample_id") , by = "sample_id") %>%
  left_join(transcriptome %>% rename_with(~ paste0("HostT_", .x)) %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  mutate(group = c(rep("fam", 10), rep("hos", 10), rep("lab", 10), rep("res", 10), rep("wild", 10))) %>% 
  mutate(habitat = c(rep("peridomestic", 20), rep("lab", 10), rep("peridomestic", 10) ,rep("wild", 10)))

multi_omic$group <- factor(multi_omic$group , levels = c("wild", "lab", "fam", "res", "hos"))

multi_omic$habitat <- factor(multi_omic$habitat, levels = c("wild", "lab", "peridomestic"))

MetaG <- paste0("MetaG_", colnames(metagenome))

MetaB <- paste0("MetaB_",colnames(metabolome))

HostT <- paste0("HostT_", colnames(transcriptome))

### scatter plot 
p1 <- 
  ggplot(multi_omic, 
         aes(x = MetaG_M00019, y = MetaB_M144) 
  )+
  geom_point(aes(fill = group, color = group),
             size = 1,
             shape = 21,
             stroke = 0.1,
             show.legend = F
  ) + 
  geom_smooth(formula = 'y ~ x',
              method = 'lm', 
              se = T,
              linewidth = 0.25
  ) + 
  stat_cor(method = "spearman",  #  "pearson" (default), "kendall", or "spearman".
           size = 2,
           label.sep = "\n",
           hjust = 0,
           vjust = 1,
           label.x.npc = 0.05,
           label.y.npc = 1,
           family = "Arial"
  ) +
  scale_color_manual(values = c("wild" = wild_color, "lab" = lab_color, "fam" = fam_color, "res" = res_color, "hos" = hos_color)) +
  scale_fill_manual(values = c("wild" = wild_color, "lab" = lab_color, "fam" = fam_color, "res" = res_color, "hos" = hos_color)) +
  # stat_regline_equation(size = 5) +
  theme_minimal() + 
  theme(
    axis.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

p1

p2 <- 
  ggplot(multi_omic, 
         aes(x = MetaG_M00535, y = MetaB_M144) 
  )+
  geom_point(aes(fill = group, color = group),
             size = 1,
             shape = 21,
             stroke = 0.1,
             show.legend = F
  ) + 
  geom_smooth(formula = 'y ~ x',
              method = 'lm', 
              se = T,
              linewidth = 0.25
  ) + 
  stat_cor(method = "spearman",  #  "pearson" (default), "kendall", or "spearman".
           size = 2,
           label.sep = "\n",
           hjust = 0,
           vjust = 1,
           label.x.npc = 0.05,
           label.y.npc = 1,
           family = "Arial"
  ) +
  scale_color_manual(values = c("wild" = wild_color, "lab" = lab_color, "fam" = fam_color, "res" = res_color, "hos" = hos_color)) +
  scale_fill_manual(values = c("wild" = wild_color, "lab" = lab_color, "fam" = fam_color, "res" = res_color, "hos" = hos_color)) +
  # stat_regline_equation(size = 5) +
  theme_minimal() + 
  theme(
    axis.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

p2

p3 <- 
  ggplot(multi_omic, 
         aes(x = MetaG_M00570, y = MetaB_M144) 
  )+
  geom_point(aes(fill = group, color = group),
             size = 1,
             shape = 21,
             stroke = 0.1,
             show.legend = F
  ) + 
  geom_smooth(formula = 'y ~ x',
              method = 'lm', 
              se = T,
              linewidth = 0.25
  ) + 
  stat_cor(method = "spearman",  #  "pearson" (default), "kendall", or "spearman".
           size = 2,
           label.sep = "\n",
           hjust = 0,
           vjust = 1,
           label.x.npc = 0.05,
           label.y.npc = 1,
           family = "Arial"
  ) +
  scale_color_manual(values = c("wild" = wild_color, "lab" = lab_color, "fam" = fam_color, "res" = res_color, "hos" = hos_color)) +
  scale_fill_manual(values = c("wild" = wild_color, "lab" = lab_color, "fam" = fam_color, "res" = res_color, "hos" = hos_color)) +
  # stat_regline_equation(size = 5) +
  theme_minimal() + 
  theme(
    axis.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

p3

p4 <- 
  ggplot(multi_omic, 
         aes(x = MetaG_M00432, y = MetaB_M144) 
  )+
  geom_point(aes(fill = group, color = group),
             size = 1,
             shape = 21,
             stroke = 0.1,
             show.legend = F
  ) + 
  geom_smooth(formula = 'y ~ x',
              method = 'lm', 
              se = T,
              linewidth = 0.25
  ) + 
  stat_cor(method = "spearman",  #  "pearson" (default), "kendall", or "spearman".
           size = 2,
           label.sep = "\n",
           hjust = 0,
           vjust = 1,
           label.x.npc = 0.05,
           label.y.npc = 1,
           family = "Arial"
  ) +
  scale_color_manual(values = c("wild" = wild_color, "lab" = lab_color, "fam" = fam_color, "res" = res_color, "hos" = hos_color)) +
  scale_fill_manual(values = c("wild" = wild_color, "lab" = lab_color, "fam" = fam_color, "res" = res_color, "hos" = hos_color)) +
  # stat_regline_equation(size = 5) +
  theme_minimal() + 
  theme(
    axis.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

p4

p5 <- 
  ggplot(multi_omic, 
         aes(y = HostT_ko04150, x = MetaB_M144) 
  )+
  geom_point(aes(fill = group, color = group),
             size = 1,
             shape = 21,
             stroke = 0.1,
             show.legend = F
  ) + 
  geom_smooth(formula = 'y ~ x',
              method = 'lm', 
              se = T,
              linewidth = 0.25
  ) + 
  stat_cor(method = "spearman",  #  "pearson" (default), "kendall", or "spearman".
           size = 2,
           label.sep = "\n",
           hjust = 0,
           vjust = 1,
           label.x.npc = 0.05,
           label.y.npc = 1,
           family = "Arial"
  ) +
  scale_color_manual(values = c("wild" = wild_color, "lab" = lab_color, "fam" = fam_color, "res" = res_color, "hos" = hos_color)) +
  scale_fill_manual(values = c("wild" = wild_color, "lab" = lab_color, "fam" = fam_color, "res" = res_color, "hos" = hos_color)) +
  # stat_regline_equation(size = 5) +
  theme_minimal() + 
  theme(
    axis.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

p5

library(patchwork)

MetaG_MetaB = p1 / p2 / p3 / p4

MetaG_MetaB

MetaB_HostT = p5

MetaB_HostT

ggplot2::ggsave("1_multi_omic_association/31_correlation.svg", plot = MetaG_MetaB, width = 4.5, height = 16, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/31_correlation.tiff", plot = MetaG_MetaB, width = 4.5, height = 16, units = "cm", dpi = 300)

ggplot2::ggsave("1_multi_omic_association/32_correlation.svg", plot = MetaB_HostT, width = 4.5, height = 16, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/32_correlation.tiff", plot = MetaB_HostT, width = 4.5, height = 16, units = "cm", dpi = 300)





#### BCAA-mTOR Axis

### Metagenome module & KO

### Metagenome BCAA module heatmap
BCAA_module_heatmap_data <- 
  read.table("1_multi_omic_association/33_BCAA_module.txt", header = T, row.names = 1) %>% 
  t() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = c(rep("Fam", 10), rep("Lab", 10), rep("Hos", 10), rep("Res", 10), rep("Wild", 10))) %>%
  mutate(group2 = ifelse(group == "Wild", "Wild", "Non-wild"))

BCAA_module_heatmap_mean_matrix <- 
  BCAA_module_heatmap_data %>%
  mutate(across(2:5, as.numeric)) %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% 
  mutate(across(where(is.numeric), ~ scales::rescale(.x, to = c(-2, 2)))) %>%
  mutate(group2 = ifelse(group == "Wild", "Wild", "Non-wild")) %>%
  group_by(group2) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

BCAA_module_heatmap_mean_matrix_long <- 
  BCAA_module_heatmap_mean_matrix %>% 
  pivot_longer(-c("group2"), names_to = "module", values_to = "enrichment")

# BCAA_module_heatmap_mean_matrix_long$group <- factor(BCAA_module_heatmap_mean_matrix_long$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))
BCAA_module_heatmap_mean_matrix_long$group2 <- factor(BCAA_module_heatmap_mean_matrix_long$group2, levels = c("Wild", "Non-wild"))

BCAA_module_heatmap_mean_matrix_long$module <- factor(BCAA_module_heatmap_mean_matrix_long$module, levels = c("M00019", "M00570", "M00535", "M00432"))

min <- min(BCAA_module_heatmap_mean_matrix_long$enrichment)
max <- max(BCAA_module_heatmap_mean_matrix_long$enrichment)

BCAA_module_heatmap <- 
  ggplot(BCAA_module_heatmap_mean_matrix_long) +
  geom_tile(aes(x = group2, y = rev(module), fill = enrichment),
            color = "white",
            height = 1,
            linewidth = 0.45,
            show.legend = T
  ) + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white",
                       high = "#FF9898FF",
                       limit = c(-2, 2),
                       breaks = c(-2, 2)
  ) +
  labs(fill = "Scaled Value") + 
  theme_void() +
  theme(
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "top",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

print(BCAA_module_heatmap)

ggplot2::ggsave("1_multi_omic_association/34_BCAA_module_heatmap.svg", plot = BCAA_module_heatmap, width = 5, height = 10, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/34_BCAA_module_heatmap.tiff", plot = BCAA_module_heatmap, width = 5, height = 10, units = "cm", dpi = 300)


### Metagenome BCAA module box plot 
key_metagenome_module <- 
  BCAA_module_heatmap_data %>%
  as.data.frame() %>%
  pivot_longer(c(2:5), names_to = "module", values_to = "enrich score") %>% 
  mutate(group = factor(rep(c(rep("Fam", 40), rep("Lab", 40), rep("Hos", 40), rep("Res", 40), rep("Wild", 40))), levels = c("Wild", "Lab", "Fam", "Res", "Hos")))

key_metagenome_module$module <- factor(key_metagenome_module$module, unique(key_metagenome_module$module))

key_metagenome_module_group <-
  key_metagenome_module %>%
  group_by(group, module) %>%
  summarise(
    mean = mean(`enrich score`),
    sd = sd(`enrich score`),
    count = n(),
    se = sd / sqrt(count),
    .groups = "drop"
  )

key_metagenome_module$top_label <- "BCAAs Biosynthetic Modules"
key_metagenome_module_group$top_label <- "BCAAs Biosynthetic Modules"

key_metagenome_module_box <- 
  ggplot(key_metagenome_module,
         aes(x = group, y = `enrich score`) 
  ) + 
  geom_point(aes(fill = group, color = group),
             stroke = 0.1,
             size = 1,
             shape = 21,
             stat = "identity",
             position = position_jitter(height = 0, width = 0.2, seed = 123),
             show.legend = F
  ) +
  geom_boxplot(aes(fill = group, color = group),
               # color = "black",
               outlier.shape = NA, 
               alpha = 0.3,
               width = 0.85,
               staplewidth = 0.4,
               linewidth = 0.1,
               show.legend = F
  ) +
  ggh4x::facet_nested_wrap(~ top_label + module, 
                           scales = "free_y",
                           ncol = 4
  ) + 
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  stat_compare_means(method = "anova",
                     # method = "kruskal.test",
                     size = 2
  ) +
  stat_compare_means(
    comparisons = list(c("Wild", "Lab"), c("Wild", "Fam"), c("Wild", "Res"), c("Wild", "Hos")),
    # method = "wilcox.test",
    method = "t.test",
    # label = "p.adj",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
    # p.adjust.method = "BH", 
    p.adjust.method = "holm", 
    # method.args = list(exact = FALSE, correct = TRUE),
    size = 3,
    vjust = 0.67,
    bracket.size = 0.25,
    step.increase = 0.075,
    tip.length = 0
  ) +
  labs(y = "Enrichemnt score") + 
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=1, angle=45, lineheight=1, color="black"),  
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = ggfun::element_roundrect(
      # fill = "#DEEFF5",
      fill = "#c9caca", 
      color = NA, 
      linewidth = 0, 
      r = unit(0.25, "snpc")
    ),
    strip.text = element_text(family = "Arial", face = "bold", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black")
  )

key_metagenome_module_box

ggplot2::ggsave("1_multi_omic_association/35_key_metagenome_module_box.svg", plot = key_metagenome_module_box , width = 10, height = 4.25, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/35_key_metagenome_module_box.tiff", plot = key_metagenome_module_box , width = 10, height = 4.25, units = "cm", dpi = 300)


### Metagenome BCAA KO heatmap
BCAA_KO_abundance <- 
  read.table("1_multi_omic_association/36_BCAA_KO_abundance.txt", header = T, row.names = 1) %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(group = c(rep("Fam", 10), rep("Lab", 10), rep("Hos", 10), rep("Res", 10), rep("Wild", 10))) %>%
  mutate(group2 = ifelse(group == "Wild", "Wild", "Non-wild"))

BCAA_KO_abundance_mean <- 
  BCAA_KO_abundance %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  mutate(across(where(is.numeric), ~ scales::rescale(.x, to = c(-2, 2)))) %>%
  mutate(group2 = ifelse(group == "Wild", "Wild", "Non-wild")) %>% 
  group_by(group2) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

BCAA_KO_abundance_mean_long <- 
  BCAA_KO_abundance_mean %>% 
  pivot_longer(-group2, names_to = "KO", values_to = "scaled TPM")

# BCAA_KO_abundance_mean_long$group <- factor(BCAA_KO_abundance_mean_long$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))
BCAA_KO_abundance_mean_long$group2 <- factor(BCAA_KO_abundance_mean_long$group2, levels = c("Wild", "Non-wild"))

BCAA_KO_abundance_mean_long$KO <- factor(BCAA_KO_abundance_mean_long$KO, levels = c("K01754",
                                                                                    "K00053",
                                                                                    "K01687",
                                                                                    "K01652",
                                                                                    "K01653",
                                                                                    "K00826",
                                                                                    "K01649",
                                                                                    "K01703",
                                                                                    "K01704",
                                                                                    "K00052",
                                                                                    "K09011",
                                                                                    "K17989"))

min <- min(BCAA_KO_abundance_mean_long$`scaled TPM`)
max <- max(BCAA_KO_abundance_mean_long$`scaled TPM`)

# Metagenomoe BCAA KO abundance Heatmap Plot 
BCAA_KO_abundance_heatmap_p <- 
  ggplot(BCAA_KO_abundance_mean_long, 
         aes(x = group2, y = KO, fill = `scaled TPM`),
  ) +
  geom_tile(color = "white",
            height = 1,
            linewidth = 0.45,
            lineend = "square",
            show.legend = T
  ) +
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF", 
                       limit = c(-2, 2),
                       breaks = c(-2, 2)
  ) +
  scale_x_discrete(position = "bottom") +
  labs(fill = "TPM Z-score", color = "TPM Z-score", size = "TPM Z-score") +
  theme_void() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

BCAA_KO_abundance_heatmap_p

ggplot2::ggsave("1_multi_omic_association/37_BCAA_KO_abundance_heatmap.svg", plot = BCAA_KO_abundance_heatmap_p, width = 4, height =8, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/37_BCAA_KO_abundance_heatmap.tiff", plot = BCAA_KO_abundance_heatmap_p, width = 4, height = 8, units = "cm", dpi = 300)


### Metagenome BCAA KO box plot
monochromeR::generate_palette("lightblue", modification = "go_lighter", n_colors = 5, view_palette = TRUE, view_labels = FALSE)

key_metagenome_KO <- 
  BCAA_KO_abundance %>% 
  rownames_to_column(var = "sample") %>% 
  as.data.frame() %>% 
  pivot_longer(-c(sample, group, group2), names_to = "KO", values_to = "TPM")

key_metagenome_KO$group  <- factor(rep(c(rep("Fam", 120), rep("Lab",120), rep("Hos", 120), rep("Res",120), rep("Wild",120))), levels = c("Wild", "Lab", "Fam", "Res", "Hos"))
key_metagenome_KO$KO <- factor(key_metagenome_KO$KO, unique(key_metagenome_KO$KO))

key_metagenome_KO$logTPM <- log2(key_metagenome_KO$TPM + 1)

key_metagenome_KO$top_label <- "BCAAs Biosynthetic Genes"

key_metagenome_KO_box <- 
  ggplot(key_metagenome_KO,
         aes(x = group, y = logTPM) 
  ) + 
  geom_point(aes(fill = group, color = group,),
             stroke = 0.1,
             size = 1,
             shape = 21,
             stat = "identity",
             position = position_jitter(height = 0, width = 0.2, seed = 123),
             show.legend = F
  ) +
  geom_boxplot(aes(fill = group, color = group),
               # color = "black",
               outlier.shape = NA, 
               alpha = 0.3,
               width = 0.85,
               staplewidth = 0.4,
               linewidth = 0.1,
               show.legend = F
  ) +
  ggh4x::facet_nested_wrap(~ top_label + KO,
                           scales = "free_y",
                           ncol = 4, 
                           nrow = 3
  ) + 
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  stat_compare_means(method = "anova",
                     # method = "kruskal.test",
                     size = 2
  ) +
  stat_compare_means(
    comparisons = list(c("Wild", "Lab"), c("Wild", "Fam"), c("Wild", "Res"), c("Wild", "Hos")),
    # method = "wilcox.test",
    method = "t.test",
    # label = "p.adj",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
    # p.adjust.method = "BH", 
    p.adjust.method = "holm", 
    # method.args = list(exact = FALSE, correct = TRUE),
    size = 3,
    vjust = 0.67,
    bracket.size = 0.25,
    step.increase = 0.075,
    tip.length = 0
  ) +
  labs(y = "log2TPM") + 
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=1, angle=45, lineheight=1, color="black"),  
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = ggfun::element_roundrect(
      # fill = "#DEEFF5",
      fill = "#c9caca", 
      color = NA, 
      linewidth = 0, 
      r = unit(0.25, "snpc")
    ),
    strip.text = element_text(family = "Arial", face = "bold", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black")
  )

key_metagenome_KO_box

ggplot2::ggsave("1_multi_omic_association/38_BCAA_KO_box.svg", plot = key_metagenome_KO_box , width = 10, height = 11.75, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/38_BCAA_KO_box.tiff", plot = key_metagenome_KO_box , width = 10, height =11.75, units = "cm", dpi = 300)


### KEGG module & KO mapping

# https://www.kegg.jp/kegg-bin/show_module?M00432
# M00432  K01702
# M00432  K01703
# M00432  K01704
# M00432  K01649
# M00432  K00052
# M00432  K00826

# https://www.kegg.jp/kegg-bin/show_module?M00535
# M00535  K01704
# M00535  K09011
# M00535  K00052
# M00535  K01703

# https://www.kegg.jp/kegg-bin/show_module?M00019
# M00019  K01653
# M00019  K00826
# M00019  K01687
# M00019  K01652
# M00019  K00053

# https://www.kegg.jp/kegg-bin/show_module?M00570
# M00570  K01754
# M00570  K01653
# M00570  K00826
# M00570  K01687
# M00570  K01652
# M00570  K00053
# M00570  K17989

### Metagenome BCAA module heatmap
BCAA_module_heatmap_mean_matrix2 <- 
  BCAA_module_heatmap_data %>%
  mutate(across(2:5, as.numeric)) %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.x))))

BCAA_module_heatmap_mean_matrix_long2 <- 
  BCAA_module_heatmap_mean_matrix2 %>% 
  pivot_longer(-group, names_to = "module", values_to = "enrichment")

BCAA_module_heatmap_mean_matrix_long2$group <- factor(BCAA_module_heatmap_mean_matrix_long2$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

BCAA_module_heatmap_mean_matrix_long2$module <- factor(BCAA_module_heatmap_mean_matrix_long2$module, levels = c("M00535", "M00019", "M00570", "M00432"))

BCAA_module_bubble_plot <- 
  ggplot(BCAA_module_heatmap_mean_matrix_long2,
         aes(x = group, y = rev(module))
  ) +
  geom_point(aes(size = enrichment, fill = group, color = group), 
             shape = 21, 
             show.legend = T
  ) +
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  labs(fill = "Scaled Value") + 
  theme_void() +
  theme(
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

print(BCAA_module_bubble_plot)


link_linewith = 0.25
linetype = 1

module <- data.frame(Module = rev(c("M00019", "M00570", "M00535", "M00432")),
                     Y = c(9.8, 7.6, 5.4, 3.2))

module$Module <- factor(module$Module, levels = c("M00019", "M00570", "M00535", "M00432"))


# plot KO bubble plot
BCAA_KO_abundance_mean2 <- 
  BCAA_KO_abundance %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  mutate(across(where(is.numeric), ~ scales::rescale(.x)))

BCAA_KO_abundance_mean_long2 <- 
  BCAA_KO_abundance_mean2 %>% 
  # pivot_longer(-c("group"), names_to = "KO", values_to = "scaled TPM")
  pivot_longer(-c("group"), names_to = "KO", values_to = "scaled TPM")

BCAA_KO_abundance_mean_long2$group <- factor(BCAA_KO_abundance_mean_long2$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

BCAA_KO_abundance_mean_long$KO <- factor(BCAA_KO_abundance_mean_long$KO, levels = c("K01754",
                                                                                    "K00053",
                                                                                    "K01687",
                                                                                    "K01652",
                                                                                    "K01653",
                                                                                    "K00826",
                                                                                    "K01649",
                                                                                    "K01703",
                                                                                    "K01704",
                                                                                    "K00052",
                                                                                    "K09011",
                                                                                    "K17989"))

BCAA_KO_abundance_bubble_p <- 
  ggplot(BCAA_KO_abundance_mean_long2, 
         aes(x = group, y = KO),
  ) +
  geom_point(aes(size = `scaled TPM`, fill = group, color = group), 
             shape = 21, 
             show.legend = T
  ) +
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_x_discrete(position = "bottom") +
  labs(fill = "Scaled Value") +
  theme_void() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

BCAA_KO_abundance_bubble_p

link_linewith = 0.25
linetype = 1

link <-
  ggplot(BCAA_KO_abundance_mean_long) +
  # KO ponit
  geom_point(aes(x = 2, y = KO), fill = "grey" , color = "grey", size = 1, shape = 21) +
  # M00432 link
  annotate("segment", x = 2, xend = 4, y = 11, yend = 7.6, color = "pink", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 10,  yend = 7.6, color = "pink", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 9,  yend = 7.6, color = "pink", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 8,  yend = 7.6, color = "pink", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 7,  yend = 7.6, color = "pink", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  # M00535 link
  annotate("segment", x = 2, xend = 4, y = 12, yend =  9.8, color = "lightblue", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 11, yend = 9.8, color = "lightblue", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 10,  yend = 9.8, color = "lightblue", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 9,  yend = 9.8, color = "lightblue", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  # M00019 link
  annotate("segment", x = 2, xend = 4, y = 7, yend = 5.4, color = "#BF281BFF", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 6, yend = 5.4, color = "#BF281BFF", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 5, yend = 5.4, color = "#BF281BFF", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 4, yend = 5.4, color = "#BF281BFF", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 3, yend = 5.4, color = "#BF281BFF", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  # M00570 link
  annotate("segment", x = 2, xend = 4, y = 7,  yend = 3.2, color = "#AC9ECEFF", linewidth = link_linewith, lineend = "square", linetype = linetype) + 
  annotate("segment", x = 2, xend = 4, y = 6,  yend = 3.2, color = "#AC9ECEFF", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 5,  yend = 3.2, color = "#AC9ECEFF", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 4,  yend = 3.2, color = "#AC9ECEFF", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 3,  yend = 3.2, color = "#AC9ECEFF", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 2,  yend = 3.2, color = "#AC9ECEFF", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  annotate("segment", x = 2, xend = 4, y = 1, yend = 3.2, color = "#AC9ECEFF", linewidth = link_linewith, lineend = "square", linetype = linetype) +
  # KO ponit
  geom_point(aes(x = 2, y = KO), fill = "grey" , color = "grey", size = 1, shape = 21) +
  # module point
  geom_point(data = module, aes(x = 4, y = Y, fill = Module, color = Module), size = 10, shape = 22, show.legend = F) + 
  scale_color_manual(values = c("M00535" = "pink", "M00432" =  "lightblue", "M00570" =  "#BF281BFF",  "M00019" = "#AC9ECEFF")) + 
  scale_fill_manual(values = c("M00535" = "pink", "M00432" =  "lightblue", "M00570" =  "#BF281BFF",  "M00019" = "#AC9ECEFF")) + 
  theme_void()

link

library(patchwork) 
Module_to_KO_link_plot <- BCAA_KO_abundance_bubble_p  + link + BCAA_module_bubble_plot + plot_layout(ncol = 3, widths = c(2, 1, 2))

Module_to_KO_link_plot

ggplot2::ggsave("1_multi_omic_association/39_BCAA_KO_abundance_heatmap_link.svg", plot = Module_to_KO_link_plot, width = 9, height = 7, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/39_BCAA_KO_abundance_heatmap_link.tiff", plot = Module_to_KO_link_plot, width =9, height = 7, units = "cm", dpi = 300)


### Key metabolome module & compuounds

### key metabolome module heatmap 
BCAA_metabolome_module <- 
  metabolome %>%
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>%
  dplyr::select(c("sample", "M144")) %>% 
  pivot_longer(-sample, names_to = "module", values_to = "enrich score") %>% 
  mutate(group = factor(c(rep("Fam", 10), rep("Hos", 10), rep("Lab", 10), rep("Res", 10), rep("Wild", 10)), levels = c("Wild", "Lab", "Fam", "Res", "Hos"))) %>%
  mutate(group2 = factor(ifelse(group == "Wild", "Wild", "Non-wild"), levels = c("Wild", "Non-wild")))

BCAA_metabolome_module_grouped <- 
  BCAA_metabolome_module %>% 
  group_by(module, group) %>%
  summarise(across(c(`enrich score`), mean, na.rm = TRUE), .groups = "drop") %>%
  group_by(module) %>%
  mutate(across(c(`enrich score`), ~ scales::rescale(.x, to = c(-2, 2)))) %>%
  mutate(group2 = factor(ifelse(group == "Wild", "Wild", "Non-wild"), levels = c("Wild", "Non-wild"))) %>%
  group_by(group2, module) %>%
  summarise(across(c(`enrich score`), mean, na.rm = TRUE), .groups = "drop")

BCAA_metabolome_module_heatmap <- 
  ggplot(BCAA_metabolome_module_grouped, 
         aes(x = group2, y = module, fill = `enrich score`),
  ) +
  geom_tile(color = "white",
            height = 1,
            linewidth = 0.45,
            lineend = "square",
            show.legend = T
  ) +
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF", 
                       limit = c(-2, 2),
                       breaks = c(-2, 2)
  ) +
  scale_x_discrete(position = "bottom") +
  labs(fill = "Scale value") +
  theme_void() +
  theme(
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 45, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

BCAA_metabolome_module_heatmap

ggplot2::ggsave("1_multi_omic_association/40_BCAA_metabolome_module_heatmap.svg", plot = BCAA_metabolome_module_heatmap , width = 2.2, height = 1, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/40_BCAA_metabolome_module_heatmap.tiff", plot = BCAA_metabolome_module_heatmap , width = 2.2, height = 1, units = "cm", dpi = 300)


### key metabolome module box
monochromeR::generate_palette("#F27127FF", modification = "go_lighter", n_colors = 5, view_palette = TRUE, view_labels = FALSE)

BCAA_metabolome_module <- 
  read.table("merge_module_Eigengenes.txt", header = T, row.names = 1) %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>%
  dplyr::select(c("sample", "MetaB_M144")) %>% 
  pivot_longer(-sample, names_to = "module", values_to = "enrich score") %>% 
  mutate(group = factor(c(rep("Wild", 10), rep("Lab", 10), rep("Fam", 10), rep("Res", 10), rep("Hos", 10)), levels = c("Wild", "Lab", "Fam", "Res", "Hos")))

BCAA_metabolome_module <- 
  BCAA_metabolome_module %>%
  mutate(module = case_when(
    module == "MetaB_M144" ~ "MetaB-M144"
  ))
factor(BCAA_metabolome_module$module, levels = unique(BCAA_metabolome_module$module))

BCAA_metabolome_module$module <- factor(BCAA_metabolome_module$module, levels = unique(BCAA_metabolome_module$module))

BCAA_metabolome_module$top_label <- "BCAAs metabolomic module"

BCAA_metabolome_module_box_plot <-
  ggplot(BCAA_metabolome_module,
         aes(x = group, y = `enrich score`)
  ) +
  geom_point(aes(fill = group, color = group),
             stroke = 0.1,
             stat = "identity",
             position = position_jitter(height = 0, width = 0.2, seed = 123),
             size = 1,
             shape = 21,
             show.legend = F
  ) +
  geom_boxplot(aes(fill = group, color = group),
               # color = "black",
               outlier.shape = NA, 
               alpha = 0.3,
               width = 0.85,
               staplewidth = 0.4,
               linewidth = 0.1,
               show.legend = F
  ) +
  facet_nested_wrap(~ top_label + module, 
             scales = "free_y",
             ncol = 3
  ) + 
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  stat_compare_means(method = "anova",
                     # method = "kruskal.test",
                     size = 2
  ) +
  stat_compare_means(
    comparisons = list(c("Wild", "Lab"), c("Wild", "Fam"), c("Wild", "Res"), c("Wild", "Hos")),
    # method = "wilcox.test",
    method = "t.test",
    # label = "p.adj",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
    # p.adjust.method = "BH", 
    p.adjust.method = "holm", 
    # method.args = list(exact = FALSE, correct = TRUE),
    size = 3,
    vjust = 0.67,
    bracket.size = 0.25,
    step.increase = 0.075,
    tip.length = 0
  ) +
  labs(y = "Enrichment score") + 
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=1, angle=45, lineheight=1, color="black"),  
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = ggfun::element_roundrect(
      # fill = "#DEEFF5",
      fill = "#c9caca", 
      color = NA, 
      linewidth = 0, 
      r = unit(0.25, "snpc")
    ),
    strip.text = element_text(family = "Arial", face = "bold", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black")
  )

BCAA_metabolome_module_box_plot

ggplot2::ggsave("1_multi_omic_association/41_BCAA_metabolome_module_box.svg", plot = BCAA_metabolome_module_box_plot, width = 5.5, height = 4.25, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/41_BCAA_metabolome_module_box.tiff", plot = BCAA_metabolome_module_box_plot, width = 5.5, height = 4.25, units = "cm", dpi = 300)


### key metabolome cpds heatmap
BCAA_non_target_Derivative_cpd <- read.table("1_multi_omic_association/42_metabolome_LOESS_PQN_missForest_matrix_noQC.tsv", header = T, row.names = 1, sep = "\t")

BCAA_non_target_Derivative_cpd_long <- 
  BCAA_non_target_Derivative_cpd %>% 
  rownames_to_column(var = "CPD_ID") %>% 
  as.data.frame() %>% 
  pivot_longer(-c("CPD_ID"), names_to = "sample", values_to = "log2intensity")

BCAA_non_target_cpd_long <- 
  BCAA_non_target_Derivative_cpd_long[BCAA_non_target_Derivative_cpd_long$CPD_ID %in% c("pos_1388", "neg_2285", "pos_24501"), ]

BCAA_non_target_cpd_long$group <- c(rep(c(rep("Wild", 10), rep("Lab", 10), rep("Fam", 10), rep("Res", 10), rep("Hos", 10)), 3))

BCAA_non_target_cpd_long$group <- factor(BCAA_non_target_cpd_long$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

BCAA_non_target_cpd_long$CPD_ID <- factor(BCAA_non_target_cpd_long$CPD_ID, levels = c("pos_1388", "neg_2285", "pos_24501"))

BCAA_non_target_cpd_long_grouped <-
  BCAA_non_target_cpd_long %>% 
  mutate(MS2_name = case_when(
    CPD_ID == "pos_1388" ~ "Isoleucine",
    CPD_ID == "neg_2285" ~ "Leucine",
    CPD_ID == "pos_24501" ~ "Valine"
  )) %>%
  group_by(group, MS2_name) %>%
  summarise(across(c(log2intensity), mean, na.rm = TRUE), .groups = "drop") %>%
  group_by(MS2_name) %>%
  mutate(across(c(log2intensity), ~ scales::rescale(.x, to = c(-2, 2)))) %>%
  mutate(group2 = factor(ifelse(group == "Wild", "Wild", "Non-wild"), levels = c("Wild", "Non-wild"))) %>%
  group_by(group2, MS2_name) %>%
  summarise(across(c(log2intensity), mean, na.rm = TRUE), .groups = "drop")

BCAA_non_target_cpd_heatmap <- 
  ggplot(BCAA_non_target_cpd_long_grouped, 
         aes(x = group2, y = MS2_name, fill = log2intensity),
  ) +
  geom_tile(color = "white",
            height = 1,
            linewidth = 0.45,
            lineend = "square",
            show.legend = T
  ) +
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF", 
                       limit = c(-2, 2),
                       breaks = c(-2, 2)
  ) +
  scale_x_discrete(position = "bottom") +
  labs(fill = "Scale value") +
  theme_void() +
  theme(
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 45, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

BCAA_non_target_cpd_heatmap

ggplot2::ggsave("1_multi_omic_association/43_BCAA_non_target_cpd_heatmap.svg", plot = BCAA_non_target_cpd_heatmap , width = 2.2, height = 3, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/43_BCAA_non_target_cpd_heatmap.tiff", plot = BCAA_non_target_cpd_heatmap , width = 2.2, height = 3, units = "cm", dpi = 300)


### key metabolome cpds box
monochromeR::generate_palette("#F27127FF", modification = "go_lighter", n_colors = 5, view_palette = TRUE, view_labels = FALSE)

BCAA_non_target_Derivative_cpd <- read.table("1_multi_omic_association/42_metabolome_LOESS_PQN_missForest_matrix_noQC.tsv", header = T, row.names = 1, sep = "\t")

BCAA_non_target_Derivative_cpd_long <- 
  BCAA_non_target_Derivative_cpd %>% 
  rownames_to_column(var = "CPD_ID") %>% 
  as.data.frame() %>% 
  pivot_longer(-c("CPD_ID"), names_to = "sample", values_to = "log2intensity")

BCAA_non_target_cpd_long <- 
  BCAA_non_target_Derivative_cpd_long[BCAA_non_target_Derivative_cpd_long$CPD_ID %in% c("pos_1388", "neg_2285", "pos_24501"), ]

BCAA_non_target_cpd_long$group <- c(rep(c(rep("Wild", 10), rep("Lab", 10), rep("Fam", 10), rep("Res", 10), rep("Hos", 10)), 3))

BCAA_non_target_cpd_long$group <- factor(BCAA_non_target_cpd_long$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

BCAA_non_target_cpd_long$CPD_ID <- factor(BCAA_non_target_cpd_long$CPD_ID, levels = c("pos_1388", "neg_2285", "pos_24501"))

BCAA_non_target_cpd_long <-
  BCAA_non_target_cpd_long %>% 
  mutate(MS2_name = case_when(
    CPD_ID == "pos_1388" ~ "L-Leucine",
    CPD_ID == "neg_2285" ~ "L-Isoleucine",
    CPD_ID == "pos_24501" ~ "L-Valine"
  ))

BCAA_non_target_cpd_long$MS2_name <- factor(BCAA_non_target_cpd_long$MS2_name, levels = c("L-Leucine", "L-Isoleucine", "L-Valine"))

BCAA_non_target_cpd_long$top_label <- "BCAA Untargeted Metabolome"

BCAA_non_target_cpd_box_plot <-
  ggplot(BCAA_non_target_cpd_long,
         aes(x = group, y = log2intensity)
  ) +
  geom_point(aes(fill = group, color = group),
             stroke = 0.1,
             size = 1,
             shape = 21,
             position = position_jitter(height = 0, width = 0.2, seed = 123),
             stat = "identity",
             show.legend = F
  ) +
  geom_boxplot(aes(fill = group, color = group),
               # color = "black",
               outlier.shape = NA, 
               alpha = 0.3,
               width = 0.85,
               staplewidth = 0.4,
               linewidth = 0.1,
               show.legend = F
  ) +
  facet_nested_wrap(~ top_label + MS2_name, 
                    scales = "free_y",
                    ncol = 3
  ) + 
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  stat_compare_means(method = "anova",
                     # method = "kruskal.test",
                     size = 2
  ) +
  stat_compare_means(
    comparisons = list(c("Wild", "Lab"), c("Wild", "Fam"), c("Wild", "Res"), c("Wild", "Hos")),
    # method = "wilcox.test",
    method = "t.test",
    # label = "p.adj",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
    # p.adjust.method = "BH", 
    p.adjust.method = "holm", 
    # method.args = list(exact = FALSE, correct = TRUE),
    size = 3,
    vjust = 0.67,
    bracket.size = 0.25,
    step.increase = 0.075,
    tip.length = 0
  ) +
  labs(y = "log10(intensity + 1)") + 
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=1, angle=45, lineheight=1, color="black"),  
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = ggfun::element_roundrect(
      # fill = "#DEEFF5",
      fill = "#c9caca", 
      color = NA, 
      linewidth = 0, 
      r = unit(0.25, "snpc")
    ),
    strip.text = element_text(family = "Arial", face = "bold", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black")
  )

BCAA_non_target_cpd_box_plot

ggplot2::ggsave("1_multi_omic_association/44_BCAA_non_target_cpd_box_plot.svg", plot = BCAA_non_target_cpd_box_plot, width = 8.25, height = 4.25, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/44_BCAA_non_target_cpd_box_plot.tiff", plot = BCAA_non_target_cpd_box_plot, width = 8.25, height = 4.25, units = "cm", dpi = 300)




### MetaB-M144 cluster and heatmap 
library(tidyverse)
library(Mfuzz)

MetaB_M144_CPD <- 
  read.table("../2.metabolome/2_WGCNA/merge_module/merge_module_gene.txt", header = T) %>%
  filter(module_ID == "MetaB_M144")

MetaB_M144_CPD_metadata <- read_xlsx("../2.metabolome/1_metabolome_matrix/metabolies_metadata.xlsx", skip = 15) %>%
  filter(ID %in% MetaB_M144_CPD$metabolome_name)

Sample_metadata <- read.table("../2.metabolome/1_metabolome_matrix/group.txt", header = T)

MetaB_M144_CPD_normalize_intensity <- 
  read.table("../2.metabolome/1_metabolome_matrix/metabolome_LOESS_PQN_missForest_matrix_noQC.tsv", header = T) %>%
  filter(ID %in% MetaB_M144_CPD$metabolome_name) %>%
  column_to_rownames("ID") %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Compound_ID") %>%
  pivot_longer(cols = -Compound_ID, names_to = "sample", values_to = "Z_Score") %>%
  left_join(Sample_metadata, by = "sample")
  
MetaB_M144_CPD_normalize_intensity$group <- factor(MetaB_M144_CPD_normalize_intensity$group , levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

MetaB_M144_CPD_cluter_Trend <- 
  ggplot(MetaB_M144_CPD_normalize_intensity, 
         aes(x = group, y = Z_Score)
  ) +
  stat_summary(aes(group = 1), 
               fun.data = mean_sdl, 
               fun.args = list(mult = 1), 
               geom = "ribbon", 
               fill = "lightblue", 
               alpha = 0.5
  ) +
  stat_summary(aes(group = 1), 
               fun = mean, 
               geom = "line", 
               color = "black", 
               linewidth = 0.5
  ) +
  stat_summary(aes(fill = group, color = group), 
               fun = mean, 
               geom = "point", 
               shape = 21,
               size = 1.5, 
               stroke = 0.5,
               show.legend = F
  ) +
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) + 
  labs(y = "Intsnsity Z-score",   
       x = ""
  ) + 
  theme_minimal() + 
  theme(
    axis.title.x = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.title.y = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=1, angle=45, lineheight=1, color="black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x = unit(2, "pt"),
    axis.ticks.length.y = unit(2, "pt"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank()
  )

print(MetaB_M144_CPD_cluter_Trend)

ggplot2::ggsave("1_multi_omic_association/45_MetaB_M144_CPD_cluter_Trend.svg", plot = MetaB_M144_CPD_cluter_Trend, width = 2.5, height = 3.5, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/45_MetaB_M144_CPD_cluter_Trend.tiff", plot = MetaB_M144_CPD_cluter_Trend, width = 2.5, height = 3.5, units = "cm", dpi = 300)


### heatmap
library(ComplexHeatmap)
library(circlize)

heatmap_mat <- 
  MetaB_M144_CPD_normalize_intensity %>%
  group_by(Compound_ID, group) %>%
  summarise(mean_Z = mean(Z_Score, na.rm = TRUE)) %>%
  mutate(zscore = scales::rescale(mean_Z, to = c(-2, 2))) %>%
  pivot_wider(id_cols = Compound_ID, names_from = group, values_from = zscore) %>%
  column_to_rownames("Compound_ID") %>%
  as.matrix()

col_fun = colorRamp2(c(-2, 0, 2), c("lightblue", "white", "pink"))

col_dist = dist(t(heatmap_mat))
col_hc = hclust(col_dist)
col_dend = as.dendrogram(col_hc)

select_cpds <- c("pos_1388", "neg_2285", "pos_24501")
row_idx <- which(rownames(heatmap_mat) %in% select_cpds)
labels <- c("L-Leucine", "L-Isoleucin", "L-Valine")

right_anno <- rowAnnotation(
  mark = anno_mark(at = row_idx, 
                   labels = labels,
                   labels_gp = gpar(fontsize = 6)))

group_colors = c(
  "Wild" = wild_color, 
  "Lab" = lab_color, 
  "Fam" = fam_color, 
  "Res" = res_color, 
  "Hos" = hos_color
)

group_info = colnames(heatmap_mat)

col_anno = columnAnnotation(
  Group = group_info,
  col = list(Group = group_colors),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  height = unit(1, "mm")
)


bottom_text_anno <- HeatmapAnnotation(
  text = anno_text(c("Wild", "Lab", "Fam", "Res", "Hos"),
                   location = 0.5,
                   just = "center",
                   rot = 45,
                   gp = gpar(fontsize = 6))
  )

MetaB_M144_CPD_heatmap <- 
  Heatmap(
  heatmap_mat, 
  name = "Z-score",
  col = col_fun, 
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  row_dend_gp = gpar(lwd = 0.1),
  column_dend_gp = gpar(lwd = 0.1),
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 6),
  bottom_annotation =  c(col_anno, bottom_text_anno),
  show_column_names = FALSE,
  column_names_side = "bottom",
  column_names_rot = 45,
  rect_gp = gpar(col = NA),
  right_annotation = right_anno,
  border = TRUE,
  border_gp = gpar(lwd = 0.1, col = "black"),
  show_heatmap_legend = TRUE
)

draw(MetaB_M144_CPD_heatmap) 
dev.off()

pdf("1_multi_omic_association/46_MetaB_M144_Heatmap.pdf", width = 3, height = 3)
draw(MetaB_M144_CPD_heatmap) 
dev.off()





### BCAA untrageted metabolomic bubble plot
BCAA_non_target_cpd_long_grouped$top_label <- "BCAA Untargeted Metabolome"

BCAA_non_target_cpd_long_grouped$sec_label <- "Scaled Intensity"

BCAA_non_target_cpd_bubble_plot <- 
  ggplot(BCAA_non_target_cpd_long_grouped,
         aes(x = group, y = MS2_name)
  ) +
  geom_point(aes(fill = group, color = group, size = log2intensity),
             stroke = 0.1,
             shape = 21,
             stat = "identity",
             show.legend = F
  ) +
  facet_nested_wrap(~ top_label + sec_label, 
                    scales = "free_y",
                    ncol = 3
  ) + 
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_size_continuous(range = c(1, 4)) +
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=1, angle=45, lineheight=1, color="black"),  
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = ggfun::element_roundrect(
      # fill = "#DEEFF5",
      fill = "#c9caca", 
      color = NA, 
      linewidth = 0, 
      r = unit(0.25, "snpc")
    ),
    strip.text = element_text(family = "Arial", face = "bold", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black")
  )

BCAA_non_target_cpd_bubble_plot

ggplot2::ggsave("1_multi_omic_association/26_BCAA_non_target_cpd_bubble_plot.svg", plot = BCAA_non_target_cpd_bubble_plot, width = 3, height = 4, units = "cm", dpi = 300)
ggplot2::ggsave("1_multi_omic_association/26_BCAA_non_target_cpd_bubble_plot.tiff", plot = BCAA_non_target_cpd_bubble_plot, width = 3, height = 4, units = "cm", dpi = 300)


### target metabolome heatmap
BCAA_target_Derivative_cpd <- read.table("1_multi_omic_association/47_BCAA_target_metabolome.txt", header = T)

BCAA_target_Derivative_cpd$Metabolite <- factor(BCAA_target_Derivative_cpd$Metabolite, levels = c("L-Leucine", "L-Isoleucine", "L-Valine"))

BCAA_target_Derivative_cpd_long <- 
  BCAA_target_Derivative_cpd %>%
  pivot_longer(-Metabolite, names_to = "sample", values_to = "concentration") %>%
  mutate(
    group = factor(rep(c(rep("Lab", 4), rep("Res", 4), rep("Fam", 4), rep("Hos", 4), rep("Wild", 4)), 3), levels = c("Wild", "Lab", "Fam", "Res", "Hos"))
  ) %>% 
  mutate(log10concentration = log10(concentration))

BCAA_target_Derivative_cpd_long_grouped <-
  BCAA_target_Derivative_cpd_long %>%
  group_by(Metabolite, group) %>%
  summarise(across(c(concentration), mean, na.rm = TRUE), .groups = "drop") %>%
  group_by(Metabolite) %>%
  mutate(across(c(concentration), ~ scales::rescale(.x, to = c(-2, 2)))) %>%
  mutate(group2 = factor(ifelse(group == "Wild", "Wild", "Non-wild"), levels = c("Wild", "Non-wild"))) %>%
  group_by(group2, Metabolite) %>%
  summarise(across(c(concentration), mean, na.rm = TRUE), .groups = "drop")

BCAA_target_cpd_heatmap <- 
  ggplot(BCAA_target_Derivative_cpd_long_grouped, 
         aes(x = group2, y = Metabolite, fill = concentration),
  ) +
  geom_tile(color = "white",
            height = 1,
            linewidth = 0.45,
            lineend = "square",
            show.legend = T
  ) +
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF", 
                       limit = c(-2, 2),
                       breaks = c(-2, 2)
  ) +
  scale_x_discrete(position = "bottom") +
  labs(fill = "Scale value") +
  theme_void() +
  theme(
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 45, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

BCAA_target_cpd_heatmap

ggplot2::ggsave("1_multi_omic_association/48_BCAA_target_cpd_heatmap.svg", plot = BCAA_target_cpd_heatmap , width = 2.2, height = 3, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/48_BCAA_target_cpd_heatmap.tiff", plot = BCAA_target_cpd_heatmap , width = 2.2, height = 3, units = "cm", dpi = 300)


### target metabolome boxplot
library(scales)

BCAA_target_Derivative_cpd_long$top_label <- "Targeted Metabolomic Profiling of Luminal BCAAs"

BCAA_target_cpd_box_plot <-
  ggplot(BCAA_target_Derivative_cpd_long,
         aes(x = group, y = concentration)
  ) +
  geom_point(aes(fill = group, color = group),
             stroke = 0.1,
             size = 1,
             shape = 21,
             position = position_jitter(height = 0, width = 0.2, seed = 123),
             stat = "identity",
             show.legend = F
  ) +
  geom_boxplot(aes(fill = group, color = group),
               outlier.shape = NA, 
               alpha = 0.3,
               width = 0.85,
               staplewidth = 0.4,
               linewidth = 0.1,
               show.legend = F
  ) +
  facet_nested_wrap(~ top_label + Metabolite, 
             scales = "free_y",
             ncol = 3
  ) + 
  scale_y_log10(
    breaks = function(x) 10^seq(floor(log10(min(x))), ceiling(log10(max(x)))),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) + 
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  stat_compare_means(method = "anova",
                     # method = "kruskal.test",
                     size = 2
  ) +
  stat_compare_means(
    comparisons = list(c("Wild", "Lab"), c("Wild", "Fam"), c("Wild", "Res"), c("Wild", "Hos")),
    # method = "wilcox.test",
    method = "t.test",
    # label = "p.adj",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
    # p.adjust.method = "BH", 
    p.adjust.method = "holm", 
    # method.args = list(exact = FALSE, correct = TRUE),
    size = 3,
    vjust = 0.67,
    bracket.size = 0.25,
    step.increase = 0.075,
    tip.length = 0
  ) +
  labs(y = "log10(Intensity + 1)") + 
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=1, angle=45, lineheight=1, color="black"),  
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = ggfun::element_roundrect(
      # fill = "#DEEFF5",
      fill = "#c9caca", 
      color = NA, 
      linewidth = 0, 
      r = unit(0.25, "snpc")
    ),
    strip.text = element_text(family = "Arial", face = "bold", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black")
  )

BCAA_target_cpd_box_plot

ggplot2::ggsave("1_multi_omic_association/49_BCAA_target_cpd_box_plot.svg", plot = BCAA_target_cpd_box_plot, width = 8.25, height = 4.25, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/49_BCAA_target_cpd_box_plot.tiff", plot = BCAA_target_cpd_box_plot, width = 8.25, height = 4.25, units = "cm", dpi = 300)



### Transcriptome key module & OG 

### transcriptome pathway heatmap
HostT_mTOR_relate_pathway <- read.table("../3.transcriptome/1_OG/Total_transcriptome_KEGG_ssGSEA_pathway.tsv", header = T, row.names = 1) %>%
  t() %>% 
  as.data.frame() %>%
  dplyr::select(c("ko04150", "ko03050", "ko04122", "ko00730", "ko04214", "ko03010", "ko00190", "ko04068", "ko04392")) %>%  
  as.matrix()

HostT_mTOR_relate_data <- 
  HostT_mTOR_relate_pathway %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  mutate(group = c(rep("Wild", 10), rep("Lab", 10), rep("Fam", 10), rep("Res", 10), rep("Hos", 10))) %>%
  dplyr::select(-sample) %>%
  group_by(group) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  mutate(across(c(-group), ~ scales::rescale(.x, to = c(-2, 2)))) %>%
  mutate(group2 = factor(ifelse(group == "Wild", "Wild", "Non-wild"), levels = c("Wild", "Non-wild"))) %>%
  group_by(group2) %>%
  summarise(across(c(-group), mean, na.rm = TRUE))

HostT_mTOR_relate_data_long <- 
  HostT_mTOR_relate_data %>%
  pivot_longer(-c("group2"), names_to = "pathway", values_to = "enrich score")

HostT_mTOR_relate_data_long$pathway <- factor(HostT_mTOR_relate_data_long$pathway, levels = c("ko04150", "ko04068", "ko04392", "ko04122", "ko00730", "ko04214" ,"ko03050", "ko03010", "ko00190"))

min <- min(HostT_mTOR_relate_data_long$`enrich score`)
max < max(HostT_mTOR_relate_data_long$`enrich score`)

HostT_mTOR_relate_pathway_heatmap <- 
  ggplot(HostT_mTOR_relate_data_long) +
  geom_tile(aes(x = group2, y = pathway, fill = `enrich score`),
            color = "white",
            height = 1,
            linewidth = 0.45,
            show.legend = T
  ) + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF",
                       limit = c(-2, 2),
                       breaks = c(-2, 2)
  ) +
  labs(fill = "Scale value") + 
  theme_void() +
  theme(
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "top",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

print(HostT_mTOR_relate_pathway_heatmap)

ggplot2::ggsave("../4.association/1_multi_omic_association/50_HostT_mTOR_relate_pathway_heatmap.svg", plot = HostT_mTOR_relate_pathway_heatmap, width = 2.2, height = 5, units = "cm", dpi = 300)
# ggplot2::ggsave("../4.association/1_multi_omic_association/50_HostT_mTOR_relate_pathway_heatmap.tiff", plot = HostT_mTOR_relate_pathway_heatmap, width = 2.2, height = 5, units = "cm", dpi = 300)


### key transcriptome pathway box plot
monochromeR::generate_palette("#6DBC90FF", modification = "go_lighter", n_colors = 5, view_palette = TRUE, view_labels = FALSE)

key_transcriptome_pathway_data <- 
  read.table("../3.transcriptome/1_OG/Total_transcriptome_KEGG_ssGSEA_pathway.tsv", header = T, row.names = 1) %>% 
  t() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  dplyr::select(c("sample", "ko04150", "ko00190", "ko03010", "ko03050", "ko04392", "ko04068", "ko04122", "ko04214", "ko00730")) %>%
  mutate(group = c(rep("Wild", 10), rep("Lab", 10), rep("Fam", 10), rep("Res", 10), rep("Hos", 10))) %>% 
  pivot_longer(-c("sample", "group"), names_to = "pathway", values_to = "enrich score")
  
key_transcriptome_pathway_data$group <- factor(key_transcriptome_pathway_data$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))
key_transcriptome_pathway_data$pathway <- factor(key_transcriptome_pathway_data$pathway, levels = c("ko04150", "ko00190", "ko04214", "ko04122", "ko00730", "ko03010", "ko03050", "ko04392", "ko04068"))

pathway_merge_name = c("ko04150" = "ko04150:mTOR", 
                       "ko00190" = "ko00190:OXPHOS",
                       "ko04122" = "ko04122:Sulfur relay system",
                       "ko00730" = "ko00730:Thiamine metabolism",
                       "ko04214" = "ko04214:Apoptosis",
                       "ko03010" = "ko03010:Ribosome", 
                       "ko03050" = "ko03050:Proteasome",
                       "ko04392" = "ko04392:Hippo",
                       "ko04068" = "ko04068:FoxO",
                       "Host Transcriptomic Pathways" = "Host Transcriptomic Pathways")

key_transcriptome_pathway_data$top_label <- "Host Transcriptomic Pathways"

key_transcriptome_pathway_data <- 
  key_transcriptome_pathway_data %>% 
  filter(pathway %in% c("ko04150", "ko00190", "ko04122", "ko04214"))

key_transcriptome_pathway_box_plot <-
  ggplot(key_transcriptome_pathway_data,
         aes(x = group, y = `enrich score`)
  ) +
  geom_point(aes(fill = group, color = group),
             stroke = 0.1,
             size = 1,
             shape = 21,
             stat = "identity",
             position = position_jitter(height = 0, width = 0.2, seed = 123),
             show.legend = F
  ) +
  geom_boxplot(aes(fill = group, color = group),
               # color = "black",
               outlier.shape = NA, 
               alpha = 0.3,
               width = 0.85,
               staplewidth = 0.4,
               linewidth = 0.1,
               show.legend = F
  ) +
  facet_nested_wrap(~ top_label+ pathway, 
             scales = "free_y",
             ncol = 2,
             nrow = 2 ,
             # labeller = as_labeller(pathway_merge_name)
  ) + 
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  stat_compare_means(method = "anova",
                     # method = "kruskal.test",
                     size = 2
  ) +
  stat_compare_means(
    comparisons = list(c("Wild", "Lab"), c("Wild", "Fam"), c("Wild", "Res"), c("Wild", "Hos")),
    # method = "wilcox.test",
    method = "t.test",
    # label = "p.adj",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
    # p.adjust.method = "BH", 
    p.adjust.method = "holm", 
    # method.args = list(exact = FALSE, correct = TRUE),
    size = 3,
    vjust = 0.67,
    bracket.size = 0.25,
    step.increase = 0.075,
    tip.length = 0
  ) +
  labs(y = "Enrichment score") + 
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=1, angle=45, lineheight=1, color="black"),  
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    panel.border =  element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = ggfun::element_roundrect(
      # fill = "#DEEFF5",
      fill = "#c9caca", 
      color = NA, 
      linewidth = 0, 
      r = unit(0.25, "snpc")
    ),
    strip.text = element_text(family = "Arial", face = "bold", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black")
  )

key_transcriptome_pathway_box_plot

ggplot2::ggsave("1_multi_omic_association/51_key_transcriptome_pathway_box.svg", plot = key_transcriptome_pathway_box_plot, width = 3, height = 4.25, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/51_key_transcriptome_pathway_box.tiff", plot = key_transcriptome_pathway_box_plot, width = 3, height = 4.25, units = "cm", dpi = 300)


### Transcriptome key gene

### mTOR relate gene expression heatmap
Gene_ID_2_KO <- read.table("../3.transcriptome/2_downstream_analysis/Gene_2_KO.txt", header = T, sep = "\t")

mTOR_relate_gene_exprission_data <- 
  read.table("../3.transcriptome/2_downstream_analysis/1_transcript_tpm_filt.tsv", header = T, sep = "\t") %>% 
  left_join(Gene_ID_2_KO, by = c("gene" = "gene")) %>%
  filter(KO %in% c("K07203", "K04688", "K03258", "K04456", "K07205", "K14995")) %>% 
  mutate(symbol = case_when(
    KO == "K07203" ~ "mTOR",
    KO == "K04688" ~ "S6K",
    KO == "K03258" ~ "eIF4B",
    KO == "K04456" ~ "AKT",
    KO == "K07205" ~ "4EBP",
    KO == "K14995" ~ "SLC38A9"
  )) %>%
  group_by(KO) %>%
  summarise(across(-c(gene, symbol), mean, na.rm = TRUE)) %>%
  dplyr::select(KO, everything()) %>%
  pivot_longer(-c(KO), names_to = "sample", values_to = "value") %>%
  mutate(symbol = case_when(
    KO == "K07203" ~ "mTOR",
    KO == "K04688" ~ "S6K",
    KO == "K03258" ~ "eIF4B",
    KO == "K04456" ~ "AKT",
    KO == "K07205" ~ "4EBP",
    KO == "K14995" ~ "SLC38A9"
  )) %>%
  filter(symbol %in% c("mTOR", "eIF4B", "SLC38A9")) %>%
  mutate(group = rep(c(rep("Wild", 10), rep("Res", 10), rep("Fam", 10), rep("Hos", 10), rep("Lab", 10)), 3)) %>%
  dplyr::select(sample, group, symbol, value)

key_transcriptome_pathway_data_new <- key_transcriptome_pathway_data %>% dplyr::select(sample, group, pathway, `enrich score`)
mTOR_relate_gene_exprission_data$group <- factor(mTOR_relate_gene_exprission_data$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))
mTOR_relate_gene_exprission_data$symbol <- factor(mTOR_relate_gene_exprission_data$symbol, levels = c("SLC38A9", "mTOR", "S6K", "eIF4B", 'AKT', "4EBP", "Survivin", "HK", "PFK", "LDH", "GLUT1"))
mTOR_relate_gene_exprission_data$top_label <- "Gut BCAAs-responsive Pathway and Gene (Transcriptome)"

mTOR_relate_gene_exprission_box_plot <-
  ggplot(mTOR_relate_gene_exprission_data,
         aes(x = group, y = value)
  ) +
  geom_point(aes(fill = group, color = group),
             stroke = 0.1,
             size = 1,
             shape = 21,
             stat = "identity",
             position = position_jitter(height = 0, width = 0.2, seed = 123),
             show.legend = F
  ) +
  geom_boxplot(aes(fill = group, color = group),
               outlier.shape = NA, 
               alpha = 0.3,
               width = 0.85,
               staplewidth = 0.4,
               linewidth = 0.1,
               show.legend = F
  ) +
  facet_nested_wrap(~ top_label + symbol, 
                    scales = "free_y",
                    ncol = 3,
                    # labeller = as_labeller(pathway_merge_name)
  ) + 
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  stat_compare_means(method = "anova",
                     # method = "kruskal.test",
                     size = 2
  ) +
  stat_compare_means(
    comparisons = list(c("Wild", "Lab"), c("Wild", "Fam"), c("Wild", "Res"), c("Wild", "Hos")),
    # method = "wilcox.test",
    method = "t.test",
    # label = "p.adj",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
    # p.adjust.method = "BH", 
    p.adjust.method = "holm", 
    # method.args = list(exact = FALSE, correct = TRUE),
    size = 3,
    vjust = 0.67,
    bracket.size = 0.25,
    step.increase = 0.075,
    tip.length = 0
  ) +
  labs(y = "Enrichment score") + 
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=1, angle=45, lineheight=1, color="black"),  
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = ggfun::element_roundrect(
      # fill = "#DEEFF5",
      fill = "#c9caca", 
      color = NA, 
      linewidth = 0, 
      r = unit(0.25, "snpc")
    ),
    strip.text = element_text(family = "Arial", face = "bold", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black")
  )

mTOR_relate_gene_exprission_box_plot

ggplot2::ggsave("1_multi_omic_association/52_mTOR_relate_gene_exprission_box_plot.svg", plot = mTOR_relate_gene_exprission_box_plot, width = 8.25, height = 4.25, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/52_mTOR_relate_gene_exprission_box_plot.tiff", plot = mTOR_relate_gene_exprission_box_plot, width = 8.25, height = 4.25, units = "cm", dpi = 300)


### mTOR expression
mTOR_exprission_data <- 
  read.table("1_multi_omic_association/53_mTOR expression data.txt", header = T, sep = "\t") 

scale_data <- 
  mTOR_exprission_data %>% 
  dplyr::select(-c(1:5)) %>%
  t() %>%
  as.data.frame()

scale_data$V1 <- scales::rescale(scale_data$V1, to = c(-2, 2))
scale_data$V2 <- scales::rescale(scale_data$V2, to = c(-2, 2))

mTOR_exprission_scale_data <- mTOR_exprission_data

mTOR_exprission_scale_data[,6:29] <- t(scale_data)

mTOR_exprission_data_long_scale <- 
  mTOR_exprission_scale_data %>% 
  pivot_longer(-c(1:5), names_to = "sample", values_to = "value") 

mTOR_exprission_data_long_scale$group <- factor(rep(c(rep("Pame", 6), rep("Pful", 6), rep("Esin", 6), rep("Bger", 6)), 2), levels = c("Esin", "Bger", "Pame", "Pful"))
mTOR_exprission_data_long_scale$assay.method <- factor(mTOR_exprission_data_long_scale$assay.method, levels = unique(mTOR_exprission_data_long_scale$assay.method))

mTOR_exprission_data_long_grouped <- 
  mTOR_exprission_data %>%
  pivot_longer(-c(1:5), names_to = "sample", values_to = "value") 

mTOR_exprission_data_long_grouped$group <- factor(rep(c(rep("Pame", 6), rep("Pful", 6), rep("Esin", 6), rep("Bger", 6)), 2), levels = c("Esin", "Bger", "Pame", "Pful"))
mTOR_exprission_data_long_grouped$group2 <- factor(ifelse(mTOR_exprission_data_long_grouped$group == "Esin", "Esin", "Other"), levels = c("Esin", "Other"))
mTOR_exprission_data_long_grouped$assay.method <- factor(mTOR_exprission_data_long_grouped$assay.method, levels = unique(mTOR_exprission_data_long_grouped$assay.method))

mTOR_exprission_data_long_grouped_scale <-
  mTOR_exprission_data_long_grouped %>%
  group_by(group, assay.method) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    n = n(),
    se_value = sd_value / sqrt(n),
    .groups = "drop"
  ) %>% 
  group_by(assay.method) %>%
  mutate(across(c(mean_value), ~ scales::rescale(.x, to = c(-2, 2)))) %>%
  mutate(group2 = ifelse(group == "Esin", "Esin", "Other")) %>%
  group_by(group2, assay.method) %>% 
  summarise(across(c(mean_value), mean, na.rm = TRUE), .groups = "drop")

### mTOR expression heatmap
mTOR_exprission_heatmap <- 
  ggplot(mTOR_exprission_data_long_grouped_scale) +
  geom_tile(aes(x = group2, y = assay.method, fill = mean_value),
            color = "white",
            height = 1,
            linewidth = 0.45,
            show.legend = T
  ) + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF",
                       limit = c(-2, 2),
                       breaks = c(-2, 2)
  ) +
  labs(fill = "Scaled value") + 
  theme_void() +
  theme(
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

mTOR_exprission_heatmap

ggplot2::ggsave("1_multi_omic_association/54_mTOR_exprission_heatmap.svg", plot = mTOR_exprission_heatmap, width = 10, height = 2, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/54_mTOR_exprission_heatmap.tiff", plot = mTOR_exprission_heatmap, width = 10, height = 2, units = "cm", dpi = 300)


###  mTOR expression bar/box plot
mTOR_exprission_data_long <- 
  mTOR_exprission_data %>%
  pivot_longer(-c(1:5), names_to = "sample", values_to = "value") %>%
  mutate(value_trans = case_when(
    assay.method == "RNA seq log2(TPM + 1)" ~ log2(value + 1),
    TRUE ~ value
  ))

mTOR_exprission_data_long$group <- factor(rep(c(rep("Pame", 6), rep("Pful", 6), rep("Esin", 6), rep("Bger", 6)), 2), levels = c("Esin", "Bger", "Pame", "Pful"))
mTOR_exprission_data_long$assay.method <- factor(mTOR_exprission_data_long$assay.method, levels = unique(mTOR_exprission_data_long$assay.method)) 

mTOR_exprission_data_long_grouped_scale$top_label <- "mTOR mRNA Expression Levels"
mTOR_exprission_data_long$top_label <- "mTOR mRNA Expression Levels"

mTOR_expression_box_plot <- 
  ggplot(mTOR_exprission_data_long,
         aes(x = group, y = value, )
  ) +
  geom_point(aes(fill = group, color = group),
             stroke = 0.1,
             size = 1,
             shape = 21,
             stat = "identity",
             position = position_jitter(height = 0, width = 0.2, seed = 123),
             show.legend = F
  ) +
  geom_boxplot(aes(fill = group, color = group),
               outlier.shape = NA, 
               alpha = 0.3,
               width = 0.85,
               staplewidth = 0.4,
               linewidth = 0.1,
               show.legend = F
  ) +
  ggh4x::facet_nested_wrap(~ top_label + assay.method, 
                           scales = "free_y",
                           ncol = 2
  ) + 
  stat_compare_means(method = "anova",
                     # method = "kruskal.test",
                     size = 2
  ) +
  stat_compare_means(
    comparisons = list(c("Esin", "Bger"), c("Esin", "Pame"), c("Esin", "Pful"), 
                       c("Bger", "Pame"), c("Bger", "Pful"), c("Pame", "Pful")),
    # method = "wilcox.test",
    method = "t.test",
    # label = "p.adj",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
    # p.adjust.method = "BH", 
    p.adjust.method = "holm", 
    # method.args = list(exact = FALSE, correct = TRUE),
    size = 3,
    vjust = 0.67,
    bracket.size = 0.25,
    step.increase = 0.075,
    tip.length = 0
  ) +
  scale_fill_manual(values = c("Bger" = Bger_color, "Esin" = Esin_color, "Pame" = Pame_color, "Pful" = Pful_color)) +
  scale_color_manual(values = c("Bger" = Bger_color, "Esin" = Esin_color, "Pame" = Pame_color, "Pful" = Pful_color)) +
  # scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=1, angle=45, lineheight=1, color="black"),  
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "square", color = "black"),
    axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "square", color = "black"),
    panel.border =  element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = ggfun::element_roundrect(
      # fill = "#DEEFF5",
      fill = "#c9caca", 
      color = NA, 
      linewidth = 0, 
      r = unit(0.25, "snpc")
    ),
    strip.text = element_text(family = "Arial", face = "bold", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black")
  )

mTOR_expression_box_plot

ggplot2::ggsave("1_multi_omic_association/55_mTOR_expression_box_plot.svg", plot = mTOR_expression_box_plot, width = 5.5, height = 4.5, units = "cm", dpi = 1200)
# ggplot2::ggsave("1_multi_omic_association/55_mTOR_expression_box_plot.tiff", plot = mTOR_expression_box_plot, width = 5.5, height = 4.5, units = "cm", dpi = 1200)


### Hemolymph and Fat Body BCAAs Levels (WST-8)
# BCAA Assay(Reagent kit for total BCAA)
library(ggplot2)
library(ggpmisc)

Hemolymph_FatBody_BCAAs <- 
  read.table("1_multi_omic_association/56_Hemolymph_FatBody_BCAA.txt", header = T, sep = "\t") %>%
  mutate(relative_concentrate = BCAA_concentrate / Volume_ul_Weight_mg)

Hemolymph_FatBody_BCAAs$group <- factor(Hemolymph_FatBody_BCAAs$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

Hemolymph_FatBody_BCAAs$Tissue <- factor(Hemolymph_FatBody_BCAAs$Tissue, levels = c("Hemolymph", "Fat_Body"))

Hemolymph_FatBody_BCAAs$top_label <- "Hemolymph and Fat Body BCAAs Levels (WST-8)"

### Hemolymp and Fat Body BCAAs box plot 

Hemolymph_FatBody_BCAA_box_plot <-
  ggplot(data = Hemolymph_FatBody_BCAAs,
         aes(x = group , y = relative_concentrate)
  ) + 
  geom_point(aes(fill = group, color = group),
             stroke = 0.1,
             size = 1,
             shape = 21,
             stat = "identity",
             position = position_jitter(height = 0, width = 0.2, seed = 123),
             show.legend = F
  ) +
  geom_boxplot(aes(fill = group, color = group),
               outlier.shape = NA, 
               alpha = 0.3,
               width = 0.85,
               staplewidth = 0.4,
               linewidth = 0.1,
               show.legend = F
  ) +
  facet_nested_wrap(~ top_label + Tissue, 
                    scales = "free_y",
                    ncol = 2,
  ) + 
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  stat_compare_means(method = "anova",
                     # method = "kruskal.test",
                     size = 2
  ) +
  stat_compare_means(
    comparisons = list(c("Wild", "Lab"), c("Wild", "Fam"), c("Wild", "Res"), c("Wild", "Hos")),
    # method = "wilcox.test",
    method = "t.test",
    # label = "p.adj",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
    # p.adjust.method = "BH", 
    p.adjust.method = "holm", 
    # method.args = list(exact = FALSE, correct = TRUE),
    size = 3,
    vjust = 0.67,
    bracket.size = 0.25,
    step.increase = 0.075,
    tip.length = 0
  ) +
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=1, angle=45, lineheight=1, color="black"),  
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = ggfun::element_roundrect(
      # fill = "#DEEFF5",
      fill = "#c9caca", 
      color = NA, 
      linewidth = 0, 
      r = unit(0.25, "snpc")
    ),
    strip.text = element_text(family = "Arial", face = "bold", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black")
  )

print(Hemolymph_FatBody_BCAA_box_plot)

ggplot2::ggsave("1_multi_omic_association/57_Hemolymph_FatBody_BCAA_box_plot.svg", plot = Hemolymph_FatBody_BCAA_box_plot, width = 5.5, height = 4.25, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/57_Hemolymph_FatBody_BCAA_box_plot.tiff", plot = Hemolymph_FatBody_BCAA_box_plot, width = 5.5, height = 4.25, units = "cm", dpi = 300)


### western bar plot(five environment of pame)
wb_result <- 
  read.table("1_multi_omic_association/58_wb.txt", header = T, row.names = 1) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample")
wb_result$`p-S6K/Tubulin` <- wb_result$`p-S6K` / wb_result$Tubulin
wb_result$`p-AKT/AKT` <- wb_result$`p-AKT` / wb_result$AKT
wb_result$`p-AKT/Tubulin` <- wb_result$`p-AKT` / wb_result$Tubulin
wb_result$group <- factor(c(rep("Wild", 3), rep("Lab", 3), rep("Fam", 3), rep("Res", 3), rep("Hos", 3)), levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

wb_result_long <- 
  wb_result %>%
  dplyr::select(-c(2,3,4,5)) %>%
  pivot_longer(-c(sample, group), names_to = "signal", values_to = "relative intensity")

wb_result_long$signal <- factor(wb_result_long$signal, levels = c("p-S6K/Tubulin", "p-AKT/AKT", "p-AKT/Tubulin")) 
wb_result_long <- 
  wb_result_long %>%
  filter(signal %in% c("p-S6K/Tubulin", "p-AKT/AKT"))

wb_result_long_group <- 
  wb_result_long %>%
  group_by(signal, group) %>%
  summarise(
    mean = mean(`relative intensity`),
    sd = sd(`relative intensity`)
  )

wb_result_long_group$top_label <- "Fat Body mTORC1 Activity Levels"

wb_result_long$top_label <- "Fat Body mTORC1 Activity Levels"

wb_result_long$group <- factor(wb_result_long$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

wb_result_long_group$group <- factor(wb_result_long_group$group, levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

wb_box_plot <- 
  ggplot(wb_result_long,
         aes(x = group, y = `relative intensity`)
  ) +
  geom_point(aes(fill = group, color = group),
             stroke = 0.1,
             size = 1,
             shape = 21,
             stat = "identity",
             position = position_jitter(height = 0, width = 0.2, seed = 123),
             show.legend = F
  ) +
  geom_boxplot(aes(fill = group, color = group),
               outlier.shape = NA, 
               alpha = 0.3,
               width = 0.85,
               staplewidth = 0.4,
               linewidth = 0.1,
               show.legend = F
  ) +
  ggh4x::facet_nested_wrap(~ top_label + signal, 
                           scales = "free_y",
                           ncol = 2
  ) + 
  stat_compare_means(method = "anova",
                     # method = "kruskal.test",
                     size = 2
  ) +
  stat_compare_means(
    comparisons = list(c("Wild", "Lab"), c("Wild", "Fam"), c("Wild", "Res"), c("Wild", "Hos")),
    # method = "wilcox.test",
    method = "t.test",
    # label = "p.adj",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
    # p.adjust.method = "BH", 
    p.adjust.method = "holm", 
    # method.args = list(exact = FALSE, correct = TRUE),
    size = 3,
    vjust = 0.67,
    bracket.size = 0.25,
    step.increase = 0.075,
    tip.length = 0
  ) +
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  labs(y = "relative intensity") + 
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=1, angle=45, lineheight=1, color="black"),  
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "square", color = "black"),
    # axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "square", color = "black"),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = ggfun::element_roundrect(
      fill = "#c9caca", 
      color = NA, 
      linewidth = 0, 
      r = unit(0.25, "snpc")
    ),
    strip.text = element_text(family = "Arial", face = "bold", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black")
  )

wb_box_plot

ggplot2::ggsave("1_multi_omic_association/59_wb_box_plot.svg", plot = wb_box_plot, width = 5.5, height = 4.25, units = "cm", dpi = 1200)
# ggplot2::ggsave("1_multi_omic_association/59_wb_box_plot.tiff", plot = wb_box_plot, width = 5.5, height = 4.25, units = "cm", dpi = 1200)


# heatmap 
wb_result <- 
  read.table("1_multi_omic_association/58_wb.txt", header = T, row.names = 1) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample")

wb_result$`p-S6K/Tubulin` <- wb_result$`p-S6K` / wb_result$Tubulin

wb_result$`p-AKT/AKT` <- wb_result$`p-AKT` / wb_result$AKT

wb_result$`p-AKT/Tubulin` <- wb_result$`p-AKT` / wb_result$Tubulin

wb_result$group <- factor(c(rep("Wild", 3), rep("Lab", 3), rep("Fam", 3), rep("Res", 3), rep("Hos", 3)), levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

wb_result_long <- 
  wb_result %>%
  dplyr::select(-c(2,3,4,5)) %>%
  pivot_longer(-c(sample, group), names_to = "signal", values_to = "relative intensity")

wb_result_long$signal <- factor(wb_result_long$signal, levels = c("p-S6K/Tubulin", "p-AKT/AKT", "p-AKT/Tubulin")) 

wb_result_long <- 
  wb_result_long %>%
  filter(signal %in% c("p-S6K/Tubulin", "p-AKT/AKT"))

wb_result_long_group_scale <- 
  wb_result_long %>%
  group_by(signal, group) %>%
  summarise(across(c(`relative intensity`), list(mean = mean, sd = sd), na.rm = TRUE)) %>%
  group_by(signal) %>% 
  mutate(across(c(`relative intensity_mean`), ~ scales::rescale(.x, to = c(-2, 2)))) %>%
  mutate(group2 = factor(ifelse(group == "Wild", "Wild", "Non-Wild"), levels = c("Wild", "Non-Wild"))) %>% 
  group_by(signal, group2) %>%
  summarise(across(`relative intensity_mean`, mean, na.rm = TRUE)) 

wb_result_long_group_scale$signal <- factor(wb_result_long_group_scale$signal, levels = c("p-AKT/AKT", "p-S6K/Tubulin"))

wb_heatmap <-
  ggplot(wb_result_long_group_scale) +
  geom_tile(aes(x = group2, y = signal, fill = `relative intensity_mean`),
            color = "white",
            height = 1,
            linewidth = 0.45,
            show.legend = T
  ) + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF",
                       limits = c(-2, 2),
                       breaks = c(-2, 0, 2)
  ) +
  labs(fill = "Scaled value") + 
  theme_void() +
  theme(
    axis.text.y = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

wb_heatmap

ggplot2::ggsave("1_multi_omic_association/60_wb_heatmap.svg", plot = wb_heatmap, width = 2.2, height = 2, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/60_wb_heatmap.tiff", plot = wb_heatmap, width = 2.2, height = 2, units = "cm", dpi = 300)


### mitochondria function assay bar plot
mito_function_data <- read.table("1_multi_omic_association/61_mito_function.txt", header = T, row.names = 1, sep = "\t") %>% 
  rownames_to_column(var = "Function") %>%
  as.data.frame() %>%
  pivot_longer(-Function, names_to = "sample", values_to = "value")

mito_function_data$group <- factor(rep(c(rep("Wild", 4), rep("Lab", 4), rep("Fam", 4), rep("Res", 4), rep("Hos", 4)), 3), levels = c("Wild", "Lab", "Fam", "Res", "Hos"))

mito_function_data$Function <- factor(mito_function_data$Function, levels = c("Mito-tracker", "TMRE", "ATP-red"))

mito_function_data_grouped <- 
  mito_function_data %>% 
  group_by(Function, group) %>%
  summarise(
    mean = mean(value),
    sd = sd(value)
  )

mito_function_data_grouped$top_label <- 'Fat Body Mitochondria Function Assay'

mito_function_data$top_label <- 'Fat Body Mitochondria Function Assay'

mito_function_box_plot <- 
  ggplot(mito_function_data,
         aes(x = group, y = value)
  ) +
  geom_boxplot(
    aes(color = group, fill = group),
    alpha = 0.3, 
    linewidth = 0.25, 
    staplewidth = 0.4,
    width = 0.6,
    outlier.shape = NA,
    show.legend = F
  ) +
  geom_point(aes(fill = group, color = group),
             shape = 21,
             stroke = 0.1,
             size = 1,
             position = position_jitter(width = 0.12, height = 0, seed = 123),
             show.legend = F
  ) +
  ggh4x::facet_nested_wrap(~ top_label + Function, 
                           scales = "free_y",
                           ncol = 3
  ) +
  stat_compare_means(method = "anova",
                     # method = "kruskal.test",
                     size = 2
  ) +
  stat_compare_means(
    comparisons = list(c("Wild", "Lab"), c("Wild", "Fam"), c("Wild", "Res"), c("Wild", "Hos")),
    # method = "wilcox.test",
    method = "t.test",
    # label = "p.adj",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
    # p.adjust.method = "BH", 
    p.adjust.method = "holm", 
    # method.args = list(exact = FALSE, correct = TRUE),
    size = 3,
    vjust = 0.67,
    bracket.size = 0.25,
    step.increase = 0.075,
    tip.length = 0
  ) +
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=1, angle=45, lineheight=1, color="black"),  
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0, angle=90, lineheight=1, color="black"), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    # axis.line.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    panel.border =  element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = ggfun::element_roundrect(
      fill = "#c9caca", 
      color = NA, 
      linewidth = 0, 
      r = unit(0.25, "snpc")
    ),
    strip.text = element_text(family = "Arial", face = "bold", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black")
  )

mito_function_box_plot

ggplot2::ggsave("1_multi_omic_association/62_mito_function_box_plot.svg", plot = mito_function_box_plot, width = 8.25, height = 4.25, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/62_mito_function_box_plot.tiff", plot = mito_function_box_plot, width = 8.25, height = 4.25, units = "cm", dpi = 300)


### mito_function_heatmap
mito_function_data_grouped_scale <- 
  mito_function_data_grouped %>% 
  group_by(Function) %>%
  mutate(across(c(mean), ~ scales::rescale(.x, to = c(-2, 2)))) %>% 
  mutate(group2 = factor(ifelse(group == "Wild", "Wild", "Non-wild"), levels = c("Wild", "Non-wild"))) %>%
  group_by(Function, group2) %>%
  summarise(across(mean, mean, na.rm = TRUE)) 

mito_function_heatmap <-
  ggplot(mito_function_data_grouped_scale) +
  geom_tile(aes(x = group2, y = Function, fill = mean),
            color = "white",
            height = 1,
            linewidth = 0.45,
            show.legend = T
  ) + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF",
                       limits = c(-2, 2),
                       breaks = c(-2, 2),
  ) +
  labs(fill = "Scaled value") + 
  theme_void() +
  theme(
    axis.text.y = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

mito_function_heatmap

ggplot2::ggsave("1_multi_omic_association/63_mito_function_heatmap.svg", plot = mito_function_heatmap, width = 2.2, height = 3, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/63_mito_function_heatmap.tiff", plot = mito_function_heatmap, width = 2.2, height = 3, units = "cm", dpi = 300)


### mito_function_bubble_plot
library(ggpubr)

mito_function_data_grouped_scale2 <- 
  mito_function_data_grouped %>% 
  group_by(Function) %>%
  mutate(across(c(mean), ~ scales::rescale(.x, to = c(-2, 2))))

mito_function_bubble_plot <-
  ggplot(mito_function_data_grouped_scale2) +
  geom_point(aes(x = group, y = Function, size = mean, fill = group, color = group)) + 
  scale_fill_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_color_manual(values = c("Wild" = wild_color, "Lab" = lab_color, "Fam" = fam_color, "Res" = res_color, "Hos" = hos_color)) +
  scale_size_continuous(range = c(3, 5)) +
  facet_wrap2(~ Function, nrow = 3, scales = "free_y") +
  theme_void() + 
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white"),
    strip.text = element_blank(),
    strip.background = element_blank()
  )

mito_function_bubble_plot

ggplot2::ggsave("1_multi_omic_association/64_mito_function_bubble_plot.svg", plot = mito_function_bubble_plot, width = 2.5, height = 5, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/64_mito_function_bubble_plot.tiff", plot = mito_function_bubble_plot, width = 2.5, height = 5, units = "cm", dpi = 300)


### OXPHOS pathway gene bubble plot
ko00190 <- read.table("../4.association/1_multi_omic_association/65_ko00190_KO_tpm.txt", header = T, sep = "\t")

scaled_tpm <- 
  ko00190 %>%
  dplyr::select(-c(1, 3:5)) %>% 
  column_to_rownames(var = "KO") %>%
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>%
  mutate(group = c(rep("wild", 10), rep("lab", 10), rep("fam", 10), rep("res", 10), rep("hos", 10))) %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames(var = "group") %>%
  mutate(across(everything(), ~ scales::rescale(.x, to = c(-2, 2)))) %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "KO")

ko00190_group <- 
  ko00190 %>% 
  dplyr::select(c(1:5)) %>%
  left_join(scaled_tpm, by = "KO")

ko00190_group$Gene <- factor(ko00190_group$Gene, levels = ko00190_group$Gene)

ko00190_group_long <- 
  ko00190_group %>%
  pivot_longer(-c(1:6), names_to = "group", values_to = "Scaled TPM")

ko00190_group_long$group <- factor(ko00190_group_long$group, levels = c("wild","lab","fam","res", "hos"))

### OXPHOX pathway gene expression box plot
ko00190 <- read.table("../4.association/1_multi_omic_association/65_ko00190_KO_tpm.txt", header = T, sep = "\t")

zscore_tpm <- ko00190 %>%
  dplyr::select(-c(1:5)) %>% 
  t() %>%
  as.data.frame() %>% 
  scale() %>%
  as.data.frame() %>%
  t()

zscore_ko00190 <- ko00190

zscore_ko00190[6:55] <-  zscore_tpm

zscore_ko00190$Symbol <- factor(zscore_ko00190$Symbol, levels = zscore_ko00190$Symbol)

zscore_ko00190_long <- 
  zscore_ko00190 %>%
  pivot_longer(-c(1:5), names_to = "sample", values_to = "tpm_z_score")

zscore_ko00190_long$sample <- factor(zscore_ko00190_long$sample, levels = unique(zscore_ko00190_long$sample))

zscore_ko00190_long$group <- factor( rep(c(rep("wild", 10), rep("lab", 10), rep("fam", 10), rep("res", 10),rep("hos", 10)) , 53),  levels = c("wild", "lab", "fam", "res", "hos")) 

zscore_ko00190_long <- 
  zscore_ko00190_long %>%
  mutate(merge_group = case_when(
    group == "wild" ~ "Wild",
    TRUE ~ "Non-wild"
  )) %>% 
  mutate(merge_group = factor(merge_group, levels = c("Wild", "Non-wild")))

zscore_ko00190_long_group <- 
  zscore_ko00190_long %>% 
  group_by(Symbol, merge_group) %>%
  summarise(
    mean_z = mean(tpm_z_score, na.rm = TRUE),
    sd = sd(tpm_z_score, na.rm = TRUE),
    se = sd / sqrt(n()),
    .groups = "drop"
  )

### flower plot
library(scales)

pval_data <- zscore_ko00190_long %>%
  group_by(Symbol) %>%
  summarise(
    esin_vals = list(tpm_z_score[merge_group == "Wild"]),
    other_vals = list(tpm_z_score[merge_group == "Non-wild"]),
    p_value = wilcox.test(unlist(esin_vals), unlist(other_vals), exact = FALSE)$p.value,
    other_higher = median(unlist(other_vals), na.rm = TRUE) > median(unlist(esin_vals), na.rm = TRUE),
    p_signif = case_when(
      p_value > 0.05 ~ "",
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*"
    ),
    .groups = "drop"
  ) %>%
  filter(p_signif != "") %>% 
  mutate(
    y_pos = max(zscore_ko00190_long$tpm_z_score, na.rm = TRUE) * 1.1
  )

ko00190_boxplot_p0 <- 
  ggplot(zscore_ko00190_long_group) +
  geom_point(
    aes(x = Symbol, y = mean_z,  color = merge_group),
    fill = "white",
    stroke = 1,
    shape = 21,
    alpha = 0.75,
    position = position_dodge(width = 0.9),
    size = 1,
    show.legend = F
  ) +
  geom_errorbar(
    aes(x = Symbol, ymin = mean_z - se, ymax = mean_z + se, color = merge_group),
    position = position_dodge(width = 0.9),
    width = 0.75,
    size = 0.5,
    show.legend = F
  ) +
  geom_vline(
    xintercept = seq(1.5, length(unique(zscore_ko00190_long$Symbol)) - 0.5, by = 1),
    linetype = "dashed", 
    color = "black", 
    alpha = 1,
    linewidth = 0.1
  ) +
  labs(x = "gene", y = "TPM Z-score") +
  geom_text(
    data = pval_data,
    aes(x = Symbol, y = y_pos, label = p_signif, color = other_higher),
    size = 6,
    angle = 90,
    hjust = 0.5,
    vjust = 0.75,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("Wild" = "#648C16FF", "Non-wild" = "#FF7200FF")) +
  scale_color_manual(
    values = c("Wild" = "#648C16FF", "Non-wild" = "#FF7200FF", "TRUE" = "#FF7200FF", "FALSE" = "#648C16FF")
  ) +
  coord_cartesian(ylim = c(-2, 3.5), clip = "off") +
  theme_void() + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "Arial", face = "plain", size = 8, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 1, color = "black"),
    axis.text.x = element_text( angle =90, vjust = 0.5, hjust = 1, size = 6, 
                                color = c(rep(paletteer_d("PrettyCols::Bright")[1], 25),
                                          rep(paletteer_d("PrettyCols::Bright")[2], 4),
                                          rep(paletteer_d("PrettyCols::Bright")[3], 4),
                                          rep(paletteer_d("PrettyCols::Bright")[4], 9),
                                          rep(paletteer_d("PrettyCols::Bright")[5], 11))),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 1, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x = unit(2, "pt"),
    axis.ticks.length.y = unit(2, "pt"),
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5, linetype = 1),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )

ko00190_boxplot_p0

ggplot2::ggsave("1_multi_omic_association/66_ko00190_boxplot.svg", plot = ko00190_boxplot_p0, width = 20, height = 10, units = "cm", dpi = 300)
# ggplot2::ggsave("1_multi_omic_association/66_ko00190_boxplot.tiff", plot = ko00190_boxplot_p0, width = 20, height = 10, units = "cm", dpi = 300)


### main heatmap
ko00190_heatmap_p0 <- 
  ggplot(ko00190_long) +
  geom_point(aes(x = Symbol, y = group, size = rescaled_mean_tpm_z_score, fill = ETC, color =  ETC),
             # color = "black",
             stroke = 0.05,
             shape = 21,
             show.legend = T
  ) + 
  annotate("rect", xmin = -Inf, xmax = 25.5, ymin = -Inf, ymax = Inf, fill = paletteer_d("PrettyCols::Bright")[1], color = NA , alpha = 0.10) + 
  annotate("rect", xmin = 0.5, xmax = 25.5, ymin = 5.6, ymax = 6, fill = paletteer_d("PrettyCols::Bright")[1], color = NA , alpha = 1) +
  annotate("text", x = 12.75, y = 5.8, label = c("ETC I"), color = "white", size = 6 / 2.845, hjust = 0.5, vjust = 0.5, fontface = "bold", family = "Arial")+
  annotate("rect", xmin = 25.5, xmax = 29.5, ymin = -Inf, ymax = Inf, fill = paletteer_d("PrettyCols::Bright")[2], color = NA , alpha = 0.15) + 
  annotate("rect", xmin = 25.5, xmax = 29.5, ymin = 5.6, ymax = 6, fill = paletteer_d("PrettyCols::Bright")[2], color = NA , alpha = 1) +
  annotate("text", x = 27.5, y = 5.8, label = c("ETC II"), color = "white", size = 6 / 2.845, hjust = 0.5, vjust = 0.5, fontface = "bold", family = "Arial")+
  annotate("rect", xmin = 29.5, xmax = 33.5, ymin = -Inf, ymax = Inf, fill = paletteer_d("PrettyCols::Bright")[3], color = NA , alpha = 0.15) + 
  annotate("rect", xmin = 29.5, xmax = 33.5, ymin = 5.6, ymax = 6, fill = paletteer_d("PrettyCols::Bright")[3], color = NA , alpha = 1) +
  annotate("text", x = 31.5, y = 5.8, label = c("ETC III"), color = "white", size = 6 / 2.845, hjust = 0.5, vjust = 0.5, fontface = "bold", family = "Arial")+
  annotate("rect", xmin = 33.5, xmax = 42.5, ymin = -Inf, ymax = Inf, fill = paletteer_d("PrettyCols::Bright")[4], color = NA , alpha = 0.15) + 
  annotate("rect", xmin = 33.5, xmax = 42.5, ymin = 5.6, ymax = 6, fill = paletteer_d("PrettyCols::Bright")[4], color = NA , alpha = 1) + 
  annotate("text", x = 38, y = 5.8, label = c("ETC IV"), color = "white", size = 6 / 2.845, hjust = 0.5, vjust = 0.5, fontface = "bold", family = "Arial")+
  annotate("rect", xmin = 42.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = paletteer_d("PrettyCols::Bright")[5], color = NA , alpha = 0.15) + 
  annotate("rect", xmin = 42.5, xmax = Inf, ymin = 5.6, ymax = 6, fill = paletteer_d("PrettyCols::Bright")[5], color = NA , alpha = 1) + 
  annotate("text", x = 48, y = 5.8, label = c("ETC V"), color = "white", size = 6 / 2.845, hjust = 0.5, vjust = 0.5, fontface = "bold", family = "Arial")+
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 4.5, ymax = Inf, fill = hos_color, color = NA , alpha = 1) +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 3.5, ymax = 4.5, fill = res_color, color = NA , alpha = 1) +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 2.5, ymax = 3.5, fill = fam_color, color = NA , alpha = 1) +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 1.5, ymax = 2.5, fill = lab_color, color = NA , alpha = 1) +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 1.5, ymax = -Inf, fill = wild_color, color = NA , alpha = 1) +
  annotate("text", x = 0.25, y = 1, label = c("Wild"), color = "black", size = 6 / 2.845, hjust = 0.5, vjust = 0.5, fontface = "bold", family = "Arial", angle = 90) +
  annotate("text", x = 0.25, y = 2, label = c("Lab"), color = "black", size = 6 / 2.845, hjust = 0.5, vjust = 0.5, fontface = "bold", family = "Arial", angle = 90) +
  annotate("text", x = 0.25, y = 3, label = c("Fam"), color = "black", size = 6 / 2.845, hjust = 0.5, vjust = 0.5, fontface = "bold", family = "Arial", angle = 90) +
  annotate("text", x = 0.25, y = 4, label = c("Res"), color = "black", size = 6 / 2.845, hjust = 0.5, vjust = 0.5, fontface = "bold", family = "Arial", angle = 90) +
  annotate("text", x = 0.25, y = 5, label = c("Hos"), color = "black", size = 6 / 2.845, hjust = 0.5, vjust = 0.5, fontface = "bold", family = "Arial", angle = 90) +
  scale_fill_manual(values = paletteer_d("PrettyCols::Bright")) +
  scale_color_manual(values = paletteer_d("PrettyCols::Bright")) +
  scale_size_continuous(range = c(1, 2)) + 
  scale_x_discrete(position = "bottom") +
  theme_void() + 
  theme(
    legend.position = "none",
    legend.key.size = unit(10, "pt"),
    legend.key.width = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "bold", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "bold", size = 6, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(2.5, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = NA),
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5),
  )

ko00190_heatmap_p0 

ggsave("1_multi_omic_association/67_HostT_ko00190_gene_heatmap.svg", plot = ko00190_heatmap_p0, width = 8.25, height = 3, units = "cm", dpi = 1200)
# ggsave("1_multi_omic_association/67_HostT_ko00190_gene_heatmap.tiff", plot = ko00190_heatmap_p0, width = 8.25, height = 3, units = "cm", dpi = 1200)

ko00190_heatmap_p0_legend <- cowplot::get_legend(ko00190_heatmap_p0 + theme(legend.position = "right"))
ko00190_heatmap_p_legend <- cowplot::plot_grid(ko00190_heatmap_p0_legend, ncol = 1)

ggsave("1_multi_omic_association/68_HostT_ko00190_gene_heatmap_legend.svg", plot = ko00190_heatmap_p_legend, width = 5, height = 5, units = "cm", dpi = 1200)
# ggsave("1_multi_omic_association/68_HostT_ko00190_gene_heatmap_legend.tiff", plot = ko00190_heatmap_p_legend, width = 5, height = 5, units = "cm", dpi = 1200)





### pathway hub-gene expression
# "ko04150", "ko00190", "ko04214", "ko04122"
Gene_ID_2_KO <- read.table("../3.transcriptome/2_downstream_analysis/Gene_2_KO.txt", header = T, sep = "\t")
pathway_2_KO <- read.table("../3.transcriptome/1_OG/KEGG_pame_gene_2_pathway.txt", header = T, sep = "\t")
mRNA_expression_tpm <- read.table("../3.transcriptome/2_downstream_analysis/1_transcript_tpm_filt.tsv", header = T, sep = "\t") 
pathway_ssGSEA <- read.table("../3.transcriptome/1_OG/Total_transcriptome_KEGG_ssGSEA_pathway.tsv", header = T, sep = "\t")

# ko04150
ko04150_gene_exprission <- 
  mRNA_expression_tpm %>% 
  left_join(pathway_2_KO, by = c("gene" = "protein")) %>%
  filter(pathway == "ko04150") %>%
  mutate(ID = paste0(gene, ":", K)) %>%
  dplyr::select(-c(gene,K,pathway)) %>%
  distinct(ID, .keep_all = TRUE) %>%
  column_to_rownames(var  = "ID") %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "sample")


ko04150_pathway_gene_corr <- 
  pathway_ssGSEA %>% dplyr::filter(gene_set == "ko04150") %>% column_to_rownames("gene_set") %>% t() %>% as.data.frame() %>% rownames_to_column(var = "sample") %>% 
  left_join(ko04150_gene_exprission, by = "sample") %>%
  column_to_rownames(var = "sample") %>%
  as.matrix() %>%
  rcorr(type = "spearman")

ko04150_corr_P <- 
  ko04150_pathway_gene_corr$r %>% as.data.frame() %>% dplyr::select(ko04150) %>% rownames_to_column(var = "id") %>% 
  left_join(ko04150_pathway_gene_corr$P %>% as.data.frame() %>% dplyr::select(ko04150) %>% rownames_to_column(var = "id"), by = "id") %>%
  filter(id != "ko04150") %>% 
  dplyr::rename(correlation = ko04150.x, p_value = ko04150.y) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  mutate(FDR_trans = -log10(FDR)) %>%
  filter(correlation > 0) %>% 
  arrange(desc(correlation))

ko04150_gene_exprission_long <- 
  ko04150_gene_exprission %>% 
  pivot_longer(-sample , names_to = "id", values_to = "tpm") %>% 
  filter(id %in% ko04150_corr_P$id) %>% 
  mutate(group = stringr::str_sub(sample, end = -3)) %>%
  group_by(id, group) %>%
  summarise(mean_tpm = mean(tpm, na.rm = TRUE)) %>%
  mutate(zscore = scales::rescale(mean_tpm, to = c(-2,2))) %>%
  mutate(id = factor(id, levels = rev(ko04150_corr_P$id))) %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos")))

ko04150_hub_gene_exprission_heatmap <- 
  ggplot(ko04150_gene_exprission_long,
         aes(x = group, y = id)
  ) +
  geom_tile(aes(fill = zscore),
            color = "white",
            height = 1,
            linewidth = 0.05
  ) + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF",
                       limits = c(-2, 2),
                       breaks = c(-2, 2),
  ) +
  labs(fill = "Scaled value") + 
  theme_void() +
  theme(
    axis.text.y = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

ko04150_hub_gene_exprission_heatmap

# ko00190
ko00190_gene_exprission <- 
  mRNA_expression_tpm %>% 
  left_join(pathway_2_KO, by = c("gene" = "protein")) %>%
  filter(pathway == "ko00190") %>%
  mutate(ID = paste0(gene, ":", K)) %>%
  dplyr::select(-c(gene,K,pathway)) %>%
  distinct(ID, .keep_all = TRUE) %>%
  column_to_rownames(var  = "ID") %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "sample")


ko00190_pathway_gene_corr <- 
  pathway_ssGSEA %>% dplyr::filter(gene_set == "ko00190") %>% column_to_rownames("gene_set") %>% t() %>% as.data.frame() %>% rownames_to_column(var = "sample") %>% 
  left_join(ko00190_gene_exprission, by = "sample") %>%
  column_to_rownames(var = "sample") %>%
  as.matrix() %>%
  rcorr(type = "spearman")

ko00190_corr_P <- 
  ko00190_pathway_gene_corr$r %>% as.data.frame() %>% dplyr::select(ko00190) %>% rownames_to_column(var = "id") %>% 
  left_join(ko00190_pathway_gene_corr$P %>% as.data.frame() %>% dplyr::select(ko00190) %>% rownames_to_column(var = "id"), by = "id") %>%
  filter(id != "ko00190") %>% 
  dplyr::rename(correlation = ko00190.x, p_value = ko00190.y) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  mutate(FDR_trans = -log10(FDR)) %>%
  filter(correlation > 0) %>% 
  arrange(desc(correlation))

ko00190_gene_exprission_long <- 
  ko00190_gene_exprission %>% 
  pivot_longer(-sample , names_to = "id", values_to = "tpm") %>% 
  filter(id %in% ko00190_corr_P$id) %>% 
  mutate(group = stringr::str_sub(sample, end = -3)) %>%
  group_by(id, group) %>%
  summarise(mean_tpm = mean(tpm, na.rm = TRUE)) %>%
  mutate(zscore = scales::rescale(mean_tpm, to = c(-2,2))) %>%
  mutate(id = factor(id, levels = rev(ko00190_corr_P$id))) %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos")))

ko00190_hub_gene_exprission_heatmap <- 
  ggplot(ko00190_gene_exprission_long,
         aes(x = group, y = id)
  ) +
  geom_tile(aes(fill = zscore),
            color = "white",
            height = 1,
            linewidth = 0.05
  ) + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF",
                       limits = c(-2, 2),
                       breaks = c(-2, 2),
  ) +
  labs(fill = "Scaled value") + 
  theme_void() +
  theme(
    axis.text.y = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

ko00190_hub_gene_exprission_heatmap


# ko04214
ko04214_gene_exprission <- 
  mRNA_expression_tpm %>% 
  left_join(pathway_2_KO, by = c("gene" = "protein")) %>%
  filter(pathway == "ko04214") %>%
  mutate(ID = paste0(gene, ":", K)) %>%
  dplyr::select(-c(gene,K,pathway)) %>%
  distinct(ID, .keep_all = TRUE) %>%
  column_to_rownames(var  = "ID") %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "sample")


ko04214_pathway_gene_corr <- 
  pathway_ssGSEA %>% dplyr::filter(gene_set == "ko04214") %>% column_to_rownames("gene_set") %>% t() %>% as.data.frame() %>% rownames_to_column(var = "sample") %>% 
  left_join(ko04214_gene_exprission, by = "sample") %>%
  column_to_rownames(var = "sample") %>%
  as.matrix() %>%
  rcorr(type = "spearman")

ko04214_corr_P <- 
  ko04214_pathway_gene_corr$r %>% as.data.frame() %>% dplyr::select(ko04214) %>% rownames_to_column(var = "id") %>% 
  left_join(ko04214_pathway_gene_corr$P %>% as.data.frame() %>% dplyr::select(ko04214) %>% rownames_to_column(var = "id"), by = "id") %>%
  filter(id != "ko04214") %>% 
  dplyr::rename(correlation = ko04214.x, p_value = ko04214.y) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  mutate(FDR_trans = -log10(FDR)) %>%
  filter(correlation > 0) %>% 
  arrange(desc(correlation))

ko04214_gene_exprission_long <- 
  ko04214_gene_exprission %>% 
  pivot_longer(-sample , names_to = "id", values_to = "tpm") %>% 
  filter(id %in% ko04214_corr_P$id) %>% 
  mutate(group = stringr::str_sub(sample, end = -3)) %>%
  group_by(id, group) %>%
  summarise(mean_tpm = mean(tpm, na.rm = TRUE)) %>%
  mutate(zscore = scales::rescale(mean_tpm, to = c(-2,2))) %>%
  mutate(id = factor(id, levels = rev(ko04214_corr_P$id))) %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos")))

ko04214_hub_gene_exprission_heatmap <- 
  ggplot(ko04214_gene_exprission_long,
         aes(x = group, y = id)
  ) +
  geom_tile(aes(fill = zscore),
            color = "white",
            height = 1,
            linewidth = 0.05
  ) + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF",
                       limits = c(-2, 2),
                       breaks = c(-2, 2),
  ) +
  labs(fill = "Scaled value") + 
  theme_void() +
  theme(
    axis.text.y = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

ko04214_hub_gene_exprission_heatmap


# ko04122
ko04122_gene_exprission <- 
  mRNA_expression_tpm %>% 
  left_join(pathway_2_KO, by = c("gene" = "protein")) %>%
  filter(pathway == "ko04122") %>%
  mutate(ID = paste0(gene, ":", K)) %>%
  dplyr::select(-c(gene,K,pathway)) %>%
  distinct(ID, .keep_all = TRUE) %>%
  column_to_rownames(var  = "ID") %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "sample")


ko04122_pathway_gene_corr <- 
  pathway_ssGSEA %>% dplyr::filter(gene_set == "ko04122") %>% column_to_rownames("gene_set") %>% t() %>% as.data.frame() %>% rownames_to_column(var = "sample") %>% 
  left_join(ko04122_gene_exprission, by = "sample") %>%
  column_to_rownames(var = "sample") %>%
  as.matrix() %>%
  rcorr(type = "spearman")

ko04122_corr_P <- 
  ko04122_pathway_gene_corr$r %>% as.data.frame() %>% dplyr::select(ko04122) %>% rownames_to_column(var = "id") %>% 
  left_join(ko04122_pathway_gene_corr$P %>% as.data.frame() %>% dplyr::select(ko04122) %>% rownames_to_column(var = "id"), by = "id") %>%
  filter(id != "ko04122") %>% 
  dplyr::rename(correlation = ko04122.x, p_value = ko04122.y) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  mutate(FDR_trans = -log10(FDR)) %>%
  filter(correlation > 0) %>% 
  arrange(desc(correlation))

ko04122_gene_exprission_long <- 
  ko04122_gene_exprission %>% 
  pivot_longer(-sample , names_to = "id", values_to = "tpm") %>% 
  filter(id %in% ko04122_corr_P$id) %>% 
  mutate(group = stringr::str_sub(sample, end = -3)) %>%
  group_by(id, group) %>%
  summarise(mean_tpm = mean(tpm, na.rm = TRUE)) %>%
  mutate(zscore = scales::rescale(mean_tpm, to = c(-2,2))) %>%
  mutate(id = factor(id, levels = rev(ko04122_corr_P$id))) %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos")))

ko04122_hub_gene_exprission_heatmap <- 
  ggplot(ko04122_gene_exprission_long,
         aes(x = group, y = id)
  ) +
  geom_tile(aes(fill = zscore),
            color = "white",
            height = 1,
            linewidth = 0.05
  ) + 
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "#FF9898FF",
                       limits = c(-2, 2),
                       breaks = c(-2, 2),
  ) +
  labs(fill = "Scaled value") + 
  theme_void() +
  theme(
    axis.text.y = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

ko04122_hub_gene_exprission_heatmap


### merge plot 
pathway_hub_gene_heatmap <-
  ko04150_hub_gene_exprission_heatmap %>%
  insert_bottom(ko00190_hub_gene_exprission_heatmap, height = 2.35) %>%
  insert_bottom(ko04122_hub_gene_exprission_heatmap, height = 0.06) %>%
  insert_bottom(ko04214_hub_gene_exprission_heatmap, height = 0.58)

pathway_hub_gene_heatmap

ggsave("1_multi_omic_association/69_HostT_pathway_hub_gene_expression_heatmap.svg", plot = pathway_hub_gene_heatmap, width = 3, height = 12, units = "cm", dpi = 1200)
# ggsave("1_multi_omic_association/69_HostT_pathway_hub_gene_expression_heatmap.tiff", plot = pathway_hub_gene_heatmap, width = 3, height = 12, units = "cm", dpi = 1200)


### pathway bubble plot
pathway_ssGSEA_long <-
  pathway_ssGSEA %>% 
  filter(gene_set %in% c("ko04150", "ko00190", "ko04122", "ko04214")) %>%
  pivot_longer(-gene_set, names_to = "sample", values_to = "ssGSEA_score") %>%
  mutate(group = stringr::str_sub(sample, end = -3)) %>%
  group_by(gene_set, group) %>%
  summarise(mean_ssGSEA_score = mean(ssGSEA_score, na.rm = TRUE)) %>%
  mutate(zscore = scales::rescale(mean_ssGSEA_score, to = c(-2,2))) %>%
  mutate(gene_set = factor(gene_set, levels = rev(c("ko04150", "ko00190", "ko04122", "ko04214")))) %>%
  mutate(group = factor(group, levels = c("wild", "lab", "fam", "res", "hos")))

# plot bubble plot
pathway_ssGSEA_bubble_plot <-
  ggplot(pathway_ssGSEA_long, 
         aes(x = group, y = gene_set),
  ) +
  geom_point(aes(size = zscore, fill = group, color = group), 
             shape = 21, 
             show.legend = T
  ) +
  scale_fill_manual(values = c("wild" = wild_color, "lab" = lab_color, "fam" = fam_color, "res" = res_color, "hos" = hos_color)) +
  scale_color_manual(values = c("wild" = wild_color, "lab" = lab_color, "fam" = fam_color, "res" = res_color, "hos" = hos_color)) +
  scale_x_discrete(position = "bottom") +
  labs(fill = "Scaled Value") +
  theme_void() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

pathway_ssGSEA_bubble_plot

ggsave("1_multi_omic_association/70_HostT_pathway_ssGSEA_bubble_plot.svg", plot = pathway_ssGSEA_bubble_plot, width = 5, height = 5, units = "cm", dpi = 1200)
# ggsave("1_multi_omic_association/70_HostT_pathway_ssGSEA_bubble_plot.tiff", plot = pathway_ssGSEA_bubble_plot, width = 5, height = 5, units = "cm", dpi = 1200)



### HostT pathway correlation plot
pathway_corr <- 
  pathway_ssGSEA %>% 
  filter(gene_set %in% c("ko04150", "ko00190", "ko04122", "ko04214")) %>%
  column_to_rownames("gene_set") %>%
  as.matrix() %>% 
  t() %>%
  rcorr(type = "spearman")

HostT_pathway_corr_plot <-
  pathway_corr$r %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "pathway_1") %>% 
  pivot_longer(-pathway_1, names_to = "pathway_2", values_to = "correlation") %>%
  left_join(pathway_corr$P %>% as.data.frame() %>% rownames_to_column(var = "pathway_1") %>% pivot_longer(-pathway_1, names_to = "pathway_2", values_to = "p_value"), by = c("pathway_1", "pathway_2")) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  mutate(significant = ifelse(FDR < 0.05, "*", "")) %>%
  mutate(pathway_1 = factor(pathway_1, levels = c("ko04150", "ko00190", "ko04122", "ko04214"))) %>%
  mutate(pathway_2 = factor(pathway_2, levels = rev(c("ko04150", "ko00190", "ko04122", "ko04214")))) %>%
  ggplot(aes(x = pathway_1, y = pathway_2)
  )+
  geom_tile(aes(fill = correlation), 
            color = "white", 
            linewidth = 0.15
  ) +
  scale_fill_gradient2(low = "lightblue", 
                       mid = "white", 
                       high = "black",
                       limits = c(-1, 1),
                       breaks = c(-1, 0, 1)
  ) +
  geom_text(aes(label = round(correlation,2)),
            color = "black",
            size = 6/2.83,
            family = "Arial",
            fontface = "bold"
  ) +
  labs(fill = "Correlation") +
  scale_x_discrete(position = "top") +
  theme_void() +
  theme(
    # axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    # axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    panel.background = element_blank(),
    legend.position = "none",
    legend.key.size = unit(5, "pt"),
    legend.key.width = unit(5, "pt"), 
    legend.key.height = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(family = "Arial", face = "plain", size = 6, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1, color = "black"),
    legend.ticks.length = unit(0, "pt"),
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "white")
  )

HostT_pathway_corr_plot

ggsave("1_multi_omic_association/71_HostT_pathway_corr_plot.svg", plot = HostT_pathway_corr_plot, width = 3, height = 3, units = "cm", dpi = 1200)
