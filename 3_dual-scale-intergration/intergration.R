### correlation between diff features of scale I and scale II
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

# setwd("D://课题进展/17_Article_Figure&Plot/1_Gut microbes produce branch-chain amino acids for controlling cockroach growth rate/dual-scale-intergration/")


#### scale I diff feature correlation

### diff_MetaG_module 
both_scale_diff_metaG <- 
  read.table("both_scale_enriched_MetaG_in_wild.txt", header = T, sep = "\t") %>%
  rbind(read.table("both_scale_depleted_MetaG_in_wild.txt", header = T, sep = "\t"))

filterd_module <- c("M00570","M00535","M00019","M00432","M00081","M00014","M00847",
                    "M00038","M00044","M00017","M00027","M00555","M00975","M00010","M00012","M00990","M00595","M00616","M00912","M00546")

scale_I_diff_MetaG_module <- 
  read.table("scale I/Esin_MetaG_up_intersection_module.tsv", header = T, sep = "\t") %>% 
  rbind(read.table("scale I/Esin_MetaG_down_intersection_module.tsv", header = T, sep = "\t")) %>%
  filter(module %in% both_scale_diff_metaG$module) %>%
  filter(module %in% filterd_module) %>%
  mutate(module = paste0("MetaG-",module)) %>%
  arrange(desc(pick(8))) %>%
  rename("module" = "feature")

scale_I_diff_MetaG_module_info <- 
  read.table("both_scale_enriched_MetaG_in_wild_info.txt", header = T, sep = "\t") %>% 
  rbind(read.table("both_scale_depleted_MetaG_in_wild_info.txt",header = T, sep = "\t")) %>%
  mutate(module = paste0("MetaG-", module)) %>%
  filter(module %in% scale_I_diff_MetaG_module$feature) %>% 
  dplyr::select(-c(3,6))

### diff_HostT_pathway
both_scale_diff_hostT <- 
  read.table("both_scale_upregulated_HostT_in_wild.txt", header = T, sep = "\t") %>%
  rbind(read.table("both_scale_downregulated_HostT_in_wild.txt", header = T, sep = "\t"))

scale_I_diff_HostT_pathway <- 
  read.table("scale I/Esin_HostT_up_intersection_pathway.tsv", header = T, sep = "\t") %>%
  rbind(read.table("scale I/Esin_HostT_down_intersection_pathway.tsv", header = T, sep = "\t")) %>%
  dplyr::select(1:9, 28:33, 22:27, 10:21) %>% 
  filter(pathway %in% both_scale_diff_hostT$pathway) %>%
  mutate(pathway = paste0("HostT-", pathway)) %>%
  arrange(desc(pick(8))) %>%
  rename("pathway" = "feature")

scale_I_diff_HostT_pathway_info <-
  read.table("both_scale_upregulated_HostT_in_wild_info.txt", header = T, sep = "\t") %>% 
  rbind(read.table("both_scale_downregulated_HostT_in_wild_info.txt",header = T, sep = "\t")) %>%
  mutate(pathway = paste0("HostT-", pathway)) %>%
  filter(pathway %in% scale_I_diff_HostT_pathway$feature)

### diff_MetaB_module 
diff_metaG_cpd <-read_xlsx("scale I/metaG_cpd.xlsx")
scale_I_metaB_module_cpd_link <- read_xlsx("scale I/merge_module_cpd.xlsx")

scale_I_key_metaB <- 
  scale_I_metaB_module_cpd_link %>%
  filter(`KEGG COMPOUND ID` %in% diff_metaG_cpd$KEGG_ID)

scale_I_diff_MetaB_module <-
  read.table("scale I/Esin_MetaB_up_intersection_module.tsv", header = T, sep = "\t") %>%
  rbind(read.table("scale I/Esin_MetaB_down_intersection_module.tsv", header = T, sep = "\t")) %>%
  rename("module" = "feature") %>%
  arrange(desc(pick(8))) %>%
  filter(feature %in% scale_I_key_metaB$MetaB_module_ID)

scale_I_diff_MetaB_module_info <-
  read.table("scale I/metabolome_module_information.txt", header = T, sep = "\t") %>% 
  filter(MetaB_module_ID %in% scale_I_diff_MetaB_module$feature)

### scale I multi-omic correlation
scale_I_MetaG_MetaB_combined_matrix <- 
  rbind(scale_I_diff_MetaB_module,
        scale_I_diff_MetaG_module) %>% 
  dplyr::select(-c(2:9)) %>% 
  column_to_rownames(var = "feature") %>%
  t() %>%
  as.matrix()

scale_I_MetaB_HostT_combined_matrix <- 
  rbind(scale_I_diff_MetaB_module,
        scale_I_diff_HostT_pathway) %>% 
  dplyr::select(-c(2:9)) %>% 
  column_to_rownames(var = "feature") %>%
  t() %>%
  as.matrix()

# Use rcorr to calculate the correlation matrix and the p-value matrix
scale_I_MetaG_MetaB_Correlation <- rcorr(scale_I_MetaG_MetaB_combined_matrix, type = "spearman")
scale_I_MetaB_HostT_Correlation <- rcorr(scale_I_MetaB_HostT_combined_matrix, type = "spearman")

# Extract the correlation matrix and the p-value matrix
scale_I_MetaG_MetaB_cor_matrix <- scale_I_MetaG_MetaB_Correlation$r[1:9, 10:29]
scale_I_MetaG_MetaB_pvalue_matrix <- scale_I_MetaG_MetaB_Correlation$P[1:9, 10:29]

scale_I_MetaB_HostT_cor_matrix <- scale_I_MetaB_HostT_Correlation$r[1:9, 10:21]
scale_I_MetaB_HostT_pvalue_matrix <- scale_I_MetaB_HostT_Correlation$P[1:9, 10:21]

# BH(Benjamini-Hochberg),BY(Benjamini-Yekutieli),Bonferroni
scale_I_MetaG_MetaB_pvalue_matrix <- 
  apply(scale_I_MetaG_MetaB_pvalue_matrix, 2, function(x) {
    p.adjust(x, method = "BH")
  })

scale_I_MetaB_HostT_pvalue_matrix <-
  apply(scale_I_MetaB_HostT_pvalue_matrix, 2, function(x) {
    p.adjust(x, method = "BH")
  })

# wide to long
scale_I_MetaG_MetaB_cor_matrix_long <- 
  scale_I_MetaG_MetaB_cor_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "MetaG", values_to = "Correlation")

scale_I_MetaG_MetaB_pvalue_matrix_long <- 
  scale_I_MetaG_MetaB_pvalue_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "MetaG", values_to = "P_value")

scale_I_MetaB_HostT_cor_matrix_long <- 
  scale_I_MetaB_HostT_cor_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "HostT", values_to = "Correlation")

scale_I_MetaB_HostT_pvalue_matrix_long <- 
  scale_I_MetaB_HostT_pvalue_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "HostT", values_to = "P_value")

# merge
scale_I_MetaG_MetaB_data <- 
  scale_I_MetaG_MetaB_cor_matrix_long %>%
  left_join(scale_I_MetaG_MetaB_pvalue_matrix_long, by = c("MetaB", "MetaG"))

scale_I_MetaB_HostT_data <- 
  scale_I_MetaB_HostT_cor_matrix_long %>%
  left_join(scale_I_MetaB_HostT_pvalue_matrix_long, by = c("MetaB", "HostT"))

scale_I_MetaG_MetaB_data$MetaG <- factor(scale_I_MetaG_MetaB_data$MetaG, levels = rev(scale_I_diff_MetaG_module$feature))
scale_I_MetaG_MetaB_data$MetaB <- factor(scale_I_MetaG_MetaB_data$MetaB, levels = scale_I_diff_MetaB_module$feature)

scale_I_MetaB_HostT_data$HostT <- factor(scale_I_MetaB_HostT_data$HostT, levels = rev(scale_I_diff_HostT_pathway$feature))
scale_I_MetaB_HostT_data$MetaB <- factor(scale_I_MetaB_HostT_data$MetaB, levels = scale_I_diff_MetaB_module$feature)

scale_I_MetaG_MetaB_data <- 
  scale_I_MetaG_MetaB_data %>%
  mutate("-log10(P_value)" = -log10(P_value)) %>%
  mutate(Significance = case_when(
    (Correlation >= 0.400 | Correlation <= -0.400) & P_value <= 0.05 ~ "*",
    TRUE ~ NA
  )) %>%
  mutate(Fill = case_when(
    (Correlation >= 0.400 | Correlation <= -0.400) & P_value < 0.05 ~ Correlation,
    TRUE ~ NA
  )) %>%
  mutate(logtrans = case_when(
    Fill != 0 ~ -log10(P_value),
    TRUE ~ NA
  ))

scale_I_MetaB_HostT_data <- 
  scale_I_MetaB_HostT_data %>%
  mutate("-log10(P_value)" = -log10(P_value)) %>%
  mutate(Significance = case_when(
    (Correlation >= 0.400 | Correlation <= -0.400) & P_value <= 0.05 ~ "*",
    TRUE ~ NA
  )) %>%
  mutate(Fill = case_when(
    (Correlation >= 0.400 | Correlation <= -0.400) & P_value < 0.05 ~ Correlation,
    TRUE ~ NA
  )) %>%
  mutate(logtrans = case_when(
    Fill != 0 ~ -log10(P_value),
    TRUE ~ NA
  ))

### MetaG_MetaB_heatmap
library(ggforce)

scale_I_MetaG_MetaB_heatmap_p0 <- 
  ggplot(scale_I_MetaG_MetaB_data, 
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
                   # labels = rev(paste0("MetaB-", metabolome_info$module))
  ) + 
  scale_y_discrete(position = "right",
                   # labels = metagenome_info$label
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

scale_I_MetaG_MetaB_heatmap_p0

### MetaG_info
scale_I_diff_MetaG_module_info <- 
  scale_I_diff_MetaG_module_info %>% 
  mutate(type = case_when(
    scale_I_diff_MetaG_module_info$category_C == "Aromatic amino acid metabolism" ~ "Aromatic amino acid metabolism",
    scale_I_diff_MetaG_module_info$category_C == "Branched-chain amino acid metabolism" ~ "Branched-chain amino acid metabolism",
    scale_I_diff_MetaG_module_info$category_C == "Cofactor and vitamin metabolism" ~ "Cofactor and vitamin metabolism",
    scale_I_diff_MetaG_module_info$category_C == "Serine and threonine metabolism" ~ "Serine and threonine metabolism",
    scale_I_diff_MetaG_module_info$category_C == "Other carbohydrate metabolism" ~ "Other carbohydrate metabolism",
    scale_I_diff_MetaG_module_info$category_C == "Sulfur metabolism" ~ "Sulfur metabolism",
    TRUE ~ "Other"
  ))

scale_I_diff_MetaG_module_info$type <- factor(scale_I_diff_MetaG_module_info$type, 
                                              levels = c("Aromatic amino acid metabolism", 
                                                         "Branched-chain amino acid metabolism", 
                                                         "Cofactor and vitamin metabolism", 
                                                         "Serine and threonine metabolism", 
                                                         "Other carbohydrate metabolism",
                                                         "Sulfur metabolism",
                                                         "Other"))
scale_I_diff_MetaG_module_info$module <- factor(scale_I_diff_MetaG_module_info$module, levels = rev(scale_I_diff_MetaG_module$feature))

scale_I_MetaG_legendry <-
  ggplot(scale_I_diff_MetaG_module_info) + 
  geom_tile(aes(x = 0, y = module, fill = type), 
            color = "white", 
            height = 1,
            linewidth = 0.5,
            show.legend = T
  ) + 
  scale_fill_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:6], "grey")) +
  scale_color_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:6], "grey")) +
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

scale_I_MetaG_legendry

scale_I_MetaG_MetaB_heatmap <- scale_I_MetaG_MetaB_heatmap_p0  %>% insert_left(scale_I_MetaG_legendry, width = 0.06)
scale_I_MetaG_MetaB_heatmap0 <- scale_I_MetaG_MetaB_heatmap_p0 + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
scale_I_MetaG_MetaB_heatmap1 <- scale_I_MetaG_MetaB_heatmap_p0 + theme(axis.text.x = element_blank())
scale_I_MetaG_MetaB_heatmap2 <- scale_I_MetaG_MetaB_heatmap_p0 + theme(axis.text.y = element_blank())

scale_I_MetaG_MetaB_heatmap
scale_I_MetaG_MetaB_heatmap0
scale_I_MetaG_MetaB_heatmap1
scale_I_MetaG_MetaB_heatmap2

ggsave("scale I/1_scale_I_MetaG_MetaB_heatmap0.svg", plot = scale_I_MetaG_MetaB_heatmap0, width = 2.75, height = 5.5, units = "cm", dpi = 300)
ggsave("scale I/1_scale_I_MetaG_MetaB_heatmap0.tiff", plot = scale_I_MetaG_MetaB_heatmap0, width = 2.75, height = 5.5, units = "cm", dpi = 300)

ggsave("scale I/1_scale_I_MetaG_MetaB_heatmap1.svg", plot = scale_I_MetaG_MetaB_heatmap1, width = 10, height = 5.5, units = "cm", dpi = 300)
ggsave("scale I/1_scale_I_MetaG_MetaB_heatmap1.tiff", plot = scale_I_MetaG_MetaB_heatmap1, width = 10, height = 5.5, units = "cm", dpi = 300)

ggsave("scale I/1_scale_I_MetaG_MetaB_heatmap2.svg", plot = scale_I_MetaG_MetaB_heatmap2, width = 2.75, height = 10, units = "cm", dpi = 300)
ggsave("scale I/1_scale_I_MetaG_MetaB_heatmap2.tiff", plot = scale_I_MetaG_MetaB_heatmap2, width = 2.75, height = 10, units = "cm", dpi = 300)

ggsave("scale I/1_scale_I_MetaG_legendry.svg", plot = scale_I_MetaG_legendry, width = 0.30, height = 5.5, units = "cm", dpi = 300)
ggsave("scale I/1_scale_I_MetaG_legendry.tiff", plot = scale_I_MetaG_legendry, width = 0.30, height = 5.5, units = "cm", dpi = 300)


### MetaB_HostT_heatmap
scale_I_MetaB_HostT_heatmap_p0 <- 
  ggplot(scale_I_MetaB_HostT_data, 
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
                   # labels = rev(paste0("MetaB-", metabolome_info$module))
  ) + 
  scale_y_discrete(position = "right",
                   # labels = paste0("HostT-", transcriptome_info$label)
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

scale_I_MetaB_HostT_heatmap_p0

# HostT_legendry
scale_I_diff_HostT_pathway_info <- 
  scale_I_diff_HostT_pathway_info %>% 
  mutate(type = case_when(
    scale_I_diff_HostT_pathway_info$Subcategory == "09101 Carbohydrate metabolism" ~ "Carbohydrate metabolism",
    scale_I_diff_HostT_pathway_info$Subcategory == "09102 Energy metabolism" ~ "Energy metabolism",
    scale_I_diff_HostT_pathway_info$Subcategory == "09108 Metabolism of cofactors and vitamins" ~ "Metabolism of cofactors and vitamins",
    scale_I_diff_HostT_pathway_info$Subcategory == "09132 Signal transduction" ~ "Signal transduction",
    scale_I_diff_HostT_pathway_info$Subcategory == "09143 Cell growth and death" ~ "Cell growth and death",
    scale_I_diff_HostT_pathway_info$Subcategory == "09157 Sensory system" ~ "Sensory system",
    TRUE ~ "Other"
  ))

scale_I_diff_HostT_pathway_info$type <- factor(scale_I_diff_HostT_pathway_info$type, 
                                               levels = c("Carbohydrate metabolism", 
                                                          "Energy metabolism", 
                                                          "Metabolism of cofactors and vitamins", 
                                                          "Signal transduction", 
                                                          "Cell growth and death", 
                                                          "Sensory system",
                                                          "Other"))
scale_I_diff_HostT_pathway_info$pathway <- factor(scale_I_diff_HostT_pathway_info$pathway, levels = rev(scale_I_diff_HostT_pathway$feature))

scale_I_HostT_legendry <-
  ggplot(scale_I_diff_HostT_pathway_info) + 
  geom_tile(aes(x = 0, y = pathway, fill = type), 
            color = "white", 
            height = 1,
            linewidth = 0.5,
            show.legend = T
  ) + 
  scale_fill_manual(values = c(paletteer::paletteer_d("wesanderson::AsteroidCity1"), "#F26386FF", "grey")) +
  scale_color_manual(values = c(paletteer::paletteer_d("wesanderson::AsteroidCity1"), "#F26386FF", "grey")) +
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

scale_I_HostT_legendry

# MetaB_legendry 
scale_I_diff_MetaB_module_info$Module.category <- factor(scale_I_diff_MetaB_module_info$Module.category, 
                                                  levels = c("Benzenoids", 
                                                             "Organic acids and derivatives", 
                                                             "Organoheterocyclic compounds", 
                                                             "Miscellaneous"))

scale_I_diff_MetaB_module_info$MetaB_module_ID <- factor(scale_I_diff_MetaB_module_info$MetaB_module_ID, levels = scale_I_diff_MetaB_module$feature)

scale_I_MetaB_legendry <-
  ggplot(scale_I_diff_MetaB_module_info) + 
  geom_tile(aes(x = MetaB_module_ID, y = 0, fill = Module.category), 
            color = "white", 
            height = 1,
            linewidth = 0.5,
            show.legend = T
  ) + 
  scale_fill_manual(values = c(paletteer::paletteer_d("LaCroixColoR::Apricot")[1:3], "grey"), drop = FALSE) +
  scale_color_manual(values = c(paletteer::paletteer_d("LaCroixColoR::Apricot")[1:3], "grey"), drop = FALSE) +
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

scale_I_MetaB_legendry

scale_I_MetaB_HostT_heatmap <- scale_I_MetaB_HostT_heatmap_p0  %>% insert_left(scale_I_HostT_legendry, width = 0.06) %>% insert_bottom(scale_I_MetaB_legendry, height = 0.06)
scale_I_MetaB_HostT_heatmap0 <- scale_I_MetaB_HostT_heatmap_p0 + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
scale_I_MetaB_HostT_heatmap1 <- scale_I_MetaB_HostT_heatmap_p0 + theme(axis.text.x = element_blank())

scale_I_MetaB_HostT_heatmap
scale_I_MetaB_HostT_heatmap0
scale_I_MetaB_HostT_heatmap1

ggsave("scale I/2_scale_I_MetaB_HostT_heatmap0.svg", plot = scale_I_MetaB_HostT_heatmap0, width = 2.75, height = 3.25, units = "cm", dpi = 300)
ggsave("scale I/2_scale_I_MetaB_HostT_heatmap0.tiff", plot = scale_I_MetaB_HostT_heatmap0, width = 2.75, height = 3.25, units = "cm", dpi = 300)

ggsave("scale I/2_scale_I_MetaB_HostT_heatmap1.svg", plot = scale_I_MetaB_HostT_heatmap1, width = 10, height = 3.25, units = "cm", dpi = 300)
ggsave("scale I/2_scale_I_MetaB_HostT_heatmap1.tiff", plot = scale_I_MetaB_HostT_heatmap1, width = 10, height = 3.25, units = "cm", dpi = 300)

ggsave("scale I/2_scale_I_HostT_legendry.svg", plot = scale_I_HostT_legendry, width = 0.30, height = 3.25, units = "cm", dpi = 300)
ggsave("scale I/2_scale_I_HostT_legendry.tiff", plot = scale_I_HostT_legendry, width = 0.30, height = 3.25, units = "cm", dpi = 300)

ggsave("scale I/2_scale_I_MetaB_legendry.svg", plot = scale_I_MetaB_legendry, width = 2.75, height = 0.30, units = "cm", dpi = 300)
ggsave("scale I/2_scale_I_MetaB_legendry.tiff", plot = scale_I_MetaB_legendry, width = 2.75, height = 0.30, units = "cm", dpi = 300)

### save all legend
scale_I_MetaG_MetaB_heatmap_p0_legend <- cowplot::get_legend(scale_I_MetaG_MetaB_heatmap_p0 + theme(legend.position = "right"))
scale_I_MetaG_legendry_legend <- cowplot::get_legend(scale_I_MetaG_legendry + theme(legend.position = "right"))
scale_I_MetaB_HostT_heatmap_p0_legend <- cowplot::get_legend(scale_I_MetaB_HostT_heatmap_p0 + theme(legend.position = "right"))
scale_I_MetaB_legendry_legend <- cowplot::get_legend(scale_I_MetaB_legendry + theme(legend.position = "right"))
scale_I_HostT_legendry_legend <- cowplot::get_legend(scale_I_HostT_legendry + theme(legend.position = "right"))

scale_I_MetaG_MetaB_HostT_legend <- cowplot::plot_grid(scale_I_MetaG_MetaB_heatmap_p0_legend,
                                                       scale_I_MetaB_HostT_heatmap_p0_legend,
                                                       scale_I_MetaG_legendry_legend, 
                                                       scale_I_MetaB_legendry_legend, 
                                                       scale_I_HostT_legendry_legend, ncol = 5)

scale_I_MetaG_MetaB_HostT_legend

ggsave("scale I/3_scale_I_MetaG_MetaB_HostT_legend.svg", plot = scale_I_MetaG_MetaB_HostT_legend, width = 20, height = 10, units = "cm", dpi = 300)
ggsave("scale I/3_scale_I_MetaG_MetaB_HostT_legend.tiff", plot = scale_I_MetaG_MetaB_HostT_legend, width = 20, height = 10, units = "cm", dpi = 300)





#### scale II diff feature correlation

### diff_MetaG_module 
both_scale_diff_metaG <- 
  read.table("both_scale_enriched_MetaG_in_wild.txt", header = T, sep = "\t") %>%
  rbind(read.table("both_scale_depleted_MetaG_in_wild.txt", header = T, sep = "\t"))

filterd_module <- c("M00570","M00535","M00019","M00432","M00081","M00014","M00847",
                    "M00038","M00044","M00017","M00027","M00555","M00975","M00010","M00012","M00990","M00595","M00616","M00912","M00546")

scale_II_diff_MetaG_module <- 
  read.table("scale II/Wild_MetaG_up_intersection_module.tsv", header = T, sep = "\t") %>% 
  rbind(read.table("scale II/Wild_MetaG_down_intersection_module.tsv", header = T, sep = "\t")) %>%
  filter(module %in% both_scale_diff_metaG$module) %>%
  filter(module %in% filterd_module) %>%
  mutate(module = paste0("MetaG-",module)) %>%
  dplyr::select(c(1:9, 50:59, 30:39, 10:19, 40:49, 20:29)) %>%
  arrange(desc(pick(8))) %>%
  rename("module" = "feature") %>%
  rename_with(~ str_remove(., "^Metagenome_"))

scale_II_diff_MetaG_module_info <- 
  read.table("both_scale_enriched_MetaG_in_wild_info.txt", header = T, sep = "\t") %>% 
  rbind(read.table("both_scale_depleted_MetaG_in_wild_info.txt",header = T, sep = "\t")) %>%
  mutate(module = paste0("MetaG-", module)) %>%
  filter(module %in% scale_II_diff_MetaG_module$feature) %>% 
  dplyr::select(-c(3,6))

### diff_HostT_pathway
both_scale_diff_hostT <- 
  read.table("both_scale_upregulated_HostT_in_wild.txt", header = T, sep = "\t") %>%
  rbind(read.table("both_scale_downregulated_HostT_in_wild.txt", header = T, sep = "\t"))

scale_II_diff_HostT_pathway <- 
  read.table("scale II/Wild_HostT_up_intersection_pathway.tsv", header = T, sep = "\t") %>%
  rbind(read.table("scale II/Wild_HostT_down_intersection_pathway.tsv", header = T, sep = "\t")) %>%
  filter(pathway %in% both_scale_diff_hostT$pathway) %>%
  mutate(pathway = paste0("HostT-", pathway)) %>%
  arrange(desc(pick(8))) %>%
  rename("pathway" = "feature")

scale_II_diff_HostT_pathway_info <-
  read.table("both_scale_upregulated_HostT_in_wild_info.txt", header = T, sep = "\t") %>% 
  rbind(read.table("both_scale_downregulated_HostT_in_wild_info.txt",header = T, sep = "\t")) %>%
  mutate(pathway = paste0("HostT-", pathway)) %>%
  filter(pathway %in% scale_II_diff_HostT_pathway$feature)

### diff_MetaB_module 
diff_metaG_cpd <-read_xlsx("scale I/metaG_cpd.xlsx")
scale_II_metaB_module_cpd_link <- read_xlsx("scale II/merge_module_cpd.xlsx")

scale_II_key_metaB <- 
  scale_II_metaB_module_cpd_link %>%
  filter(`KEGG Compound ID` %in% diff_metaG_cpd$KEGG_ID)

scale_II_diff_MetaB_module <-
  read.table("scale II/Wild_MetaB_up_intersection_module.tsv", header = T, sep = "\t") %>%
  rbind(read.table("scale II/Wild_MetaB_down_intersection_module.tsv", header = T, sep = "\t")) %>%
  rename("module" = "feature") %>%
  arrange(desc(pick(8))) %>%
  filter(feature %in% scale_II_key_metaB$module_ID)

scale_II_diff_MetaB_module_info <-
  read.table("scale II/metabolome_module_information.txt", header = T, sep = "\t") %>% 
  filter(module_ID %in% scale_II_diff_MetaB_module$feature)


### scale II multi-omic correlation
scale_II_MetaG_MetaB_combined_matrix <- 
  rbind(scale_II_diff_MetaB_module,
        scale_II_diff_MetaG_module) %>% 
  dplyr::select(-c(2:9)) %>% 
  column_to_rownames(var = "feature") %>%
  t() %>%
  as.matrix()

scale_II_MetaB_HostT_combined_matrix <- 
  rbind(scale_II_diff_MetaB_module,
        scale_II_diff_HostT_pathway) %>% 
  dplyr::select(-c(2:9)) %>% 
  column_to_rownames(var = "feature") %>%
  t() %>%
  as.matrix()

# Use rcorr to calculate the correlation matrix and the p-value matrix
scale_II_MetaG_MetaB_Correlation <- rcorr(scale_II_MetaG_MetaB_combined_matrix, type = "spearman")
scale_II_MetaB_HostT_Correlation <- rcorr(scale_II_MetaB_HostT_combined_matrix, type = "spearman")

# Extract the correlation matrix and the p-value matrix
scale_II_MetaG_MetaB_cor_matrix <- scale_II_MetaG_MetaB_Correlation$r[1:13, 14:33]
scale_II_MetaG_MetaB_pvalue_matrix <- scale_II_MetaG_MetaB_Correlation$P[1:13, 14:33]

scale_II_MetaB_HostT_cor_matrix <- scale_II_MetaB_HostT_Correlation$r[1:13, 14:25]
scale_II_MetaB_HostT_pvalue_matrix <- scale_II_MetaB_HostT_Correlation$P[1:13, 14:25]

# BH(Benjamini-Hochberg),BY(Benjamini-Yekutieli),Bonferroni
scale_II_MetaG_MetaB_pvalue_matrix <- 
  apply(scale_II_MetaG_MetaB_pvalue_matrix, 2, function(x) {
    p.adjust(x, method = "BH")
  })

scale_II_MetaB_HostT_pvalue_matrix <-
  apply(scale_II_MetaB_HostT_pvalue_matrix, 2, function(x) {
    p.adjust(x, method = "BH")
  })

# wide to long
scale_II_MetaG_MetaB_cor_matrix_long <- 
  scale_II_MetaG_MetaB_cor_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "MetaG", values_to = "Correlation")

scale_II_MetaG_MetaB_pvalue_matrix_long <- 
  scale_II_MetaG_MetaB_pvalue_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "MetaG", values_to = "P_value")

scale_II_MetaB_HostT_cor_matrix_long <- 
  scale_II_MetaB_HostT_cor_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "HostT", values_to = "Correlation")

scale_II_MetaB_HostT_pvalue_matrix_long <- 
  scale_II_MetaB_HostT_pvalue_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "MetaB") %>%
  pivot_longer(cols = -MetaB, names_to = "HostT", values_to = "P_value")

# merge
scale_II_MetaG_MetaB_data <- 
  scale_II_MetaG_MetaB_cor_matrix_long %>%
  left_join(scale_II_MetaG_MetaB_pvalue_matrix_long, by = c("MetaB", "MetaG"))

scale_II_MetaB_HostT_data <- 
  scale_II_MetaB_HostT_cor_matrix_long %>%
  left_join(scale_II_MetaB_HostT_pvalue_matrix_long, by = c("MetaB", "HostT"))

scale_II_MetaG_MetaB_data$MetaG <- factor(scale_II_MetaG_MetaB_data$MetaG, levels = rev(scale_I_diff_MetaG_module$feature))
scale_II_MetaG_MetaB_data$MetaB <- factor(scale_II_MetaG_MetaB_data$MetaB, levels = scale_II_diff_MetaB_module$feature)

scale_II_MetaB_HostT_data$HostT <- factor(scale_II_MetaB_HostT_data$HostT, levels = rev(scale_I_diff_HostT_pathway$feature))
scale_II_MetaB_HostT_data$MetaB <- factor(scale_II_MetaB_HostT_data$MetaB, levels = scale_II_diff_MetaB_module$feature)

scale_II_MetaG_MetaB_data <- 
  scale_II_MetaG_MetaB_data %>%
  mutate("-log10(P_value)" = -log10(P_value)) %>%
  mutate(Significance = case_when(
    (Correlation >= 0.400 | Correlation <= -0.400) & P_value <= 0.05 ~ "*",
    TRUE ~ NA
  )) %>%
  mutate(Fill = case_when(
    (Correlation >= 0.400 | Correlation <= -0.400) & P_value < 0.05 ~ Correlation,
    TRUE ~ NA
  )) %>%
  mutate(logtrans = case_when(
    Fill != 0 ~ -log10(P_value),
    TRUE ~ NA
  ))

scale_II_MetaB_HostT_data <- 
  scale_II_MetaB_HostT_data %>%
  mutate("-log10(P_value)" = -log10(P_value)) %>%
  mutate(Significance = case_when(
    (Correlation >= 0.400 | Correlation <= -0.400) & P_value <= 0.05 ~ "*",
    TRUE ~ NA
  )) %>%
  mutate(Fill = case_when(
    (Correlation >= 0.400 | Correlation <= -0.400) & P_value < 0.05 ~ Correlation,
    TRUE ~ NA
  )) %>%
  mutate(logtrans = case_when(
    Fill != 0 ~ -log10(P_value),
    TRUE ~ NA
  ))

### MetaG_MetaB_heatmap
library(ggforce)

scale_II_MetaG_MetaB_heatmap_p0 <- 
  ggplot(scale_II_MetaG_MetaB_data, 
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
                   # labels = rev(paste0("MetaB-", metabolome_info$module))
  ) + 
  scale_y_discrete(position = "right",
                   # labels = metagenome_info$label
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

scale_II_MetaG_MetaB_heatmap_p0

### MetaG_info
scale_II_diff_MetaG_module_info <- 
  scale_II_diff_MetaG_module_info %>% 
  mutate(type = case_when(
    scale_II_diff_MetaG_module_info$category_C == "Aromatic amino acid metabolism" ~ "Aromatic amino acid metabolism",
    scale_II_diff_MetaG_module_info$category_C == "Branched-chain amino acid metabolism" ~ "Branched-chain amino acid metabolism",
    scale_II_diff_MetaG_module_info$category_C == "Cofactor and vitamin metabolism" ~ "Cofactor and vitamin metabolism",
    scale_II_diff_MetaG_module_info$category_C == "Serine and threonine metabolism" ~ "Serine and threonine metabolism",
    scale_II_diff_MetaG_module_info$category_C == "Other carbohydrate metabolism" ~ "Other carbohydrate metabolism",
    scale_II_diff_MetaG_module_info$category_C == "Sulfur metabolism" ~ "Sulfur metabolism",
    TRUE ~ "Other"
  ))

scale_II_diff_MetaG_module_info$type <- factor(scale_II_diff_MetaG_module_info$type, 
                                              levels = c("Aromatic amino acid metabolism", 
                                                         "Branched-chain amino acid metabolism", 
                                                         "Cofactor and vitamin metabolism", 
                                                         "Serine and threonine metabolism", 
                                                         "Other carbohydrate metabolism",
                                                         "Sulfur metabolism",
                                                         "Other"))
scale_II_diff_MetaG_module_info$module <- factor(scale_II_diff_MetaG_module_info$module, levels = rev(scale_I_diff_MetaG_module$feature))

scale_II_MetaG_legendry <-
  ggplot(scale_II_diff_MetaG_module_info) + 
  geom_tile(aes(x = 0, y = module, fill = type), 
            color = "white", 
            height = 1,
            linewidth = 0.5,
            show.legend = T
  ) + 
  scale_fill_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:6], "grey")) +
  scale_color_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")[1:6], "grey")) +
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

scale_II_MetaG_legendry

scale_II_MetaG_MetaB_heatmap <- scale_II_MetaG_MetaB_heatmap_p0  %>% insert_left(scale_II_MetaG_legendry, width = 0.06)
scale_II_MetaG_MetaB_heatmap0 <- scale_II_MetaG_MetaB_heatmap_p0 + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
scale_II_MetaG_MetaB_heatmap1 <- scale_II_MetaG_MetaB_heatmap_p0 + theme(axis.text.x = element_blank())
scale_II_MetaG_MetaB_heatmap2 <- scale_II_MetaG_MetaB_heatmap_p0 + theme(axis.text.y = element_blank())

scale_II_MetaG_MetaB_heatmap
scale_II_MetaG_MetaB_heatmap0
scale_II_MetaG_MetaB_heatmap1
scale_II_MetaG_MetaB_heatmap2

ggsave("scale II/1_scale_II_MetaG_MetaB_heatmap0.svg", plot = scale_II_MetaG_MetaB_heatmap0, width = 3.75, height = 5.5, units = "cm", dpi = 300)
ggsave("scale II/1_scale_II_MetaG_MetaB_heatmap0.tiff", plot = scale_II_MetaG_MetaB_heatmap0, width = 3.75, height = 5.5, units = "cm", dpi = 300)

ggsave("scale II/1_scale_II_MetaG_MetaB_heatmap1.svg", plot = scale_II_MetaG_MetaB_heatmap1, width = 10, height = 5.5, units = "cm", dpi = 300)
ggsave("scale II/1_scale_II_MetaG_MetaB_heatmap1.tiff", plot = scale_II_MetaG_MetaB_heatmap1, width = 10, height = 5.5, units = "cm", dpi = 300)

ggsave("scale II/1_scale_II_MetaG_MetaB_heatmap2.svg", plot = scale_II_MetaG_MetaB_heatmap2, width = 3.75, height = 10, units = "cm", dpi = 300)
ggsave("scale II/1_scale_II_MetaG_MetaB_heatmap2.tiff", plot = scale_II_MetaG_MetaB_heatmap2, width = 3.75, height = 10, units = "cm", dpi = 300)

ggsave("scale II/1_scale_II_MetaG_legendry.svg", plot = scale_II_MetaG_legendry, width = 0.30, height = 5.5, units = "cm", dpi = 300)
ggsave("scale II/1_scale_II_MetaG_legendry.tiff", plot = scale_II_MetaG_legendry, width = 0.30, height = 5.5, units = "cm", dpi = 300)


### MetaB_HostT_heatmap
scale_II_MetaB_HostT_heatmap_p0 <- 
  ggplot(scale_II_MetaB_HostT_data, 
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
                   # labels = rev(paste0("MetaB-", metabolome_info$module))
  ) + 
  scale_y_discrete(position = "right",
                   # labels = paste0("HostT-", transcriptome_info$label)
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

scale_II_MetaB_HostT_heatmap_p0

# HostT_legendry
scale_II_diff_HostT_pathway_info <- 
  scale_II_diff_HostT_pathway_info %>% 
  mutate(type = case_when(
    scale_II_diff_HostT_pathway_info$Subcategory == "09101 Carbohydrate metabolism" ~ "Carbohydrate metabolism",
    scale_II_diff_HostT_pathway_info$Subcategory == "09102 Energy metabolism" ~ "Energy metabolism",
    scale_II_diff_HostT_pathway_info$Subcategory == "09108 Metabolism of cofactors and vitamins" ~ "Metabolism of cofactors and vitamins",
    scale_II_diff_HostT_pathway_info$Subcategory == "09132 Signal transduction" ~ "Signal transduction",
    scale_II_diff_HostT_pathway_info$Subcategory == "09143 Cell growth and death" ~ "Cell growth and death",
    scale_II_diff_HostT_pathway_info$Subcategory == "09157 Sensory system" ~ "Sensory system",
    TRUE ~ "Other"
  ))

scale_II_diff_HostT_pathway_info$type <- factor(scale_II_diff_HostT_pathway_info$type, 
                                               levels = c("Carbohydrate metabolism", 
                                                          "Energy metabolism", 
                                                          "Metabolism of cofactors and vitamins", 
                                                          "Signal transduction", 
                                                          "Cell growth and death", 
                                                          "Sensory system",
                                                          "Other"))
scale_II_diff_HostT_pathway_info$pathway <- factor(scale_II_diff_HostT_pathway_info$pathway, levels = rev(scale_I_diff_HostT_pathway$feature))

scale_II_HostT_legendry <-
  ggplot(scale_II_diff_HostT_pathway_info) + 
  geom_tile(aes(x = 0, y = pathway, fill = type), 
            color = "white", 
            height = 1,
            linewidth = 0.5,
            show.legend = T
  ) + 
  scale_fill_manual(values = c(paletteer::paletteer_d("wesanderson::AsteroidCity1"), "#F26386FF", "grey")) +
  scale_color_manual(values = c(paletteer::paletteer_d("wesanderson::AsteroidCity1"), "#F26386FF", "grey")) +
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

scale_II_HostT_legendry

# MetaB_legendry 
scale_II_diff_MetaB_module_info$Module.Category <- factor(scale_II_diff_MetaB_module_info$Module.Category, 
                                                         levels = c("Organic acids and derivatives", 
                                                                    "Organoheterocyclic compounds", 
                                                                    "Lipids and lipid-like molecules",
                                                                    "Phenylpropanoids and polyketides",
                                                                    "Organic oxygen compounds", 
                                                                    "Miscellaneous"))

scale_II_diff_MetaB_module_info$module_ID <- factor(scale_II_diff_MetaB_module_info$module_ID, levels = scale_II_diff_MetaB_module$feature)

scale_II_MetaB_legendry <-
  ggplot(scale_II_diff_MetaB_module_info) + 
  geom_tile(aes(x = module_ID, y = 0, fill = Module.Category), 
            color = "white", 
            height = 1,
            linewidth = 0.5,
            show.legend = T
  ) + 
  scale_fill_manual(values = c(paletteer::paletteer_d("LaCroixColoR::Apricot")[2:6], "grey"), drop = FALSE) +
  scale_color_manual(values = c(paletteer::paletteer_d("LaCroixColoR::Apricot")[2:6], "grey"), drop = FALSE) +
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

scale_II_MetaB_legendry

scale_II_MetaB_HostT_heatmap <- scale_II_MetaB_HostT_heatmap_p0  %>% insert_left(scale_II_HostT_legendry, width = 0.06) %>% insert_bottom(scale_II_MetaB_legendry, height = 0.06)
scale_II_MetaB_HostT_heatmap0 <- scale_II_MetaB_HostT_heatmap_p0 + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
scale_II_MetaB_HostT_heatmap1 <- scale_II_MetaB_HostT_heatmap_p0 + theme(axis.text.x = element_blank())

scale_II_MetaB_HostT_heatmap
scale_II_MetaB_HostT_heatmap0
scale_II_MetaB_HostT_heatmap1

ggsave("scale II/2_scale_II_MetaB_HostT_heatmap0.svg", plot = scale_II_MetaB_HostT_heatmap0, width = 3.75, height = 3.25, units = "cm", dpi = 300)
ggsave("scale II/2_scale_II_MetaB_HostT_heatmap0.tiff", plot = scale_II_MetaB_HostT_heatmap0, width = 3.75, height = 3.25, units = "cm", dpi = 300)

ggsave("scale II/2_scale_II_MetaB_HostT_heatmap1.svg", plot = scale_II_MetaB_HostT_heatmap1, width = 10, height = 3.25, units = "cm", dpi = 300)
ggsave("scale II/2_scale_II_MetaB_HostT_heatmap1.tiff", plot = scale_II_MetaB_HostT_heatmap1, width = 10, height = 3.25, units = "cm", dpi = 300)

ggsave("scale II/2_scale_II_HostT_legendry.svg", plot = scale_II_HostT_legendry, width = 0.30, height = 3.25, units = "cm", dpi = 300)
ggsave("scale II/2_scale_II_HostT_legendry.tiff", plot = scale_II_HostT_legendry, width = 0.30, height = 3.25, units = "cm", dpi = 300)

ggsave("scale II/2_scale_II_MetaB_legendry.svg", plot = scale_II_MetaB_legendry, width = 3.75, height = 0.30, units = "cm", dpi = 300)
ggsave("scale II/2_scale_II_MetaB_legendry.tiff", plot = scale_II_MetaB_legendry, width = 3.75, height = 0.30, units = "cm", dpi = 300)

### save all legend
scale_II_MetaG_MetaB_heatmap_p0_legend <- cowplot::get_legend(scale_II_MetaG_MetaB_heatmap_p0 + theme(legend.position = "right"))
scale_II_MetaG_legendry_legend <- cowplot::get_legend(scale_II_MetaG_legendry + theme(legend.position = "right"))
scale_II_MetaB_HostT_heatmap_p0_legend <- cowplot::get_legend(scale_II_MetaB_HostT_heatmap_p0 + theme(legend.position = "right"))
scale_II_MetaB_legendry_legend <- cowplot::get_legend(scale_II_MetaB_legendry + theme(legend.position = "right"))
scale_II_HostT_legendry_legend <- cowplot::get_legend(scale_II_HostT_legendry + theme(legend.position = "right"))

scale_II_MetaG_MetaB_HostT_legend <- cowplot::plot_grid(scale_II_MetaG_MetaB_heatmap_p0_legend,
                                                       scale_II_MetaB_HostT_heatmap_p0_legend,
                                                       scale_II_MetaG_legendry_legend, 
                                                       scale_II_MetaB_legendry_legend, 
                                                       scale_II_HostT_legendry_legend, ncol = 5)

scale_II_MetaG_MetaB_HostT_legend

ggsave("scale II/3_scale_II_MetaG_MetaB_HostT_legend.svg", plot = scale_II_MetaG_MetaB_HostT_legend, width = 20, height = 10, units = "cm", dpi = 300)
ggsave("scale II/3_scale_II_MetaG_MetaB_HostT_legend.tiff", plot = scale_II_MetaG_MetaB_HostT_legend, width = 20, height = 10, units = "cm", dpi = 300)




























### scale_I_network <- 
scale_I_total_diff_multi_omic_feature <- 
  rbind(scale_I_diff_MetaG_module,
        scale_I_diff_MetaB_module,
        scale_I_diff_HostT_pathway
  )

### correlation between diff multi-omic feature
scale_I_correlation <- 
  scale_I_total_diff_multi_omic_feature %>%
  dplyr::select(1,10:33) %>% 
  column_to_rownames(var = "feature") %>%
  as.matrix()

library(Hmisc)
library(reshape2)
library(igraph)
library(ggraph)
library(dplyr)

scale_I_cor_res <- rcorr(t(scale_I_correlation), type = "spearman")

scale_I_edge_r <- melt(scale_I_cor_res$r, value.name = "R_value")
scale_I_edge_p <- melt(scale_I_cor_res$P, value.name = "P_value")

scale_I_edge_list <- cbind(scale_I_edge_r, P_value = scale_I_edge_p$P_value)
colnames(scale_I_edge_list)[1:2] <- c("from", "to")

scale_I_edge_list_clean <- 
  scale_I_edge_list %>%
  filter(!is.na(P_value)) %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  filter(abs(R_value) > 0.40 & FDR < 0.05) %>%
  filter(
    (grepl("^MetaG-", from) & grepl("^MetaB_", to)) |
      (grepl("^MetaB_", from) & grepl("^HostT-", to))
  )

scale_I_valid_metaB <- 
  scale_I_edge_list_clean %>% 
  filter(grepl("^HostT-", to)) %>% 
  pull(from) %>% 
  unique()

scale_I_valid_metaG <- 
  scale_I_edge_list_clean %>% 
  filter(to %in% scale_I_valid_metaB & grepl("^MetaG-", from)) %>% 
  pull(from) %>% 
  unique()

scale_I_edge_list_cascade <- 
  scale_I_edge_list_clean %>%
  filter(
    (from %in% scale_I_valid_metaG & to %in% scale_I_valid_metaB) | 
      (from %in% scale_I_valid_metaB & grepl("^HostT-", to))
  )

scale_I_g <- graph_from_data_frame(scale_I_edge_list_cascade, directed = TRUE)

V(scale_I_g)$Omics <- case_when(
  grepl("^MetaG-", V(scale_I_g)$name) ~ "MetaG",
  grepl("^MetaB_", V(scale_I_g)$name) ~ "MetaB",
  grepl("^HostT-", V(scale_I_g)$name) ~ "HostT"
)

library(igraph)
library(dplyr)

BCAA_mTOR_nodes <- c("MetaG-M00535", "MetaG-M00019", "MetaG-M00570", "MetaG-M00432", 
                     "MetaB_M183", "MetaB_M134", 
                     "HostT-ko04150", "HostT-ko04214", "HostT-ko00190")

existing_targets <- intersect(BCAA_mTOR_nodes, V(scale_I_g)$name)

other_nodes <- V(scale_I_g)$name[!(V(scale_I_g)$name %in% existing_targets)]

g_others <- delete_vertices(scale_I_g, existing_targets) %>%  as_undirected(mode = "collapse")

set.seed(123)
others_leiden_res <- 
  cluster_leiden(
  g_others, 
  objective_function = "CPM", 
  weights = E(g_others)$weight, 
  resolution_parameter = 0.3
)

final_clusters <- rep(NA, vcount(scale_I_g))
names(final_clusters) <- V(scale_I_g)$name

final_clusters[existing_targets] <- "Cluster_BCAA_mTOR"

final_clusters[V(g_others)$name] <- paste0("Cluster_", others_leiden_res$membership)

V(scale_I_g)$Cluster <- final_clusters

cluster_counts <- table(V(scale_I_g)$Cluster)
keep_clusters <- names(cluster_counts[cluster_counts > 1])
nodes_to_keep <- V(scale_I_g)$name[V(scale_I_g)$Cluster %in% keep_clusters]

scale_I_g_filtered <- subgraph(scale_I_g, nodes_to_keep)

cluster_stats <- table(V(scale_I_g_filtered)$Cluster)
print(cluster_stats)



library(ggraph)
library(ggforce)
library(dplyr)

nodes <- as_data_frame(scale_I_g_filtered, what = "vertices")
edges <- as_data_frame(scale_I_g_filtered, what = "edges")

edges_intra <- edges %>%
  left_join(nodes %>% dplyr::select(name, cluster_from = Cluster), by = c("from" = "name")) %>%
  left_join(nodes %>% dplyr::select(name, cluster_to = Cluster), by = c("to" = "name")) %>%
  filter(cluster_from == cluster_to)

g_plot <- graph_from_data_frame(d = edges_intra, directed = TRUE, vertices = nodes)

set.seed(20260101)
lay_apear <- 
  create_layout(g_plot, layout = 'fr', niter = 1500) %>%
  left_join(scale_I_total_diff_multi_omic_feature, by = c("name" = "feature")) %>%
  mutate(point_size = -log10(adj.P.Val)) %>%
  mutate(name = gsub("^MetaB_|^MetaG-|^HostT-", "", name))

final_apear_p <- 
  ggraph(lay_apear) +
  #geom_mark_hull(
  #  aes(x = x, y = y, fill = Cluster),
  #  concavity = 0.5,
  #  expand = unit(2, "mm"),
  #  alpha = 0.1,
  #  linewidth = 0,
  #  color = NA,
  #  show.legend = FALSE
  #) +
  geom_edge_link(
    aes(edge_width = abs(R_value), color = R_value), 
    alpha = 0.5,
    end_cap = circle(0, 'mm') 
  ) + 
  scale_edge_color_gradient2(
    low = "lightblue", 
    mid = "white", 
    high = "pink", 
    midpoint = 0, 
    name = "Spearman R"
  ) +
  scale_edge_width(range = c(0.2, 0.8), guide = "none") +
  geom_node_point(aes(fill = Omics, size = point_size, color = Omics),
                  shape = 21, 
                  stroke = 0,
                  
  ) +
  scale_fill_manual(values = c("MetaG" = "#D12E6CFF", "MetaB" = "#B5C8E2FF", "HostT" = "#E9C39BFF")) +
  scale_color_manual(values = c("MetaG" = "#D12E6CFF", "MetaB" = "#B5C8E2FF", "HostT" = "#E9C39BFF")) +
  scale_size_continuous(range = c(1, 4)) +
  geom_node_text(
    aes(label = name), 
    repel = TRUE, 
    point.padding = unit(0.5, "lines"),
    size = 5/2.845, 
    fontface = "plain",
    bg.color = NA, 
    bg.r = 0
  ) +
  labs(title = "Scale I Multi-omic Modular Network") +
  theme_void() +
  coord_cartesian(clip = "off")

print(final_apear_p)

ggsave("scale I/4_scale_I_Multi-omic_Modular_Network.svg", plot = final_apear_p, width = 10, height = 10, units = "cm", dpi = 300)
ggsave("scale I/4_scale_I_Multi-omic_Modular_Network.tiff", plot = final_apear_p, width = 10, height = 10, units = "cm", dpi = 300)






















### scale_II_network <- 
scale_II_total_diff_multi_omic_feature <- 
  rbind(scale_II_diff_MetaG_module,
        scale_II_diff_MetaB_module,
        scale_II_diff_HostT_pathway
  )

### correlation between diff multi-omic feature
scale_II_correlation <- 
  scale_II_total_diff_multi_omic_feature %>%
  dplyr::select(1,10:59) %>% 
  column_to_rownames(var = "feature") %>%
  as.matrix()

library(Hmisc)
library(reshape2)
library(igraph)
library(ggraph)
library(dplyr)

scale_II_cor_res <- rcorr(t(scale_II_correlation), type = "spearman")

scale_II_edge_r <- melt(scale_II_cor_res$r, value.name = "R_value")
scale_II_edge_p <- melt(scale_II_cor_res$P, value.name = "P_value")

scale_II_edge_list <- cbind(scale_II_edge_r, P_value = scale_II_edge_p$P_value)
colnames(scale_II_edge_list)[1:2] <- c("from", "to")

scale_II_edge_list_clean <- 
  scale_II_edge_list %>%
  filter(!is.na(P_value)) %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  filter(abs(R_value) > 0.40 & FDR < 0.05) %>%
  filter(
    (grepl("^MetaG-", from) & grepl("^MetaB_", to)) |
      (grepl("^MetaB_", from) & grepl("^HostT-", to))
  )

scale_II_valid_metaB <- 
  scale_II_edge_list_clean %>% 
  filter(grepl("^HostT-", to)) %>% 
  pull(from) %>% 
  unique()

scale_II_valid_metaG <- 
  scale_II_edge_list_clean %>% 
  filter(to %in% scale_II_valid_metaB & grepl("^MetaG-", from)) %>% 
  pull(from) %>% 
  unique()

scale_II_edge_list_cascade <- 
  scale_II_edge_list_clean %>%
  filter(
    (from %in% scale_II_valid_metaG & to %in% scale_II_valid_metaB) | 
      (from %in% scale_II_valid_metaB & grepl("^HostT-", to))
  )

scale_II_g <- graph_from_data_frame(scale_II_edge_list_cascade, directed = TRUE)

V(scale_II_g)$Omics <- case_when(
  grepl("^MetaG-", V(scale_II_g)$name) ~ "MetaG",
  grepl("^MetaB_", V(scale_II_g)$name) ~ "MetaB",
  grepl("^HostT-", V(scale_II_g)$name) ~ "HostT"
)

library(igraph)
library(dplyr)

BCAA_mTOR_nodes <- c("MetaG-M00535", "MetaG-M00019", "MetaG-M00570", "MetaG-M00432", 
                     "MetaB_M144",
                     "HostT-ko04150", "HostT-ko04214", "HostT-ko00190")

existing_targets <- intersect(BCAA_mTOR_nodes, V(scale_II_g)$name)

other_nodes <- V(scale_II_g)$name[!(V(scale_II_g)$name %in% existing_targets)]

g_others <- delete_vertices(scale_II_g, existing_targets) %>%  as_undirected(mode = "collapse")

set.seed(123)
others_leiden_res <- 
  cluster_leiden(
    g_others, 
    objective_function = "CPM", 
    weights = E(g_others)$weight, 
    resolution_parameter = 0.3
  )

final_clusters <- rep(NA, vcount(scale_II_g))
names(final_clusters) <- V(scale_II_g)$name

final_clusters[existing_targets] <- "Cluster_BCAA_mTOR"

final_clusters[V(g_others)$name] <- paste0("Cluster_", others_leiden_res$membership)

V(scale_II_g)$Cluster <- final_clusters

cluster_counts <- table(V(scale_II_g)$Cluster)
keep_clusters <- names(cluster_counts[cluster_counts > 1])
nodes_to_keep <- V(scale_II_g)$name[V(scale_II_g)$Cluster %in% keep_clusters]

scale_II_g_filtered <- subgraph(scale_II_g, nodes_to_keep)

cluster_stats <- table(V(scale_II_g_filtered)$Cluster)
print(cluster_stats)



library(ggraph)
library(ggforce)
library(dplyr)

nodes <- as_data_frame(scale_II_g_filtered, what = "vertices")
edges <- as_data_frame(scale_II_g_filtered, what = "edges")

edges_intra <- edges %>%
  left_join(nodes %>% dplyr::select(name, cluster_from = Cluster), by = c("from" = "name")) %>%
  left_join(nodes %>% dplyr::select(name, cluster_to = Cluster), by = c("to" = "name")) %>%
  filter(cluster_from == cluster_to)

g_plot <- graph_from_data_frame(d = edges_intra, directed = TRUE, vertices = nodes)

set.seed(20260202)
lay_apear <- 
  create_layout(g_plot, layout = 'fr', niter = 1500) %>%
  left_join(scale_II_total_diff_multi_omic_feature, by = c("name" = "feature")) %>%
  mutate(point_size = -log10(adj.P.Val)) %>%
  mutate(name = gsub("^MetaB_|^MetaG-|^HostT-", "", name))

final_apear_p <- 
  ggraph(lay_apear) + 
  #geom_mark_hull(
  #  aes(x = x, y = y, fill = Cluster),
  #  concavity = 0.5,
  #  expand = unit(2, "mm"),
  #  alpha = 0.1,
  #  size = 0,
  #  color = NA,
  #  show.legend = FALSE
  #) +
  geom_edge_link(
    aes(edge_width = abs(R_value), color = R_value), 
    alpha = 0.5,
    end_cap = circle(0, 'mm') 
  ) + 
  scale_edge_color_gradient2(
    low = "lightblue", 
    mid = "white", 
    high = "pink", 
    midpoint = 0, 
    name = "Spearman R"
  ) +
  scale_edge_width(range = c(0.2, 0.8)) +
  geom_node_point(aes(fill = Omics, size = point_size, color = Omics),
                  shape = 21, 
                  stroke = 0,
                  
  ) +
  scale_fill_manual(values = c("MetaG" = "#D12E6CFF", "MetaB" = "#B5C8E2FF", "HostT" = "#E9C39BFF")) +
  scale_color_manual(values = c("MetaG" = "#D12E6CFF", "MetaB" = "#B5C8E2FF", "HostT" = "#E9C39BFF")) +
  scale_size_continuous(range = c(1, 4)) +
  geom_node_text(
    aes(label = name), 
    repel = TRUE, 
    point.padding = unit(0.5, "lines"),
    size = 5/2.845, 
    fontface = "plain",
    bg.color = NA, 
    bg.r = 0
  ) +
  labs(title = "Scale I Multi-omic Modular Network") +
  theme_void() +
  coord_cartesian(clip = "off")

print(final_apear_p)

ggsave("scale II/4_scale_II_Multi-omic_Modular_Network.svg", plot = final_apear_p, width = 10, height = 10, units = "cm", dpi = 300)
ggsave("scale II/4_scale_II_Multi-omic_Modular_Network.tiff", plot = final_apear_p, width = 10, height = 10, units = "cm", dpi = 300)





### Process of biological signal molecule transmission
# remotes::install_github("hrbrmstr/ggchicklet") # this package can plot round rect 
library(ggchicklet)
library(ggplot2)
library(ggforce)

# Set the attributes of all elements
Substrate_point_fill = "#FFD966"
Substrate_point_color = "orange"
Substrate_point_linewidth = 0.75
Substrate_point_r = 0.4
Substrate_text_color = "orange"
Substrate_text_size = 6
Substrate_text_hjust = 0
Substrate_text_vjust = 0

Final_point_fill = "#A9D18E"
Final_point_color = "#548235"
Final_point_linewidth = 0.75
Final_point_r = 0.4
Final_text_color = "#548235"
Final_text_size = 6
Final_text_hjust = 0.5
Final_text_vjust = 0

mid_point_fill = "#A5A5A5"
mid_point_color = "#787878"
mid_point_linewidth = 0.75
mid_point_r = 0.4
mid_text_color = "black"
mid_text_size = 4
mid_text_hjust = 0.5
mid_text_vjust = 0

rect_fill = "#4472C4"
rect_color = "#2F528F"
rect_radius = unit(0.2, "npc")
rect_linetype = "solid"
rect_linewidth = 0.75
rect_alpha = 0.7

KO_text_color = "white"
KO_text_size = 6
KO_text_hjust = 0.5
KO_text_vjust = 0.5

arrow <- arrow(length = unit(0.015, "npc"), angle = 15, type = "closed", ends = "last")
arrow2 <- arrow(length = unit(0.03, "npc"), angle = 15, type = "closed", ends = "last")
point_arrow_dis = 0.15
rect_arrow_dis = 0.15

### plot module background block
BCAA_pathway_plot0 <- 
  ggplot() +
  xlim(c(-1.5, 64)) +
  ylim(c(-1.5, 34.5)) +       
  geom_rect(aes(xmin = -1.5, xmax = 30, ymin = -Inf, ymax = Inf), color = "NA", alpha = 0.2, fill = "lightblue") + 
  geom_rect(aes(xmin = 30, xmax = 42.5, ymin = -Inf, ymax = Inf), color = "NA", alpha = 0.2, fill = "pink") + 
  geom_rect(aes(xmin = 42.5, xmax = Inf, ymin = -Inf, ymax = Inf), color = "NA", alpha = 0.2, fill = "lightyellow") + 
  geom_rect(aes(xmin = 42.5, xmax = 58, ymin = -Inf, ymax = 12), color = "NA", alpha = 0.1, fill = "orange") + 
  theme_void()

BCAA_pathway_plot0

BCAA_pathway_plot <- 
  BCAA_pathway_plot0 +
  # M00570 pathway plot(Isoleucine biosynthesis)
  geom_circle(aes(x0 = 1, y0 = 1 + 16, r = Substrate_point_r), fill = Substrate_point_fill, color = Substrate_point_color, linewidth = Substrate_point_linewidth) + 
  geom_text(aes(x = -1, y = 1.7 + 16), label = "L-Threonine", size = Substrate_text_size, color = Substrate_text_color, hjust = Substrate_text_hjust, vjust = Substrate_text_vjust) + 
  geom_segment(aes(x = 1 + Substrate_point_r + point_arrow_dis, xend = 4 - rect_arrow_dis, y = 1 + 16, yend = 1 + 16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 4, xmax = 8, ymin = 0.5 + 16, ymax = 1.5 + 16), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(x = 6, y = 1 + 16), label = "K01754", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust)  +
  geom_segment(aes(x = 6, xend = 6, y = 1.5 + rect_arrow_dis + 16, yend = 3.5 - mid_point_r - point_arrow_dis + 16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 6, y0 = 3.5 + 16, r = mid_point_r), fill = mid_point_fill, color = mid_point_color, linewidth = mid_point_linewidth) + 
  geom_text(aes(x = 6, y = 4.2 + 16), label = "2-Oxobutanoate", size = mid_text_size, color = mid_text_color, hjust = mid_text_hjust, vjust = mid_text_vjust) +
  geom_segment(aes(x = 6 + mid_point_r + point_arrow_dis, xend = 9 - rect_arrow_dis, y = 3.5 + 16, yend = 3.5 + 16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 9, xmax = 13, ymin = 3.6 + 16, ymax = 4.6 + 16), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_rrect(aes(xmin = 9, xmax = 13, ymin = 2.4 + 16, ymax = 3.4 + 16), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(y = 4.2 + 16, x = 11), label = "K01652", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust) +
  geom_text(aes(y = 2.9 + 16, x = 11), label = "K01653", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust) +
  geom_segment(aes(x = 11 , xend = 11 , y = 2.4 - rect_arrow_dis + 16, yend = 1 + mid_point_r + point_arrow_dis + 16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 11, y0 = 1 + 16, r = mid_point_r), fill = mid_point_fill, color = mid_point_color, linewidth = mid_point_linewidth) + 
  geom_text(aes(x = 11, y = 0 + 16), label = "(S)-2-Aceto-2-hydroxybutanoate", size = mid_text_size, color = mid_text_color, hjust = mid_text_hjust, vjust = mid_text_vjust) +
  geom_segment(aes(x = 11 + mid_point_r + point_arrow_dis, xend = 14 - rect_arrow_dis, y = 1 + 16, yend = 1 + 16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 14, xmax = 18, ymin = 0.5 + 16, ymax = 1.5 + 16), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(y = 1 + 16, x = 16), label = "K00053", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust) +
  geom_segment(aes(x = 16, xend = 16, y = 1.5 + rect_arrow_dis + 16, yend = 3.5 - mid_point_r - point_arrow_dis +16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 16, y0 = 3.5 + 16, r = mid_point_r), fill = mid_point_fill, color = mid_point_color, linewidth = mid_point_linewidth) +
  geom_text(aes(x = 16, y = 4.2 + 16), label = "(R)-2,3-Dihydroxy-3-\nmethylvalerate", size = mid_text_size, color = mid_text_color, hjust = mid_text_hjust, vjust = mid_text_vjust) +
  geom_segment(aes(x = 16 + mid_point_r + point_arrow_dis, xend = 19 - rect_arrow_dis, y = 3.5 + 16, yend = 3.5 + 16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 19, xmax = 23, ymin = 3 + 16, ymax = 4 + 16), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(y = 3.5 + 16, x = 21), label = "K01687", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust) +
  geom_segment(aes(x = 21 , xend = 21, y = 3 - rect_arrow_dis + 16, yend = 1 + mid_point_r + point_arrow_dis + 16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 21, y0 = 1 + 16, r = mid_point_r), fill = mid_point_fill, color = mid_point_color, linewidth = mid_point_linewidth) +
  geom_text(aes(x = 21, y = 0 + 16), label = "(S)-3-Methyl-2-oxopentanoate", size = mid_text_size, color = mid_text_color, hjust = mid_text_hjust, vjust = mid_text_vjust) +
  geom_segment(aes(x = 21 + mid_point_r + point_arrow_dis, xend = 24 - rect_arrow_dis, y = 1 + 16, yend = 1 + 16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 24, xmax = 28, ymin = 0.5 + 16, ymax = 1.5 + 16), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(y = 1 + 16, x = 26), label = "K00826", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust) +
  geom_segment(aes(x = 26 , xend = 26, y = 1.5 + rect_arrow_dis + 16, yend = 3.5 - Final_point_r - point_arrow_dis + 16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 26, y0 = 3.5 + 16, r = Final_point_r), fill = Final_point_fill, color = Final_point_color, linewidth = Final_point_linewidth) +
  geom_text(aes(x = 26, y = 4.2 + 16), label = "Isoleucine", size = Final_text_size, color = Final_text_color, hjust = Final_text_hjust, vjust = Final_text_vjust) +
  geom_text(aes(x = 15, y = 6 + 16), label = "M00570:Isoleucine biosynthesis", size = 8, color = "#005050", hjust = 0.5, vjust = 0.5) +
  # M00019 pathway plot(Valine/isoleucine biosynthesis)
  geom_circle(aes(x0 = 1, y0 = 9, r = Substrate_point_r), fill = Substrate_point_fill, color = Substrate_point_color, linewidth = Substrate_point_linewidth) +
  geom_text(aes(x = -1, y = 9.7), label = "2-Oxobutanoate", size = Substrate_text_size, color = Substrate_text_color, hjust = Substrate_text_hjust, vjust = Substrate_text_vjust) + 
  geom_segment(aes(x = 1 + Substrate_point_r + point_arrow_dis, xend = 4 - rect_arrow_dis, y = 9, yend = 9), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 4, xmax = 8, ymin = 9.1, ymax = 10.1), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_rrect(aes(xmin = 4, xmax = 8, ymin = 7.9, ymax = 8.9), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(x = 6, y = 9.6), label = "K01652", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust)  +
  geom_text(aes(x = 6, y = 8.4), label = "K01653", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust)  +
  geom_segment(aes(x = 6, xend = 6, y = 10.1 + rect_arrow_dis, yend = 11.5 - mid_point_r - point_arrow_dis), linewidth = 0.75, linetype = "solid", arrow = arrow) + 
  geom_circle(aes(x0 = 6, y0 = 11.5, r = mid_point_r), fill = mid_point_fill, color = mid_point_color, linewidth = mid_point_linewidth) + 
  geom_text(aes(x = 6, y = 12.2), label = "(S)-2-Aceto-2-hydroxybutanoate", size = mid_text_size, color = mid_text_color, hjust = mid_text_hjust, vjust = mid_text_vjust) +
  geom_segment(aes(x = 6 + mid_point_r + point_arrow_dis, xend = 9 - rect_arrow_dis, y = 11.5, yend = 11.5), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 9, xmax = 13, ymin = 11, ymax = 12), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(y = 11.5, x = 11), label = "K00053", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust) +
  geom_segment(aes(x = 11, xend = 11, y = 11 - rect_arrow_dis, yend = 9 + mid_point_r + point_arrow_dis), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 11, y0 = 9, r = mid_point_r), fill = mid_point_fill, color = mid_point_color, linewidth = mid_point_linewidth) + 
  geom_text(aes(x = 11, y = 7), label = "(R)-2,3-Dihydroxy-\n3-methylvalerate", size = mid_text_size, color = mid_text_color, hjust = mid_text_hjust, vjust = mid_text_vjust) +
  geom_segment(aes(x = 11 + mid_point_r + point_arrow_dis, xend = 14 - rect_arrow_dis, y = 9, yend = 9), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 14, xmax = 18, ymin = 8.5, ymax = 9.5), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(y = 9, x = 16), label = "K01687", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust) +
  geom_segment(aes(x = 16, xend = 16, y = 9.5 + rect_arrow_dis, yend = 11.5 - mid_point_r - point_arrow_dis), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 16, y0 = 11.5, r = mid_point_r), fill = mid_point_fill, color = mid_point_color, linewidth = mid_point_linewidth) +
  geom_text(aes(x = 16, y = 12.2), label = "(S)-3-Methyl-2-oxopentanoate", size = mid_text_size, color = mid_text_color, hjust = mid_text_hjust, vjust = mid_text_vjust)  +
  geom_segment(aes(x = 16 + mid_point_r + point_arrow_dis, xend = 19 - rect_arrow_dis, y = 11.5, yend = 11.5), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 19, xmax = 23, ymin = 11, ymax = 12), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(y = 11.5, x = 21), label = "K00826", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust) +
  geom_segment(aes(x = 23 + rect_arrow_dis, xend = 26 - Final_point_r - point_arrow_dis, y = 11.5, yend = 11.5), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 26, y0 = 11.5, r = Final_point_r), fill = Final_point_fill, color = Final_point_color, linewidth = Final_point_linewidth) +
  geom_text(aes(x = 26, y = 12.2), label = "Isoleucine", size = Final_text_size, color = Final_text_color, hjust = Final_text_hjust, vjust = Final_text_vjust) +
  geom_segment(aes(x = 24, xend = 24, y = 11.5, yend = 9), linewidth = 0.75, linetype = "solid") +
  geom_segment(aes(x = 24, xend = 26 - Final_point_r - point_arrow_dis, y = 9, yend = 9), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 26, y0 = 9, r = Final_point_r), fill = Final_point_fill, color = Final_point_color, linewidth = Final_point_linewidth) +
  geom_text(aes(x = 26, y = 9.7), label = "Valine", size = Final_text_size, color = Final_text_color, hjust = Final_text_hjust, vjust = Final_text_vjust) +
  geom_text(aes(x = 15, y = 14), label = "M00019:Valine/isoleucine biosynthesis", size = 8, color = "#005050", hjust = 0.5, vjust = 0.5) + 
  # M00535 pathway plot(Isoleucine biosynthesis)
  geom_circle(aes(x0 = 1, y0 = 17 -16, r = Substrate_point_r), fill = Substrate_point_fill, color = Substrate_point_color, linewidth = Substrate_point_linewidth) +
  geom_text(aes(x = -1, y = 17.7 - 16), label = "Pyruvate", size = Substrate_text_size, color = Substrate_text_color, hjust = Substrate_text_hjust, vjust = Substrate_text_vjust) + 
  geom_segment(aes(x = 1 + Substrate_point_r + point_arrow_dis, xend = 4 - rect_arrow_dis, y = 17 -16, yend = 17 -16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 4, xmax = 8, ymin = 16.5 -16, ymax = 17.5 -16), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(x = 6, y = 17 -16), label = "K09011", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust)  +
  geom_segment(aes(x = 6, xend = 6, y = 17.5 + rect_arrow_dis -16, yend = 19.5 - mid_point_r - point_arrow_dis -16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 6, y0 = 19.5 -16, r = mid_point_r), fill = mid_point_fill, color = mid_point_color, linewidth = mid_point_linewidth) + 
  geom_text(aes(x = 6, y = 20.2 -16), label = "D-Citramalic acid", size = mid_text_size, color = mid_text_color, hjust = mid_text_hjust, vjust = mid_text_vjust)  +
  geom_segment(aes(x = 6 + mid_point_r + point_arrow_dis, xend = 9 - rect_arrow_dis, y = 19.5 -16, yend = 19.5 -16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 9, xmax = 13, ymin = 19.6 -16, ymax = 20.6 -16), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_rrect(aes(xmin = 9, xmax = 13, ymin = 18.4 -16, ymax = 19.4 -16), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(y = 20.1 -16, x = 11), label = "K01703", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust) +
  geom_text(aes(y = 18.9 -16, x = 11), label = "K01704", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust) +
  geom_segment(aes(x = 13 + rect_arrow_dis, xend = 16 - mid_point_r - point_arrow_dis, y = 19.5 -16, yend = 19.5 -16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 16, y0 = 19.5 -16, r = mid_point_r), fill = mid_point_fill, color = mid_point_color, linewidth = mid_point_linewidth) + 
  geom_text(aes(x = 16, y = 20.2 -16), label = "(2R,3S)-3-Methylmalate", size = mid_text_size, color = mid_text_color, hjust = mid_text_hjust, vjust = mid_text_vjust) +
  geom_segment(aes(x = 16 + mid_point_r + point_arrow_dis, xend = 19 - rect_arrow_dis, y = 19.5 -16, yend = 19.5 -16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 19, xmax = 23, ymin = 19 -16, ymax = 20 -16), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(y = 19.5 -16, x = 21), label = "K00052", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust) +
  geom_segment(aes(x = 23 + rect_arrow_dis, xend = 26 - Final_point_r - point_arrow_dis, y = 19.5 -16, yend = 19.5 -16), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 26, y0 = 19.5 -16, r = Final_point_r), fill = Final_point_fill, color = Final_point_color, linewidth = Final_point_linewidth) +
  geom_text(aes(x = 26, y = 20.2 -16), label = "2-Oxobutanoate", size = Final_text_size, color = Final_text_color, hjust = Final_text_hjust, vjust = Final_text_vjust) +
  geom_text(aes(x = 15, y = 22 -16), label = "M00535:Isoleucine biosynthesis", size = 8, color = "#005050", hjust = 0.5, vjust = 0.5) + 
  # link between M00535 and M00019
  #geom_segment(aes(x = 26, xend = 26 , y = 19.5 - Final_point_r - point_arrow_dis - 16, yend = 15.5 - 16), linewidth = 0.75, linetype = "dashed", lineend = "square") +
  #geom_segment(aes(x = 26, xend = 1 , y = 15.5 - 16, yend = 15.5 -16), linewidth = 0.75, linetype = "dashed", lineend = "square") + 
  #geom_segment(aes(x = 1, xend = 1 , y = 15.5 - 16, yend = 9.7 + Substrate_point_r + point_arrow_dis - 16), linewidth = 0.75, linetype = "dashed", lineend = "square", arrow = arrow) + 
  # M00432 pathway plot(Leucine biosynthesis)
  geom_circle(aes(x0 = 1, y0 = 25, r = Substrate_point_r), fill = Substrate_point_fill, color = Substrate_point_color, linewidth = Substrate_point_linewidth) +
  geom_text(aes(x = -1, y = 25.7), label = "2-Oxoisovalerate", size = Substrate_text_size, color = Substrate_text_color, hjust = Substrate_text_hjust, vjust = Substrate_text_vjust) + 
  geom_segment(aes(x = 1 + Substrate_point_r + point_arrow_dis, xend = 4 - rect_arrow_dis, y = 25, yend = 25), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 4, xmax = 8, ymin = 24.5, ymax = 25.5), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(x = 6, y = 25), label = "K01649", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust)  +
  geom_segment(aes(x = 6, xend = 6, y = 25.5 + rect_arrow_dis, yend = 27.5 - mid_point_r - point_arrow_dis), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 6, y0 = 27.5, r = mid_point_r), fill = mid_point_fill, color = mid_point_color, linewidth = mid_point_linewidth) + 
  geom_text(aes(x = 6, y = 28.2), label = "α-Isopropylmalate", size = mid_text_size, color = mid_text_color, hjust = mid_text_hjust, vjust = mid_text_vjust) + 
  geom_segment(aes(x = 6 + mid_point_r + point_arrow_dis, xend = 9 - rect_arrow_dis, y = 27.5, yend = 27.5), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 9, xmax = 13, ymin = 27.6, ymax = 28.6), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_rrect(aes(xmin = 9, xmax = 13, ymin = 26.4, ymax = 27.4), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(y = 28.1, x = 11), label = "K01703", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust) +
  geom_text(aes(y = 26.9, x = 11), label = "K01704", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust) +
  geom_segment(aes(x = 13 + rect_arrow_dis, xend = 16 - mid_point_r - point_arrow_dis, y = 27.5, yend = 27.5), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 16, y0 = 27.5, r = mid_point_r), fill = mid_point_fill, color = mid_point_color, linewidth = mid_point_linewidth) + 
  geom_text(aes(x = 16, y = 28.2), label = "3-Isopropylmalate", size = mid_text_size, color = mid_text_color, hjust = mid_text_hjust, vjust = mid_text_vjust) +
  geom_segment(aes(x = 16, xend = 16, y = 27.5 - mid_point_r - point_arrow_dis, yend = 25.5 + rect_arrow_dis), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 14, xmax = 18, ymin = 24.5, ymax = 25.5), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(y = 25, x = 16), label = "K00052", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust) +
  geom_segment(aes(x = 18 + rect_arrow_dis, xend = 21 - mid_point_r - point_arrow_dis, y = 25, yend = 25), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 21, y0 = 25, r = mid_point_r), fill = mid_point_fill, color = mid_point_color, linewidth = mid_point_linewidth) +
  geom_text(aes(x = 21, y = 24), label = "2-Oxoisocaproate", size = mid_text_size, color = mid_text_color, hjust = mid_text_hjust, vjust = mid_text_vjust) +
  geom_segment(aes(x = 21, xend = 21, y = 25 + mid_point_r + point_arrow_dis, yend = 27 - rect_arrow_dis), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_rrect(aes(xmin = 19, xmax = 23, ymin = 27, ymax = 28), radius = rect_radius, fill = rect_fill, color = rect_color, linetype = rect_linetype, alpha = rect_alpha, size = rect_linewidth) +
  geom_text(aes(y = 27.5, x = 21), label = "K00826", size = KO_text_size, color = KO_text_color, hjust = KO_text_hjust, vjust = KO_text_vjust) +
  geom_segment(aes(x = 23 + rect_arrow_dis, xend = 26 - Final_point_r - point_arrow_dis, y = 27.5, yend = 27.5), linewidth = 0.75, linetype = "solid", arrow = arrow) +
  geom_circle(aes(x0 = 26, y0 = 27.5, r = Final_point_r), fill = Final_point_fill, color = Final_point_color, linewidth = Final_point_linewidth) +
  geom_text(aes(x = 26, y = 28.2), label = "Leucine", size = Final_text_size, color = Final_text_color, hjust = Final_text_hjust, vjust = Final_text_vjust) + 
  geom_text(aes(x = 15, y = 30), label = "M00432:Leucine biosynthesis", size = 8, color = "#005050", hjust = 0.5, vjust = 0.5) +
  # plot legend
  geom_circle(aes(x0 = 1, y0 = 34, r = Substrate_point_r), fill = Substrate_point_fill, color = Substrate_point_color, linewidth = Substrate_point_linewidth) +
  geom_text(aes(x = 2, y = 34), label = "Substrate", size = Substrate_text_size, color = Substrate_text_color, hjust = 0, vjust = 0.5) + 
  geom_circle(aes(x0 = 13.5, y0 = 34, r = mid_point_r), fill = mid_point_fill, color = mid_point_color, linewidth = mid_point_linewidth) +
  geom_text(aes(x = 14.5, y = 34), label = "Intermediate", size = 6, color = mid_text_color, hjust = 0, vjust = 0.5) +
  geom_circle(aes(x0 = 26, y0 = 34, r = Final_point_r), fill = Final_point_fill, color = Final_point_color, linewidth = Final_point_linewidth) +
  geom_text(aes(x = 27, y = 34), label = "Product", size = Final_text_size, color = Final_text_color, hjust = 0, vjust = 0.5) + 
  coord_fixed(ratio = 1) + 
  # plot link to metabolome & transcriptome
  geom_segment(aes(x = 26 + Final_point_r + rect_arrow_dis, xend = 30, y = 27.5, yend = 27.5), linewidth = 0.75, linetype = "dashed", lineend = "square") + 
  geom_segment(aes(x = 26 + Final_point_r + rect_arrow_dis, xend = 33 - rect_arrow_dis, y = 27.5 - 8 , yend = 27.5 - 8), linewidth = 0.75, linetype = "dashed", lineend = "square") + 
  geom_segment(aes(x = 26 + Final_point_r + rect_arrow_dis, xend = 30, y = 27.5 - 16 , yend = 27.5 - 16), linewidth = 0.75, linetype = "dashed", lineend = "square") + 
  geom_segment(aes(x = 30, xend = 30, y = 27.5, yend = 27.5 - 16), linewidth = 0.75, linetype = "dashed", lineend = "square") + 
  geom_segment(aes(x = 26 + Final_point_r + rect_arrow_dis, xend = 33 - rect_arrow_dis, y = 27.5 - 18.5 , yend = 27.5 - 18.5), linewidth = 0.75, linetype = "dashed", lineend = "square") +
  geom_rect(aes(xmin = 33, xmax = 39, ymin = 27.5 - 7.5 + 0.25, ymax = 27.5 - 8.5 - 0.25), fill = NA, color = "grey", linewidth = 1, lineend = "square") +
  geom_rect(aes(xmin = 33, xmax = 39, ymin = 27.5 - 7.5 - 10.5 + 0.25, ymax = 27.5 - 8.5 - 10.5 - 0.25), fill = NA, color = "grey", linewidth = 1, lineend = "square") + 
  geom_segment(aes(x = 39 + rect_arrow_dis, xend = 47 - rect_arrow_dis, y = 27.5 - 8 , yend = 27.5 - 8), linewidth = 0.75, linetype = "solid", lineend = "square", arrow = arrow2) +
  geom_segment(aes(x = 39 + rect_arrow_dis, xend = 47 - rect_arrow_dis, y = 27.5 - 8 -10.5 , yend = 27.5 - 8 - 10.5), linewidth = 0.75, linetype = "solid", lineend = "square", arrow = arrow2) + 
  geom_rect(aes(xmin = 47, xmax = 55, ymin = 27.5 - 7.5 + 0.25, ymax = 27.5 - 8.5 - 0.25), fill = NA, color = "grey", linewidth = 1, lineend = "square") +
  geom_rect(aes(xmin = 47, xmax = 55, ymin = 27.5 - 7.5 - 10.5 + 0.25, ymax = 27.5 - 8.5 - 10.5 - 0.25), fill = NA, color = "grey", linewidth = 1, lineend = "square") + 
  geom_text(aes(x = 51, y = 19.5), label = "mTOR mRNA", size = 6, color = mid_text_color, hjust = 0.5, vjust = 0.5) +
  geom_text(aes(x = 51, y = 9 ), label = "mTOR activity", size = 6, color = mid_text_color, hjust = 0.5, vjust = 0.5) + 
  geom_segment(aes(x = 39 + rect_arrow_dis, xend = 47 - rect_arrow_dis, y = 27.5 - 8.5 - 0.25 - rect_arrow_dis, yend = 27.5 - 7.5 - 10.5 + 0.25 + rect_arrow_dis), linewidth = 0.75, linetype = "solid", lineend = "square", arrow = arrow2) + 
  geom_segment(aes(x = 55 + rect_arrow_dis, xend = 57, y = 27.5 - 8, yend = 27.5 - 8), linewidth = 0.75, linetype = "dashed", lineend = "square") +
  geom_segment(aes(x = 55 + rect_arrow_dis, xend = 57, y = 27.5 - 8 - 10.5, yend = 27.5 - 8 - 10.5), linewidth = 0.75, linetype = "dashed", lineend = "square") + 
  geom_segment(aes(x = 57, xend = 57, y = 27.5 - 8, yend = 27.5 - 8 - 10.5), linewidth = 0.75, linetype = "dashed", lineend = "square") + 
  geom_segment(aes(x = 57, xend = 60, y = 14.25, yend = 14.25), linewidth = 0.75, linetype = "solid", lineend = "square") + 
  geom_segment(aes(x = 51, xend = 51, y = 27.5 - 8 + 0.75 + rect_arrow_dis, yend = 27.5 - 8 + 0.75 + rect_arrow_dis + 3), linewidth = 0.75, linetype = "solid", lineend = "square", arrow = arrow2) + 
  geom_segment(aes(x = 51, xend = 51, y = 27.5 - 8 - 10.5 - 0.75 - rect_arrow_dis, yend = 27.5 - 8 - 10.5 - 0.75 - rect_arrow_dis - 3), linewidth = 0.75, linetype = "solid", lineend = "square", arrow = arrow2) +
  geom_segment(aes(x = 51, xend = 51, y = 27.5 - 8 - 10.5 - 0.75 - rect_arrow_dis, yend = 27.5 - 8 - 10.5 - 0.75 - rect_arrow_dis - 3), linewidth = 0.75, linetype = "solid", lineend = "square", arrow = arrow2) +
  geom_segment(aes(x = 60, xend = 64, y = 14.25, yend = 17.25), color = "pink", alpha = 0.8791304, linewidth = 0.75, linetype = "solid", lineend = "square", arrow = arrow) + 
  geom_segment(aes(x = 60, xend = 64, y = 14.25, yend = 15.75), color = "pink", alpha = 0.6000000, linewidth = 0.75, linetype = "solid", lineend = "square", arrow = arrow) +
  geom_segment(aes(x = 60, xend = 64, y = 14.25, yend = 14.25), color = "pink", alpha = 0.4747826, linewidth = 0.75, linetype = "solid", lineend = "square", arrow = arrow) + 
  geom_segment(aes(x = 60, xend = 64, y = 14.25, yend = 12.75), color = "lightblue", alpha = 0.2930435, linewidth = 0.75, linetype = "solid", lineend = "square", arrow = arrow) + 
  geom_segment(aes(x = 60, xend = 64, y = 14.25, yend = 11.25), color = "lightblue", alpha = 0.6295652, linewidth = 0.75, linetype = "solid", lineend = "square", arrow = arrow) + 
  geom_circle(aes(x0 = 60, y0 = 14.25, r = Final_point_r), fill = Final_point_fill, color = Final_point_color, linewidth = Final_point_linewidth)

print(BCAA_pathway_plot)

ggplot2::ggsave("1_multi_omic_association/46_Process of biological signal molecule transmission.svg", plot = BCAA_pathway_plot, width = 40, height = 100, units = "cm", dpi = 300)
ggplot2::ggsave("1_multi_omic_association/46_Process of biological signal molecule transmission.tiff", plot = BCAA_pathway_plot, width = 40, height = 100, units = "cm", dpi = 60)





































