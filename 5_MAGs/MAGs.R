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
library(ANCOMBC)
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






### 

set.seed(123)
mags_data <- read.table("I://1_Article/four cockroach/Figure1/Figure1a/bac_tree/bac_bin_qa.txt", header = TRUE, row.names = 1)


mags_data$Quality <- ifelse(
  mags_data$Completeness >= 90 & mags_data$Contamination <= 5,
  "Near-complete",
  ifelse(
    mags_data$Completeness >= 90 & mags_data$Contamination <= 10,
    "High-quality",
    "Medium-quality"
  )
)

mags_data$Quality <- factor(mags_data$Quality, levels = c("Near-complete", "High-quality", "Medium-quality"))


legend_labels <- c(
  "Near-complete" = "Near-complete\nCompleteness >= 90%\nContamination <= 10%\n5S + 16S + 23S\n>=18 tRNA",
  "High-quality" = "High-quality\nCompleteness >= 90%\nContamination <= 10%",
  "Medium-quality" = "Medium-quality\nCompleteness >= 50%\nContamination <= 10%"
)

mags_data$Quality <- factor(mags_data$Quality,
                            levels = names(legend_labels),
                            labels = legend_labels)


g1 <- 
  ggplot(mags_data, 
         aes(x = Completeness, y = Contamination, color = Quality)
  ) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Near-complete\nCompleteness >= 90%\nContamination <= 10%\n5S + 16S + 23S\n>=18 tRNA" = "#80B1D3",
                                "High-quality\nCompleteness >= 90%\nContamination <= 10%" = "#FDB462",
                                "Medium-quality\nCompleteness >= 50%\nContamination <= 10%" = "#8DD3C7"),
                     labels = legend_labels
  ) +
  xlim(50, 100) +
  ylim(-0, 10) +
  labs(x = "Completeness (%)", 
       y = "Contamination (%)", 
       color = "Quality Category"
  ) +
  guides(color = guide_legend(override.aes = list(size = 8), position = "right"))+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 2, linetype = "solid"),
    plot.title = element_text(family = "Arial", face = "bold", size = 20, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1.5, color = "black"),
    axis.title = element_text(family = "Arial", face = "bold", size = 20, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
    axis.text.x = element_text(family = "Arial", face = "bold", size = 20, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
    axis.text.y = element_text(family = "Arial", face = "bold", size = 20, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
    axis.ticks = element_line(linewidth = 2.0, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "right",
    legend.key.size = unit(15, "pt"),
    legend.text = element_text(family = "Arial", face = "bold", size = 20, hjust = 0, vjust = 0, angle = 0, lineheight = 1.2, color = "black"),
    legend.text.position = "bottom",
    legend.title = element_text(family = "Arial", face = "bold", size = 20, hjust = 0, vjust = 0, angle = 0, lineheight = 1.2, color = "black"),
    legend.title.position = "top",
    legend.key.spacing.y = unit(20, "pt")
  )

g1

library(ggExtra)

g1_with_hist <- 
  ggMarginal(
  g1,
  type        = "densigram",
  bins        = 40,
  groupFill   = TRUE,
  size        = 5
)

g1_with_hist


ggsave("I://1_Article/four cockroach/Figure1/Figure1a/genome_quanilty.svg", plot = g1_with_hist, width = 30, height = 30, units = "cm", dpi = 1200)
ggsave("I://1_Article/four cockroach/Figure1/Figure1a/genome_quanilty.tiff", plot = g1_with_hist, width = 30, height = 30, units = "cm", dpi = 1200)









quality_dist <- mags_data %>%
  count(Quality) %>%
  mutate(Percentage = n / sum(n) * 100,
         Label = paste0(Quality, "\n", round(Percentage, 1), "% (n=", n, ")"))

quality_dist$Quality <- factor(
  quality_dist$Quality,
  levels = c(
    "Near-complete\nCompleteness >= 90%\nContamination <= 10%\n5S + 16S + 23S\n>=18 tRNA",
    "High-quality\nCompleteness >= 90%\nContamination <= 10%",
    "Medium-quality\nCompleteness >= 50%\nContamination <= 10%"
  )
)


pie_chart <- 
  ggplot(quality_dist, 
         aes(x = 2, y = n, fill = Quality)
  ) +
  geom_bar(stat = "identity", 
           width = 1, 
           color = "white"
  ) +
  scale_fill_manual(values = c("Near-complete\nCompleteness >= 90%\nContamination <= 10%\n5S + 16S + 23S\n>=18 tRNA" = "#80B1D3",
                               "High-quality\nCompleteness >= 90%\nContamination <= 10%" = "#FDB462",
                               "Medium-quality\nCompleteness >= 50%\nContamination <= 10%" = "#8DD3C7")) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = sprintf("%.2f%%", Percentage), x = 2 + 0.5),  
            position = position_stack(vjust = 0.5),
            color = "black",
            fontface = "bold",
            size = 20
  ) +
  scale_x_continuous(limits = c(0.1, 3)) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
    legend.position = "none",
  )

pie_chart


ggsave("Figure1a/pie_chart.svg", plot = pie_chart, width = 30, height = 30, units = "cm", dpi = 1200)
ggsave("Figure1a/pie_chart.tiff", plot = pie_chart, width = 30, height = 30, units = "cm", dpi = 1200)







library(ggbeeswarm)
library(ggnewscale) 
library(ggh4x)

plot_data <- mags_data %>%
  dplyr::select(Completeness, Contamination) %>%
  pivot_longer(everything(), names_to = "Metric", values_to = "Value")


p <- 
  ggplot(plot_data,
       aes(x = Value, y = "")
  ) +
  new_scale_colour() +
  geom_jitter(aes(colour = Metric),
              width  = 0.05,
              height = 0.2,
              size   = 2,
              alpha  = 0.7) +
  scale_colour_manual(values = c(
    Completeness  = "#007D82FF",
    Contamination = "#C38961FF"
  )) +
  guides(colour = "none") +
  new_scale_colour() +
  geom_violin(aes(fill = Metric), alpha = 0.5) +
  scale_fill_manual(values = c(
    Completeness  = "#004042FF",
    Contamination = "#9F5630FF"
  )) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_boxplot(aes(fill = Metric, alpha = 1),
               width = 0.1,
               outlier.shape = NA) +
  scale_fill_manual(values = c(
    Completeness  = "grey",
    Contamination = "grey"
  )) +
  guides(fill = "none") +
  facet_wrap(~ Metric, nrow = 2, scales = "free_x", strip.position = "top") +
  labs(title = "",
       x = NULL,
       y = "Value (%)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.border =  element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"), 
        strip.text = element_text(family = "Arial", face = "bold", size = 20, hjust = 0, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(family = "Arial", face = "bold", size = 20, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1.2, color = "black"),
        axis.ticks.x = element_line(linewidth = 2.0, linetype = "solid", lineend = "butt", color = "black"),
        axis.ticks.x.length = unit(0.2, "cm")
        )

p

p2 <- p + facetted_pos_scales(
  x = list(
    Metric == "Completeness" ~ scale_x_continuous(breaks = c(0,50,100),limits = c(0, 100)),
    Metric == "Contamination" ~ scale_x_continuous(breaks = c(0,50,100),limits = c(0, 100))
  )
)
  
p2

ggsave("Figure1a/box.svg", plot = p2, width = 30, height = 15, units = "cm", dpi = 1200)
ggsave("Figure1a/box.tiff", plot = p2, width = 30, height = 15, units = "cm", dpi = 1200)









library(tidyr)
library(dplyr)


BCAA_module <- c("M00019", "M00432", "M00535", "M00570")

BCAA_KO <- c("K01754")

BCAA_TR_protein <- c("K26605", "K26606", "K11250")

KO_data <- read.table("G://bin_2_ko.txt")

KO_data_wide <- 
  KO_data %>%
  count(V1, V2) %>%
  pivot_wider(
    names_from = V2,
    values_from = n,
    values_fill = 0
  ) %>% 
  dplyr::select(1, BCAA_KO, BCAA_TR_protein)

KEGG_module_completeness <- read.table("G://pathway_coverage.tsv", header = T) %>% dplyr::select(1, BCAA_module)

KEGG_BCAA <- 
  KEGG_module_completeness %>% left_join(KO_data_wide, by = c("id_genome" = "V1")) %>%
  mutate(M00535_modify = case_when(
    as.numeric(M00535) >= 0.75 ~ as.numeric(M00535),
    TRUE ~ 0)) %>%
  mutate(M00019_modify = case_when(
    as.numeric(M00019) >= 0.75 ~ as.numeric(M00019),
    TRUE ~ 0)) %>%
  mutate(M00570_modify = case_when(
    as.numeric(M00570) >= 0.75 & as.numeric(K01754) != 0 ~ as.numeric(M00570),
    TRUE ~ 0)) %>%
  mutate(M00432_modify = case_when(
    as.numeric(M00432) >= 0.75 ~ as.numeric(M00432),
    TRUE ~ 0))
  

write.table(KEGG_BCAA, "G://KEGG_BCAA.txt", sep = "\t", quote = F, row.names = F)








### MAGs_module_coverage_plot

MAGs_module_completeness_data <- read.table("Figure1a/bac_tree/pathway_coverage_complete_info.tsv", header = T, sep = "\t")

MAGs_module_completeness_data_clean <- 
  MAGs_module_completeness_data %>% 
  mutate(bin_quanlity = case_when(
    Completeness >= 90 & Contamination <= 10 ~ "HQ",
    Completeness >= 50 & Completeness < 90 & Contamination <= 10 ~ "MQ"
  )) %>%
  separate(classification, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  separate(Bin.Id, into = c("group", "bin", "number"), sep = "_", remove = F) %>%
  dplyr::select(-c("bin", "number")) %>%
  # filter(bin_quanlity == "HQ") %>%
  arrange(phylum, class, order, family, genus, species, desc(M00019), desc(M00535), desc(M00570), desc(M00432))

MAGs_module_completeness_data_clean$Bin.Id <-factor(MAGs_module_completeness_data_clean$Bin.Id, levels = MAGs_module_completeness_data_clean$Bin.Id )
MAGs_module_completeness_data_clean$group <- factor(MAGs_module_completeness_data_clean$group, levels = c("Esin", "Bger", "Pame", "Pful"))
MAGs_module_completeness_data_clean$phylum <- factor(MAGs_module_completeness_data_clean$phylum, levels = unique(MAGs_module_completeness_data_clean$phylum))

Module_info <- 
  read.table("Figure1a/bac_tree/3_KEGG_module_information.tsv", header = T, sep = "\t") %>%
  arrange(category_B, category_C, module_ID)
  
Module_info$module_ID <- factor(Module_info$module_ID, levels = Module_info$module_ID)
Module_info$category_B <- factor(Module_info$category_B, levels = unique(Module_info$category_B))
Module_info$category_C <- factor(Module_info$category_C, levels = unique(Module_info$category_C))

MAGs_module_completeness_data_clean_long <- 
  MAGs_module_completeness_data_clean %>%
  pivot_longer(-c(1:11, "bin_quanlity"), names_to = "module_ID", values_to = "module_compltetness") %>%
  left_join(Module_info, by = "module_ID")
  
Esin_MAGs_module_completeness_data_clean_long <-
  MAGs_module_completeness_data_clean_long %>%
  filter(group == "Esin" & module_ID %in% c("M00535", "M00019", "M00570", "M00432")) 

Bger_MAGs_module_completeness_data_clean_long <-
  MAGs_module_completeness_data_clean_long %>%
  filter(group == "Bger" & module_ID %in% c("M00535", "M00019", "M00570", "M00432")) 

Pame_MAGs_module_completeness_data_clean_long <-
  MAGs_module_completeness_data_clean_long %>%
  filter(group == "Pame" & module_ID %in% c("M00535", "M00019", "M00570", "M00432")) 

Pful_MAGs_module_completeness_data_clean_long <-
  MAGs_module_completeness_data_clean_long %>%
  filter(group == "Pful" & module_ID %in% c("M00535", "M00019", "M00570", "M00432")) 


p_Esin <- 
  ggplot(Esin_MAGs_module_completeness_data_clean_long) +
  geom_tile(aes(x = module_ID, y = Bin.Id, fill = module_compltetness)) +
  scale_fill_gradientn(
    colours = c("white", "white", "white", "white", "lightblue", "yellow"),
    values  = c(0, 0.25, 0.5, 0.749, 0.75, 1) 
  )

p_Esin

                      p_Bger <- 
  ggplot(Bger_MAGs_module_completeness_data_clean_long) +
  geom_tile(aes(x = module_ID, y = Bin.Id, fill = module_compltetness)) +
  scale_fill_gradientn(
    colours = c("white", "white", "white", "white", "lightblue", "yellow"),
    values  = c(0, 0.25, 0.5, 0.749, 0.75, 1) 
  )

p_Bger


p_Pame <- 
  ggplot(Pame_MAGs_module_completeness_data_clean_long) +
  geom_tile(aes(x = module_ID, y = Bin.Id, fill = module_compltetness)) +
  scale_fill_gradientn(
    colours = c("white", "white", "white", "white", "lightblue", "yellow"),
    values  = c(0, 0.25, 0.5, 0.749, 0.75, 1) 
  )

p_Pame


p_Pful <- 
  ggplot(Pful_MAGs_module_completeness_data_clean_long) +
  geom_tile(aes(x = module_ID, y = Bin.Id, fill = module_compltetness)) +
  scale_fill_gradientn(
    colours = c("white", "white", "white", "white", "lightblue", "yellow"),
    values  = c(0, 0.25, 0.5, 0.749, 0.75, 1) 
  )

p_Pful





p_phylum <- 
  ggplot(MAGs_module_completeness_data_clean_long) +
  geom_tile(aes(x = 0, y = Bin.Id, fill = phylum)) +
  theme_void()

p_phylum

p_class <- 
  ggplot(MAGs_module_completeness_data_clean_long) +
  geom_tile(aes(x = 0, y = Bin.Id, fill = class)) +
  theme_void()

p_class

p_group <- 
  ggplot(MAGs_module_completeness_data_clean_long) +
  geom_tile(aes(x = 0, y = Bin.Id, fill = group)) +
  scale_fill_manual(values = c("Esin" = Esin_color, "Bger" = Bger_color, "Pame" = Pame_color, "Pful" = Pful_color)) + 
  theme_void()

p_group

p_module_category_B <- 
  ggplot(Module_info) +
  geom_tile(aes(x = module_ID, y = 0, fill = category_B)) +
  theme_void()

p_module_category_B

p_module_category_C <- 
  ggplot(Module_info) +
  geom_tile(aes(x = module_ID, y = 0, fill = category_C)) +
  theme_void()

p_module_category_C

p_module <- 
  p_Esin %>% 
  # insert_left(p_class, width = 0.01) %>%
  insert_left(p_phylum, width = 0.01) %>%
  insert_left(p_group, width = 0.01) %>%
  insert_top(p_module_category_C, height = 0.01) %>%
  insert_top(p_module_category_B, height = 0.01)

p_module


ggsave()






