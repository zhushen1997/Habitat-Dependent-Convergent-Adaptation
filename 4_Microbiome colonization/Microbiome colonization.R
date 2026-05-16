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
library(MetaNet)
library(igraph)
library(dplyr)
library(pcutils)
library(ggfun)

# load Arial font
font_import(pattern = "arial")
loadfonts()
windowsFonts()

### growth curve
paletteer_d("MoMAColors::Connors")
paletteer_d("nationalparkcolors::Arches")

grouth_curve <- 
  read.table("grouth curve.txt", header = T) %>%
  mutate(exp = c(rep("Bacteria-BCAAs Axis", 24), rep("mTOR-Grouth Axis", 18)))

grouth_curve_long <- 
  grouth_curve %>% 
  pivot_longer(-c("group", "time", "exp"), names_to = "sample", values_to = "weight") %>%
  dplyr::select(-c(sample))

grouth_curve_long$weight_single <- grouth_curve_long$weight / 20 * 1000
grouth_curve_long$time <-factor(grouth_curve_long$time, levels = unique(grouth_curve_long$time))
grouth_curve_long$group <-factor(grouth_curve_long$group, levels = unique(grouth_curve_long$group))
grouth_curve_long$top_label <- "Grouth Curve"


grouth_color <- c("ABX50" = "#729ECEFF", "ABX20" = "#67BF5CFF", "ABX20+FMT" = "#AD8BC9FF", "ABX50+BCAAs" = "#FF9E4AFF", 
                  "LY294002" = "#CDCC5DFF", "Rapamycin" = "#ED665DFF","CR" = "#A2A2A2FF")

grouth_curve_plot <- 
  ggplot(grouth_curve_long %>% filter(!time %in% c("0", "10")),
         aes(x = group, y = weight_single) 
  )+ 
  geom_point(aes(fill = group, color = group),
             # color = "black",
             stroke = 0.1,
             size = 1,
             shape = 21,
             position = position_jitter(height = 0, width = 0.2, seed = 123),
             show.legend = F
  ) +
  geom_boxplot(aes(fill = group, color = group),
               # color = "black",
               outlier.shape = NA, 
               alpha = 0.3,
               width = 0.85,
               staplewidth = 0.4,
               linewidth = 0.25,
               show.legend = F
  ) +
  # annotate("segment", x = 4.5, xend = 4.5 , y = -Inf, yend = Inf, linetype = 2, linewidth = 0.25 , color = "black") + 
  ggh4x::facet_nested_wrap(~ top_label + time, 
                           scales = "free_y",
                           ncol = 6
  ) + 
  # ylim(c(0, 55)) + 
  stat_compare_means(
    comparisons = rev(list(c("LY294002", "CR"),
                           c("Rapamycin", "CR"), c("LY294002", "Rapamycin"))),
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
  stat_compare_means(
    comparisons = rev(list( 
                       c("ABX50", "ABX50+BCAAs"),
                       c("ABX20", "ABX50+BCAAs"), c("ABX50", "ABX20+FMT"), 
                       c("ABX20+FMT", "ABX50+BCAAs"), c("ABX20", "ABX20+FMT"), c("ABX50", "ABX20"))),
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
  labs(y = "Body weight(mg)") + 
  scale_fill_manual(values = grouth_color) + 
  scale_color_manual(values = grouth_color) + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill = NA, color = "black", linetype = 1, linewidth = 0.25),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.key.size = unit(10, "pt"), 
    legend.key.width = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.text = element_blank(),
    strip.background = ggfun::element_roundrect(
      # fill = "#DEEFF5",
      fill = "#c9caca", 
      color = NA, 
      linewidth = 0, 
      r = unit(0.25, "snpc")
    ),
    strip.text = element_text(family = "Arial", face = "plain", size=8, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black"),
  )

grouth_curve_plot

ggsave("1_grouth_curve_plot.svg", plot = grouth_curve_plot, width = 13, height = 4, units = "cm", dpi = 300)
ggsave("1_grouth_curve_plot.tiff", plot = grouth_curve_plot, width = 13, height = 4, units = "cm", dpi = 300)














strip_colors <- c("#94c7c0", "#94c7c0")

paletteer::paletteer_d("MetBrewer::Benedictus")

p_growth <- 
  ggplot(grouth_curve_long, 
         aes(x = time, y = weight_single, color = group, group = group)
  ) +
  stat_summary(fun = mean, 
               geom = "line", 
               linewidth = 0.15,
               show.legend = F
  ) +
  stat_summary(fun = mean, 
               geom = "point", 
               size = 0.5,
               show.legend = F
  ) +
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               width = 0.1,
               linewidth = 0.15,
               show.legend = F
  ) +
  facet_wrap2(~ exp, 
              scales = "free_y", 
              ncol = 2,
              strip = strip_themed(background_x = elem_list_rect(fill = strip_colors))
  ) +
  scale_fill_manual(values = grouth_color) + 
  scale_color_manual(values = grouth_color) +
  labs(
    x = "Developmental Time (Days)", 
    y = "Weight per Individual (mg)"
  ) +
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black"), 
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"),
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill = NA, color = "black", linetype = 1, linewidth = 0.25),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.key.size = unit(10, "pt"), 
    legend.key.width = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.text = element_blank(),
    strip.background = element_rect(fill = NA, color = "black", linetype = 1, linewidth = 0.25),
    strip.text = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black")
  )

print(p_growth)

ggsave("1_grouth_curve_plot_2.svg", plot = p_growth, width = 12.88, height = 2, units = "cm", dpi = 300)
ggsave("1_grouth_curve_plot_2.tiff", plot = p_growth, width = 12.88, height = 2, units = "cm", dpi = 300)



### bubble plot

grouth_summary <- 
  grouth_curve_long %>%
  group_by(exp, time, group) %>%
  summarise(mean_weight = mean(weight_single, na.rm = TRUE), .groups = "drop")

grouth_summary$group <- factor(grouth_summary$group, levels = rev(levels(grouth_curve_long$group)))


p_growth_bubble <- 
  ggplot(grouth_summary, 
         aes(x = time, y = group, size = mean_weight, fill = group, color = group)
  ) +
  geom_point(shape = 21, 
             stroke = 0.25, 
             alpha = 0.8,
             show.legend = T
  ) +
  scale_size_continuous(range = c(0.25, 4), 
                        name = "Mean Weight (mg)"
  ) +
  facet_wrap2(~ exp, 
              scales = "free_y", 
              ncol = 2,
              strip = strip_themed(background_x = elem_list_rect(fill = strip_colors))
  ) + 
  scale_fill_manual(values = grouth_color) + 
  scale_color_manual(values = grouth_color) + 
  labs(x = "", y = "") +
  theme_void() +
  theme(
    axis.text.x = element_text(family = "Arial", size = 6, color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x = unit(2, "pt"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(family = "Arial", face = "plain", size = 6, color = "black", lineheight = 1),
    legend.position = "none"
  )

print(p_growth_bubble)

ggsave("1_grouth_curve_plot_3.svg", plot = p_growth_bubble, width = 8, height = 1.5, units = "cm", dpi = 300)
ggsave("1_grouth_curve_plot_3.tiff", plot = p_growth_bubble, width = 8, height = 1.5, units = "cm", dpi = 300)


library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)

exp1 <- grouth_curve_long %>% filter(exp == "Bacteria-BCAAs Axis")
exp1$group <- factor(exp1$group, levels = c("ABX50", "ABX20", "ABX20+FMT", "ABX50+BCAAs"))
exp2 <- grouth_curve_long %>% filter(exp == "mTOR-Grouth Axis")

m1 <- lm(weight ~ group * time, 
           data = exp1)

summary(m1)
anova(m1)

emm1 <- emmeans(m1, ~ group | time)
pairs(emm1, adjust = "tukey")
cld(emm1, Letters = letters, adjust = "tukey")

mean_df1 <- 
  exp1 %>%
  group_by(group, time) %>%
  summarise(mean_weight = mean(weight), .groups = "drop")

start_end1 <- 
  mean_df1 %>%
  filter(time %in% c("0", "50")) %>%
  tidyr::pivot_wider(names_from = time, values_from = mean_weight)

day0_val1 <- mean_df1 %>% filter(time == "0") %>% pull(mean_weight)
day50_val1 <- mean_df1 %>% filter(time == "50") %>% pull(mean_weight)
group_names1 <- mean_df1 %>% filter(time == "0") %>% pull(group)

growth_rate1 <- (day50_val1 - day0_val1) / 50
rate_df1 <- data.frame(group = group_names1, rate = growth_rate1)

grouth_rate1 <- 
ggplot(rate_df1, 
       aes(x = rev(group), y = rate, fill = group)
  ) +
  geom_col(aes(color = group , fill = group),
           width = 0.6
  ) +
  scale_fill_manual(values = grouth_color) + 
  scale_color_manual(values = grouth_color) + 
  coord_flip() + 
  theme_void() +
  theme(
    axis.text.x = element_text(family = "Arial", size = 6, color = "black", hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x = unit(2, "pt"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(family = "Arial", face = "plain", size = 6, color = "black", lineheight = 1),
    legend.position = "none"
  )

grouth_rate1

ggsave("1_grouth_curve_plot_4.svg", plot = grouth_rate1 , width = 1.5, height = 1.5, units = "cm", dpi = 300)
ggsave("1_grouth_curve_plot_4.tiff", plot = grouth_rate1 , width = 1.5, height = 1.5, units = "cm", dpi = 300)




m2 <- lm(weight ~ group * time, 
         data = exp2)

summary(m2)
anova(m2)

emm2 <- emmeans(m2, ~ group | time)
pairs(emm2, adjust = "tukey")
cld(emm2, Letters = letters, adjust = "sidak")

mean_df2 <- 
  exp2 %>%
  group_by(group, time) %>%
  summarise(mean_weight = mean(weight), .groups = "drop")

start_end2 <- 
  mean_df2 %>%
  filter(time %in% c("0", "50")) %>%
  tidyr::pivot_wider(names_from = time, values_from = mean_weight)

day0_val2 <- mean_df2 %>% filter(time == "0") %>% pull(mean_weight)
day50_val2 <- mean_df2 %>% filter(time == "50") %>% pull(mean_weight)
group_names2 <- mean_df2 %>% filter(time == "0") %>% pull(group)

growth_rate2 <- (day50_val2 - day0_val2) / 50
rate_df2 <- data.frame(group = group_names2, rate = growth_rate2)

grouth_rate2 <- 
  ggplot(rate_df2, 
         aes(x = rev(group), y = rate, fill = group)
  ) +
  geom_col(aes(color = group , fill = group),
           width = 0.6
  ) +
  scale_fill_manual(values = grouth_color) + 
  scale_color_manual(values = grouth_color) + 
  coord_flip() + 
  theme_void() +
  theme(
    axis.text.x = element_text(family = "Arial", size = 6, color = "black", hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1), 
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.length.x = unit(2, "pt"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(family = "Arial", face = "plain", size = 6, color = "black", lineheight = 1),
    legend.position = "none"
  )

grouth_rate2

ggsave("1_grouth_curve_plot_5.svg", plot = grouth_rate2 , width = 1.5, height = 1.5, units = "cm", dpi = 300)
ggsave("1_grouth_curve_plot_5.tiff", plot = grouth_rate2 , width = 1.5, height = 1.5, units = "cm", dpi = 300)








### western bar plot(different treatment pf Pame)
wb_result <- 
  read.table("wb/wb.txt", header = T, row.names = 1) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample")
wb_result$`p-S6K/Tubulin` <- wb_result$`p-S6K` / wb_result$Tubulin
wb_result$`p-AKT/AKT` <- wb_result$`p-AKT` / wb_result$AKT
wb_result$`p-AKT/Tubulin` <- wb_result$`p-AKT` / wb_result$Tubulin
wb_result$group <- factor(c(rep("CR", 3), rep("ABX50", 3), rep("ABX20", 3), rep("ABX20+FMT", 3), rep("ABX50+BCAAs", 3), rep("LY294002", 3), rep("Rapamycin", 3)), 
                          levels = c("ABX50", "ABX20", "ABX20+FMT", "ABX50+BCAAs", "LY294002", "Rapamycin", "CR"))

wb_result_long <- 
  wb_result %>%
  dplyr::select(-c(2,3,4,5)) %>%
  pivot_longer(-c(sample, group), names_to = "singnal", values_to = "relativa intensity")

wb_result_long$singnal <- factor(wb_result_long$singnal , levels = c("p-S6K/Tubulin", "p-AKT/AKT", "p-AKT/Tubulin"))

wb_result_long_group <- 
  wb_result_long %>%
  group_by(singnal, group) %>%
  summarise(
    mean = mean(`relativa intensity`),
    sd = sd(`relativa intensity`)
  )

wb_plot <- 
  ggplot() +
  #geom_boxplot(data = wb_result_long,
  #             aes(x = group, y =`relativa intensity`, fill = group, color = group),
  #             outlier.shape = NA, 
  #             alpha = 0.3,
  #             width = 0.85,
  #             staplewidth = 0.4,
  #             linewidth = 0.25,
  #             show.legend = F
  #) +
  geom_point(data = wb_result_long,
             aes(x = group, y =`relativa intensity`, fill = group,  color = group),
             shape = 21,
             size = 1,
             position = position_jitter(height = 0, width = 0.2, seed = 123),
             show.legend = F
  ) +
  geom_col(data = wb_result_long_group,
           aes(x = group, y = mean, color = group, fill = group),
           alpha = 0.3, 
           linewidth = 0.25,
           width = 0.75,
           show.legend = F
  ) +
  geom_errorbar(data = wb_result_long_group,
                aes(x = group, ymin = mean - sd, ymax = mean + sd, color = group),
                width = 0.3,
                linewidth = 0.25,
                linetype = 1,
                show.legend = F
  ) +
  facet_wrap(~ singnal, 
             scale = "free_y",
             nrow = 3
  ) + 
  geom_signif(data = wb_result_long,
              aes(x = group, y = `relativa intensity`), 
              comparisons = list(c("LY294002", "Rapamycin"), c("Rapamycin", "CR"), c("LY294002", "CR")),
              map_signif_level = function(p) {
                ifelse(p < 0.001, "***", 
                       ifelse(p < 0.01, "**", 
                              ifelse(p < 0.05, "*", "ns")))
              },
              test = "t.test",
              test.args = list(exact = FALSE), 
              textsize = 3,
              vjust = 0.6,
              tip_length = 0,
              size = 0.25,
              step_increase = 0.3
  ) +
  geom_signif(data = wb_result_long,
              aes(x = group, y = `relativa intensity`), 
              comparisons = list(c("ABX50", "ABX20"), c("ABX20", "ABX20+FMT"), c("ABX20+FMT", "ABX50+BCAAs"),
                                 c("ABX50", "ABX20+FMT"), c("ABX20", "ABX50+BCAAs"), 
                                 c("ABX50", "ABX50+BCAAs")),
              map_signif_level = function(p) {
                ifelse(p < 0.001, "***", 
                       ifelse(p < 0.01, "**", 
                              ifelse(p < 0.05, "*", "ns")))
              },
              test = "t.test",
              test.args = list(exact = FALSE), 
              textsize = 3,
              vjust = 0.6,
              tip_length = 0,
              size = 0.25,
              step_increase = 0.3
  ) +
  coord_cartesian(ylim = c(0, NA)) + 
  scale_fill_manual(values = c("ABX50" = "#729ECEFF", "ABX20" = "#67BF5CFF", "ABX20+FMT" = "#AD8BC9FF", "ABX50+BCAAs" = "#FF9E4AFF", "LY294002" = "#CDCC5DFF", "Rapamycin" = "#ED665DFF","CR" = "#A2A2A2FF")) +
  scale_color_manual(values = c("ABX50" = "#729ECEFF", "ABX20" = "#67BF5CFF", "ABX20+FMT" = "#AD8BC9FF", "ABX50+BCAAs" = "#FF9E4AFF", "LY294002" = "#CDCC5DFF", "Rapamycin" = "#ED665DFF","CR" = "#A2A2A2FF")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=1, vjust=1, angle=45, lineheight=1, color="black"),  
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"), 
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

wb_plot

ggplot2::ggsave("2_wb_plot.svg", plot = wb_plot, width = 3, height = 12, units = "cm", dpi = 300)
ggplot2::ggsave("2_wb_plot.tiff", plot = wb_plot, width = 3, height = 12, units = "cm", dpi = 300)




# BCAA Assay(Reagent kit for total BCAA)

library(ggplot2)
library(ggpmisc)

ST_curve <- 
  read.table("ST_curve.txt", header = T, sep = "\t") %>%
  pivot_longer(-c(BCAA.concentrate), names_to = "repetition", values_to = "OD450")


Measured_data <- 
  read.table("Measured_data.txt", header = T, sep = "\t") %>%
  mutate(exp = rep("Hemolymph BCAA", 21))

Measured_data$group <- factor(Measured_data$group, levels = c("ABX50", "ABX20", "ABX20+FMT", "ABX50+BCAAs", "LY294002", "Rapamycin", "CR"))

ST_curve_plot <-
  ggplot(data = ST_curve, aes(y = BCAA.concentrate, x = OD450)) + 
  geom_point(size = 0.5, show.legend = F) + 
  geom_smooth(method = 'lm', 
              formula = 'y ~ x', 
              se = TRUE, 
              size = 0.25
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*")),
    formula = y ~ x, 
    parse = TRUE,
    rr.digits = 4,
    coef.digits = 4,
    label.x = "left",
    label.y = "top",
    size = 2
  ) + 
  labs(y = "BCAA concentrate(μM)", x = "OD450") + 
  geom_point(data = Measured_data, 
             aes(y = BCAA_concentrate, x = OD450, color = group, fill = group, shape = group), 
             size = 1.5, 
             alpha = 0.75,
             show.legend = T
  ) + 
  scale_fill_manual(values = c("ABX50" = "#729ECEFF", "ABX20" = "#67BF5CFF", "ABX20+FMT" = "#AD8BC9FF", "ABX50+BCAAs" = "#FF9E4AFF", "LY294002" = "#CDCC5DFF", "Rapamycin" = "#ED665DFF", "CR" ="#A2A2A2FF")) + 
  scale_color_manual(values = c("ABX50" = "#729ECEFF", "ABX20" = "#67BF5CFF", "ABX20+FMT" = "#AD8BC9FF", "ABX50+BCAAs" = "#FF9E4AFF", "LY294002" = "#CDCC5DFF", "Rapamycin" = "#ED665DFF", "CR" ="#A2A2A2FF")) + 
  scale_shape_manual(values = c("ABX50" = 16, "ABX20" = 17, "ABX20+FMT" = 18, "ABX50+BCAAs" = 19, "LY294002" = 20, "Rapamycin" = 21, "CR" = 22)) +
  theme(
    axis.title.x = element_text(family = "Arial", face = "plain", size = 8, hjust=0.5, vjust=0.5, angle= 0, lineheight=1, color="black"),
    axis.title.y = element_text(family = "Arial", face = "plain", size = 8, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"),
    axis.text.x = element_text(family = "Arial", face = "plain", size = 6, hjust=0.5, vjust=0.5, angle= 0, lineheight=1, color="black"),
    axis.text.y = element_text(family = "Arial", face = "plain", size = 6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"),
    axis.ticks.x = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.ticks.y = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = "black"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(fill = NA, color = "black", linetype = 1, linewidth = 0.25),
    panel.background = element_blank(),
    panel.grid = element_line(color = "grey90", linewidth = 0.10),
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.35), 
    legend.background = element_blank(),
    legend.text = element_text(family = "Arial", face = "plain", size = 4, hjust= 0 , vjust=0.5, angle= 0, lineheight=1, color="black"),
    legend.title = element_blank(),
    legend.justification = c(0.5, 0.5),
    legend.key.size = unit(10, "pt"), 
    legend.key.width = unit(10, "pt"), 
    legend.key.height = unit(10, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    strip.background = element_rect(fill = "#F6955EFF", color = "black", linewidth = 0.25, linetype = "solid"),
    strip.text = element_text(family = "Arial", face = "plain", size=8, hjust=0.5, vjust=0.5, angle=0, lineheight=1, color="black")
  )

print(ST_curve_plot)

ggplot2::ggsave("3_Hemolymph_BCAA_ST_curve.svg", plot = ST_curve_plot, width = 6, height = 6, units = "cm", dpi = 300)
ggplot2::ggsave("3_Hemolymph_BCAA_ST_curve.tiff", plot = ST_curve_plot, width = 6, height = 6, units = "cm", dpi = 300)




### bar lot 

Measured_data$log_concentrate <- log10(Measured_data$BCAA_concentrate)

Measured_data_grouped <- 
  Measured_data %>% 
  group_by(group) %>% 
  summarise(
    n = n(), 
    mean_BCAA = mean(log_concentrate, na.rm = TRUE), 
    sd_BCAA = sd(log_concentrate, na.rm = TRUE), 
    se_BCAA = sd_BCAA / sqrt(n)
  ) 

Hemolymph_BCAA_bar_plot <-
  ggplot() + 
  geom_point(data = Measured_data,
             aes(x = group, y = log_concentrate, fill = group, color = group),
             stroke = 0.1,
             size = 1,
             shape = 21,
             position = position_jitter(height = 0, width = 0.2, seed = 123),
             show.legend = F
  ) +
  geom_col(data = Measured_data_grouped,
           aes(x = group, y = mean_BCAA, color = group, fill = group),
           alpha = 0.3, 
           linewidth = 0.25,
           width = 0.75,
           show.legend = F
  ) +
  geom_errorbar(data = Measured_data_grouped,
                aes(x = group, ymin = mean_BCAA - sd_BCAA, ymax = mean_BCAA + sd_BCAA, color = group),
                width = 0.3,
                linewidth = 0.25,
                linetype = 1,
                show.legend = F
  ) +
  facet_wrap(~ exp,
             scale = "free_y",
             nrow = 1
  ) + 
  annotate("segment", x = 4.5, xend = 4.5 , y = -Inf, yend = Inf, linetype = 2, linewidth = 0.25 , color = "black") + 
  geom_signif(data = Measured_data,
              aes(x = group, y = log_concentrate), 
              comparisons = list(c("LY294002", "Rapamycin"), c("Rapamycin", "CR"), c("LY294002", "CR")),
              map_signif_level = function(p) {
                ifelse(p < 0.001, "***", 
                       ifelse(p < 0.01, "**", 
                              ifelse(p < 0.05, "*", "ns")))
              },
              test = "t.test",
              test.args = list(exact = FALSE), 
              textsize = 3,
              vjust = 0.7,
              tip_length = 0,
              size = 0.25,
              step_increase = 0.25
  ) +
  geom_signif(data = Measured_data,
              aes(x = group, y = log_concentrate), 
              comparisons = list(c("ABX50", "ABX20"), c("ABX20", "ABX20+FMT"), c("ABX20+FMT", "ABX50+BCAAs"),
                                 c("ABX50", "ABX20+FMT"), c("ABX20", "ABX50+BCAAs"), 
                                 c("ABX50", "ABX50+BCAAs")),
              map_signif_level = function(p) {
                ifelse(p < 0.001, "***", 
                       ifelse(p < 0.01, "**", 
                              ifelse(p < 0.05, "*", "ns")))
              },
              test = "t.test",
              test.args = list(exact = FALSE), 
              textsize = 3,
              vjust = 0.7,
              tip_length = 0,
              size = 0.25,
              step_increase = 0.25
  ) +
  coord_cartesian(ylim = c(0, NA)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     labels = scales::label_math(10^.x)
  ) + 
  scale_fill_manual(values = c("ABX50" = "#729ECEFF", "ABX20" = "#67BF5CFF", "ABX20+FMT" = "#AD8BC9FF", "ABX50+BCAAs" = "#FF9E4AFF", "LY294002" = "#CDCC5DFF", "Rapamycin" = "#ED665DFF", "CR" ="#A2A2A2FF")) + 
  scale_color_manual(values = c("ABX50" = "#729ECEFF", "ABX20" = "#67BF5CFF", "ABX20+FMT" = "#AD8BC9FF", "ABX50+BCAAs" = "#FF9E4AFF", "LY294002" = "#CDCC5DFF", "Rapamycin" = "#ED665DFF", "CR" ="#A2A2A2FF")) + 
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial", face = "plain", size=6, hjust=1, vjust=1, angle=45, lineheight=1, color="black"),  
    axis.text.y = element_text(family = "Arial", face = "plain", size=6, hjust=0.5, vjust=0.5, angle=90, lineheight=1, color="black"), 
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

print(Hemolymph_BCAA_bar_plot)

ggplot2::ggsave("3_Hemolymph_BCAA_concentrate.svg", plot = Hemolymph_BCAA_bar_plot, width = 3, height = 5, units = "cm", dpi = 300)
ggplot2::ggsave("3_Hemolymph_BCAA_concentrate.tiff", plot = Hemolymph_BCAA_bar_plot, width = 3, height = 5, units = "cm", dpi = 300)


### bubble plot
Hemolymph_BCAA_bubble_plot <- 
  ggplot(data = Measured_data_grouped) + 
  geom_point(aes(x = group, y = 1, size = mean_BCAA, fill = group, color = group),
             stroke = 0.05,
             shape = 21,
             show.legend = T
  ) + 
  scale_fill_manual(values = c("ABX50" = "#729ECEFF", "ABX20" = "#67BF5CFF", "ABX20+FMT" = "#AD8BC9FF", "ABX50+BCAAs" = "#FF9E4AFF", "LY294002" = "#CDCC5DFF", "Rapamycin" = "#ED665DFF", "CR" ="#A2A2A2FF")) + 
  scale_color_manual(values = c("ABX50" = "#729ECEFF", "ABX20" = "#67BF5CFF", "ABX20+FMT" = "#AD8BC9FF", "ABX50+BCAAs" = "#FF9E4AFF", "LY294002" = "#CDCC5DFF", "Rapamycin" = "#ED665DFF", "CR" ="#A2A2A2FF")) +
  scale_size_continuous(range = c(1, 3)) + 
  scale_x_discrete(position = "top") + 
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
    legend.ticks = element_line(linewidth = 0.25, linetype = "solid", lineend = "butt", color = NA)
  )

Hemolymph_BCAA_bubble_plot















