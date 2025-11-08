# Weiling Li (wel4007@med.cornell.edu)
# Fig. 7

library(readr)
library(forcats)
library(dplyr)
library(ggplot2)
library(ggpubr)

library(tidyverse)
library(forcats)
library(readr)

# Define the pathways of interest
paths <- c(
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS",
  "TCGA_ANDROGEN_RESPONSE",
  "HALLMARK_MYC_TSARGETS_V2",
  "BELTRAN_NEPC_UP",
  "HALLMARK_PANCREAS_BETA_CELLS",
  "HALLMARK_KRAS_SINGALING_DN",
  "SMID_BREAST_CANCER_LUMINAL_A_UP"
)

setwd('Users/wel4007/gsea_home/output/')
#####################
# Read data
enriched_pos <- read_tsv('gsea_home/output/jul23/ar_allrest26_more.Gsea.1753305747559/gsea_report_for_CRPC-AR_1753305747559.tsv')
enriched_neg <- read_tsv('gsea_home/output/jul23/ar_allrest26_more.Gsea.1753305747559/gsea_report_for_REST_1753305747559.tsv')

# Filter paths
pos_sel <- enriched_pos %>% filter(NAME %in% paths)
neg_sel <- enriched_neg %>% filter(NAME %in% paths)

# Columns to ensure numeric for p-values
cols_to_numeric <- c("FDR q-val", "FWER p-val", "NOM p-val")

# Convert relevant columns to numeric with suppressWarnings to avoid warnings
pos_sel <- pos_sel %>%
  mutate(
    NES_num = suppressWarnings(as.numeric(NES)),
    NES = if_else(is.na(NES_num), 0, NES_num),
    across(all_of(cols_to_numeric), ~ suppressWarnings(as.numeric(.)))
  ) %>%
  select(-NES_num)

neg_sel <- neg_sel %>%
  mutate(
    NES_num = suppressWarnings(as.numeric(NES)),
    NES = if_else(is.na(NES_num), 0, NES_num),
    across(all_of(cols_to_numeric), ~ suppressWarnings(as.numeric(.)))
  ) %>%
  select(-NES_num)

# Add significance and group labels
pos_sel <- pos_sel %>%
  mutate(
    significant = if_else(
      !is.na(NES) & `FDR q-val` <= 0.05 & `FWER p-val` <= 0.05 & `NOM p-val` <= 0.05,
      TRUE,
      FALSE
    ),
    group = if_else(significant, "pos_sig", "nonsig"),
    dataset = "pos"
  )

neg_sel <- neg_sel %>%
  mutate(
    significant = if_else(
      !is.na(NES) & `FDR q-val` <= 0.05 & `FWER p-val` <= 0.05 & `NOM p-val` <= 0.05,
      TRUE,
      FALSE
    ),
    group = if_else(significant, "neg_sig", "nonsig"),
    dataset = "neg"
  )

# Combine datasets
combined <- bind_rows(pos_sel, neg_sel)

# Prepare plot data
plot_data <- combined %>%
  mutate(
    name = NAME,
    val = NES,
    fill_group = group
  ) %>%
  arrange(val) %>%
  mutate(name = factor(name, levels = unique(name)))

n_bars <- nrow(plot_data)
# Plot
ar <- ggplot(plot_data, aes(x = name, y = val, fill = fill_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.7, width = 0.7) +
  scale_fill_manual(
    values = c(
      pos_sig = "black",
      neg_sig = "black",
      nonsig = "grey"
    ),
    labels = c(
      pos_sig = "Positive Significant (≤ 0.05)",
      neg_sig = "Negative Significant (≤ 0.05)",
      nonsig = "Not Significant"
    )
  ) +
  #coord_flip() +
  ylab("Normalized Enrichment Score") +
  xlab(NULL) +
  coord_flip(ylim = c(min(plot_data$val, na.rm = TRUE), max(plot_data$val, na.rm = TRUE) * 1.1)) +
  ggtitle("CRPC-Adeno") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )
#ggsave("Downloads/AR_vs_rest_plot2.pdf", plot = ar, width = 8, height =  n_bars, limitsize = FALSE)
print(ar)

#######
# Read data
enriched_pos <- read_tsv('gsea_home/output/jul23/ne_allrest26_more.Gsea.1753305817769/gsea_report_for_CRPC-NE_1753305817769.tsv')
enriched_neg <- read_tsv('gsea_home/output/jul23/ne_allrest26_more.Gsea.1753305817769/gsea_report_for_REST_1753305817769.tsv')

# Filter paths
pos_sel <- enriched_pos %>% filter(NAME %in% paths)
neg_sel <- enriched_neg %>% filter(NAME %in% paths)

# Columns to ensure numeric for p-values
cols_to_numeric <- c("FDR q-val", "FWER p-val", "NOM p-val")

# Convert relevant columns to numeric with suppressWarnings to avoid warnings
pos_sel <- pos_sel %>%
  mutate(
    NES_num = suppressWarnings(as.numeric(NES)),
    NES = if_else(is.na(NES_num), 0, NES_num),
    across(all_of(cols_to_numeric), ~ suppressWarnings(as.numeric(.)))
  ) %>%
  select(-NES_num)

neg_sel <- neg_sel %>%
  mutate(
    NES_num = suppressWarnings(as.numeric(NES)),
    NES = if_else(is.na(NES_num), 0, NES_num),
    across(all_of(cols_to_numeric), ~ suppressWarnings(as.numeric(.)))
  ) %>%
  select(-NES_num)

# Add significance and group labels
pos_sel <- pos_sel %>%
  mutate(
    significant = if_else(
      !is.na(NES) & `FDR q-val` <= 0.05 & `FWER p-val` <= 0.05 & `NOM p-val` <= 0.05,
      TRUE,
      FALSE
    ),
    group = if_else(significant, "pos_sig", "nonsig"),
    dataset = "pos"
  )

neg_sel <- neg_sel %>%
  mutate(
    significant = if_else(
      !is.na(NES) & `FDR q-val` <= 0.05 & `FWER p-val` <= 0.05 & `NOM p-val` <= 0.05,
      TRUE,
      FALSE
    ),
    group = if_else(significant, "neg_sig", "nonsig"),
    dataset = "neg"
  )

# Combine datasets
combined <- bind_rows(pos_sel, neg_sel)

# Prepare plot data
plot_data <- combined %>%
  mutate(
    name = NAME,
    val = NES,
    fill_group = group
  ) %>%
  arrange(val) %>%
  mutate(name = factor(name, levels = unique(name)))

n_bars <- nrow(plot_data)
# Plot
ne <- ggplot(plot_data, aes(x = name, y = val, fill = fill_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.7, width = 0.7) +
  scale_fill_manual(
    values = c(
      pos_sig ="black",# "blue",
      neg_sig = "black",#"lightblue",
      nonsig = "grey"
    ),
    labels = c(
      pos_sig = "Positive Significant (≤ 0.05)",
      neg_sig = "Negative Significant (≤ 0.05)",
      nonsig = "Not Significant"
    )
  ) +
  #coord_flip() +
  ylab("Normalized Enrichment Score") +
  xlab(NULL) +
  coord_flip(ylim = c(min(plot_data$val, na.rm = TRUE), max(plot_data$val, na.rm = TRUE) * 1.1)) +
  ggtitle("CRPC-NE") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )
#ggsave("Downloads/NE_vs_rest_plot2.pdf", plot = ne, width = 8, height =  n_bars, limitsize = FALSE)
print(ne)
##########
# Read data
enriched_pos <- read_tsv('gsea_home/output/jul23/mbc_allrest26_more.Gsea.1753305858072/gsea_report_for_MBC_1753305858072.tsv')
enriched_neg <- read_tsv('gsea_home/output/jul23/mbc_allrest26_more.Gsea.1753305858072/gsea_report_for_REST_1753305858072.tsv')

# Filter paths
pos_sel <- enriched_pos %>% filter(NAME %in% paths)
neg_sel <- enriched_neg %>% filter(NAME %in% paths)

# Columns to ensure numeric for p-values
cols_to_numeric <- c("FDR q-val", "FWER p-val", "NOM p-val")

# Convert relevant columns to numeric with suppressWarnings to avoid warnings
pos_sel <- pos_sel %>%
  mutate(
    NES_num = suppressWarnings(as.numeric(NES)),
    NES = if_else(is.na(NES_num), 0, NES_num),
    across(all_of(cols_to_numeric), ~ suppressWarnings(as.numeric(.)))
  ) %>%
  select(-NES_num)

neg_sel <- neg_sel %>%
  mutate(
    NES_num = suppressWarnings(as.numeric(NES)),
    NES = if_else(is.na(NES_num), 0, NES_num),
    across(all_of(cols_to_numeric), ~ suppressWarnings(as.numeric(.)))
  ) %>%
  select(-NES_num)

# Add significance and group labels
pos_sel <- pos_sel %>%
  mutate(
    significant = if_else(
      !is.na(NES) & `FDR q-val` <= 0.05 & `FWER p-val` <= 0.05 & `NOM p-val` <= 0.05,
      TRUE,
      FALSE
    ),
    group = if_else(significant, "pos_sig", "nonsig"),
    dataset = "pos"
  )

neg_sel <- neg_sel %>%
  mutate(
    significant = if_else(
      !is.na(NES) & `FDR q-val` <= 0.05 & `FWER p-val` <= 0.05 & `NOM p-val` <= 0.05,
      TRUE,
      FALSE
    ),
    group = if_else(significant, "neg_sig", "nonsig"),
    dataset = "neg"
  )

# Combine datasets
combined <- bind_rows(pos_sel, neg_sel)

# Prepare plot data
plot_data <- combined %>%
  mutate(
    name = NAME,
    val = NES,
    fill_group = group
  ) %>%
  arrange(val) %>%
  mutate(name = factor(name, levels = unique(name)))

n_bars <- nrow(plot_data)
# Plot
mbc <- ggplot(plot_data, aes(x = name, y = val, fill = fill_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.7, width = 0.7) +
  scale_fill_manual(
    values = c(
      pos_sig = "black",
      neg_sig = "black", 
      nonsig = "grey"
    ),
    labels = c(
      pos_sig = "Positive Significant (≤ 0.05)",
      neg_sig = "Negative Significant (≤ 0.05)",
      nonsig = "Not Significant"
    )
  ) +
  coord_flip() +
  ylab("Normalized Enrichment Score") +
  xlab(NULL) +
  ylim(min(plot_data$val, na.rm = TRUE), max(plot_data$val, na.rm = TRUE) * 1.1) + # dynamic y limits
  ggtitle("MBC") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )
#ggsave("Downloads/MBC_vs_rest_plot2.pdf", plot = mbc, width = 8, height =  n_bars, limitsize = FALSE)
print(mbc)

#######
# oct
enriched_pos <- read_tsv('gsea_home/output/jun18/mbc27.2controls2pdx.correction.Gsea.1750361963846/gsea_report_for_1_1750361963846.tsv')
enriched_neg <- read_tsv('gsea_home/output/jun18/mbc27.2controls2pdx.correction.Gsea.1750361963846/gsea_report_for_0_1750361963846.tsv')

# Filter paths
pos_sel <- enriched_pos %>% filter(NAME %in% paths)
neg_sel <- enriched_neg %>% filter(NAME %in% paths)

# Columns to ensure numeric for p-values
cols_to_numeric <- c("FDR q-val", "FWER p-val", "NOM p-val")

# Convert relevant columns to numeric with suppressWarnings to avoid warnings
pos_sel <- pos_sel %>%
  mutate(
    NES_num = suppressWarnings(as.numeric(NES)),
    NES = if_else(is.na(NES_num), 0, NES_num),
    across(all_of(cols_to_numeric), ~ suppressWarnings(as.numeric(.)))
  ) %>%
  select(-NES_num)

neg_sel <- neg_sel %>%
  mutate(
    NES_num = suppressWarnings(as.numeric(NES)),
    NES = if_else(is.na(NES_num), 0, NES_num),
    across(all_of(cols_to_numeric), ~ suppressWarnings(as.numeric(.)))
  ) %>%
  select(-NES_num)

# Add significance and group labels
pos_sel <- pos_sel %>%
  mutate(
    significant = if_else(
      !is.na(NES) & `FDR q-val` <= 0.05 & `FWER p-val` <= 0.05 & `NOM p-val` <= 0.05,
      TRUE,
      FALSE
    ),
    group = if_else(significant, "pos_sig", "nonsig"),
    dataset = "pos"
  )

neg_sel <- neg_sel %>%
  mutate(
    significant = if_else(
      !is.na(NES) & `FDR q-val` <= 0.05 & `FWER p-val` <= 0.05 & `NOM p-val` <= 0.05,
      TRUE,
      FALSE
    ),
    group = if_else(significant, "neg_sig", "nonsig"),
    dataset = "neg"
  )

# Combine datasets
combined <- bind_rows(pos_sel, neg_sel)

# Prepare plot data
plot_data <- combined %>%
  mutate(
    name = NAME,
    val = NES,
    fill_group = group
  ) %>%
  arrange(val) %>%
  mutate(name = factor(name, levels = unique(name)))

n_bars <- nrow(plot_data)
# Plot
mbc <- ggplot(plot_data, aes(x = name, y = val, fill = fill_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.7, width = 0.7) +
  scale_fill_manual(
    values = c(
      pos_sig = "black",
      neg_sig = "black",
      nonsig = "grey"
    ),
    labels = c(
      pos_sig = "Positive Significant (≤ 0.05)",
      neg_sig = "Negative Significant (≤ 0.05)",
      nonsig = "Not Significant"
    )
  ) +
  #coord_flip() +
  ylab("Normalized Enrichment Score") +
  xlab(NULL) +
  coord_flip(ylim = c(min(plot_data$val, na.rm = TRUE), max(plot_data$val, na.rm = TRUE) * 1.1)) +
  ggtitle("MBC") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )

#ggsave("Downloads/MBC_bypdxpred.pdf", plot = mbc, width = 8, height =  n_bars, limitsize = FALSE)
print(mbc)
