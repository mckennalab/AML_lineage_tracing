##### packages #####
setwd("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/other/pAML_geo_data_042825")
set.seed(2468)
`%nin%` = Negate(`%in%`)

library(BiocManager)
library(tibble)
library(ggplot2)
library(stringr)
library(tidyr)
library(Seurat)
library(dplyr)
library(cowplot)
library(rstatix)
library(ggpubr)


##### violin plots (figure 6A) #####
#load sketched seurat object (integrated following sketch workflow)
merged_pt <- readRDS("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/other/pAML_geo_data_042825/GSE235063_data/all_samples_merged_sketched_rpcaint_proj_110825.rds")

dx_subset <- subset(merged_pt, timepoint %in% "diagnosis")
data <- FetchData(dx_subset, vars = c("resist_lin1", "Lambo_et_al_ID", "Subgroup", "Malignant"))

#all subtypes together 
pairwise.test <- data %>% wilcox_test(resist_lin1 ~ Malignant)
ggplot(data, aes(x=Malignant, y=resist_lin1)) + 
  geom_violin() +
  geom_boxplot(width = 0.25) + theme_cowplot() +
  stat_pvalue_manual(
    pairwise.test, label = "p", 
    y.position = 0.5
  ) + scale_y_continuous(limits = c(-0.2, 0.5)) + 
  ggtitle("All patients")


#run for each subtype (figure 6A) or cell type (figure S8)
pairwise.test <- data[data$Subgroup %in% "KMT2A",] %>% wilcox_test(resist_lin1 ~ Malignant)
ggplot(data[data$Subgroup %in% "KMT2A",], aes(x=Malignant, y=resist_lin1)) + 
  geom_violin() +
  geom_boxplot(width = 0.25) + theme_cowplot() +
  stat_pvalue_manual(
    pairwise.test, label = "p", 
    y.position = 0.5
  ) + scale_y_continuous(limits = c(-0.2, 0.5)) + 
  ggtitle("KMT2A")

#adjust pvals (after running for all subtpes)
pvals <- data.frame(pvals = c(0,0,0,0,5.56E-07))
pvals$padj <- p.adjust(pvals$pvals, method="BH")

#same process for Mumme et al data (harmony integrated instead of sketched)

##### scatter plots (figure 6B, S9A-f) #####
data <- FetchData(malig_pt, vars = c("orig.ident", "Lambo_et_al_ID", "timepoint", "Overall_survival_days", "Disease_free_days", "Progression_free_survival_days", "Risk_Group", "aml_blasttype", "Phase", "Classified_Celltype", "Subgroup", "resist_lin1", "sig_wo_kmt2a_genes1"))

pt_map <- data[, c(1:7,11)] %>% distinct()
pt_prop <- data %>% group_by(orig.ident, Classified_Celltype) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n), hscs_per_100 = (100 * n)/sum(n))
pt_prop <- pt_prop %>% inner_join(pt_map, by = "orig.ident")
pt_prop_min <- pt_prop[pt_prop$Lambo_et_al_ID %nin% c("AML3", "AML11", "AML14", "AML20"), c(2,4:7)] #remove patients without matched dx relapse samples

pt_prop_min_wide <- pt_prop_min %>% pivot_wider(names_from = timepoint, values_from = c(freq, hscs_per_100)) 
pt_prop_min_wide[is.na(pt_prop_min_wide)] <- 0
pt_prop_min_wide <- pt_prop_min_wide %>% mutate(diff = hscs_per_100_relapse - hscs_per_100_diagnosis, code = paste0(Lambo_et_al_ID, "_", Classified_Celltype))
pt_prop_min_wide <- pt_prop_min_wide[, -c(1,2)]

sig_by_blasttype <- data %>% group_by(Lambo_et_al_ID, timepoint, Classified_Celltype) %>% 
  summarize(mean = mean(resist_lin1), max = max(resist_lin1))
sig_by_blasttype_wide <- sig_by_blasttype %>% pivot_wider(names_from = timepoint, values_from = c(mean, max))
sig_by_blasttype_wide <- sig_by_blasttype_wide %>% mutate(code = paste0(Lambo_et_al_ID, "_", Classified_Celltype))

full_df <- pt_prop_min_wide %>% right_join(sig_by_blasttype_wide, by = "code")
pt_orig_map <- data[, c(2,4:7,11)] %>% distinct()
full_df <- full_df %>% left_join(pt_orig_map, by = "Lambo_et_al_ID")
full_df_hsc <- full_df[full_df$Classified_Celltype == "HSC", ]

ggplot(full_df_hsc, aes(x = mean_diagnosis, y = diff)) + 
  geom_point(aes(colour = Subgroup)) +
  geom_smooth(method = "lm", se = TRUE) + stat_cor(method="pearson")
