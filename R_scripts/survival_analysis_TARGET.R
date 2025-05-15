# Survival analysis 

##### packages #####
library(ggplot2)
library(survival)
library(corrplot)
library(devtools)
library(ShortRead) 
library(data.table)
library(dplyr)
library(stringr)
library(BiocManager)
library(bayesbio)
library(grid)
library(cowplot)
library(org.Hs.eg.db)
library(biomaRt)



##### load modules, get human orthologs for C1498 modules #####
hl60_mods <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS30_RS31_humancellEP/GEX/hotspot/hl60_megamega_modlist.csv")

# split modules into individual lists, clean up (depending on module csv format)
mega_mod1 <- c(hl60_mods$Module.1)
mega_mod1 <- mega_mod1[mega_mod1 != ""]

mega_mod2 <- c(hl60_mods$Module.2)
mega_mod2 <- mega_mod2[mega_mod2 != ""]

mega_mod3 <- c(hl60_mods$Module.3)
mega_mod3 <- mega_mod3[mega_mod3 != ""]

mega_mod4 <- c(hl60_mods$Module.4)
mega_mod4 <- mega_mod4[mega_mod4 != ""]

mega_mod5 <- c(hl60_mods$Module.5)
mega_mod5 <- mega_mod5[mega_mod5 != ""]

mega_mod6 <- c(hl60_mods$Module.6)
mega_mod6 <- mega_mod6[mega_mod6 != ""]

mega_mod7 <- c(hl60_mods$Module.7)
mega_mod7 <- mega_mod7[mega_mod7 != ""]

mega_mod8 <- c(hl60_mods$Module.8)
mega_mod8 <- mega_mod8[mega_mod8 != ""]

mega_mod9 <- c(hl60_mods$Module.9)
mega_mod9 <- mega_mod9[mega_mod9 != ""]

mega_mod10 <- c(hl60_mods$Module.10)
mega_mod10 <- mega_mod10[mega_mod10 != ""]

mega_mod11 <- c(hl60_mods$Module.11)
mega_mod11 <- mega_mod11[mega_mod11 != ""]

mega_mod12 <- c(hl60_mods$Module.12)
mega_mod12 <- mega_mod12[mega_mod12 != ""]

mega_mod13 <- c(hl60_mods$Module.13)
mega_mod13 <- mega_mod13[mega_mod13 != ""]

mega_mod14 <- c(hl60_mods$Module.14)
mega_mod14 <- mega_mod14[mega_mod14 != ""]

mega_mod15 <- c(hl60_mods$Module.15)
mega_mod15 <- mega_mod15[mega_mod15 != ""]

mega_mod16 <- c(hl60_mods$Module.16)
mega_mod16 <- mega_mod16[mega_mod16 != ""]

mega_mod17 <- c(hl60_mods$Module.17)
mega_mod17 <- mega_mod17[mega_mod17 != ""]

mega_mod18 <- c(hl60_mods$Module.18)
mega_mod18 <- mega_mod18[mega_mod18 != ""]

# load C1498 modules and convert to human gene orthologs
c1498_mods <- read.csv("/dartfs/rc/lab/M/McKennaLab/projects/hannah/aml/analysis/c1498_lineage_NEW/rs22_rs28/trees_res1/HotSpot/241113_rs22_rs28_combo_joined_rm_MT-RB/24114_rs22_rs28_combo_joined_rm_MT_RB_hotspot_modules_annotated.csv")
f1_mega_mods <- c1498_mods[c1498_mods$clone == "F1_Vanilla_rmWT_megatree", ]
f6_mega_mods <- c1498_mods[c1498_mods$clone == "F6_Vanilla_rmWT_megatree", ]
orthologs <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
mouse <- split.data.frame(orthologs,orthologs$Common.Organism.Name)[[2]]
human <- split.data.frame(orthologs,orthologs$Common.Organism.Name)[[1]]
mouse <- mouse[,c(1,4)]
human <- human[,c(1,4)]
mh_data <- merge.data.frame(mouse,human,by = "DB.Class.Key",all.y = TRUE) 
mh_data <- mh_data[,-1]

make_hs_list = function(df) {
  names(df) <- "gene"
  df <- df %>% inner_join(mh_data, by = c('gene' = 'Symbol.x'))
  df <- c(df$Symbol.y)
}
f1m1 <- make_hs_list(as.data.frame(f1_mega_mods$Gene[f1_mega_mods$Module == 1]))
f1m2 <- make_hs_list(as.data.frame(f1_mega_mods$Gene[f1_mega_mods$Module == 2]))
f1m3 <- make_hs_list(as.data.frame(f1_mega_mods$Gene[f1_mega_mods$Module == 3]))
f1m4 <- make_hs_list(as.data.frame(f1_mega_mods$Gene[f1_mega_mods$Module == 4]))
f6m1 <- make_hs_list(as.data.frame(f6_mega_mods$Gene[f6_mega_mods$Module == 1]))
f6m2 <- make_hs_list(as.data.frame(f6_mega_mods$Gene[f6_mega_mods$Module == 2]))
f6m3 <- make_hs_list(as.data.frame(f6_mega_mods$Gene[f6_mega_mods$Module == 3]))
f6m4 <- make_hs_list(as.data.frame(f6_mega_mods$Gene[f6_mega_mods$Module == 4]))
f6m5 <- make_hs_list(as.data.frame(f6_mega_mods$Gene[f6_mega_mods$Module == 5]))

# load pre-resistant signature genes
pre_resist_sig <- c("CLSTN2", "SPATA6", "SGCZ", "PRICKLE1", "IGF2BP2", "CAMK2D", "MEIS1", "KCND2", "MKX", "BMP2", "TRHDE", "APBA2", "PPFIA2", "CPNE8", "PON2", "LAMP5", "LNCAROD", "CPLANE1", "HOXB-AS3", "NAV3", "SDK1", "IRF8", "CES1", "PTPRN2", "HIST1H1D", "NCAM1", "LHX2", "SASH1", "AK1", "ELANE", "PPDPF", "UNCX", "CFD", "AC107223.1", "EXT1", "BEX3", "CLEC11A", "LINC02169", "TNS3", "TTC28", "KLHL29", "CKAP4", "AL163541.1", "PRDX2", "BEX1", "MAP1LC3A", "HOXA9", "PSD3", "AC090796.1", "SERPINB6", "WDR49", "ANXA2", "CD36", "MEF2C", "SLC35F1", "HNMT", "AL713998.1", "MRAS", "DAPK1", "DTNA", "ACAA2", "HLA-B", "MGMT", "ASGR2", "STARD13")


##### prepping TARGET-AML data #####
fn <- "https://cbioportal-datahub.s3.amazonaws.com/aml_target_gdc.tar.gz"
download.file(fn,destfile="tmp.tar.gz")
untar("tmp.tar.gz",list=TRUE)  ## check contents

untar("tmp.tar.gz",files="aml_target_gdc/data_mrna_seq_fpkm_zscores_ref_all_samples.txt")
target_mrna_fpkm_zscores <- read.csv("aml_target_gdc/data_mrna_seq_fpkm_zscores_ref_all_samples.txt", sep = "\t")

untar("tmp.tar.gz",files="aml_target_gdc/data_clinical_sample.txt")
clin_sample_data <- read.csv("aml_target_gdc/data_clinical_sample.txt", sep = "\t")
names(clin_sample_data) <- clin_sample_data[4,]
clin_sample_data <- clin_sample_data[5:2486,]

untar("tmp.tar.gz",files="aml_target_gdc/data_timeline_status.txt")
timeline <- read.csv("aml_target_gdc/data_timeline_status.txt", sep = "\t")

untar("tmp.tar.gz",files="aml_target_gdc/data_clinical_patient.txt")
clin_data <- read.csv("aml_target_gdc/data_clinical_patient.txt", sep = "\t")
names(clin_data) <- clin_data[4,]
clin_data <- clin_data[5:2350,]
clin_data_full <- clin_data %>% right_join(clin_sample_data, by = ("PATIENT_ID"))


# adding extra clin info from GDC
clin_val_info <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/other/c1498_hl60_module_analysis/TARGET_AML_clinvalidation_clindata.csv")
clin_val_info_min <- clin_val_info[, c(1,2,5,6,7,8,9,14,15,18,49,50,51,52,53,54,55,57,58,59,64)]
discovery_info <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/other/c1498_hl60_module_analysis/TARGET_AML_discovery_clindata.csv")
discovery_info_min <- discovery_info[, c(1,2,5,6,7,8,9,14,15,18,49,50,51,52,53,54,55,57,58,59,64)]
#get rid of duplicate rows, clean things up, then merge
clin_val_info_min$MRD...at.end.of.course.1[clin_val_info_min$MRD...at.end.of.course.1 == ""] <- NA
discovery_info_min$MRD...at.end.of.course.1[discovery_info_min$MRD...at.end.of.course.1 == ""] <- NA
clin_val_info_min$MRD...at.end.of.course.2[clin_val_info_min$MRD...at.end.of.course.2 == ""] <- NA
discovery_info_min$MRD...at.end.of.course.2[discovery_info_min$MRD...at.end.of.course.2 == ""] <- NA
clin_val_info_min$MRD...at.end.of.course.1[clin_val_info_min$MRD...at.end.of.course.1 == "N/A"] <- NA
discovery_info_min$MRD...at.end.of.course.1[discovery_info_min$MRD...at.end.of.course.1 == "N/A"] <- NA
clin_val_info_min$MRD...at.end.of.course.2[clin_val_info_min$MRD...at.end.of.course.2 == "N/A"] <- NA
discovery_info_min$MRD...at.end.of.course.2[discovery_info_min$MRD...at.end.of.course.2 == "N/A"] <- NA
discovery_info_min$Risk.group[discovery_info_min$Risk.group %in% "Unknown"] <- ""
clin_val_info_min$Risk.group[clin_val_info_min$Risk.group %in% "Unknown"] <- ""
new_clin_data <- rbind(clin_val_info_min, discovery_info_min) %>% distinct()
clin_data_full <- clin_data_full %>% left_join(new_clin_data, by = c("PATIENT_ID" = "TARGET.USI"))

#censorship data 
censor_data <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/other/c1498_hl60_module_analysis/clinical.project-target-aml.2024-12-25/follow_up.tsv", sep = "\t")
censor_data_min <- censor_data[censor_data$days_to_follow_up != "'--", c(2,19)]
censor_data_event <- censor_data[censor_data$first_event != "'--", c(2, 37)]
censor_data_top <- censor_data_min %>% group_by(case_submitter_id) %>% top_n(1, days_to_follow_up) %>% distinct()
censor_data_top <- censor_data_top %>% full_join(censor_data_event, by = "case_submitter_id")
clin_data_full <- clin_data_full %>% inner_join(censor_data_top, by = c("PATIENT_ID" = "case_submitter_id"))
clin_data_full$days_to_follow_up <- as.numeric(clin_data_full$days_to_follow_up)
clin_data_full$AGE <- as_numeric(as.character(clin_data_full$AGE))
clin_data_full$status <- 0
clin_data_full$status[clin_data_full$VITAL_STATUS == "Dead"] <- 1
clin_data_full$status[clin_data_full$first_event %in% c("Censored", "Second Malignant Neoplasm")] <- 0
clin_data_full <- clin_data_full %>% 
  mutate(months_to_follow_up = round(days_to_follow_up/30.437, digit=3))

# map gene names - sometimes the site is down
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), mart = ensembl)
mapping <- mapping %>% distinct()
write.csv(mapping, "/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/references/hs_entrez_hgnc_mapping.csv") #for when the site is down
target_zscores_matched <- target_mrna_fpkm_zscores %>% inner_join(mapping, by = c("Entrez_Gene_Id" = "entrezgene_id"))
target_zscores_matched <- target_zscores_matched[, -1]


##### score patients on expression #####
score_df <- as.data.frame(names(target_zscores_matched[, -2203]))
names(score_df) <- "patient"

# add patient scores for modules and gene sets (not a great method, but it works)
pre_resist_data <- target_zscores_matched[target_zscores_matched$hgnc_symbol %in% pre_resist_sig, ]
rownames(pre_resist_data) <- pre_resist_data$hgnc_symbol
pre_resist_data <- pre_resist_data[, -2203]
pre_resist_t <- t(pre_resist_data)
score_list <- c()
for (i in 1:nrow(pre_resist_t)) {
  x <- as.data.frame(t(pre_resist_t[i,]))
  x <- as.numeric(x)
  score <- sum(x, na.rm = TRUE)/sum(!is.na(x))
  score_list <- c(score_list, score)
}
pre_resist_min <- as.data.frame(rownames(pre_resist_t))
names(pre_resist_min) <- "patient"
pre_resist_min$pre_resist <- score_list
score_df <- score_df %>% left_join(pre_resist_min, by = "patient") 

prot_11s_data <- target_zscores_matched[target_zscores_matched$hgnc_symbol %in% c("PSME1", "PSME2"), ]
rownames(prot_11s_data) <- prot_11s_data$hgnc_symbol
prot_11s_data <- prot_11s_data[, -2203]
prot_11s_t <- t(prot_11s_data)
score_list <- c()
for (i in 1:nrow(prot_11s_t)) {
  x <- as.data.frame(t(prot_11s_t[i,]))
  x <- as.numeric(x)
  score <- sum(x, na.rm = TRUE)/sum(!is.na(x))
  score_list <- c(score_list, score)
}
prot_11s_min <- as.data.frame(rownames(prot_11s_t))
names(prot_11s_min) <- "patient"
prot_11s_min$prot_11s <- score_list
score_df <- score_df %>% left_join(prot_11s_min, by = "patient") #repeat for all gene sets of interest


#merge score_df with clinical info 
score_df$Sample_ID <- score_df$patient
score_df$Sample_ID <- gsub('\\.', '-', score_df$Sample_ID)
score_df <- score_df %>% inner_join(clin_data_full, by = c("Sample_ID" = "SAMPLE_ID"))

#remove the few samples with sorted samples (unclear how they fit in with the rest of the patient data)
score_df_nosort <- score_df[!grepl("Sort", score_df$PATIENT_ID),]
score_df_nosort <- score_df_nosort[!grepl("sort", score_df_nosort$PATIENT_ID),]
score_df <- score_df_nosort

# split into dfs based on sample type
post_tx_bm <- score_df[score_df$SAMPLE_TYPE %in% "Blood Derived Cancer - Bone Marrow, Post-treatment", ]
post_tx_blood <- score_df[score_df$SAMPLE_TYPE %in% "Blood Derived Cancer - Peripheral Blood, Post-treatment", ]
dx_bm <- score_df[score_df$SAMPLE_TYPE %in% "Primary Blood Derived Cancer - Bone Marrow", ] 
dx_blood <- score_df[score_df$SAMPLE_TYPE %in% "Primary Blood Derived Cancer - Peripheral Blood", ] 
all_bm <- score_df[score_df$SAMPLE_TYPE %in% c("Blood Derived Cancer - Bone Marrow, Post-treatment", "Primary Blood Derived Cancer - Bone Marrow"), ] 
all_blood <- score_df[score_df$SAMPLE_TYPE %in% c("Primary Blood Derived Cancer - Peripheral Blood", "Recurrent Blood Derived Cancer - Peripheral Blood"), ] 
all_dx <- score_df[score_df$SAMPLE_TYPE %in% c("Primary Blood Derived Cancer - Peripheral Blood", "Primary Blood Derived Cancer - Bone Marrow"), ] 


##### coxph, forest plots #####
model <- coxph(Surv(months_to_follow_up, status) ~ AGE + SEX + prot_11s + pre_resist, data = dx_bm)
ggforest(model, data = dx_bm) + theme_cowplot()

# with risk group 
dx_bm_w_riskcat <- dx_bm[!is.na(dx_bm$Risk.group), ]
# remove unclear risk group categories
dx_bm_w_riskcat <- dx_bm_w_riskcat[dx_bm_w_riskcat$Risk.group %nin% c("", "10", "30"), ]
dx_bm_w_riskcat$Risk.group <- factor(dx_bm_w_riskcat$Risk.group, levels = c("Standard Risk", "Low Risk", "High Risk"))

model <- coxph(Surv(months_to_follow_up, status) ~ AGE + SEX + Risk.group + prot_11s + pre_resist, data = dx_bm_w_riskcat)
ggforest(model, data = dx_bm_w_riskcat) + theme_cowplot()


##### kaplan meiers #####
sd <- sd(dx_bm$pre_resist)
mean <- mean(dx_bm$pre_resist)
highExpr <- mean + sd
lowExpr <- mean - sd
dx_bm$level <- ifelse(dx_bm$pre_resist >= highExpr, 'High',
                      ifelse(dx_bm$pre_resist <= lowExpr, 'Low', 'Mid'))
# relevel the factors to have mid as the ref level
dx_bm$level <- factor(dx_bm$level,
                      levels = c('Mid', 'Low', 'High'))

ggsurvplot(survfit(Surv(months_to_follow_up, status) ~ level,
                   data = dx_bm),
           data = dx_bm,
           risk.table = TRUE,
           pval = TRUE,
           break.time.by = 12,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE, 
           censor.shape=124, title = "f6_mod1")



##### expression level violin plots #####
# by outcome at 5 years post dx, need to pull from follow_up data
censor_data <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/other/c1498_hl60_module_analysis/clinical.project-target-aml.2024-12-25/follow_up.tsv", sep = "\t")
censor_data <- censor_data[censor_data$case_submitter_id %in% score_df$PATIENT_ID, ]
censor_data_timetoevent <- censor_data[censor_data$days_to_first_event != "'--", c(2,18)]
censor_data_event <- censor_data[censor_data$first_event != "'--", c(2, 37)]
first_event_info <- censor_data_timetoevent %>% full_join(censor_data_event, by = "case_submitter_id")
first_event_info$days_to_first_event <- as.numeric(first_event_info$days_to_first_event)
no_event_within5 <- first_event_info[first_event_info$days_to_first_event > 1826.25, ]
no_event_within5$yr5_status <- "DF"
first_event_within5 <- first_event_info[first_event_info$days_to_first_event < 1826.25, ]
first_event_within5 <- first_event_within5[first_event_within5$first_event %nin% "Censored", ]
first_event_within5$yr5_status <- first_event_within5$first_event
yr5_status <- rbind(no_event_within5, first_event_within5)
yr5_status$yr5_status[yr5_status$yr5_status %in% c("Death", "Death without Remission")] <- "Death"
score_df_5yr <- score_df %>% inner_join(yr5_status, by = c("PATIENT_ID" = "case_submitter_id"))
dx_bm_5yr <- score_df_5yr[score_df_5yr$SAMPLE_TYPE %in% "Primary Blood Derived Cancer - Bone Marrow", ] #check patients are all unique
dx_blood_5yr <- score_df_5yr[score_df_5yr$SAMPLE_TYPE %in% "Primary Blood Derived Cancer - Peripheral Blood", ] #check patients are all unique

library(rstatix)
library(ggbeeswarm)
pairwise.test <- dx_blood_5yr[dx_blood_5yr$yr5_status %nin% "Second Malignant Neoplasm", ] %>% wilcox_test(pre_resist ~ yr5_status, p.adjust.method = "bonferroni") 
pairwise.test <- pairwise.test[pairwise.test$p.adj.signif != "ns", ]
ggplot(dx_blood_5yr[dx_blood_5yr$yr5_status %nin% "Second Malignant Neoplasm", ], aes(x=factor(yr5_status, levels = c("DF", "Death", "Induction Failure", "Relapse")), y=pre_resist)) + 
  geom_violin(aes(fill = yr5_status), trim=TRUE, scale = "width") +
  geom_quasirandom(size = 1) +
  geom_boxplot(width=.1) + theme_cowplot() +
  stat_pvalue_manual(
    pairwise.test, label = "p.adj.signif", 
    y.position = c(1.3, 1.4, 1.5)
  )

pairwise.test <- dx_bm_5yr[dx_bm_5yr$yr5_status %nin% "Second Malignant Neoplasm", ] %>% wilcox_test(pre_resist ~ yr5_status, p.adjust.method = "bonferroni")
pairwise.test <- pairwise.test[pairwise.test$p.adj.signif != "ns", ]
ggplot(dx_bm_5yr[dx_bm_5yr$yr5_status %nin% "Second Malignant Neoplasm", ], aes(x=factor(yr5_status, levels = c("DF", "Death", "Induction Failure", "Relapse")), y=pre_resist)) + 
  geom_violin(aes(fill = yr5_status), trim=TRUE, scale = "width") +
  geom_quasirandom(size = 0.5) +
  geom_boxplot(width=.1) + theme_cowplot() +
  stat_pvalue_manual(
    pairwise.test, label = "p.adj.signif", 
    y.position = c(2, 2.1)
  )


##### calculate correlation between gene set expression and clinical factors #####
corr_df_bm <- dx_bm_subset_w_pcts[, c(2,37,47,86,87)] #get columns of interest
corr_df_blood <- dx_blood_subset_w_pcts[, c(2,37,47,86,87)]

corr_num_factors <- cor(corr_df_blood, method = "kendall")
corrplot(
  corr_num_factors,
  method = "color")
