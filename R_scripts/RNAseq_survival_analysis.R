##### packages #####
library(ggplot2)
library(survival)
library(corrplot)
library(ShortRead) 
library(data.table)
library(dplyr)
library(stringr)
library(bayesbio)
library(grid)
library(cowplot)
library(org.Hs.eg.db)
library(biomaRt)
library(rstatix)
library(ggbeeswarm)
`%nin%` = Negate(`%in%`)



##### load modules, get human orthologs for C1498 modules #####
hl60_mods <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS30_RS31_humancellEP/GEX/hotspot/hl60_megamega_modlist.csv")

# split modules into individual lists, clean up (saved in a weird csv format)
mega_mod1 <- hl60_mods$Module.1[hl60_mods$Module.1 != ""]
mega_mod2 <- hl60_mods$Module.2[hl60_mods$Module.2 != ""]
mega_mod3 <- hl60_mods$Module.3[hl60_mods$Module.3 != ""]
mega_mod4 <- hl60_mods$Module.4[hl60_mods$Module.4 != ""]
mega_mod5 <- hl60_mods$Module.5[hl60_mods$Module.5 != ""]
# repeat for all modules of interest


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
f6m1 <- make_hs_list(as.data.frame(f6_mega_mods$Gene[f6_mega_mods$Module == 1]))
f6m2 <- make_hs_list(as.data.frame(f6_mega_mods$Gene[f6_mega_mods$Module == 2]))
# repeat for all modules of interest

resist_lin <- c("CLSTN2", "SPATA6", "SGCZ", "PRICKLE1", "IGF2BP2", "CAMK2D", "MEIS1", "KCND2", "MKX", "BMP2", "TRHDE", "APBA2", "PPFIA2", "CPNE8", "PON2", "LAMP5", "LNCAROD", "CPLANE1", "HOXB-AS3", "NAV3", "SDK1", "IRF8", "CES1", "PTPRN2", "H1-3", "NCAM1", "LHX2", "SASH1", "AK1", "ELANE", "PPDPF", "UNCX", "CFD", "AC107223.1", "EXT1", "BEX3", "CLEC11A", "LINC02169", "TNS3", "TTC28", "KLHL29", "CKAP4", "AL163541.1", "PRDX2", "BEX1", "MAP1LC3A", "HOXA9", "PSD3", "AC090796.1", "SERPINB6", "WDR49", "ANXA2", "CD36", "MEF2C", "SLC35F1", "HNMT", "AL713998.1", "MRAS", "DAPK1", "DTNA", "ACAA2", "HLA-B", "MGMT", "ASGR2", "STARD13")
resist_lin_subset <- resist_lin[resist_lin %nin% c("MEIS1", "IRF8", "MEF2C", "HOXA9")]


##### loading and cleaning TARGET-AML data #####
#same process for BEAT AML and TCGA LAML datasets
fn <- "https://cbioportal-datahub.s3.amazonaws.com/aml_target_gdc.tar.gz"
download.file(fn,destfile="tmp.tar.gz")
untar("tmp.tar.gz",list=TRUE)  ## check contents

untar("tmp.tar.gz",files="aml_target_gdc/data_mrna_seq_fpkm_zscores_ref_all_samples.txt")
target_mrna_fpkm_zscores <- read.csv("aml_target_gdc/data_mrna_seq_fpkm_zscores_ref_all_samples.txt", sep = "\t")
genes <- read.csv("/dartfs/rc/lab/M/McKennaLab/projects/saxe/data/other/c1498_hl60_module_analysis/ensembl_biomart_canonical_transcripts_per_hgnc.txt", sep = "\t")
new_mapping <- genes %>% dplyr::select(hgnc_symbol, entrez_gene_id)
target_zscores_matched <- new_mapping %>% right_join(target_mrna_fpkm_zscores, by = c("entrez_gene_id" = "Entrez_Gene_Id"))

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


# adding extra clinical info from GDC, cleaning it up 
clin_val_info <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/other/c1498_hl60_module_analysis/TARGET_AML_clinvalidation_clindata.csv")
clin_val_info_min <- clin_val_info[, c(1,2,3,5,6,7,8,9,13,14,15,18:36,38,39,40,41,42,44,45,46,49,50,51,52,53,54,55,57,58,59,62,64)]
discovery_info <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/other/c1498_hl60_module_analysis/TARGET_AML_discovery_clindata.csv")
discovery_info_min <- discovery_info[, c(1,2,3,5,6,7,8,9,13,14,15,18:36,38,39,40,41,42,44,45,46,49,50,51,52,53,54,55,57,58,59,62,64)]
other_samples <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/other/c1498_hl60_module_analysis/target_AML_clindata_AML1031_20230720.csv")
other_samples_min <- other_samples[, c(1,2,3,5,6,7,8,9,13,14,15,18:36,38,39,40,41,42,44,45,46,49,50,51,52,53,54,55,57,58,59,62,64)]


#merge them, get rid of duplicate rows, cleaning data
clin_val_info_min$MRD...at.end.of.course.1[clin_val_info_min$MRD...at.end.of.course.1 == ""] <- NA
discovery_info_min$MRD...at.end.of.course.1[discovery_info_min$MRD...at.end.of.course.1 == ""] <- NA
other_samples_min$MRD...at.end.of.course.1[other_samples_min$MRD...at.end.of.course.1 == ""] <- NA
clin_val_info_min$MRD...at.end.of.course.2[clin_val_info_min$MRD...at.end.of.course.2 == ""] <- NA
discovery_info_min$MRD...at.end.of.course.2[discovery_info_min$MRD...at.end.of.course.2 == ""] <- NA
other_samples_min$MRD...at.end.of.course.2[other_samples_min$MRD...at.end.of.course.2 == ""] <- NA
clin_val_info_min$MRD...at.end.of.course.1[clin_val_info_min$MRD...at.end.of.course.1 == "N/A"] <- NA
discovery_info_min$MRD...at.end.of.course.1[discovery_info_min$MRD...at.end.of.course.1 == "N/A"] <- NA
other_samples_min$MRD...at.end.of.course.1[other_samples_min$MRD...at.end.of.course.1 == "N/A"] <- NA
clin_val_info_min$MRD...at.end.of.course.2[clin_val_info_min$MRD...at.end.of.course.2 == "N/A"] <- NA
discovery_info_min$MRD...at.end.of.course.2[discovery_info_min$MRD...at.end.of.course.2 == "N/A"] <- NA
other_samples_min$MRD...at.end.of.course.2[other_samples_min$MRD...at.end.of.course.2 == "N/A"] <- NA
discovery_info_min$Risk.group[discovery_info_min$Risk.group %in% "Unknown"] <- ""
other_samples_min$Risk.group[other_samples_min$Risk.group %in% "Unknown"] <- ""
clin_val_info_min$Risk.group[clin_val_info_min$Risk.group %in% "Unknown"] <- ""
clin_val_info_min$Primary.Cytogenetic.Code[clin_val_info_min$Primary.Cytogenetic.Code == ""] <- NA
discovery_info_min$Primary.Cytogenetic.Code[discovery_info_min$Primary.Cytogenetic.Code == ""] <- NA
other_samples_min$Primary.Cytogenetic.Code[other_samples_min$Primary.Cytogenetic.Code == ""] <- NA
clin_val_info_min$Cytogenetic.Complexity[clin_val_info_min$Cytogenetic.Complexity == "N/A"] <- NA
discovery_info_min$Cytogenetic.Complexity[discovery_info_min$Cytogenetic.Complexity == "N/A"] <- NA
other_samples_min$Cytogenetic.Complexity[other_samples_min$Cytogenetic.Complexity == "N/A"] <- NA
clin_val_info_min$FLT3.ITD.positive.[clin_val_info_min$FLT3.ITD.positive. == ""] <- "Unknown"
discovery_info_min$FLT3.ITD.positive.[discovery_info_min$FLT3.ITD.positive. == ""] <- "Unknown"
other_samples_min$FLT3.ITD.positive.[other_samples_min$FLT3.ITD.positive. == ""] <- "Unknown"
clin_val_info_min$FLT3.ITD.positive.[clin_val_info_min$FLT3.ITD.positive. == "NO"] <- "No"
discovery_info_min$FLT3.ITD.positive.[discovery_info_min$FLT3.ITD.positive. == "NO"] <- "No"
other_samples_min$FLT3.ITD.positive.[other_samples_min$FLT3.ITD.positive. == "NO"] <- "No"
clin_val_info_min$FLT3.ITD.positive.[clin_val_info_min$FLT3.ITD.positive. == "YES"] <- "Yes"
discovery_info_min$FLT3.ITD.positive.[discovery_info_min$FLT3.ITD.positive. == "YES"] <- "Yes"
other_samples_min$FLT3.ITD.positive.[other_samples_min$FLT3.ITD.positive. == "YES"] <- "Yes"
clin_val_info_min$FLT3.ITD.allelic.ratio[clin_val_info_min$FLT3.ITD.allelic.ratio == "N/A"] <- NA
discovery_info_min$FLT3.ITD.allelic.ratio[discovery_info_min$FLT3.ITD.allelic.ratio == "N/A"] <- NA
other_samples_min$FLT3.ITD.allelic.ratio[other_samples_min$FLT3.ITD.allelic.ratio == "N/A"] <- NA
clin_val_info_min$CEBPA.mutation[clin_val_info_min$CEBPA.mutation %in% c("", "Unknown")] <- NA
discovery_info_min$CEBPA.mutation[discovery_info_min$CEBPA.mutation %in% c("", "Unknown")] <- NA
other_samples_min$CEBPA.mutation[other_samples_min$CEBPA.mutation %in% c("", "Unknown")] <- NA
clin_val_info_min$CEBPA.mutation[clin_val_info_min$CEBPA.mutation == "NO"] <- "No"
discovery_info_min$CEBPA.mutation[discovery_info_min$CEBPA.mutation == "NO"] <- "No"
other_samples_min$CEBPA.mutation[other_samples_min$CEBPA.mutation == "NO"] <- "No"
clin_val_info_min$CEBPA.mutation[clin_val_info_min$CEBPA.mutation == "YES"] <- "Yes"
discovery_info_min$CEBPA.mutation[discovery_info_min$CEBPA.mutation == "YES"] <- "Yes"
other_samples_min$CEBPA.mutation[other_samples_min$CEBPA.mutation == "YES"] <- "Yes"
clin_val_info_min$NPM.mutation[clin_val_info_min$NPM.mutation == "NO"] <- "No"
discovery_info_min$NPM.mutation[discovery_info_min$NPM.mutation == "NO"] <- "No"
other_samples_min$NPM.mutation[other_samples_min$NPM.mutation == "NO"] <- "No"
clin_val_info_min$NPM.mutation[clin_val_info_min$NPM.mutation == "YES"] <- "Yes"
discovery_info_min$NPM.mutation[discovery_info_min$NPM.mutation == "YES"] <- "Yes"
other_samples_min$NPM.mutation[other_samples_min$NPM.mutation == "YES"] <- "Yes"
clin_val_info_min$NPM.mutation[clin_val_info_min$NPM.mutation %in% c("", "Unknown")] <- NA
discovery_info_min$NPM.mutation[discovery_info_min$NPM.mutation %in% c("", "Unknown")] <- NA
other_samples_min$NPM.mutation[other_samples_min$NPM.mutation %in% c("", "Unknown")] <- NA
clin_val_info_min$WT1.mutation[clin_val_info_min$WT1.mutation == "NO"] <- "No"
discovery_info_min$WT1.mutation[discovery_info_min$WT1.mutation == "NO"] <- "No"
other_samples_min$WT1.mutation[other_samples_min$WT1.mutation == "NO"] <- "No"
clin_val_info_min$WT1.mutation[clin_val_info_min$WT1.mutation == "YES"] <- "Yes"
discovery_info_min$WT1.mutation[discovery_info_min$WT1.mutation == "YES"] <- "Yes"
other_samples_min$WT1.mutation[other_samples_min$WT1.mutation == "YES"] <- "Yes"
clin_val_info_min$WT1.mutation[clin_val_info_min$WT1.mutation %in% c("", "Unknown")] <- NA
discovery_info_min$WT1.mutation[discovery_info_min$WT1.mutation %in% c("", "Unknown")] <- NA
other_samples_min$WT1.mutation[other_samples_min$WT1.mutation %in% c("", "Unknown")] <- NA
clin_val_info_min$Risk.group[clin_val_info_min$Risk.group %in% c("", 10, 30)] <- NA
discovery_info_min$Risk.group[discovery_info_min$Risk.group %in% c("", 10, 30)] <- NA
other_samples_min$Risk.group[other_samples_min$Risk.group %in% c("", 10, 30)] <- NA
clin_val_info_min$Refractory.Timepoint.sent.for.Induction.Failure.Project[clin_val_info_min$Refractory.Timepoint.sent.for.Induction.Failure.Project == ""] <- NA
discovery_info_min$Refractory.Timepoint.sent.for.Induction.Failure.Project[discovery_info_min$Refractory.Timepoint.sent.for.Induction.Failure.Project == ""] <- NA
other_samples_min$Refractory.Timepoint.sent.for.Induction.Failure.Project[other_samples_min$Refractory.Timepoint.sent.for.Induction.Failure.Project == ""] <- NA
other_samples_min$CNS.Site.of.Relapse.Induction.Failure[other_samples_min$CNS.Site.of.Relapse.Induction.Failure == "Not done"] <- "Not Done"
other_samples_min$Chloroma.Site.of.Relapse.Induction.Failure[other_samples_min$Chloroma.Site.of.Relapse.Induction.Failure == "Not done"] <- "Not Done"
other_samples_min$Bone.Marrow.Site.of.Relapse.Induction.Failure[other_samples_min$Bone.Marrow.Site.of.Relapse.Induction.Failure == "Not done"] <- "Not Done"

new_clin_data <- rbind(clin_val_info_min, discovery_info_min, other_samples_min) %>% distinct() 
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

#found separate metadata, adds annotations for a couple patients
gdc_md_sample <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/other/c1498_hl60_module_analysis/clinical_data_GDCtarget_tcgabiolinks.csv", row.names = 1)
pts_with_sample_info <- unique(gdc_md_sample$sample)
sample_info <- gdc_md_sample[gdc_md_sample$definition %nin% c("Blood Derived Normal", "Bone Marrow Normal"), ]
sample_info_new_min <- sample_info[, c(3,8)] %>% distinct()

sample_type_info <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/other/c1498_hl60_module_analysis/080625_removedsort_TARGET_score_df_best_option_use_this_fixed_gene.csv", row.names = 1) 
sample_type_min <- sample_type_info[, c(1,78)]
sample_type_min$sample <- gsub("\\.", "-", sample_type_min$patient)

new <- sample_info_new_min$sample[sample_info_new_min$sample %nin% sample_type_min$sample] #388 new!
sample_type_min <- sample_type_min[, 2:3]
sample_type_min <- sample_type_min %>% mutate(patient = substring(sample, 1, 16)) %>% distinct()
sample_info_new_min <- sample_info_new_min %>% mutate(patient = substring(sample, 1, 16)) %>% distinct()
names(sample_type_min) <- c("sample_type", "sample", "patient")
names(sample_info_new_min) <- c("sample", "sample_type", "patient")

all_sample_meta <- rbind(sample_type_min, sample_info_new_min) %>% distinct() #2372 total
all_the_metadata_icanfind <- all_sample_meta %>% inner_join(clin_data_full, by = c("sample" = "SAMPLE_ID"))




##### score patients on expression #####
score_df <- as.data.frame(names(target_zscores_matched[, -c(1:2)]))
names(score_df) <- "patient"

#resistant lin
resist_lin_data <- target_zscores_matched[target_zscores_matched$hgnc_symbol %in% resist_lin, ]
rownames(resist_lin_data) <- resist_lin_data$hgnc_symbol
resist_lin_data <- resist_lin_data[, -c(1:2)]
resist_lin_t <- t(resist_lin_data)

score_list <- c()
for (i in 1:nrow(resist_lin_t)) {
  x <- as.data.frame(t(resist_lin_t[i,]))
  x <- as.numeric(x)
  score <- sum(x, na.rm = TRUE)/sum(!is.na(x))
  score_list <- c(score_list, score)
}
resist_lin_min <- as.data.frame(rownames(resist_lin_t))
names(resist_lin_min) <- "patient"
resist_lin_min$resist_lin <- score_list
score_df <- score_df %>% left_join(resist_lin_min, by = "patient")

#resistant lin subset
resist_lin_sub_data <- target_zscores_matched[target_zscores_matched$hgnc_symbol %in% resist_lin_subset, ]
rownames(resist_lin_sub_data) <- resist_lin_sub_data$hgnc_symbol
resist_lin_sub_data <- resist_lin_sub_data[, -c(1:2)]
resist_lin_sub_t <- t(resist_lin_sub_data)

score_list <- c()
for (i in 1:nrow(resist_lin_sub_t)) {
  x <- as.data.frame(t(resist_lin_sub_t[i,]))
  x <- as.numeric(x)
  score <- sum(x, na.rm = TRUE)/sum(!is.na(x))
  score_list <- c(score_list, score)
}
resist_lin_sub_min <- as.data.frame(rownames(resist_lin_sub_t))
names(resist_lin_sub_min) <- "patient"
resist_lin_sub_min$resist_lin_sub <- score_list
score_df <- score_df %>% left_join(resist_lin_sub_min, by = "patient")

#immunoprot 11s
prot_11s_data <- target_zscores_matched[target_zscores_matched$hgnc_symbol %in% c("PSME1", "PSME2"), ]
rownames(prot_11s_data) <- prot_11s_data$hgnc_symbol
prot_11s_data <- prot_11s_data[, -c(1:2)]
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
score_df <- score_df %>% left_join(prot_11s_min, by = "patient")
#repeat for any signature or module of interest

#merging with clin inf0
score_df$Sample_ID <- score_df$patient
score_df$Sample_ID <- gsub('\\.', '-', score_df$Sample_ID)
score_df <- score_df %>% inner_join(all_the_metadata_icanfind, by = c("Sample_ID" = "sample")) #2332 samples

#get one sample per patient at diagnosis. prioritize bone marrow but use blood if no bone marrow
dx_bm <- score_df[score_df$sample_type %in% "Primary Blood Derived Cancer - Bone Marrow", ] #1594, all unique
dx_blood <- score_df[score_df$sample_type %in% "Primary Blood Derived Cancer - Peripheral Blood", ] #312, all unique
all_dx_unique <- rbind(dx_bm, dx_blood[dx_blood$PATIENT_ID %nin% intersect(dx_blood$PATIENT_ID, dx_bm$PATIENT_ID), ]) #1890

all_dx_unique$Risk.group <- factor(all_dx_unique$Risk.group, levels = c("Standard Risk", "Low Risk", "High Risk"))
all_dx_unique$Primary.Cytogenetic.Code[all_dx_unique$Primary.Cytogenetic.Code %in% "Unknown"] <- NA
all_dx_unique$Primary.Cytogenetic.Code <- factor(all_dx_unique$Primary.Cytogenetic.Code, levels = c("Normal", "MLL", "inv(16)", "t(8;21)", "Other"))
all_dx_unique$FLT3.ITD.positive.[all_dx_unique$FLT3.ITD.positive. == "Unknown"] <- NA

all_dx_unique_subset <- all_dx_unique[!is.na(all_dx_unique$Primary.Cytogenetic.Code), ]
all_dx_unique_subset <- all_dx_unique_subset[!is.na(all_dx_unique_subset$FLT3.ITD.positive.), ]
all_dx_unique_subset <- all_dx_unique_subset[!is.na(all_dx_unique_subset$WBC.at.Diagnosis), ]


##### coxph, forest plots #####
model <- coxph(Surv(months_to_follow_up, status) ~ AGE + SEX + Primary.Cytogenetic.Code + FLT3.ITD.positive. + resist_lin, data = all_dx_unique_subset)
ggforest(model, data = all_dx_unique_subset) + theme_cowplot()


##### kaplan meiers #####
sd <- sd(all_dx_unique$resist_lin)
mean <- mean(all_dx_unique$resist_lin)
highExpr <- mean + sd
lowExpr <- mean - sd
all_dx_unique$level <- ifelse(all_dx_unique$resist_lin >= highExpr, 'High',
                      ifelse(all_dx_unique$resist_lin <= lowExpr, 'Low', 'Mid'))
# relevel the factors to have mid as the ref level
all_dx_unique$level <- factor(all_dx_unique$level,
                      levels = c('Mid', 'Low', 'High'))

ggsurvplot(survfit(Surv(months_to_follow_up, status) ~ level,
                   data = all_dx_unique),
           data = all_dx_unique,
           risk.table = TRUE,
           pval = TRUE,
           break.time.by = 12,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE, 
           censor.shape=124, title = "resist_lin") 


##### comparing to random gene sets (figure S9H)
clean <- target_zscores_matched[!is.na(target_zscores_matched$hgnc_symbol), ]
clean <- clean[rowSums(is.na(clean[,c(3:2615)])) != 2613, ]

model <- coxph(Surv(months_to_follow_up, status) ~ AGE + SEX + resist_lin, data = all_dx_unique)
ggforest(model, data = all_dx_unique) + theme_cowplot() #not sig
real_hr <- summary(model)$coefficients[3,2]
real_pval <- summary(model)$coefficients[3,5]
real_conf_low <- summary(model)$conf.int[3,3]
real_conf_high <- summary(model)$conf.int[3,4]

hrs <- c()
pvals <- c()
conf_low <- c()
conf_high <- c()
nperm <- 1000
for (n in 1:nperm) {
  score_df_temp <- as.data.frame(names(clean[, -c(1:2)]))
  names(score_df_temp) <- "patient"
  randos_int <- sample(1:nrow(clean), 61) #run for random iterations
  rando_sig_data <- clean[randos_int, ]
  rownames(rando_sig_data) <- rando_sig_data$hgnc_symbol
  rando_sig_data <- rando_sig_data[, -c(1:2)]
  rando_sig_t <- t(rando_sig_data)
  
  rando_score_list <- c()
  
  for (i in 1:nrow(rando_sig_t)) {
    x <- as.data.frame(t(rando_sig_t[i,]))
    x <- as.numeric(x)
    score <- sum(x, na.rm = TRUE)/sum(!is.na(x))
    rando_score_list <- c(rando_score_list, score)
  }
  
  rando_sig_min <- as.data.frame(rownames(rando_sig_t))
  names(rando_sig_min) <- "patient"
  rando_sig_min$rando_score <- rando_score_list
  score_df_temp <- score_df_temp %>% left_join(rando_sig_min, by = "patient")
  
  score_df_temp$Sample_ID <- score_df_temp$patient
  score_df_temp$Sample_ID <- gsub('\\.', '-', score_df_temp$Sample_ID)
  
  score_df_temp <- score_df_temp %>% inner_join(all_the_metadata_icanfind, by = c("Sample_ID" = "sample"))
  
  dx_bm_temp <- score_df_temp[score_df_temp$sample_type %in% "Primary Blood Derived Cancer - Bone Marrow", ] #all unique
  dx_blood_temp <- score_df_temp[score_df_temp$sample_type %in% "Primary Blood Derived Cancer - Peripheral Blood", ] #all unique
  all_dx_unique_temp <- rbind(dx_bm_temp, dx_blood_temp[dx_blood_temp$PATIENT_ID %nin% intersect(dx_blood_temp$PATIENT_ID, dx_bm_temp$PATIENT_ID), ])
  
  model <- coxph(Surv(months_to_follow_up, status) ~ AGE + SEX + rando_score, data = all_dx_unique_temp)
  
  hrs <- c(hrs, summary(model)$coefficients[3,2])
  pvals <- c(pvals, summary(model)$coefficients[3,5])
  conf_low <- c(conf_low, summary(model)$conf.int[3,3])
  conf_high <- c(conf_high, summary(model)$conf.int[3,4])
}

result_df <- data.frame(hazard = hrs, pval = pvals, low_conf = conf_low, high_conf = conf_high)
hr_pval_df$padj <- p.adjust(hr_pval_df$pval, method="BH") 

ggplot(result_df[-1,], aes(x = hazard)) +
  geom_density() + 
  geom_vline(aes(xintercept = real_hr), 
             color="blue", linetype="solid") + 
  geom_vline(aes(xintercept = real_conf_high), 
  color="red", linetype="dashed") + 
  geom_vline(aes(xintercept = real_conf_low), 
  color="red", linetype="dashed") + xlim(0,4) + theme_cowplot()

perm_pval <- length(result_df[result_df$hazard > real_hr])

##### correlation of signature genes (figure S9G) #####
sig_gene_set <- target_zscores_matched[target_zscores_matched$hgnc_symbol %in% resist_lin, ]
t_df <- as.data.frame(t(sig_gene_set))
names(t_df) <- t_df[1,]
sig_df <- t_df[-c(1:2),] %>% mutate_all(function(x) as.numeric(as.character(x))) %>% rownames_to_column(var = "patient")

sig_df$Sample_ID <- sig_df$patient
sig_df$Sample_ID <- gsub('\\.', '-', sig_df$Sample_ID)
sig_df <- sig_df %>% inner_join(all_the_metadata_ever, by = c("Sample_ID" = "sample"))

dx_bm <- sig_df[sig_df$sample_type %in% "Primary Blood Derived Cancer - Bone Marrow", ] #all unique
dx_blood <- sig_df[sig_df$sample_type %in% "Primary Blood Derived Cancer - Peripheral Blood", ] #all unique
all_dx_unique <- rbind(dx_bm, dx_blood[dx_blood$PATIENT_ID %nin% intersect(dx_blood$PATIENT_ID, dx_bm$PATIENT_ID), ])

corr_num_factors <- cor(all_dx_unique[, 2:62], method = "pearson")
corrplot(
  corr_num_factors,
  method = "color",
  order = "hclust")
