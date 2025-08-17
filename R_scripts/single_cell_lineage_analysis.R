# adding FLARE lineage information to single cell objects, downstream analysis

##### packages #####
library(ggplot2)
library(devtools)
library(ShortRead) 
library(data.table)
library(Seurat)
library(dplyr)
library(stringr)
library(BiocManager)
library(radiant.data)
library(ggnewscale)
library(grid)
library(scCustomize)
library(cowplot)
library(rstatix)
library(ggpubr)
`%nin%` = Negate(`%in%`)


##### adding FLARE lineage information #####
lin_assignments <- read.csv("/dartfs/rc/lab/M/McKennaLab/projects/hannah/aml/analysis/c1498_lineage_NEW/rs22_rs28/final_tables_res1/rs22_rs28_combined_cells_founder_and_casclone_annotation.csv", row.names = 1)
ftags_min <- as.data.frame(lin_assignments[,3])
row.names(ftags_min) <- lin_assignments$GEX_bc
sep_demux_combo_joined <- AddMetaData(
  object = sep_demux_combo_joined,
  metadata = ftags_min,
  col.name = "founder"
)

castags_min <- as.data.frame(lin_assignments[,6])
row.names(castags_min) <- lin_assignments$GEX_bc
sep_demux_combo_joined <- AddMetaData(
  object = sep_demux_combo_joined,
  metadata = castags_min,
  col.name = "cas_clone"
)

# remove cells without a founder assignment
strict_founder_cells <- subset(sep_demux_combo_joined, founder != "NA")


##### DEGs by founder, control cells only #####
controls <- strict_founder_cells[, strict_founder_cells$group == "Control"]
controls <- SetIdent(controls, value = controls$founder)
f1_degs <- FindMarkers(controls, ident.1 = 1)
f6_degs <- FindMarkers(controls, ident.1 = 6)


##### bar plot - founders by exp group  (figure 2c) #####
founder_splits <- as.data.frame(strict_founder_cells$group)
names(founder_splits) <- "group"
founder_splits$founder <- strict_founder_cells$founder
prop_df <- founder_splits %>% count(group, founder)
prop_df <- prop_df[!(is.na(prop_df$founder)), ]
extra_row <- c("Resistant", 11, 0)
extra_row1 <- c("Persistent", 11, 0)
prop_df <- rbind(prop_df, extra_row, extra_row1)
prop_df$founder <- as.factor(prop_df$founder)
prop_df$n <- as.numeric(prop_df$n)

prop_df %>%
  group_by(group) %>%
  mutate(per =  100 *n/sum(n)) %>% 
  ungroup

ggplot(prop_df, aes(fill = factor(founder, levels=c(11,6,1)), y = n, x = factor(group, levels=c("Control", "Persistent", "Resistant")))) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + ylab("cell count") + xlab(NULL) +
  scale_fill_manual(values=c("#2a731f", "#5a9835", "#cfeac7")) + theme_cowplot() + theme(legend.position = "none")



##### dot plot of clone proportions by exp group (figure 2d) #####
library(ggprism)
# isolate cells with clone assignemnts
strict_founder_cells_cas <- strict_founder_cells[, strict_founder_cells$cas_clone != "F1.unassigned"]
strict_founder_cells_cas <- strict_founder_cells_cas[, strict_founder_cells_cas$cas_clone != "F6.unassigned"]
strict_founder_cells_cas <- strict_founder_cells_cas[, strict_founder_cells_cas$cas_clone != "F11.unassigned"]

clone_df <- as.data.frame(strict_founder_cells_cas$cas_clone)
names(clone_df) <- "clone"
clone_df$group <- strict_founder_cells_cas$group
clone_df$founder <- strict_founder_cells_cas$founder

resist_clones <- as.data.frame(table(strict_founder_cells_cas$cas_clone[strict_founder_cells_cas$group == "Resistant"]))
resist_clones$prop <- resist_clones$Freq/sum(resist_clones$Freq)
resist_clones$group <- "Resistant"

persist_clones <- as.data.frame(table(strict_founder_cells_cas$cas_clone[strict_founder_cells_cas$group == "Persistent"]))
persist_clones$prop <- persist_clones$Freq/sum(persist_clones$Freq)
persist_clones$group <- "Persistent"

ct_clones <- as.data.frame(table(strict_founder_cells_cas$cas_clone[strict_founder_cells_cas$group == "Control"]))
ct_clones$prop <- ct_clones$Freq/sum(ct_clones$Freq)
ct_clones$group <- "Control"

clone_summ <- rbind(resist_clones, persist_clones, ct_clones)
founder_ids <- clone_df[, c(1,3)] %>% distinct()
clone_summ <- clone_summ %>% inner_join(founder_ids, by = c("Var1" = "clone"))

f1_split <- clone_summ[clone_summ$founder == 1, ]
f6_split <- clone_summ[clone_summ$founder == 6, ]
f11_split <- clone_summ[clone_summ$founder == 11, ]

clone_summ_byf <- rbind(f1_split, f6_split, f11_split)
clone_summ_byf$Var1 <- as.character(clone_summ_byf$Var1)
clone_summ_byf$founder <- as.factor(clone_summ_byf$founder)
clone_summ_byf$Var1 <- factor(clone_summ_byf$Var1, levels = unique(clone_summ_byf$Var1))
clone_summ_byf %>% 
  ggplot(aes(x=group, y = factor(Var1, levels=c("F1.1", "F1.2", "F1.3", "F1.4", "F1.5", "F1.6", "F1.7", "F1.8", "F1.9", "F1.10", "F1.11", "F6.1", "F6.2", "F6.3", "F6.4", "F6.5", "F6.6", "F6.7", "F6.8", "F6.9", "F6.10", "F6.11", "F6.12", "F11.1")), color = founder, size = prop + 1)) + 
  geom_point() + coord_flip() + theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10)) + 
  ylab('Clone') + xlab(NULL) + 
  scale_color_manual(values = c("#9cd379", "#5a9835", "#2a731f"))


##### violin plots of Hotspot module expression #####
vln_df = data.frame(module = strict_founder_cells$hotspot_f1_mod1_1, cluster = strict_founder_cells$group)
pairwise.test <- vln_df %>% wilcox_test(module ~ cluster, p.adjust.method = "bonferroni")
pairwise.test <- pairwise.test[pairwise.test$p.adj.signif != "ns", ]
ggplot(vln_df, aes(x=cluster, y=module)) + 
  geom_violin(aes(fill = cluster), trim=TRUE, scale = "width") +
  geom_boxplot(width=.1) + theme_cowplot() +
  stat_pvalue_manual(
    pairwise.test, label = "p.adj.signif", 
    y.position = c(1.2, 1.3)
  )

##### fGSEA #####
library(msigdbr)
library(presto)
library(fgsea)
library(tidyverse)

controls <- strict_founder_cells[, strict_founder_cells$group == "Control"]
controls <- NormalizeData(controls)
DefaultAssay(controls) <- "RNA"
controls <- SetIdent(controls, value = controls$founder)
f_genes <- wilcoxauc(controls, "founder", layer = 'data')

f1_genes<- f_genes %>%
  dplyr::filter(group == 1) %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
f6_genes<- f_genes %>%
  dplyr::filter(group == 6) %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

ranks_f1 <- deframe(f1_genes)
ranks_f6 <- deframe(f6_genes)

m_df<- msigdbr(species = "Mus musculus", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

fgseaRes_f1 <- fgsea(fgsea_sets, stats = ranks_f1)
fgseaRes_f6 <- fgsea(fgsea_sets, stats = ranks_f6)

fgseaResTidy_f1 <- fgseaRes_f1 %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy_f6 <- fgseaRes_f6 %>%
  as_tibble() %>%
  arrange(desc(NES))

# only plot the top 20 pathways
ggplot(fgseaResTidy_f1 %>% filter(padj < 0.001) %>% head(n= 30), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES - f1") + 
  theme_minimal()
ggplot(fgseaResTidy_f6 %>% filter(padj < 0.001) %>% head(n= 30), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES - f6") + 
  theme_minimal()



##### fish/alluvial plots for HL60s founder proportion over time (figure 3c) #####
start_table <- as.data.frame(table(strict_founder_hl60s$founder[strict_founder_hl60s$sample == "Pre_exp_pool"]))
names(start_table) <- c("Founder", 0)
amid_table <- as.data.frame(table(strict_founder_hl60s$founder[strict_founder_hl60s$sample == "A_midpoint"]))
names(amid_table) <- c("Founder", 36)
bmid_sub_table <- as.data.frame(table(strict_founder_hl60s$founder[strict_founder_hl60s$sample == "B_midpoint"]))
names(bmid_sub_table) <- c("Founder", 36)
cmid_sub_table <- as.data.frame(table(strict_founder_hl60s$founder[strict_founder_hl60s$sample == "Control_midpoint"]))
names(cmid_sub_table) <- c("Founder", 43)
aep_sub_table <- as.data.frame(table(strict_founder_hl60s$founder[strict_founder_hl60s$sample == "A"]))
names(aep_sub_table) <- c("Founder", 89)
bep_sub_table <- as.data.frame(table(strict_founder_hl60s$founder[strict_founder_hl60s$sample == "B"]))
names(bep_sub_table) <- c("Founder", 89)
cep_sub_table <- as.data.frame(table(strict_founder_hl60s$founder[strict_founder_hl60s$sample == "Control"]))
names(cep_sub_table) <- c("Founder", 89)

a_counts <- start_table %>% full_join(amid_table, by = "Founder")
a_counts <- a_counts %>% full_join(aep_sub_table, by = "Founder")
a_counts[is.na(a_counts)] <- 0
a_counts[, -1] <- lapply(a_counts[, -1], function(x) x/sum(x, na.rm=TRUE) )
a_counts_prop_df <- a_counts %>%
  pivot_longer(!Founder, names_to = "timepoint", values_to = "prop")
a_counts_prop_df$timepoint <- as.numeric(a_counts_prop_df$timepoint)
ggplot(a_counts_prop_df, aes(y = prop, x = timepoint, fill = factor(Founder, levels = c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "H13", "H14", "H15")))) +
  geom_alluvium(aes(alluvium = factor(Founder, levels = c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "H13", "H14", "H15"))), alpha= .9, color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() + theme_cowplot() + ylab("Proportion of reads") + guides(fill = guide_legend(title = "Founder")) + scale_fill_manual(values = rev(as.vector(stepped(15)))) 

b_counts <- start_table %>% full_join(bmid_sub_table, by = "Founder")
b_counts <- b_counts %>% full_join(bep_sub_table, by = "Founder")
b_counts[is.na(b_counts)] <- 0
b_counts[, -1] <- lapply(b_counts[, -1], function(x) x/sum(x, na.rm=TRUE) )
b_counts_prop_df <- b_counts %>%
  pivot_longer(!Founder, names_to = "timepoint", values_to = "prop")
b_counts_prop_df$timepoint <- as.numeric(b_counts_prop_df$timepoint)
ggplot(b_counts_prop_df, aes(y = prop, x = timepoint, fill = factor(Founder, levels = c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "H13", "H14", "H15")))) +
  geom_alluvium(aes(alluvium = factor(Founder, levels = c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "H13", "H14", "H15"))), alpha= .9, color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() + theme_cowplot() + ylab("Proportion of reads") + guides(fill = guide_legend(title = "Founder")) + scale_fill_manual(values = rev(as.vector(stepped(15)))) 

c_counts <- start_table %>% full_join(cmid_sub_table, by = "Founder")
c_counts <- c_counts %>% full_join(cep_sub_table, by = "Founder")
c_counts[is.na(c_counts)] <- 0
c_counts[, -1] <- lapply(c_counts[, -1], function(x) x/sum(x, na.rm=TRUE) )
c_counts_prop_df <- c_counts %>%
  pivot_longer(!Founder, names_to = "timepoint", values_to = "prop")
c_counts_prop_df$timepoint <- as.numeric(c_counts_prop_df$timepoint)
ggplot(c_counts_prop_df, aes(y = prop, x = timepoint, fill = factor(Founder, levels = c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "H13", "H14", "H15")))) +
  geom_alluvium(aes(alluvium = factor(Founder, levels = c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "H13", "H14", "H15"))), alpha= .9, color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() + theme_cowplot() + ylab("Proportion of reads") + guides(fill = guide_legend(title = "Founder")) + scale_fill_manual(values = rev(as.vector(stepped(15)))) 


