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



##### PATH analysis #####
library(ape)
library(PATH)
library(expm)
library(tidytree)
library(Matrix)
library(patchwork)
library(parallel)
library(msigdbr)
library(magrittr)
library(tidyverse)
library(fgsea)
library(ComplexHeatmap)

rosetta <- read.table("/dartfs/rc/lab/M/McKennaLab/resources/cellranger_versions/cellranger-7.2.0/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz") #use updated whitelist for GEM-X 10x runs


DefaultAssay(rs22_rs28_combo_joined) <- "RNA"
var_features <- (Loadings(object = rs22_rs28_combo_joined[["pca"]])[, 1:20])
rs22_new_md <- FetchData(rs22_rs28_combo_joined, vars=c("integrated_snn_res.0.2", "group", rownames(var_features)[1:3000]), layer='data') # fetch data 
rs22_new_md$bc_stripped <- gsub('_\\d', '', rownames(rs22_new_md)) # remove any numbers etc. from the cellID
rs22_new_md$rows <- rownames(rs22_new_md)
rs22_new_md <- merge(rs22_new_md, rosetta, by.x='bc_stripped', by.y = 'V1') # merge with whitelist 
rs22_new_md$exp <- rs22_rs28_combo_joined@meta.data[rs22_new_md$rows, ]$exp
rs22_new_md$lin_cellID <- paste0(rs22_new_md$exp, '_' , rs22_new_md$V2) # add the experimental prefix (i.e. rs22_barcodestring) if applicable
rownames(rs22_new_md) <- rs22_new_md$lin_cellID # assign barcodes that are compatible with tree labels to row names 


#founder 1 tree
mega_tree <- ape::read.tree("/dartfs/rc/lab/M/McKennaLab/projects/hannah/aml/analysis/c1498_lineage_NEW/rs22_rs28/trees_res1/F1_newick/rs22_rs28_F1_Vanilla_rmWT_megatree.newick")
#founder 6 tree
mega_tree_6 <- ape::read.tree("/dartfs/rc/lab/M/McKennaLab/projects/hannah/aml/analysis/c1498_lineage_NEW/rs22_rs28/trees_res1/F6_newick/rs22_rs28_F6_Vanilla_rmWT_megatree.newick")


# function to get autocorrelations 
path_autocorr <- function(tree, metadata){ #prepping tree and running path_inf
  X <-  apply(metadata[tree$tip.label, -1], 2, as.numeric)
  X <- as(X, 'sparseMatrix')
  X <- cbind(PATH::catMat(metadata[tree$tip.label, ]$group), X)
  Winv <- inv_tree_dist(tree, node = TRUE, norm = FALSE)
  modxcor <- xcor(X, Winv)
  Idf <- reshape2::melt(modxcor$phy_cor, 
                      value.name = "I")
  Zdf <- reshape2::melt(modxcor$Z.score, 
                      value.name = "Z")
  Pdf <- reshape2::melt(modxcor$one.sided.pvalue, value.name='p.val')
  Zp_df <- full_join(Zdf, Pdf, by=c('Var1', 'Var2'))
  df <- full_join(Idf, Zp_df, by=c("Var1", "Var2"))
  df <- df %>% mutate(Var1 = as.factor(Var1), 
                      Var2 = as.factor(Var2))
  df$adj_p.val <-p.adjust(df$p.val, method="BH")
  auto_df <- df %>%
    filter(Var1 == Var2) %>%
    filter(adj_p.val < 0.05)
  auto_df <- auto_df[, -1]

  return(auto_df)
}

# does the same thing as the function, for when you'd rather do it step by step or want the cross correlations
X1 <-  apply(rs22_new_md[mega_tree$tip.label, 4:3003], 2, as.numeric)
X1 <- as(X1, 'sparseMatrix')
X1 <- cbind(PATH::catMat(rs22_new_md[mega_tree$tip.label, ]$group), X1) 
Winv1 <- inv_tree_dist(mega_tree, node = TRUE, norm = FALSE)
modxcor1 <- xcor(X1, Winv1)
Idf1 <- reshape2::melt(modxcor1$phy_cor, 
                       value.name = "I")
Zdf1 <- reshape2::melt(modxcor1$Z.score, 
                       value.name = "Z")
Pdf1 <- reshape2::melt(modxcor1$one.sided.pvalue, value.name='p.val')
Zp_df1 <- full_join(Zdf1, Pdf1, by=c('Var1', 'Var2'))
df1 <- full_join(Idf1, Zp_df1, by=c("Var1", "Var2"))
df1 <- df1 %>% mutate(Var1 = as.factor(Var1), 
                      Var2 = as.factor(Var2))
df1$adj_p.val <-p.adjust(df1$p.val, method="BH")

auto_df1 <- df1 %>%
  filter(Var1 == Var2) %>%
  filter(adj_p.val < 0.05)
auto_df1 <- auto_df1[, -1]

#get the top z scores
auto_df1 <- auto_df1[-grep("Rpl", auto_df1$Var2), ]
auto_df1 <- auto_df1[-grep("Rps", auto_df1$Var2), ]
auto_df1 <- auto_df1[-grep("mt-", auto_df1$Var2), ]
auto_df1 <- auto_df1 %>% arrange(desc(Z))
f1_genes <- auto_df1$Var2[1:15] #for the main figure


X6 <-  apply(rs22_new_md[mega_tree_6$tip.label, 4:3003], 2, as.numeric)
X6 <- as(X6, 'sparseMatrix') 
X6 <- cbind(PATH::catMat(rs22_new_md[mega_tree_6$tip.label, ]$group), X6) 
Winv6 <- inv_tree_dist(mega_tree_6, node = TRUE, norm = FALSE)
modxcor6 <- xcor(X6, Winv6)
Idf6 <- reshape2::melt(modxcor6$phy_cor, 
                       value.name = "I")
Zdf6 <- reshape2::melt(modxcor6$Z.score, 
                       value.name = "Z")
Pdf6 <- reshape2::melt(modxcor6$one.sided.pvalue, value.name='p.val')
Zp_df6 <- full_join(Zdf6, Pdf6, by=c('Var1', 'Var2'))
df6 <- full_join(Idf6, Zp_df6, by=c("Var1", "Var2"))
df6 <- df6 %>% mutate(Var1 = as.factor(Var1), 
                      Var2 = as.factor(Var2))
df6$adj_p.val <-p.adjust(df6$p.val, method="BH")

auto_df6 <- df6 %>%
  filter(Var1 == Var2) %>%
  filter(adj_p.val < 0.05)
auto_df6 <- auto_df6[, -1]

# get the top z scores
auto_df6 <- auto_df6[-grep("Rpl", auto_df6$Var2), ]
auto_df6 <- auto_df6[-grep("Rps", auto_df6$Var2), ]
auto_df6 <- auto_df6[-grep("mt-", auto_df6$Var2), ]
auto_df6 <- auto_df6 %>% arrange(desc(Z))
f6_genes <- auto_df6$Var2[1:15]

# make joined gene list, get rid of duplicates
gene_list <- c(f1_genes, f6_genes)
# set order if applicable
gene_list_order <- c("Control", "Persistent", "Resistant", "Mouse", "Adprh", "Mctp1", "Itsn1", "Cd7", "Slc9a9", "Vim", "Lpcat2", "Gng5", "Grap2", "Olfm3", "Sugct", "Rapsn", "Hba-a1", "Fcer1g", "Npm1", "Kcnk13", "Cpa3", "Ybx1", "Eef1a1", "Eif4a1", "Gm10076", "Hsp90aa1", "Scd2", "Malat1")

f1_gene_df <- df1[df1$Var1 %in% gene_list_order, ]
f1_gene_df <- f1_gene_df %>% filter(Var1 == Var2)
f1_gene_df$founder <- "F1"

f6_gene_df <- df6[df6$Var1 %in% gene_list_order, ]
f6_gene_df <- f6_gene_df %>% filter(Var1 == Var2)
f6_gene_df$founder <- "F6"

big_df <- rbind(f1_gene_df, f6_gene_df)
extra_row <- c("Mouse", "Mouse", 0.00, 0.00, 0.00, 0.00, "F1")
big_df <- rbind(big_df, extra_row)
big_df$I <- as.numeric(big_df$I)
big_df$Z <- as.numeric(big_df$Z)

# Phylogenetic auto-correlation bar plot (figure 2h)
maxz <- max(abs(big_df$Z))
big_df %>%
  ggplot(aes(x = factor(Var1, levels = gene_list_order), y = Z, fill = founder)) +
  geom_bar(stat="identity", position = "dodge", width = 0.5) +  
  labs(fill="Gene") +
  ylab("Phylogenetic auto-correlation\n(z score)") + 
  xlab("Gene") +
  geom_hline(yintercept = qnorm(0.05, lower.tail = F), 
             col="black", lty=2, linewidth=1) +
  ggtitle("Variable feature heritability", 
          subtitle = "z score") + 
  scale_fill_manual(values = c("#bee0b4", "#72a645")) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 12), legend.position = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, maxz))


# cross correlation heatmap
res_mat <- modxcor$Z.score
res_mat[is.na(res_mat)] <- 0 
Heatmap(res_mat, 
        #show_row_name=FALSE,
        #right_annotation = ha,
        row_names_side = "left", 
        #row_names_gp = gpar(fontsize = 4), 
        name='Phylogenetic Correlation Z-score')


##### State transitions + permutation analysis #####
new_md <- FetchData(f1, vars=c("clust_region", "clust_region_fine"), layer='data') # fetch data 
new_md$bc_stripped <- gsub('_\\d', '', rownames(new_md)) #remove any numbers etc. from the cellID
new_md$rows <- rownames(new_md)
new_md <- merge(new_md, rosetta, by.x='bc_stripped', by.y = 'V1') # merge with whitelist 
new_md$exp <- "rs22"
new_md$lin_cellID <- paste0(new_md$exp, '_' , new_md$V2) # this specific dataset has a prefix, so add the experimental prefix (i.e. rs22_barcodestring) 
rownames(new_md) <- new_md$lin_cellID

# prune the tree to only include leaves of interest (control cells only for fig. S5A)
f1.1_todrop <- f1.1_tree$tip.label[f1.1_tree$tip.label %nin% rownames(new_md)]
ctrl_f1.1_tree <- ape::drop.tip(f1.1_tree, f1.1_todrop)

transition_matrix <- function(tree, metadata){ #prepping tree and running path_inf
  tree_md <- tree
  clust_region <- metadata[tree_md$tip.label, ]$clust_region
  tree_md$clust_region <- clust_region
  tree_md$edge.length <- rep(1, length(tree_md$edge))
  Pinf <- PATH_inf(tree = tree_md, cell_states = "clust_region", impute_branches = TRUE, sample_rate_est = 10^-6)
  return(Pinf$Pt)
}

perm_tree_test <- function(tree, metadata, nperm){
  true_transitions <- transition_matrix(tree, metadata) #get actual transition matrix
  list_set <- list(true_transitions)
  
  for(i in 1:nperm){
    tree_random <- tree
    tree_random$tip.label <- sample(tree_random$tip.label)
    list_set[[i+1]] <- transition_matrix(tree_random, metadata)
  }
  return(list_set)
}

perm_list <- perm_tree_test(f1.2_tree, new_md, 1000)

# prepping results for plotting
perm_df <- data.frame(matrix(ncol=16,nrow=1, dimnames=list(NA, c("adprh2adprh", "adprh2ctrl", "adprh2persist", "adprh2resist", "ctrl2adprh", "ctrl2ctrl", "ctrl2persist", "ctrl2resist", "persist2adprh", "persist2ctrl", "persist2persist", "persist2resist", "resist2adprh", "resist2ctrl", "resist2persist", "resist2resist"))))
for(perm in 1:length(perm_list)){
  temp_row <- c()
  mt <- perm_list[[perm]]
  for(i in 1:nrow(mt)){
    for(j in 1:ncol(mt)){
      val <- mt[i,j]
      temp_row <- c(temp_row, val)
    }
  }
  if (sum(is.na(temp_row)) == 0) {
    perm_df <- rbind(perm_df, temp_row)
  }
}
perm_df <- perm_df[-1,] %>% rownames_to_column() # get rid of extra row with NAs
perm_df <- perm_df[1:7501,] # isolate just the first 7500 permutations

perm_df_longer <- perm_df %>% pivot_longer(!rowname, names_to = "transition", values_to = "prob")

# for plotting
perm_df_boxplot <- perm_df_longer[perm_df_longer$rowname != 2, ] # isolate permutations
perm_df_true <- perm_df_longer[perm_df_longer$rowname == 2, ] #isolate true transition probabilities

# get p vals
perm_df_test <- perm_df[-1,] # remove true transitions

pvals <- c()
for (i in 2:17){
  test_mean <- mean(perm_df_test[,i]) # find permutation mean
  test_diff <- abs(test_mean - perm_df[1, i]) # calculate difference between permutation mean and true transition prob
  pval <- (sum(perm_df_test[,i] >= (test_mean + test_diff) | perm_df_test[,i] <= (test_mean - test_diff)))/length(perm_df_test[,i]) # calculate p vals
  pvals <- c(pvals, pval)
}

pval_df <- data.frame(matrix(unique(perm_df_boxplot$transition), nrow=length(unique(perm_df_boxplot$transition)), byrow=TRUE))
pval_df$pval <- pvals
pval_df$pvals_adj <- p.adjust(pval_df$pval, method = "fdr", n = length(pval_df$pval))
names(pval_df) <- c("transition", "pval", "pval_adj")
# need to add group1 and group2 for stat_pvalue_manual to be happy
pval_df$group1 <- 1
pval_df$group2 <- 1

pval_df <- add_significance(pval_df, "pval_adj", "sig_level")
pval_df <- pval_df %>%
    filter(sig_level != "ns")

#plot the spread of tree permutations, true transitions in red
ggplot(perm_df_boxplot, aes(x=transition, y=prob)) + 
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) + 
  ggplot2::geom_point(perm_df_true, mapping = aes(x=transition, y=prob), color = "red", size = 4) +
  theme_cowplot() +
  stat_pvalue_manual(
  pval_df, x = "transition", y.position = 1.05,
  label = "sig_level", size = 8) 


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


