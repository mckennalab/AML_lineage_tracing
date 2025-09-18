# blood draw lineage analysis from bulk sequencing

##### packages #####
library(stringi)
library(cowplot)
library(data.table)
library(stringr)
library(zeallot)
library(R.utils) 
library(rlang)
library(dplyr)
library(ggplot2)
library(janitor)
library(ggalluvial)
library(tidyr)
library(vegan)
library(Biostrings)
`%nin%` = Negate(`%in%`)

##### load and prep data #####
# output from processing with lineage_seq_preprocessing.R
np_28 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS36_blooddraws/np_d28_castags.csv")
np_43 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS36_blooddraws/np_d43_castags.csv")
np_57 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS36_blooddraws/np_d57_castags.csv")
np_67 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS36_blooddraws/np_d67_castags.csv")
rp_28 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS36_blooddraws/rp_d28_castags.csv")
rp_43 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS36_blooddraws/rp_d43_castags.csv")
rp_50 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS36_blooddraws/rp_d50_castags.csv")
lp_28 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS36_blooddraws/lp_d28_castags.csv")
lp_43 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS36_blooddraws/lp_d43_castags.csv")
lp_67 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS36_blooddraws/lp_d67_castags.csv")
rlp_28 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS36_blooddraws/rlp_d28_castags.csv")
rlp_43 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS36_blooddraws/rlp_d43_castags.csv")
rlp_57 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS36_blooddraws/rlp_d57_castags.csv")
rlp_67 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS36_blooddraws/rlp_d67_castags.csv")

castag_masterlist <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS36_blooddraws/rs24_rs19_tagmasterlist_w_founders_120924.csv", row.names = 1) #founder-associated clone tags from bulk sequencing
castag_masterlist_min <- castag_masterlist[, 1:2]

# filter reads
np28 <- as.data.frame(table(np_28$rc_castag)) %>%
  filter(Freq > 2) %>%
  filter(!grepl("-{2,14}", Var1))
np43 <- as.data.frame(table(np_43$rc_castag)) %>%
  filter(Freq > 2) %>%
  filter(!grepl("-{2,14}", Var1))
np57 <- as.data.frame(table(np_57$rc_castag)) %>%
  filter(Freq > 2) %>%
  filter(!grepl("-{2,14}", Var1))
np67 <- as.data.frame(table(np_67$rc_castag)) %>%
  filter(Freq > 2) %>%
  filter(!grepl("-{2,14}", Var1))
rp28 <- as.data.frame(table(rp_28$rc_castag)) %>%
  filter(Freq > 2) %>%
  filter(!grepl("-{2,14}", Var1))
rp43 <- as.data.frame(table(rp_43$rc_castag)) %>%
  filter(Freq > 2) %>%
  filter(!grepl("-{2,14}", Var1))
rp50 <- as.data.frame(table(rp_50$rc_castag)) %>%
  filter(Freq > 2) %>%
  filter(!grepl("-{2,14}", Var1))
lp28 <- as.data.frame(table(lp_28$rc_castag)) %>%
  filter(Freq > 2) %>%
  filter(!grepl("-{2,14}", Var1))
lp43 <- as.data.frame(table(lp_43$rc_castag)) %>%
  filter(Freq > 2) %>%
  filter(!grepl("-{2,14}", Var1))
lp67 <- as.data.frame(table(lp_67$rc_castag)) %>%
  filter(Freq > 2) %>%
  filter(!grepl("-{2,14}", Var1))
rlp28 <- as.data.frame(table(rlp_28$rc_castag)) %>%
  filter(Freq > 2) %>%
  filter(!grepl("-{2,14}", Var1))
rlp43 <- as.data.frame(table(rlp_43$rc_castag)) %>%
  filter(Freq > 2) %>%
  filter(!grepl("-{2,14}", Var1))
rlp57 <- as.data.frame(table(rlp_57$rc_castag)) %>%
  filter(Freq > 2) %>%
  filter(!grepl("-{2,14}", Var1))
rlp67 <- as.data.frame(table(rlp_67$rc_castag)) %>%
  filter(Freq > 2) %>%
  filter(!grepl("-{2,14}", Var1))

# join with castag_masterlist 
np28_sub <- np_28 %>% inner_join(castag_masterlist_min, by = c("rc_castag" = "tag")) %>% filter(rc_castag %in% np28$Var1)
np43_sub <- np_43 %>% inner_join(castag_masterlist_min, by = c("rc_castag" = "tag")) %>% filter(rc_castag %in% np43$Var1)
np57_sub <- np_57 %>% inner_join(castag_masterlist_min, by = c("rc_castag" = "tag")) %>% filter(rc_castag %in% np57$Var1)
np67_sub <- np_67 %>% inner_join(castag_masterlist_min, by = c("rc_castag" = "tag")) %>% filter(rc_castag %in% np67$Var1)
rp28_sub <- rp_28 %>% inner_join(castag_masterlist_min, by = c("rc_castag" = "tag")) %>% filter(rc_castag %in% rp28$Var1)
rp43_sub <- rp_43 %>% inner_join(castag_masterlist_min, by = c("rc_castag" = "tag")) %>% filter(rc_castag %in% rp43$Var1)
rp50_sub <- rp_50 %>% inner_join(castag_masterlist_min, by = c("rc_castag" = "tag")) %>% filter(rc_castag %in% rp50$Var1)
lp28_sub <- lp_28 %>% inner_join(castag_masterlist_min, by = c("rc_castag" = "tag")) %>% filter(rc_castag %in% lp28$Var1)
lp43_sub <- lp_43 %>% inner_join(castag_masterlist_min, by = c("rc_castag" = "tag")) %>% filter(rc_castag %in% lp43$Var1)
lp67_sub <- lp_67 %>% inner_join(castag_masterlist_min, by = c("rc_castag" = "tag")) %>% filter(rc_castag %in% lp67$Var1)
rlp28_sub <- rlp_28 %>% inner_join(castag_masterlist_min, by = c("rc_castag" = "tag")) %>% filter(rc_castag %in% rlp28$Var1)
rlp57_sub <- rlp_57 %>% inner_join(castag_masterlist_min, by = c("rc_castag" = "tag")) %>% filter(rc_castag %in% rlp57$Var1)
rlp67_sub <- rlp_67 %>% inner_join(castag_masterlist_min, by = c("rc_castag" = "tag")) %>% filter(rc_castag %in% rlp67$Var1)


##### f1 v f6 line plot #####
#mouse matching
#np = M1
#lp = M2
#rp = M3
#rlp = M4

# day 13 proportions come from founder tag sequencing
NP_13 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS20_reanalysis/NP_13_founders.csv", row.names = 1)
RP_13 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS20_reanalysis/RP_13_founders.csv", row.names = 1)
LP_13 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS20_reanalysis/LP_13_founders.csv", row.names = 1)
RLP_13 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS20_reanalysis/RLP_13_founders.csv", row.names = 1)

NP13 <- NP_13[, c(1,2,19,20)]
RP13 <- RP_13[, c(1,2,19,20)]
LP13 <- LP_13[, c(1,2,19,20)]
RLP13 <- RLP_13[, c(1,2,19,20)]

clone_ids <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/other/neoclone_IDs.csv", header = TRUE, row.names = 1) #pre-identified founder tags
clone_ids <- clone_ids[clone_ids$clone %in% c(1,5,6,11), ]

# prep data for each mouse/timepoint, add place holders for missing data, not the best approach but it works
np13_sub_table <- as.data.frame(table(NP13$ftag))
np13_sub_table <- np13_sub_table %>% inner_join(clone_ids, by = c("Var1" = "ID"))
np13_counts <- np13_sub_table %>%
  group_by(clone) %>%
  summarize(Mean = mean(Freq))
np13_counts_t <- as.data.frame(t(np13_counts))
np13_counts_t[1,1] <- "F1"
np13_counts_t[1,2] <- "F6"
np13_counts_t[2,2] <- 0
np13_sub_table <- np13_counts_t %>%
  row_to_names(row_number = 1) %>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
np13_sub_table$timepoint <- 13

np28_sub_table <- as.data.frame(table(np28_sub$rc_castag))
np28_sub_table <- np28_sub_table %>% inner_join(castag_masterlist_min, by = c("Var1" = "tag"))
np28_counts <- np28_sub_table %>%
  group_by(founder) %>%
  summarize(Mean = mean(Freq))
np28_counts_t <- as.data.frame(t(np28_counts))
np28_sub_table <- np28_counts_t %>%
  row_to_names(row_number = 1) %>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
np28_sub_table$timepoint <- 28

np43_sub_table <- as.data.frame(table(np43_sub$rc_castag))
np43_sub_table <- np43_sub_table %>% inner_join(castag_masterlist_min, by = c("Var1" = "tag"))
np43_counts <- np43_sub_table %>%
  group_by(founder) %>%
  summarize(Mean = mean(Freq))
np43_counts_t <- as.data.frame(t(np43_counts))
np43_sub_table <- np43_counts_t %>%
  row_to_names(row_number = 1) 
np43_sub_table$F1 <- 0
np43_sub_table <- np43_sub_table%>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
np43_sub_table$timepoint <- 43

np57_sub_table <- as.data.frame(table(np57_sub$rc_castag))
np57_sub_table <- np57_sub_table %>% inner_join(castag_masterlist_min, by = c("Var1" = "tag"))
np57_counts <- np57_sub_table %>%
  group_by(founder) %>%
  summarize(Mean = mean(Freq))
np57_counts_t <- as.data.frame(t(np57_counts))
np57_sub_table <- np57_counts_t %>%
  row_to_names(row_number = 1) 
np57_sub_table <- np57_sub_table%>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
np57_sub_table$timepoint <- 57

np67_sub_table <- as.data.frame(table(np67_sub$rc_castag))
np67_sub_table <- np67_sub_table %>% inner_join(castag_masterlist_min, by = c("Var1" = "tag"))
np67_counts <- np67_sub_table %>%
  group_by(founder) %>%
  summarize(Mean = mean(Freq))
np67_counts_t <- as.data.frame(t(np67_counts))
np67_sub_table <- np67_counts_t %>%
  row_to_names(row_number = 1) 
np67_sub_table <- np67_sub_table%>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
np67_sub_table$timepoint <- 67

# individual mouse line plot
np_lineplot <- rbind(np13_sub_table, np28_sub_table, np43_sub_table, np57_sub_table, np67_sub_table)


rp13_sub_table <- as.data.frame(table(RP13$ftag))
rp13_sub_table <- rp13_sub_table %>% inner_join(clone_ids, by = c("Var1" = "ID"))
rp13_counts <- rp13_sub_table %>%
  group_by(clone) %>%
  summarize(Mean = mean(Freq))
rp13_counts_t <- as.data.frame(t(rp13_counts))
rp13_counts_t[1,1] <- "F1"
rp13_counts_t[1,2] <- "F6"
rp13_counts_t[2,2] <- 0
rp13_sub_table <- rp13_counts_t %>%
  row_to_names(row_number = 1) %>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
rp13_sub_table$timepoint <- 13

rp28_sub_table <- as.data.frame(table(rp28_sub$rc_castag))
rp28_sub_table <- rp28_sub_table %>% inner_join(castag_masterlist_min, by = c("Var1" = "tag"))
rp28_counts <- rp28_sub_table %>%
  group_by(founder) %>%
  summarize(Mean = mean(Freq))
rp28_counts_t <- as.data.frame(t(rp28_counts))
rp28_sub_table <- rp28_counts_t %>%
  row_to_names(row_number = 1) 
rp28_sub_table <- rp28_sub_table%>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
rp28_sub_table$timepoint <- 28

rp43_sub_table <- as.data.frame(table(rp43_sub$rc_castag))
rp43_sub_table <- rp43_sub_table %>% inner_join(castag_masterlist_min, by = c("Var1" = "tag"))
rp43_counts <- rp43_sub_table %>%
  group_by(founder) %>%
  summarize(Mean = mean(Freq))
rp43_counts_t <- as.data.frame(t(rp43_counts))
rp43_sub_table <- rp43_counts_t %>%
  row_to_names(row_number = 1) 
rp43_sub_table <- rp43_sub_table%>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
rp43_sub_table$timepoint <- 43

rp50_sub_table <- as.data.frame(table(rp50_sub$rc_castag))
rp50_sub_table <- rp50_sub_table %>% inner_join(castag_masterlist_min, by = c("Var1" = "tag"))
rp50_counts <- rp50_sub_table %>%
  group_by(founder) %>%
  summarize(Mean = mean(Freq))
rp50_counts_t <- as.data.frame(t(rp50_counts))
rp50_sub_table <- rp50_counts_t %>%
  row_to_names(row_number = 1) 
rp50_sub_table <- rp50_sub_table%>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
rp50_sub_table$timepoint <- 50

# individual mouse line plot
rp_lineplot <- rbind(rp13_sub_table, rp28_sub_table, rp43_sub_table, rp50_sub_table)


lp13_sub_table <- as.data.frame(table(LP13$ftag))
lp13_sub_table <- lp13_sub_table %>% inner_join(clone_ids, by = c("Var1" = "ID"))
lp13_counts <- lp13_sub_table %>%
  group_by(clone) %>%
  summarize(Mean = mean(Freq))
lp13_counts_t <- as.data.frame(t(lp13_counts))
lp13_counts_t[1,1] <- "F1"
lp13_counts_t[1,2] <- "F5"
lp13_counts_t[1,3] <- "F6"
lp13_counts_t[1,4] <- "F11"
lp13_sub_table <- lp13_counts_t %>%
  row_to_names(row_number = 1) %>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
lp13_sub_table$timepoint <- 13
lp13_sub_table <- lp13_sub_table[, -c(2,4)]

lp28_sub_table <- as.data.frame(table(lp28_sub$rc_castag))
lp28_sub_table <- lp28_sub_table %>% inner_join(castag_masterlist_min, by = c("Var1" = "tag"))
lp28_counts <- lp28_sub_table %>%
  group_by(founder) %>%
  summarize(Mean = mean(Freq))
lp28_counts_t <- as.data.frame(t(lp28_counts))
lp28_sub_table <- lp28_counts_t %>%
  row_to_names(row_number = 1) 
lp28_sub_table <- lp28_sub_table%>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
lp28_sub_table$timepoint <- 28

lp67_sub_table <- as.data.frame(table(lp67_sub$rc_castag))
lp67_sub_table <- lp67_sub_table %>% inner_join(castag_masterlist_min, by = c("Var1" = "tag"))
lp67_counts <- lp67_sub_table %>%
  group_by(founder) %>%
  summarize(Mean = mean(Freq))
lp67_counts_t <- as.data.frame(t(lp67_counts))
lp67_sub_table <- lp67_counts_t %>%
  row_to_names(row_number = 1) 
lp67_sub_table$F1 <- 0
lp67_sub_table <- lp67_sub_table%>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
lp67_sub_table$timepoint <- 67

# individual mouse line plot
lp_lineplot <- rbind(lp13_sub_table, lp28_sub_table, lp67_sub_table)


rlp13_sub_table <- as.data.frame(table(RLP13$ftag))
rlp13_sub_table <- rlp13_sub_table %>% inner_join(clone_ids, by = c("Var1" = "ID"))
rlp13_counts <- rlp13_sub_table %>%
  group_by(clone) %>%
  summarize(Mean = mean(Freq))
rlp13_counts_t <- as.data.frame(t(rlp13_counts))
rlp13_counts_t[1,1] <- "F1"
rlp13_counts_t[1,2] <- "F6"
rlp13_counts_t[2,2] <- 0
rlp13_sub_table <- rlp13_counts_t %>%
  row_to_names(row_number = 1) %>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
rlp13_sub_table$timepoint <- 13

rlp28_sub_table <- as.data.frame(table(rlp28_sub$rc_castag))
rlp28_sub_table <- rlp28_sub_table %>% inner_join(castag_masterlist_min, by = c("Var1" = "tag"))
rlp28_counts <- rlp28_sub_table %>%
  group_by(founder) %>%
  summarize(Mean = mean(Freq))
rlp28_counts_t <- as.data.frame(t(rlp28_counts))
rlp28_sub_table <- rlp28_counts_t %>%
  row_to_names(row_number = 1) 
rlp28_sub_table <- rlp28_sub_table%>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
rlp28_sub_table$timepoint <- 28

rlp57_sub_table <- as.data.frame(table(rlp57_sub$rc_castag))
rlp57_sub_table <- rlp57_sub_table %>% inner_join(castag_masterlist_min, by = c("Var1" = "tag"))
rlp57_counts <- rlp57_sub_table %>%
  group_by(founder) %>%
  summarize(Mean = mean(Freq))
rlp57_counts_t <- as.data.frame(t(rlp57_counts))
rlp57_sub_table <- rlp57_counts_t %>%
  row_to_names(row_number = 1) 
rlp57_sub_table <- rlp57_sub_table%>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
rlp57_sub_table$timepoint <- 57

rlp67_sub_table <- as.data.frame(table(rlp67_sub$rc_castag))
rlp67_sub_table <- rlp67_sub_table %>% inner_join(castag_masterlist_min, by = c("Var1" = "tag"))
rlp67_counts <- rlp67_sub_table %>%
  group_by(founder) %>%
  summarize(Mean = mean(Freq))
rlp67_counts_t <- as.data.frame(t(rlp67_counts))
rlp67_sub_table <- rlp67_counts_t %>%
  row_to_names(row_number = 1) 
rlp67_sub_table$F1 <- 0
rlp67_sub_table <- rlp67_sub_table%>%
  mutate(ratio = as.numeric(F1) / (as.numeric(F6)+as.numeric(F1)))
rlp67_sub_table$timepoint <- 67

# individual mouse line plot
rlp_lineplot <- rbind(rlp13_sub_table, rlp28_sub_table, rlp57_sub_table, rlp67_sub_table)

# combine mice for 1 plot
lineplot_full <- rbind(np_lineplot, rp_lineplot, lp_lineplot, rlp_lineplot)
ggplot(lineplot_full, aes(x = timepoint, y = ratio, group = mouse)) +
  geom_line(aes(color = mouse), linewidth = 1) +
  geom_point(aes(shape = mouse, color = mouse), size = 4) + 
  scale_shape_manual(values = 15:18) + 
  ylab("F1 proportion") +
  xlab("Days post engraftment") +
  theme_cowplot()


##### fish/alluvial type plots for clone proportion over time #####
# prep data, again not the best method but it works
np28_sub_table <- as.data.frame(table(np28_sub$rc_castag))
np28_sub_table <- np28_sub_table %>% inner_join(castag_masterlist_min, by = c("Var1" = "tag"))
np28_counts <- np28_sub_table %>%
  group_by(cas_clone.1) %>%
  summarize(Mean = mean(Freq))
names(np28_counts) <- c("Clone", 28)

np43_sub_table <- as.data.frame(table(np43_sub$rc_castag))
np43_sub_table <- np43_sub_table %>% inner_join(cell_castags, by = c("Var1" = "CasTag"))
np43_counts <- np43_sub_table %>%
  group_by(cas_clone.1) %>%
  summarize(Mean = mean(Freq))
names(np43_counts) <- c("Clone", 43)

np57_sub_table <- as.data.frame(table(np57_sub$rc_castag))
np57_sub_table <- np57_sub_table %>% inner_join(cell_castags, by = c("Var1" = "CasTag"))
np57_counts <- np57_sub_table %>%
  group_by(cas_clone.1) %>%
  summarize(Mean = mean(Freq))
names(np57_counts) <- c("Clone", 57)

np67_sub_table <- as.data.frame(table(np67_sub$rc_castag))
np67_sub_table <- np67_sub_table %>% inner_join(cell_castags, by = c("Var1" = "CasTag"))
np67_counts <- np67_sub_table %>%
  group_by(cas_clone.1) %>%
  summarize(Mean = mean(Freq))
names(np67_counts) <- c("Clone", 67)

# merge the different timepoints
np_counts <- np28_counts %>% full_join(np43_counts, by = "Clone")
np_counts <- np_counts %>% full_join(np57_counts, by = "Clone")
np_counts <- np_counts %>% full_join(np67_counts, by = "Clone")
np_counts[is.na(np_counts)] <- 0
np_counts[, -1] <- lapply(np_counts[, -1], function(x) x/sum(x, na.rm=TRUE) )
np_plot_prop_df <- np_counts %>%
  pivot_longer(!Clone, names_to = "timepoint", values_to = "prop")
np_plot_prop_df$timepoint <- as.numeric(np_plot_prop_df$timepoint)

# plot
ggplot(np_plot_prop_df, aes(y = prop, x = timepoint, fill = Clone)) +
  geom_alluvium(aes(alluvium = Clone), alpha= .9, color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() + theme_cowplot() + ylab("Proportion of reads")

# remaining mice
rp28_sub_table <- as.data.frame(table(rp28_sub$rc_castag))
rp28_sub_table <- rp28_sub_table %>% inner_join(cell_castags, by = c("Var1" = "CasTag"))
rp28_counts <- rp28_sub_table %>%
  group_by(cas_clone.1) %>%
  summarize(Mean = mean(Freq))
names(rp28_counts) <- c("Clone", 28)

rp43_sub_table <- as.data.frame(table(rp43_sub$rc_castag))
rp43_sub_table <- rp43_sub_table %>% inner_join(cell_castags, by = c("Var1" = "CasTag"))
rp43_counts <- rp43_sub_table %>%
  group_by(cas_clone.1) %>%
  summarize(Mean = mean(Freq))
names(rp43_counts) <- c("Clone", 43)

rp50_sub_table <- as.data.frame(table(rp50_sub$rc_castag))
rp50_sub_table <- rp50_sub_table %>% inner_join(cell_castags, by = c("Var1" = "CasTag"))
rp50_counts <- rp50_sub_table %>%
  group_by(cas_clone.1) %>%
  summarize(Mean = mean(Freq))
names(rp50_counts) <- c("Clone", 50)

rp_counts <- rp28_counts %>% full_join(rp43_counts, by = "Clone")
rp_counts <- rp_counts %>% full_join(rp50_counts, by = "Clone")
rp_counts[is.na(rp_counts)] <- 0
rp_counts[, -1] <- lapply(rp_counts[, -1], function(x) x/sum(x, na.rm=TRUE) )
rp_plot_prop_df <- rp_counts %>%
  pivot_longer(!Clone, names_to = "timepoint", values_to = "prop")
rp_plot_prop_df$timepoint <- as.numeric(rp_plot_prop_df$timepoint)

ggplot(rp_plot_prop_df, aes(y = prop, x = timepoint, fill = Clone)) +
  geom_alluvium(aes(alluvium = Clone), alpha= .9, color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() + theme_cowplot() + ylab("Proportion of reads")


lp67_sub_table <- as.data.frame(table(lp67_sub$rc_castag))
lp67_sub_table <- lp67_sub_table %>% inner_join(cell_castags, by = c("Var1" = "CasTag"))
lp67_counts <- lp67_sub_table %>%
  group_by(cas_clone.1) %>%
  summarize(Mean = mean(Freq))
names(lp67_counts) <- c("Clone", 67)

lp67_counts[, -1] <- lapply(lp67_counts[, -1], function(x) x/sum(x, na.rm=TRUE) )
lp_plot_prop_df <- lp67_counts %>%
  pivot_longer(!Clone, names_to = "timepoint", values_to = "prop")

ggplot(lp_plot_prop_df, aes(y = prop, x = timepoint, fill = Clone)) +
  geom_alluvium(aes(alluvium = Clone), alpha= .9, color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() + theme_cowplot() + ylab("Proportion of reads")


rlp57_sub_table <- as.data.frame(table(rlp57_sub$rc_castag))
rlp57_sub_table <- rlp57_sub_table %>% inner_join(cell_castags, by = c("Var1" = "CasTag"))
rlp57_counts <- rlp57_sub_table %>%
  group_by(cas_clone.1) %>%
  summarize(Mean = mean(Freq))
names(rlp57_counts) <- c("Clone", 57)

rlp_plot_prop_df <- rlp57_counts %>%
  pivot_longer(!Clone, names_to = "timepoint", values_to = "prop")

ggplot(rlp_plot_prop_df, aes(y = prop, x = timepoint, fill = Clone)) +
  geom_alluvium(aes(alluvium = Clone), alpha= .9, color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() + theme_cowplot() + ylab("Proportion of reads")
