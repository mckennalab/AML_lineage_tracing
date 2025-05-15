# plotting editing timeline and patterns (figure 1)

##### packages #####
library(ggplot2)
library(cowplot)
library(devtools)
library(ShortRead) 
library(data.table)
library(dplyr)
library(stringr)
library(BiocManager)
library(radiant.data)
library(bayesbio)
library(grid)
`%nin%` = Negate(`%in%`)



# load files (processed via lineage_seq_preprocessing.R)
ctrl <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS33_c1498_timeline/control.csv")
d13_p <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS33_c1498_timeline/day13.csv")
d44_p <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS33_c1498_timeline/day44_presubset.csv")
d64_p_c1 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS33_c1498_timeline/day64_realcell_subset.csv")
d87_p_c2 <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS33_c1498_timeline/day87_c2_subset.csv")





##### editing over time (figure 1c) #####

# extracting columns with target site info
cntrl_sites <- ctrl[, c(4:11)]
d13_sites <- d13_p[, c(4:11)]
d44_sites <- d44_p[, c(4:11)]
d64_c1_sites <- d64_p_c1[, c(4:11)]
d87_c2_sites <- d87_p_c2[, c(4:11)]

# count the number of edited target sites for each read
cntrl_sites <- cntrl_sites %>%  mutate(edit_count = rowSums(.[1:8] != "NONE"))
d13_sites <- d13_sites %>%  mutate(edit_count = rowSums(.[1:8] != "NONE"))
d44_sites <- d44_sites %>%  mutate(edit_count = rowSums(.[1:8] != "NONE"))
d64_c1_sites <- d64_c1_sites %>%  mutate(edit_count = rowSums(.[1:8] != "NONE"))
d87_c2_sites <- d87_c2_sites %>%  mutate(edit_count = rowSums(.[1:8] != "NONE"))

# get it all into 1 df for plotting (certainly not the best method)
sample <- c(replicate(133088, 0), replicate(85711, 13), replicate(1170303, 44), replicate(210253, 64), replicate(1738606, 87))
prop_df <- as.data.frame(sample)
prop_df$count <- c(cntrl_sites$edit_count, d13_sites$edit_count, d44_sites$edit_count, d64_c1_sites$edit_count, d87_c2_sites$edit_count)

prop_df$count <- as_numeric(prop_df$count)
summary_df <- prop_df %>% dplyr::group_by(sample) %>%
  dplyr::summarise(mean = mean(count), 
                   sd = sd(count), 
                   n = n(),
                   se = sd / sqrt(n))


ggplot(summary_df, aes(x=sample, y=mean, fill = "test")) + 
  geom_bar(stat="identity", color="#b3cef2", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin = ifelse(mean - sd < 0, 0, mean - sd), ymax = ifelse(mean + sd > 8, 8, mean + sd)), width=.2,
                position=position_dodge(0.05)) +
  scale_fill_manual(values = "#b3cef2") + 
  theme_cowplot() +
  theme(legend.position = "none") +
  ylab("Edited targets per lineage tracer") +
  xlab("Day")



##### plot edit type (figure 1d) #####
clone_ids <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/other/neoclone_IDs.csv", header = TRUE, row.names = 1) #pre-identified founder tags
clone_ids <- clone_ids[clone_ids$clone %in% c(1,6), ] #just interested in F1 and F6
d87_p_c2 <- d87_p_c2[d87_p_c2$ftag %in% clone_ids$ID, ] 

# determine if edits are unique to one target or span multiple targets
edit_types <- as.data.frame(d87_p_c2$readName)
edit_types <- edit_types %>% 
  mutate(t1 = ifelse(d87_p_c2$target1 == "NONE", "unedited", ifelse(d87_p_c2$target1 == d87_p_c2$target2, "multi", "single"))) %>% 
  mutate(t2 = ifelse(d87_p_c2$target2 == "NONE", "unedited", ifelse(d87_p_c2$target2 == d87_p_c2$target1, "multi", ifelse(d87_p_c2$target2 == d87_p_c2$target3, "multi", "single")))) %>% 
  mutate(t3 = ifelse(d87_p_c2$target3 == "NONE", "unedited", ifelse(d87_p_c2$target3 == d87_p_c2$target2, "multi", ifelse(d87_p_c2$target3 == d87_p_c2$target4, "multi", "single")))) %>% 
  mutate(t4 = ifelse(d87_p_c2$target4 == "NONE", "unedited", ifelse(d87_p_c2$target4 == d87_p_c2$target3, "multi", ifelse(d87_p_c2$target4 == d87_p_c2$target5, "multi", "single")))) %>% 
  mutate(t5 = ifelse(d87_p_c2$target5 == "NONE", "unedited", ifelse(d87_p_c2$target5 == d87_p_c2$target4, "multi", ifelse(d87_p_c2$target5 == d87_p_c2$target6, "multi", "single")))) %>% 
  mutate(t6 = ifelse(d87_p_c2$target6 == "NONE", "unedited", ifelse(d87_p_c2$target6 == d87_p_c2$target5, "multi", ifelse(d87_p_c2$target6 == d87_p_c2$target7, "multi", "single")))) %>% 
  mutate(t7 = ifelse(d87_p_c2$target7 == "NONE", "unedited", ifelse(d87_p_c2$target7 == d87_p_c2$target6, "multi", ifelse(d87_p_c2$target7 == d87_p_c2$target8, "multi", "single")))) %>% 
  mutate(t8 = ifelse(d87_p_c2$target8 == "NONE", "unedited", ifelse(d87_p_c2$target8 == d87_p_c2$target7, "multi", "single")))

plot_df <- rbind(as.data.frame(table(edit_types$t1)), as.data.frame(table(edit_types$t2)), as.data.frame(table(edit_types$t3)), as.data.frame(table(edit_types$t4)), as.data.frame(table(edit_types$t5)), as.data.frame(table(edit_types$t6)), as.data.frame(table(edit_types$t7)), as.data.frame(table(edit_types$t8)))
plot_df$target <- c(replicate(3, 1), replicate(3, 2), replicate(3, 3), replicate(3, 4), replicate(3, 5), replicate(3, 6), replicate(3, 7), replicate(3, 8))

ggplot(plot_df, aes(fill = factor(Var1, levels=c("unedited", "single", "multi")), y = Freq, x = target)) + 
  geom_bar(position="fill", stat="identity") +
  theme_cowplot() +
  scale_fill_manual(values=c("#d6e5cb", "#628D56", "#335a30")) +
  scale_x_continuous(breaks = 1:8) +
  ylab("Proportion") + 
  xlab("Lineage recorder target site") + 
  guides(fill=guide_legend(title="Edit type"))



##### editing by integration site in founder 1 (figure 1e) #####
d87_p_c2$summary <- paste(d87_p_c2$target1, d87_p_c2$target2, d87_p_c2$target3, d87_p_c2$target4, d87_p_c2$target5, d87_p_c2$target6, d87_p_c2$target7, d87_p_c2$target8, sep = " ")

clone_ids <- clone_ids[clone_ids$clone == 1, ] #just interested in F1 for this plot
d87_p_c2 <- d87_p_c2[d87_p_c2$ftag %in% clone_ids$ID, ]

# identify fully unedited target arrays
d87_p_c2$status <- "edited"
d87_p_c2$status[d87_p_c2$summary == "NONE NONE NONE NONE NONE NONE NONE NONE"] <- "unedited"

edit_by_int <- d87_p_c2 %>% count(ftag, status)

ggplot(edit_by_int, aes(fill = factor(status, levels=c("unedited","edited" )), y = n, x = ftag)) + 
  geom_bar(position="fill", stat="identity") +
  theme_cowplot() +
  scale_fill_manual(values=c("#b3cef2", "#3f5d9b")) +
  theme(axis.text.x=element_blank()) + 
  ylab("Proportion edited") + 
  xlab("Founder 1 lineage recorder integration") + 
  guides(fill=guide_legend(title="Status"))

