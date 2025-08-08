library(data.table)
library(tidyverse)
aml_genes = fread("aml_CRISPR_screen_genes.tsv",header=T)
mouse = fread("broadgpp-brie-library-contents.txt")
human = fread("broadgpp-brunello-library-contents.txt")


# Load dplyr library
library(dplyr)

# names that have to get remapped - change name for human
# unique(aml_genes$mouse_id[!is.element(aml_genes$mouse_id,mouse$`Target Gene Symbol`)])


# Set of names to remap and their new names
target_names <- c("SYLT3", "LINC01089", "MIR142HG", "MORRBID", "CARD19",  "DNMT3a", "DNMT3b", "mTOR", "PDGFRa", "PDGFRb", "TOP-1", "MACROH2A2", "RIPOR2", "AL163541.1", "AL713998.1", "CPLANE1", "BEX3",    "AC008568.1", "ARMH1",    "CASC11", "OVAAL", "DUXAP8")
new_names <- c(   "SYLT3", "------",    "------",    "------", "C9orf89", "DNMT3A", "DNMT3B", "MTOR", "PDGFRA", "PDGFRB", "TOP1",  "H2AFY2",    "FAM65B", "------",     "------",     "C5orf42", "NGFRAP1", "------",     "C1orf228", "------", "------", "------")

# Create a named vector for remapping
name_map <- setNames(new_names, target_names)

aml_genes$remapped_human = aml_genes$human_id

# Remap names in 'name_col' if they are in 'target_names'
aml_genes <- aml_genes %>%
  mutate(remapped_human = case_when(
    remapped_human %in% names(name_map) ~ name_map[remapped_human],
    TRUE ~ remapped_human  # Keep original name if not in target_names
  ))

aml_genes = aml_genes[aml_genes$human_id != "",]
aml_genes = aml_genes[!duplicated(aml_genes$human_id),]
human_full_aml_genes = left_join(aml_genes,human,by=c("remapped_human"="Target Gene Symbol"))


# mouse
# surprizes -
# Frmd4b -- can't seem to find it in the Brie library?
# Csff1r -> Csf1r
# Smi -> SMO was the mapping for your human mapping 
# Ptprd -> not sure about, some options
# Rnf213 - not sure why I couldn't match it
target_names = c("Trbc2", "Cd244a", "Ighv1-77", "AI480526", "Mir142hg", "Morrbid", "Immp2l", "Ighv1-82", "Gm20627", "Card19",        "Gm19331", "Frmd4b", "Csff1r",   "Smi",    "Trp52",  "Muc16", "Ptprd",    "Rnf213",   "Robo2",  "Stat2", "Bcl211", "Cfb", "H2-Dma", "Mcx1", "Inpp4b", "Dshs1", "Apc6", "Macroh2a2", "Ripor2", "Mef2c", "A100a10", "Cplane1",      "Bex3", "Armh1", "Gm40318", "Sycp2l", "Igfr1", "Rnf25", "Grp160", "Aig1")
new_names =    c("------", "Cd244",  "--------", "--------", "--------", "-------", "------", "--------", "-------", "1110007C09Rik", "-------", "------", "------",  "------", "------", "------", "--------", "--------", "--------", "---", "---",    "---", "---",    "---", "---",    "---",    "---", "H2afy2",    "Fam65b", "---",    "---",    "2410089E03Rik", "---", "---",    "---",     "---",    "---",   "---",   "---",    "---")


# Create a named vector for remapping
name_map <- setNames(new_names, target_names)
aml_genes = fread("aml_CRISPR_screen_genes.tsv",header=T)

aml_genes$remapped_mouse = aml_genes$mouse_id

# Remap names in 'name_col' if they are in 'target_names'
aml_genes <- aml_genes %>%
  mutate(remapped_mouse = case_when(
    remapped_mouse %in% names(name_map) ~ name_map[remapped_mouse],
    TRUE ~ remapped_mouse  # Keep original name if not in target_names
  ))

aml_genes = aml_genes[aml_genes$mouse_id != "",]
aml_genes = aml_genes[!duplicated(aml_genes$mouse_id),]

mouse_full_aml_genes = left_join(aml_genes,mouse,by=c("remapped_mouse"="Target Gene Symbol"))
mouse_full_aml_genes$dataset = "mouse"
human_full_aml_genes$dataset = "human"

mouse_full_aml_genes = mouse_full_aml_genes[!is.na(mouse_full_aml_genes$`sgRNA Target Sequence`),]
human_full_aml_genes = human_full_aml_genes[!is.na(human_full_aml_genes$`sgRNA Target Sequence`),]

write.table(mouse_full_aml_genes,file="mouse_full_aml_genes_guides.txt",sep="\t",quote=F,row.names = F)
write.table(human_full_aml_genes,file="human_full_aml_genes_guides.txt",sep="\t",quote=F,row.names = F)