library("orthogene")

# ----------------------------------------------------------------------------------------------------------------
# COAD and PDAC gene lists
# ----------------------------------------------------------------------------------------------------------------
coad = fread("COAD_mutated_genes_no_PCT.txt")
pdac = fread("PDAC_mutated_genes_no_PCT.txt")
colnames(pdac) = c("gene","mutsig","mut_count","saMPles_with_mut","saMPle_count","freq","is_cancer_gene")
colnames(coad) = c("gene","mutsig","mut_count","saMPles_with_mut","saMPle_count","freq","is_cancer_gene")

coad_genes = coad[coad$freq >= 10]
pdac_genes = pdac[pdac$freq >= 5]

# add notch, MYC, FGFR1, CDKN2B, and CXCL1 to the list
new_data <- data.table(
  gene = c("NOTCH1", "MYC", "FGFR1", "CDKN2B", "CXCL1"),
  mutsig = c(NA, NA, NA, NA, NA),
  mut_count = c(0, 0, 0, 0, 0),  # Made-up values
  saMPles_with_mut = c(0, 0, 0, 0, 0),  # Made-up values
  saMPle_count = c(150, 150, 150, 150, 150),
  freq = c(0.0, 0.0, 0.0, 0.0, 0.0),  # Made-up values
  is_cancer_gene = c("Yes", "Yes", "Yes", "Yes", "No")  # Made-up values
)

# Combine existing and new data
pdac_genes <- rbind(pdac_genes, new_data)

# ----------------------------------------------------------------------------------------------------------------
# some stem cell lists
# ----------------------------------------------------------------------------------------------------------------

stem_cell = as.data.table(c("LGR5", "MYC", "FGFR1", "CDKN2B", "CXCL1", "ZNF800","ATOH1","SOX9","GFI1",
                            "KLF4","SPDEF","GFI1","POU2F3","SOX4","SPIB","HES1","HNF4A","HNF4G",
                            "SPIB","IL-4","IL-13","HES1","ATOH1","ARX","PAX4","SPIB","POU2F3","SOX4",
                            "SPDEF","GFI1","NEUROG3",
                            "NOTCH1","NOTCH2","NOTCH2NL","NOTCH3","NOTCH4"))
stem_cell$mutsig = 0

stem_cell$mut_count = 0
stem_cell$saMPles_with_mut = 0
stem_cell$saMPle_count = 0
stem_cell$freq = 0
stem_cell$is_cancer_gene = "No"

# ----------------------------------------------------------------------------------------------------------------
# master list
# ----------------------------------------------------------------------------------------------------------------
full_list = stem_cell[,1]
full_list$target_lib = "stem_cell"
colnames(full_list) = c("gene","target_lib")

pdac_temp = pdac_genes[,1]
pdac_temp$target_lib = "pdac_list"
colnames(pdac_temp) = c("gene","target_lib")

coad_temp = coad_genes[,1]
coad_temp$target_lib = "coad_list"
colnames(coad_temp) = c("gene","target_lib")

full_list = rbind(full_list,coad_temp,pdac_temp)
# ----------------------------------------------------------------------------------------------------------------
# pathway libraries
# ----------------------------------------------------------------------------------------------------------------

gene_table = fread("all_new_genes_and_targets_output.txt")
colnames(gene_table) = c("gene","desc","target_lib")
gene_table_subset = gene_table[,c("gene","target_lib")]
full_list = rbind(full_list,gene_table_subset)
gene_table_unique = full_list[!duplicated(full_list$gene),]


gene_df <- orthogene::convert_orthologs(gene_df = gene_table_unique,
                                        gene_input = "gene", 
                                        gene_output = "columns", 
                                        input_species = "human",
                                        output_species = "mouse",
                                        non121_strategy = "drop_both_species") 

matched_genes = left_join(full_list,gene_df,by=c("gene"="input_gene"))

# ----------------------------------------------------------------------------------------------------------------
# combine to make a human gene list
# ----------------------------------------------------------------------------------------------------------------
brunello = fread("broadgpp-brunello-library-contents.txt")
human_guides = left_join(matched_genes,brunello,by=c("gene"="Target Gene Symbol"))
human_guides$species = "human"
human_guides_filtered = human_guides[!is.element(human_guides$target_lib,c("bmp_library","wnt_library","notch_library")),]
write.table(human_guides_filtered,"human_guides_filtered.txt",quote=F,row.names = F,sep="\t")

brie = fread("broadgpp-brie-library-contents.txt")
mouse_guides = left_join(matched_genes[!is.na(matched_genes$ortholog_gene),],brie,by=c("ortholog_gene"="Target Gene Symbol"))
mouse_guides$species = "mouse"
write.table(mouse_guides,"mouse_guides_filtered.txt",quote=F,row.names = F,sep="\t")

# ----------------------------------------------------------------------------------------------------------------
# make a master table for the twist order
# ----------------------------------------------------------------------------------------------------------------

mouse_aml_genes = fread("mouse_full_aml_genes_guides.txt")
human_aml_genes = fread("human_full_aml_genes_guides.txt")

total_guide_table = mouse_aml_genes[,c("mouse_id","sgRNA Target Sequence")]
colnames(total_guide_table) = c("gene","guide")
total_guide_table$target_lib = "AML"
total_guide_table$species = "mouse"

tmp_add = house_aml_genes[,c("human_id","sgRNA Target Sequence")]
colnames(tmp_add) = c("gene","guide")
tmp_add$target_lib = "AML"
tmp_add$species = "human"
total_guide_table = rbind(total_guide_table,tmp_add)

tmp_add = human_guides_filtered[,c("gene","sgRNA Target Sequence","target_lib")]
colnames(tmp_add) = c("gene","guide","target_lib")
tmp_add$species = "human"
total_guide_table = rbind(total_guide_table,tmp_add)

tmp_add = mouse_guides[,c("gene","sgRNA Target Sequence","target_lib")]
colnames(tmp_add) = c("gene","guide","target_lib")
tmp_add$species = "mouse"
total_guide_table = rbind(total_guide_table,tmp_add)

write.table(total_guide_table,file="2024_11_02_total_guide_table.txt",sep="\t",quote=F,row.names=F)

# drop the human notch genes and add controls for each mouse and human
total_guide_table2 = total_guide_table[total_guide_table$target_lib != "notch_library" & total_guide_table$species != "human",]