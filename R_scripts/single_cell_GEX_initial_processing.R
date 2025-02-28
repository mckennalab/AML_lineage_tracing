##### packages #####
library(ggplot2)
library(devtools)
library(ShortRead) 
library(data.table)
library(Seurat)
library(dplyr)
library(sctransform)
library(glmGamPoi)
library(stringr)
library(BiocManager)
library(viridis)
library(EnhancedVolcano)
library(radiant.data)
library(RColorBrewer)
library(ggnewscale)
library(grid)
library(scCustomize)
library(cowplot)
library(pals)
library(forcats)
library(ggsignif)
`%nin%` = Negate(`%in%`)


##### creating seurat object, demultiplexing MULTIseq, initial analysis #####

# creating object
lane_a_data <- Read10X(data.dir = "/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS22_invitro_EP/GEX/cellranger_multi_062624/laneA/outs/per_sample_outs/laneA/count/sample_filtered_feature_bc_matrix")
lane_a <- CreateSeuratObject(counts = lane_a_data$`Gene Expression`)
lane_a[['Sample']] <- CreateAssayObject(counts = lane_a_data$`Antibody Capture`) 
lane_a[["pct_mt"]] <- PercentageFeatureSet(lane_a, pattern = "^mt-") #change to "^MT" for human cells
VlnPlot(lane_a, features = c("nFeature_RNA", "nCount_RNA", "pct_mt"), ncol = 3) #to decide appropriate filtering cut offs
lane_a <- subset(lane_a, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & pct_mt < 8)
lane_a <- SCTransform(lane_a, vst.flavor = "v2", verbose = TRUE)

# demultiplexing
DefaultAssay(lane_a) <- "Sample"
lane_a <- NormalizeData(lane_a, assay = "Sample", normalization.method = "CLR")

lane_a <- MULTIseqDemux(lane_a,
                   assay = "Sample",
                   autoThresh = TRUE,
                   maxiter = 10,
                   qrange = seq(0.01, 0.99, by=0.02),
                   verbose = TRUE)

lane_a <- subset(lane_a, MULTI_ID == "Doublet", invert = TRUE)
lane_a <- subset(lane_a, MULTI_ID == "Negative", invert = TRUE) 

# initial analysis
DefaultAssay(lane_a) <- "SCT"
lane_a <- RunPCA(lane_a, verbose = FALSE)
ElbowPlot(lane_a, ndims = 50, reduction = "pca")
lane_a <- RunUMAP(lane_a, dims = 1:30, verbose = FALSE)
lane_a <- FindNeighbors(lane_a, dims = 1:30, verbose = FALSE)
lane_a <- FindClusters(lane_a, verbose = FALSE, resolution = 0.2)


##### integrating seurat objects (when applicable) #####

data_list <- c(lane_a, lane_b)
features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3000)
data_list <- PrepSCTIntegration(object.list = data_list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = data_list, normalization.method = "SCT",
                                  anchor.features = features)
sep_demux_combo <- IntegrateData(anchorset = anchors, normalization.method = "SCT") 

DefaultAssay(sep_demux_combo) <- "integrated"
sep_demux_combo <- RunPCA(sep_demux_combo, verbose = FALSE)
ElbowPlot(sep_demux_combo, ndims = 50, reduction = "pca")
sep_demux_combo <- RunUMAP(sep_demux_combo, dims = 1:30, verbose = FALSE)
sep_demux_combo <- FindNeighbors(sep_demux_combo, dims = 1:30, verbose = FALSE)
sep_demux_combo <- FindClusters(sep_demux_combo, verbose = FALSE, resolution = 0.4)

DefaultAssay(sep_demux_combo) <- "RNA"
sep_demux_combo_joined <- JoinLayers(sep_demux_combo)


##### cell type annotation for blood + bone marrow samples #####
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(dbplyr)

ls('package:celldex')
ref.set <- celldex::MouseRNAseqData()

blood_celldex <- as.SingleCellExperiment(blood)
pred.cnts <- SingleR::SingleR(test = blood_celldex, ref = ref.set, labels = ref.set$label.main)
lbls.keep <- table(pred.cnts$labels)>10
blood$SingleR.labels <- ifelse(lbls.keep[pred.cnts$labels], pred.cnts$labels, 'Other')
DimPlot(blood, reduction='umap', group.by='SingleR.labels')

pred.cnts_fine <- SingleR::SingleR(test = blood_celldex, ref = ref.set, labels = ref.set$label.fine)
lbls.keep_fine <- table(pred.cnts_fine$labels)>10
blood$SingleR.labels_fine <- ifelse(lbls.keep_fine[pred.cnts_fine$labels], pred.cnts_fine$labels, 'Other')
DimPlot(blood, reduction='umap', group.by='SingleR.labels_fine')



##### doublet removal for non multiplexed samples #####
library(irlba)
library(scDblFinder)

blood_sce <- as.SingleCellExperiment(blood)
blood_sce <- scDblFinder(blood_sce) 
blood_calls <- as.data.frame(blood_sce$scDblFinder.class)
rownames(blood_calls) <- rownames(blood_sce@colData)
blood <- AddMetaData(
  object = blood,
  metadata = blood_calls,
  col.name = "scDblFinder_class"
)
blood <- subset(blood, subset = scDblFinder_class  == "singlet")



##### cell cycle annotations #####

# mouse cells 
s.genes <- (c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Mlf1ip","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8"))
g2m.genes <- (c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2","AurkB","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa"))

strict_founder_cells <- CellCycleScoring(
  object = strict_founder_cells,
  g2m.features = g2m.genes,
  s.features = s.genes
)

# human cells 
hl60s_no_mid <- CellCycleScoring(
  object = hl60s_no_mid,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)