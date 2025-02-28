# Hotspot module overlap analysis
library(GeneOverlap)

# load modules, get human orthologs for C1498 modules
hl60_mods <- read.csv("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/RS30_RS31_humancellEP/GEX/hotspot/hl60_megamega_modlist.csv")
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


##### calculate jaccard indices #####
c1498_list <- list(f1m1, f1m2, f1m3, f1m4, f6m1, f6m2, f6m3, f6m4, f6m5)
names(c1498_list) <- c('f1m1', 'f1m2', 'f1m3', 'f1m4', 'f6m1', 'f6m2', 'f6m3', 'f6m4', 'f6m5')
hl60_list <- list(mega_mod1, mega_mod2, mega_mod3, mega_mod4, mega_mod5, mega_mod6, mega_mod7, mega_mod8, mega_mod9, mega_mod10, mega_mod11, mega_mod12, mega_mod13, mega_mod14, mega_mod15, mega_mod16, mega_mod17, mega_mod18)
names(hl60_list) <- c('hl60_m1', 'hl60_m2', 'hl60_m3', 'hl60_m4', 'hl60_m5', 'hl60_m6', 'hl60_m7', 'hl60_m8', 'hl60_m9', 'hl60_m10', 'hl60_m11', 'hl60_m12', 'hl60_m13', 'hl60_m14', 'hl60_m15', 'hl60_m16', 'hl60_m17', 'hl60_m18')

gom.obj <- newGOM(c1498_list, hl60_list, genome.size=NULL, 
                  spec='hg19.gene')
drawHeatmap(gom.obj, what = "Jaccard", adj.p = TRUE, grid.col="Blues", note.col="black", ncolused = 9) 



##### expression correlation between modules and EMT and cell cycle gene sets #####
DefaultAssay(strict_founder_cells) <- "RNA"
DefaultAssay(hl60s_no_mid) <- "RNA"


#mouse EMT sets 
reactome_tgfb_emt <- list(c('Tgfb1','Arhgef18','Pard6a','Tgfbr1','Rhoa','Ubc','Ubb','Rps27a','Pard3','Prkcz','Tgfbr2','F11r','Smurf1','Uba52rt','Cgn','Uba52'))
gotzmann_emt_up <- list(c('Gata2','Mpdz','Cenpa','Chek1','Kifc5b','Pros1','Abcb1b','Fn1','Pdgfa','Rhox5','Tgfb3','Serpine1','Pla2g4a','Map4','Ppic','Inhba','Dck','Chrnb1','Col5a2','Ctsl','Irgm1','Acta2','Ctla2a','Gadd45b','Mcm3','Oat','Inhbb','Fbln2','Tpp2','Mki67','Glg1','Ncam1','Pkp1','Cct6a','Top1','Col18a1','Edn1','Tead2','Rpsa','Col4a1','Col4a2','Clu','Furin','Plk4','Igfbp7','Ly6a','Ccl2','Csf1','Il4ra','Cdh2','Zfx','Ccn2','Snai1','Odc1','Cyp1b1','Col1a1','Rbl1','Rbbp6','Traf3','Brca1','Smarcd1','Smad1','Pafah1b1','Timp1','Ctsl','Jun','Col3a1','Tnc','Nr2f1','Itgb1','Hmgb1','Ccl7','Timp3'))
gotzmann_emt_dn <- list(c('Pigq','Elp5','Exoc7','Aqp8','Hs3st1','Cnih1','Slc22a18','Ech1','Atp5f1b','Anxa8','Tcea3','Cd82','Fdft1','Xrn2','Fnta','Pigf','Sdf2','Trim25','Polr1d','Stac','Plxna2','Ddost','Enpp1','Pfkl','Mcpt2','Or1j1','Itgb4','Arhgdib','Btc','Anxa7','Raly','Atp5mc1','Eps8','Eps15','Ddx19a','Ahcyl','Gdi2','Il1r1','Tmem165','Dynlt1b','Ppp1cb','Calm2','Srp14','Zfp36l1','Pthlh','Fkbp2','Abcd3','Bax','Wt1','Krt15','Nfkbia','Gas8','Pa2g4','Atp6v1a','Cd53','Slc3a2','Srsf1','Clcn3','Cbx1','Ptprm','Sec23a','Gnai1','Clns1a','Tle1','Lbp','Tob1','Hells','Hccs','Lgals3','Gjb4','Prkg2','Gpd1','Tfam','Ivl','Rab18','Man1a','Acadm','Klra2','Pdcd2','Atp6v1e1','Capza2','Ube2h','Acadl','Cpe','Serpinb6a','Cxcl5','Stip1','Stk3','Dbp','Vdac1','Gss','Npepps','Klf3','Dmbt1','Psmc1','Sephs2','Tpd52','Ddx24','Fgfbp1','Rock1','Ostf1','Cyp2j6','Ufd1','Anxa11','Noct','Cops5','Fxyd5','Nherf1','Fhl1','Dok1','Ppfibp2','Pam','Gata5','Aip','Bnc1','Sec22b','Hsd17b10','Alad','Pfn1','Tpm3','Edn2','Prkacb','Cyp3a13','Mmp13','Mcpt1','Serpine2','Itih2','Bcap31','Sprr1a','Selenop','S100a13','Rnf2','Hnrnph1','Sos1','Mrc1','Cd53','Psmc5','Slc35a1','H3f3a','Snap23','Fosl1','Phyh','Gm14270','Top2a','Gstp2','Ptprr','Cd55','Ywhaz','Polr2c','Ywhah','Ywhaq','Arf4','Psme1','Sdc4','Cxcl1','Mthfd2','Gstm1','Trp53','Myc','Ugt1a2','Foxa2','Tm4sf1','Gnb5','Slc1a5','Tnnt2','Cdkn2a','H2-T23','Penk','Tcea1','Egr1','Hmbs','Krt19','Abcb1a','Loricrin','H2-T10','Ier2','Sod1','Cftr','Il1rn','Atf2','S100a8','Ctnnb1','Scp2','Cpt2','Pura','Tkt','Pparg','Cdkn1a','Mapk14','Stom','Klf4','Slc12a1','Mxi1','Gna11','H2-K1','Oaz1','Casp3','Lgals9','Ntan1','Casp4','Tfam','Ucp2','Eya2','Mrpl23','Cx3cl1','Fos','Cdh1','Ap2a2','Gja1','Hmmr','Kcnab1','Hmga2','Amd1','Cct3'))
jechlinger_emt_up <- list(c('Ifitm3','Htra1','B2m','Rras','Ctsz','Ackr3','Serpinh1','Fmo1','Upp1','Cxcl1','Mthfd2','Prl2c3','Srm','Phgdh','Ifit3','Procr','Csn3','Ada','Cfh','Lamb1','C4b','Ccl2','Vim','Cdh2','Pmp22','Ptgs1','Ppic','Adss1','Cdh15','Mmp12','S100a8','Mmp2','Pdgfra','Ptpn22','Bcl3','Snai1','Sdc2','Cyp1b1','Il11','Vldlr','Stat1','Dab2','Slc3a2','Cxcl5','Pla2g7','Asns','Ifit1','Gbp3','Slpi','Irf7','Rnaset2b','Galk1','Sparc','Pdgfrb','Col3a1','Dcn','Tnc','Isg15','Pcolce','Cck','Gas1','Col6a2','Col6a1','Mmp13','Cd68','Inhba','Tnxb','Ddr2','Hif1a','Sdc1'))
jechlinger_emt_dn <- list(c("Flna","Actn4","Vamp8","Atp1a1","Itgb5","Chka","Sgk1","Padi2","Fbp2","Nrp1","Car2","Pcx","Sat1","Bmp4","Krt14","Egr1","F3","Id1","Zfp239","Tgfb3","Tgm2","Zfp36","Btg2","Id2","Ccn2","Epcam","Cyp2f2","Thbs1","Inmt","Jup","Grb7","Prkcz","Plk2","Tiam1","Ctsh","Atf3","Stat5a","Klf2","Bcl6","Fzd8","Kitl","Serpinb5","Arhgef1","Numb","Irf6","Fos","Myh9","Cdh1","Egr2","Itpr1","Nr4a1","Dusp1","Tsc22d1","Hmmr","Id4","Hmga2","Pkp1","Amd1","Ctnnd1","Timp3","Klf10","Nnt"))
s.genes <- list(c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Mlf1ip","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8"))
g2m.genes <- list(c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2","AurkB","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa"))

# correlation in strict founder cells (C1498 in vitro)
strict_founder_cells <- AddModuleScore(strict_founder_cells, features = reactome_tgfb_emt, name="reactome_tgfb_emt") 
strict_founder_cells <- AddModuleScore(strict_founder_cells, features = gotzmann_emt_up, name="gotzmann_emt_up") 
strict_founder_cells <- AddModuleScore(strict_founder_cells, features = gotzmann_emt_dn, name="gotzmann_emt_dn") 
strict_founder_cells <- AddModuleScore(strict_founder_cells, features = jechlinger_emt_up, name="jechlinger_emt_up") 
strict_founder_cells <- AddModuleScore(strict_founder_cells, features = jechlinger_emt_dn, name="jechlinger_emt_dn")
strict_founder_cells <- AddModuleScore(strict_founder_cells, features = s.genes, name="s.genes")
strict_founder_cells <- AddModuleScore(strict_founder_cells, features = g2m.genes, name="g2m.genes")

mouse_emt_corr <- cbind(as.data.frame(strict_founder_cells$hotspot_f1_mod1_1), as.data.frame(strict_founder_cells$hl60_mod2_1), as.data.frame(strict_founder_cells$reactome_tgfb_emt1), as.data.frame(strict_founder_cells$gotzmann_emt_up1), as.data.frame(strict_founder_cells$gotzmann_emt_dn1), as.data.frame(strict_founder_cells$jechlinger_emt_up1), as.data.frame(strict_founder_cells$jechlinger_emt_dn1), as.data.frame(strict_founder_cells$s.genes1), as.data.frame(strict_founder_cells$g2m.genes1))
names(mouse_emt_corr) <- c("f1m1", "mega_mod2", "reactome_tgfb_emt", "gotzmann_emt_up", "gotzmann_emt_dn", "jechlinger_emt_up", "jechlinger_emt_dn", "s.genes", "g2m.genes")
sfc_m1_test_cormat_kendall <- cor(mouse_emt_corr, method = "kendall")
sfc_m1_test_cormat_pearson <- cor(mouse_emt_corr, method = "pearson")
sfc_m1_test_cormat_spearman <- cor(mouse_emt_corr, method = "spearman")
corrplot(
  sfc_m1_test_cormat_kendall)



#human EMT sets
gotzmann_emt_up <- list(c('GATA2','MPDZ','CENPA','CHEK1','KIFC1','PROS1','ABCB1','FN1','PDGFA','TGFB3','SERPINE1','PLA2G4A','MAP4','PPIC','INHBA','DCK','CHRNB1','COL5A2','CTSV','IRGM','ACTA2','GADD45B','MCM3','OAT','INHBB','FBLN2','TPP2','MKI67','GLG1','NCAM1','PKP1','CCT6A','TOP1','COL18A1','EDN1','TEAD2','RPSA','COL4A1','COL4A2','CLU','FURIN','PLK4','IGFBP7','CSF1','IL4R','CDH2','ZFX','CCN2','SNAI1','ODC1','CYP1B1','COL1A1','RBL1','RBBP6','TRAF3','BRCA1','SMARCD1','SMAD1','PAFAH1B1','TIMP1','CTSV','JUN','COL3A1','TNC','NR2F1','ITGB1','HMGB1','CCL7','TIMP3'))
gotzmann_emt_dn <- list(c("IGQ","ELP5","EXOC7","AQP8","HS3ST1","CNIH1","SLC22A18","ECH1","ATP5F1B","ANXA8L1","TCEA3","CD82","FDFT1","XRN2","FNTA","PIGF","SDF2","TRIM25","POLR1D","STAC","PLXNA2","DDOST","ENPP1","PFKL","RPL18","ITGB4","ARHGDIB","BTC","ANXA7","RALY","ATP5MC1","EPS8","EPS15","DDX19A","AHCY","GDI2","IL1R1","TMEM165","DYNLT1","PPP1CB","CALM2","SRP14","ZFP36L1","PTHLH","FKBP2","ABCD3","BAX","WT1","KRT15","NFKBIA","GAS8","PA2G4","ATP6V1A","CD53","SLC3A2","SRSF1","CLCN3","CBX1","PTPRM","SEC23A","GNAI1","CLNS1A","TLE1","LBP","TOB1","HELLS","HCCS","LGALS3","GJB4","PRKG2","GPD1","TFAM","IVL","RAB18","MAN1A1","ACADM","PDCD2","ATP6V1E1","CAPZA2","UBE2H","ACADL","CPE","SERPINB6","CXCL6","STIP1","STK3","DBP","VDAC1","GSS","NPEPPS","KLF3","DMBT1","PSMC1","SEPHS2","TPD52","DDX24","FGFBP1","ROCK1","OSTF1","CYP2J2","UFD1","ANXA11","NOCT","COPS5","FXYD5","NHERF1","FHL1","DOK1","PPFIBP2","PAM","GATA5","AIP","BNC1","SEC22B","HSD17B10","ALAD","PFN1","TPM3","EDN2","PRKACB","CYP3A7","MMP13","SERPINE2","ITIH2","BCAP31","SPRR1A","SELENOP","S100A13","RNF2","HNRNPH1","SOS1","MRC1","PSMC5","SLC35A1","H3-3A","SNAP23","FOSL1","PHYH","HAX1","TOP2A","GSTP1","PTPRR","CD55","YWHAZ","POLR2C","YWHAH","YWHAQ","ARF4","PSME1","SDC4","CXCL3","MTHFD2","GSTM1","TP53","MYC","UGT1A10","FOXA2","TM4SF1","GNB5","SLC1A5","TNNT2","CDKN2A","HLA-E","PENK","TCEA1","EGR1","HMBS","KRT19","ABCB1","LORICRIN","IER2","SOD1","CFTR","IL1RN","ATF2","S100A8","CTNNB1","SCP2","CPT2","PURA","TKT","PPARG","CDKN1A","MAPK14","STOM","KLF4","SLC12A1","MXI1","GNA11","HLA-A","OAZ1","CASP3","LGALS9","NTAN1","CASP4","UCP2","EYA2","MRPL23","CX3CL1","FOS","CDH1","AP2A2","GJA1","HMMR","KCNAB1","HMGA2","AMD1","CCT3"))
jechlinger_emt_dn <- list(c('FLNA','VAMP8','ACTN4','ATP1A1','ITGB5','CHKA','SGK1','PADI2','FBP2','NRP1','CA2','PC','SAT1','BMP4','KRT14','EGR1','F3','ID1','ZNF239','TGFB3','TGM2','ZFP36','BTG2','ID2','CCN2','EPCAM','CYP2F1','THBS1','INMT','JUP','GRB7','PRKCZ','PLK2','TIAM1','CTSH','ATF3','STAT5A','KLF2','BCL6','FZD8','KITLG','SERPINB5','ARHGEF1','NUMB','IRF6','FOS','ATP1A1','ITGB5','MYH9','CDH1','EGR2','ITPR1','NR4A1','DUSP1','TSC22D1','HMMR','ID4','HMGA2','PKP1','AMD1','CTNND1','TIMP3','KLF10','NNT'))
jechlinger_emt_up <- list(c('IFITM3','HTRA1','B2M','RRAS','CTSZ','ACKR3','SERPINH1','FMO1','UPP1','CXCL3','MTHFD2','SRM','PHGDH','IFIT3','PROCR','CSN3','ADA','CFH','LAMB1','C4B','VIM','CDH2','PMP22','PTGS1','PPIC','ADSS1','CDH15','MMP12','S100A8','MMP2','PDGFRA','PTPN22','BCL3','SNAI1','SDC2','CYP1B1','IL11','VLDLR','STAT1','DAB2','SLC3A2','CXCL6','PLA2G7','ASNS','IFIT1B','GBP4','SLPI','IRF7','RNASET2','GALK1','SPARC','PDGFRB','COL3A1','DCN','TNC','ISG15','PCOLCE','CCK','GAS1','COL6A2','COL6A1','MMP13','CD68','INHBA','TNXB','DDR2','HIF1A','SDC1'))
reactome_tgfb_emt <- list(c('RHOA','PRKCZ','FKBP1A','PARD6A','ARHGEF18','TGFB1','TGFBR1','CGN','RPS27A','PARD3','UBC','F11R','TGFBR2','UBB','SMURF1','UBA52'))
s.genes <- list(c(cc.genes$s.genes))
g2m.genes <- list(c(cc.genes$g2m.genes))

# correlation in HL60s 
hl60s_no_mid <- AddModuleScore(hl60s_no_mid, features = reactome_tgfb_emt, name="reactome_tgfb_emt") 
hl60s_no_mid <- AddModuleScore(hl60s_no_mid, features = gotzmann_emt_up, name="gotzmann_emt_up") 
hl60s_no_mid <- AddModuleScore(hl60s_no_mid, features = gotzmann_emt_dn, name="gotzmann_emt_dn") 
hl60s_no_mid <- AddModuleScore(hl60s_no_mid, features = jechlinger_emt_up, name="jechlinger_emt_up") 
hl60s_no_mid <- AddModuleScore(hl60s_no_mid, features = jechlinger_emt_dn, name="jechlinger_emt_dn")
hl60s_no_mid <- AddModuleScore(hl60s_no_mid, features = s.genes, name="s.genes")
hl60s_no_mid <- AddModuleScore(hl60s_no_mid, features = g2m.genes, name="g2m.genes")

human_emt_corr <- cbind(as.data.frame(hl60s_no_mid$f1m11), as.data.frame(hl60s_no_mid$mega_mod2_1), as.data.frame(hl60s_no_mid$reactome_tgfb_emt1), as.data.frame(hl60s_no_mid$gotzmann_emt_up1), as.data.frame(hl60s_no_mid$gotzmann_emt_dn1), as.data.frame(hl60s_no_mid$jechlinger_emt_up1), as.data.frame(hl60s_no_mid$jechlinger_emt_dn1), as.data.frame(hl60s_no_mid$s.genes1), as.data.frame(hl60s_no_mid$g2m.genes1))
names(human_emt_corr) <- c("f1m1", "mega_mod2", "reactome_tgfb_emt", "gotzmann_emt_up", "gotzmann_emt_dn", "jechlinger_emt_up", "jechlinger_emt_dn", "s.genes", "g2m.genes")
cormat_kendall <- cor(human_emt_corr, method = "kendall")
cormat_pearson <- cor(human_emt_corr, method = "pearson")
cormat_spearman <- cor(human_emt_corr, method = "spearman")
corrplot(
  cormat_kendall)