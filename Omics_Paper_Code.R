library(rhdf5)
library(hdf5r)
library(Seurat)
library(scCustomize)

Donor <-as.data.frame(h5read("~/HumanFetalBrainPool.h5", "/shoji/Donor", ))
Age <-as.data.frame(h5read("~/HumanFetalBrainPool.h5", "/shoji/Age", ))
Tissue <-as.data.frame(h5read("~/HumanFetalBrainPool.h5", "/shoji/Tissue", ))

a <- as.numeric(rownames(subset(Age, `h5read("~/HumanFetalBrainPool.h5", "/shoji/Age", ` == "5.5")))
b <- as.numeric(rownames(subset(Tissue, `h5read("~/HumanFetalBrainPool.h5", "/shoji/Tissue", ` == "Forebrain")))
c <- intersect(a, b)

Gene <-as.data.frame(h5read("~/HumanFetalBrainPool.h5", "/shoji/Gene", ))
Gene$Gene <- Gene$`h5read("~/HumanFetalBrainPool.h5", "/shoji/Gene", `
Gene$`h5read("~/HumanFetalBrainPool.h5", "/shoji/Gene", ` <- NULL
cnts_GW7.5 <- as.data.frame(h5read("~/HumanFetalBrainPool.h5", "/shoji/Expression", index = list(1:59480, c)))
cnts_GW7.5$Gene <- Gene$Gene
cnts_GW7.5 <- as.data.frame(cnts_GW7.5[!duplicated(cnts_GW7.5$Gene),])
rownames(cnts_GW7.5) <- cnts_GW7.5$Gene
cnts_GW7.5$Gene <- NULL

CellID <-as.data.frame(h5read("~/HumanFetalBrainPool.h5", "/shoji/CellID", ))
CellID$CellID <- CellID$`h5read("~/HumanFetalBrainPool.h5", "/shoji/CellID", `
CellID$`h5read("~/HumanFetalBrainPool.h5", "/shoji/CellID", ` <- NULL
CellID <- as.data.frame(CellID[c,])
CellID$keep <- lapply(strsplit(as.character(CellID$CellID), ":"), "[", 2)
colnames(cnts_GW7.5) <- CellID$keep

GW7.5_forebrain <- CreateAssayObject(counts = cnts_GW7.5)
GW7.5_forebrain <- CreateSeuratObject(counts = GW7.5_forebrain@counts, project = "GW7.5_Forebrain")
GW7.5_forebrain <- NormalizeData(GW7.5_forebrain)
GW7.5_forebrain <- FindVariableFeatures(GW7.5_forebrain, method = "vst", nfeatures = 2000)
GW7.5_forebrain <- ScaleData(GW7.5_forebrain, verbose = FALSE)
GW7.5_forebrain <- RunPCA(GW7.5_forebrain, verbose = FALSE)
ElbowPlot(GW7.5_forebrain, ndims = 50)
GW7.5_forebrain <- RunUMAP(GW7.5_forebrain, dims = 1:30, verbose = FALSE)
GW7.5_forebrain <- FindNeighbors(GW7.5_forebrain)
GW7.5_forebrain <- FindClusters(GW7.5_forebrain)

#FeaturePlots
FeaturePlot_scCustom(GW7.5_forebrain, features = c("HES1", "EOMES", "NEUROD2","NR4A2", "LHX5", "GAD2", "AIF1", "COL1A1", "ESAM", "HBA2", "PDGFRA"))


#Cell Type Annotation
GW7.5_forebrain <- RenameIdents(GW7.5_forebrain,`1` = "Radial Glia", `1` = "Fibroblast", , `2` = "Radial Glia", `3` = "Radial Glia", `4` = "Fibroblast", `5` = "Radial Glia", `6` = "Fibroblast", `7` = "Interneuron", 
                                `8` = "Radial Glia", , `9` = "Fibroblast", `10` = "Radial Glia", `11` = "Cajal-Retzius Neuron", `12` = "Interneuron", `13` = "Interneuron", `14` = "Radial Glia",
                                `15` = "Intermediate Progenitor", , `16` = "Subplate Neuron", `17` = "Fibroblast", `18` = "Radial Glia", `19` = "Radial Glia", `20` = "Erythrocyte", `21` = "Cajal-Retzius Neuron",
                                `22` = "Interneuron", , `23` = "Microglia", `24` = "Vascular Cell")
GW7.5_forebrain$Celltype <- Idents(GW7.5_forebrain)


#DEG Lists
cell_type_markers <- FindAllMarkers(a, only.pos = TRUE, logfc.threshold = 0.4, min.pct = 0.10, min.cells.group = 10)
write.csv(cell_type_markers, "~/Gw7.5_forebrain_celltypes.csv")


#Make Pseudobulk Object
GW7.5_forebrain <- SetIdent(GW14_Cortex, value = "Celltype")
Pial_Cells <- subset(GW7.5_forebrain, idents = c("Vascular Cell", "Fibroblast"))
aggregate.Pial_Cells <- as.data.frame(AggregateExpression(Pial_Cells, assays = "RNA", slot = 'counts', group.by = "orig.ident"))
write.csv(aggregate.Pial_Cells, "~/GW7.5_forebrain_Pial_Cells.psuedobulk.csv")
