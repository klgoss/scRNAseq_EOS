library(Seurat)
library(dplyr)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(DESeq2)
library(tibble)

set.seed(123)

#read patient A01 data
a01.data <- Read10X("A01_EOS/filtered_feature_bc_matrix/")
# create seurat object
a01 <- CreateSeuratObject(counts = a01.data, project = "A01_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("A01_EOS/A01_EOS_keep.csv")
a01 <- subset(a01, cells = keep$Barcode)

# quality control
a01[["percent.mt"]] <- PercentageFeatureSet(a01, pattern = "^MT-", assay = "RNA")
VlnPlot(a01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
a01 <- subset(a01, subset = nFeature_RNA > 100 & nFeature_RNA < 1500 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 2000)

# read patient A02 data
a02.data <- Read10X("A02_EOS/filtered_feature_bc_matrix/")
# create seurat object
a02 <- CreateSeuratObject(counts = a02.data, project = "A02_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("A02_EOS/A02_EOS_keep.csv")
a02 <- subset(a02, cells = keep$Barcode)

# quality control
a02[["percent.mt"]] <- PercentageFeatureSet(a02, pattern = "^MT-")
VlnPlot(a02, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(a02, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
a02 <- subset(a02, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 3000)

# laod patient A03 data
a03.data <- Read10X("A03_EOS/filtered_feature_bc_matrix/")
# create seurat object
a03 <- CreateSeuratObject(counts = a03.data, project = "A03_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("A03_EOS/A03_EOS_keep.csv")
a03 <- subset(a03, cells = keep$Barcode)

# quality control
a03[["percent.mt"]] <- PercentageFeatureSet(a03, pattern = "^MT-")
VlnPlot(a03, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
a03 <- subset(a03, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 3000)

#  laod patient A04 data
a04.data <- Read10X("A04_EOS/filtered_feature_bc_matrix/")
# create seurat object
a04 <- CreateSeuratObject(counts = a04.data, project = "A04_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("A04_EOS/A04_EOS_keep.csv")
a04 <- subset(a04, cells = keep$Barcode)

# quality control
a04[["percent.mt"]] <- PercentageFeatureSet(a04, pattern = "^MT-")
VlnPlot(a04, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
a04 <- subset(a04, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 3000)

#  laod patient A05 data
a05.data <- Read10X("A05_EOS/filtered_feature_bc_matrix/")
# create seurat object
a05 <- CreateSeuratObject(counts = a05.data, project = "A05_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("A05_EOS/A05_EOS_keep.csv")
a05 <- subset(a05, cells = keep$Barcode)

# quality control
a05[["percent.mt"]] <- PercentageFeatureSet(a05, pattern = "^MT-")
VlnPlot(a05, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
a05 <- subset(a05, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 3000)

#  laod patient A06 data
a06.data <- Read10X("A06_EOS/filtered_feature_bc_matrix/")
# create seurat object
a06 <- CreateSeuratObject(counts = a06.data, project = "A06_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("A06_EOS/A06_EOS_keep.csv")
a06 <- subset(a06, cells = keep$Barcode)

# quality control
a06[["percent.mt"]] <- PercentageFeatureSet(a06, pattern = "^MT-")
VlnPlot(a06, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
a06 <- subset(a06, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 2000)

# read patient b01 data
b01.data <- Read10X("B01_EOS/filtered_feature_bc_matrix/")
# create seurat object
b01 <- CreateSeuratObject(counts = b01.data, project = "B01_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("B01_EOS/B01-EOS-keep.csv")
b01 <- subset(b01, cells = keep$Barcode)

# quality control
b01[["percent.mt"]] <- PercentageFeatureSet(b01, pattern = "^MT-")
VlnPlot(b01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
b01 <- subset(b01, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 4000)

# read patient b02 data
b02.data <- Read10X("B02_EOS/filtered_feature_bc_matrix/")
# create seurat object
b02 <- CreateSeuratObject(counts = b02.data, project = "B02_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("B02_EOS/B02_EOS_keep.csv")
b02 <- subset(b02, cells = keep$Barcode)

# quality control
b02[["percent.mt"]] <- PercentageFeatureSet(b02, pattern = "^MT-")
VlnPlot(b02, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
b02 <- subset(b02, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 3000)

# load patient B03 data
b03.data <- Read10X("B03_EOS/filtered_feature_bc_matrix/")
# create seurat object
b03 <- CreateSeuratObject(counts = b03.data, project = "B03_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("B03_EOS/B03_EOS_keep.csv")
b03 <- subset(b03, cells = keep$Barcode)

# quality control
b03[["percent.mt"]] <- PercentageFeatureSet(b03, pattern = "^MT-")
VlnPlot(b03, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
b03 <- subset(b03, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 3000)

# load patient B05 data and create Seurat object
b05.data <- Read10X("B05_EOS/filtered_feature_bc_matrix/")
b05 <- CreateSeuratObject(counts = b05.data, project = "B05_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("B05_EOS/B05_EOS_keep.csv")
b05 <- subset(b05, cells = keep$Barcode)

# quality control
b05[["percent.mt"]] <- PercentageFeatureSet(b05, pattern = "^MT-")
VlnPlot(b05, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
b05 <- subset(b05, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 3000)

# load patient C01 data and create Seurat object
c01.data <- Read10X("C01_EOS/filtered_feature_bc_matrix/")
c01 <- CreateSeuratObject(counts=c01.data, project = 'C01_EOS', min.cells = 3, min.features = 0)
keep <- read.csv("C01_EOS/C01_EOS_keep.csv")
c01 <- subset(c01, cells = keep$Barcode)
c01[["percent.mt"]] <- PercentageFeatureSet(c01, pattern = "^MT-")
VlnPlot(c01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
c01 <- subset(c01, subset = nFeature_RNA > 100 & nFeature_RNA < 1200 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 2000)

# load patient C02 data and create Seurat object
c02.data <- Read10X("C02_EOS/filtered_feature_bc_matrix/")
c02 <- CreateSeuratObject(counts = c02.data, project = "C02_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("C02_EOS/C02_EOS_keep.csv")
c02 <- subset(c02, cells = keep$Barcode)
c02[["percent.mt"]] <- PercentageFeatureSet(c02, pattern = "^MT-")
VlnPlot(c02, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
c02 <- subset(c02, subset = nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 3000)

# load patient C03 data and create Seurat object
c03.data <- Read10X("C03_EOS/filtered_feature_bc_matrix/")
c03 <- CreateSeuratObject(counts = c03.data, project = "C03_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("C03_EOS/C03_EOS_keep.csv")
c03 <- subset(c03, cells = keep$Barcode)
c03[["percent.mt"]] <- PercentageFeatureSet(c03, pattern = "^MT-")
VlnPlot(c03, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
c03 <- subset(c03, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 2000)

# load patient C04 data and create Seurat object
c04.data <- Read10X("C04_EOS/filtered_feature_bc_matrix/")
c04 <- CreateSeuratObject(counts = c04.data, project = "C04_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("C04_EOS/C04_EOS_keep.csv")
c04 <- subset(c04, cells = keep$Barcode)
c04[["percent.mt"]] <- PercentageFeatureSet(c04, pattern = "^MT-")
VlnPlot(c04, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
c04 <- subset(c04, subset = nFeature_RNA > 100 & nFeature_RNA < 1500 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 2000)

a01$severity <- "mild"
a02$severity <- "mild"
a03$severity <- "mild"
a04$severity <- "mild"
a05$severity <- "mild"
a06$severity <- "mild"
b01$severity <- "severe"
b02$severity <- "severe"
b03$severity <- "severe"
b05$severity <- "severe"
c01$severity <- "healthy"
c02$severity <- "healthy"
c03$severity <- "healthy"
c04$severity <- "healthy"

a01$patient <- "A01"
a02$patient <- "A02"
a03$patient <- "A03"
a04$patient <- "A04"
a05$patient <- "A05"
a06$patient <- "A06"
b01$patient <- "B01"
b02$patient <- "B02"
b03$patient <- "B03"
b05$patient <- "B05"
c01$patient <- "C01"
c02$patient <- "C02"
c03$patient <- "C03"
c04$patient <- "C04"

# analyze combined datasets
obj <- merge(a01, y = c(a02, a03, a04, a05, a06, b01, b02, b03, b05, c01, c02, c03, c04), add.cell.ids = c("a01", "a02", "a03", "a04", "a05", "a06", "b01", "b02", "b03", "b05", "c01", "c02", "c03", "c04"), project = "EOS")

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$patient)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, vars.to.regress = "percent.mt")
obj <- RunPCA(obj)

ElbowPlot(obj, ndims=30) #ndims to use = 20

obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.1)
obj <- RunUMAP(obj, dims = 1:20)

DimPlot(obj)
DimPlot(obj, group.by = "patient")
DimPlot(obj, group.by = "severity")

# Integrate dataset
obj <- IntegrateLayers(object = obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca")

# re-join layers after integration
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:20, reduction = "integrated.cca")

#saveRDS(obj, file = "EOS_integrated.rds")

DimPlot(obj, group.by = "patient")
DimPlot(obj, group.by = "severity")
FeaturePlot(obj, features = "percent.mt", cols=c("grey", "red"))
FeaturePlot(obj, features = "nCount_RNA",cols=c("grey", "red"))
FeaturePlot(obj, features = "nFeature_RNA", cols=c("grey", "red"))

clusters <- FindAllMarkers(obj, only.pos = T, logfc.threshold = 0.25, min.pct = 0.1)
clusters <- clusters[clusters$p_val_adj < 0.05,]

top20 <- clusters %>%
  group_by(cluster) %>%
  top_n(20, avg_log2FC)

FeaturePlot(obj, features = "CSF3R", order=T)
FeaturePlot(obj, features = "FCGR3B", order=T)
FeaturePlot(obj, features = "CAMP", order=T)
FeaturePlot(obj, features = "MPO", order=T)
FeaturePlot(obj, features = "ELANE", order=T)
FeaturePlot(obj, features = "CXCR1", order=T)

FeaturePlot(obj, features = "FCGR3A", order=T)
FeaturePlot(obj, features = "ADGRE1", order=T)
FeaturePlot(obj, features = "IL5RA", order=T)
FeaturePlot(obj, features = "CLC", order=T)
FeaturePlot(obj, features = "CCR3", order=T)
FeaturePlot(obj, features = "ITGAM", order=T)
FeaturePlot(obj, features = "ADGRE1", order=T)
FeaturePlot(obj, features = "CPA3", order=T)

FeaturePlot(obj, features = "DYSF", order=T)
FeaturePlot(obj, features = "PLBD1", order=T)
FeaturePlot(obj, features = "ARG1", order=T)
FeaturePlot(obj, features = "ABCA13", order=T)
FeaturePlot(obj, features = "CRISP3", order=T)
FeaturePlot(obj, features = "MMP8", order=T)
FeaturePlot(obj, features = "S100A12", order=T)
FeaturePlot(obj, features = "FCAR", order=T)
FeaturePlot(obj, features = "MMP9", order=T)
FeaturePlot(obj, features = "CLEC4D", order=T)
FeaturePlot(obj, features = "PTGDR2", order=T)


obj <- FindSubCluster(obj, cluster = 6, subcluster.name = "sub6", graph.name = "RNA_snn", resolution = 0.3)
DimPlot(obj, group.by = "sub6")

# clusteres 6_1 and 7 do not express EOS markers at high level. Let's remove those
Idents(obj) <- obj$sub6
obj <- subset(obj, idents = c("6_1", "7"), invert=T)

DimPlot(obj, label = T, group.by = "sub6")

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$patient)

obj <- NormalizeData(obj)

# look at classic EOS markers again
FeaturePlot(obj, features = "ADGRE1", order=T)
FeaturePlot(obj, features = "IL5RA", order=T)
FeaturePlot(obj, features = "CLC", order=T)
FeaturePlot(obj, features = "SIGLEC8", order=T)
FeaturePlot(obj, features = "RNASE2", order=T)
FeaturePlot(obj, features = "RNASE3", order=T)

# re-join layers 
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

#saveRDS(obj, file = "EOS_integrated_updated.rds")

obj <- readRDS("EOS_integrated_updated.rds")

########################################
# QC plots for supplementary figure 1 #
#######################################
qc_1 <- VlnPlot(obj, features = "nCount_RNA", group.by = "patient") + theme(legend.position = "none")
qc_2 <- VlnPlot(obj, features = "nFeature_RNA", group.by = "patient") + theme(legend.position = "none")
qc_3 <- VlnPlot(obj, features = "percent.mt", group.by = "patient") + theme(legend.position = "none")

ggsave(qc_1, filename="../AsthmaSeq-EOS-paper/figures/QC_VLN_1.svg")
ggsave(qc_2, filename="../AsthmaSeq-EOS-paper/figures/QC_VLN_2.svg")
ggsave(qc_3, filename="../AsthmaSeq-EOS-paper/figures/QC_VLN_3.svg")

mycolors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4" ,"#B2DF8A", "#33A02C","#FB9A99", "orchid")

qc_4 <- FeaturePlot(obj, features = 'nCount_RNA', order=T)
qc_4_1 <- DimPlot(obj, group.by = "patient", cols = mycolors, shuffle=T) + ggtitle("")
qc_5 <- FeaturePlot(obj, features = 'nFeature_RNA', order=T)
qc_6 <- FeaturePlot(obj, features = 'percent.mt', order=T)

ggsave(qc_4, filename="../AsthmaSeq-EOS-paper/figures/QC_VLN_4.png")
ggsave(qc_4_1, filename = "../AsthmaSeq-EOS-paper/figures/QC_VLN4_1.png")
ggsave(qc_5, filename="../AsthmaSeq-EOS-paper/figures/QC_VLN_5.png")
ggsave(qc_6, filename="../AsthmaSeq-EOS-paper/figures/QC_VLN_6.png")

##########################
# HEALTHY VS ASTHMA EOS #
#########################

obj$severity1 <- ifelse(obj$severity == "healthy", yes = "healthy", no = "asthma")

# downsample asthma cells so that we have a closer number to healthy. There are 4214 healthy cells
#Idents(obj) <- obj$severity1
#downsample <- subset(obj, downsample = 4200)

#downsample <- downsample %>%
#  NormalizeData() %>%
#  ScaleData() %>% RunPCA()

#ElbowPlot(downsample)

#downsample <- FindNeighbors(downsample, reduction = "integrated.cca", dims = 1:15)
#downsample <- FindClusters(downsample, resolution = 0.3)
#downsample <- RunUMAP(downsample, dims = 1:15, reduction = "integrated.cca")

#downsample <- readRDS("../AsthmaSeq-EOS-paper/REVISION/EOS_downsampled.rds")


###################
### FIGURE 3A ####
##################

umap1 <- DimPlot(downsample, pt.size=1) + theme_void() + theme(legend.position="none")
ggsave(umap1, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthyVasthma_umap1.png")

umap2 <- DimPlot(downsample, group.by = "severity1", cols = c("red", "#008ECE"), shuffle=T, pt.size=1) + theme_void() + theme(legend.position="none") + ggtitle("")
ggsave(umap2, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthyVasthma_umap2.png")

mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))

mycolors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4" ,"#B2DF8A", "#33A02C","#FB9A99", "orchid")

umap3 <- DimPlot(downsample, group.by = "patient", shuffle=T, pt.size=1, cols = mycolors) + theme_void() + theme(legend.position="none") + ggtitle("")
ggsave(umap3, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthyVasthma_umap3.png")

###################
### FIGURE 3B ####
##################
bar <- dittoSeq::dittoBarPlot(downsample, var = "severity1", group.by = "RNA_snn_res.0.3", color.panel = c("red", "#008ECE"))
ggsave(bar, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthyVasthma_bar.svg")


###################
### FIGURE 3C ####
##################
umap_HealthInf_highlight <- DimPlot(downsample, cells.highlight = healthy_inflammatory, pt.size=0.75, cols = "lightgray")  + theme_void() + theme(legend.position="none")
ggsave(umap_HealthInf_highlight, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/HealthyvAsthmaEos_HealthInflamHighlight.png")


###################
### FIGURE 3D ####
##################

healthy_vs_asthma_eos_downsample <- FindMarkers(downsample, ident.1 = "healthy", ident.2 = "asthma", group.by = "severity1", logfc.threshold = 0.25, min.pct = 0.1)
healthy_vs_asthma_eos_downsample <- healthy_vs_asthma_eos_downsample[healthy_vs_asthma_eos_downsample$p_val_adj < 0.05,]
write.table(healthy_vs_asthma_eos_downsample, file = "HealthyVAsthma_DEG.tsv", sep = "\t", quote=F)

keyvals <- ifelse(
  healthy_vs_asthma_eos_downsample$avg_log2FC < 0, 'red',
  ifelse(healthy_vs_asthma_eos_downsample$avg_log2FC > 0, '#008ECE',
         'black'))

names(keyvals)[keyvals == 'red'] <- 'Asthma'
names(keyvals)[keyvals == '#008ECE'] <- 'Healthy'

volcano <- EnhancedVolcano::EnhancedVolcano(healthy_vs_asthma_eos_downsample, lab = rownames(healthy_vs_asthma_eos_downsample), x='avg_log2FC',y = 'p_val_adj',FCcutoff = 0.25,pCutoff = 0.05, title = 'Healthy vs Asthma', subtitle = "", legendPosition = 'bottom', colCustom = keyvals, gridlines.major = F, gridlines.minor = F)
ggsave(volcano, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/HealthyVAsthma_volcano.svg")


###################
### FIGURE 3E ####
##################
rnase2 <- VlnPlot(downsample, "RNASE2", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(RNASE2~ident, data = rnase2$data)
ggsave(rnase2, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/RNASE2.svg")

rnase3 <- VlnPlot(downsample, "RNASE3", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(RNASE3~ident, data = rnase3$data)
ggsave(rnase3, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/rnase3.svg")

ccr3 <- VlnPlot(downsample, "CCR3", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(CCR3~ident, data = ccr3$data)

b2m <- VlnPlot(downsample, "B2M", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(B2M~ident, data = b2m$data)

tapbp <- VlnPlot(downsample, "TAPBP", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(TAPBP~ident, data = tapbp$data)

hla_a <- VlnPlot(downsample, "HLA-A", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(`HLA-A`~ident, data = hla_a$data)

hla_b <- VlnPlot(downsample, "HLA-B", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(`HLA-B`~ident, data = hla_b$data)

hla_c <- VlnPlot(downsample, "HLA-C", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(`HLA-C`~ident, data = hla_c$data)

cd74 <- VlnPlot(downsample, "CD74", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")

vln <- VlnPlot(downsample, features = c("HLA-A", "HLA-C", "TAPBP", "RNASE2", "RNASE3", "CCR3"), group.by = "severity1", flip=T, stack=T) + theme(legend.position="none")
ggsave(vln, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/HealthyVAsthma_violin.svg")


###################
#### FIGURE 4A ###
##################
healthy_vs_asthma_eos_downsample <- FindMarkers(downsample, ident.1 = "healthy", ident.2 = "asthma", group.by = "severity1", logfc.threshold = 0.25, min.pct = 0.1)
healthy_vs_asthma_eos_downsample <- healthy_vs_asthma_eos_downsample[order(-healthy_vs_asthma_eos_downsample$avg_log2FC),]

gene_list <- healthy_vs_asthma_eos_downsample$avg_log2FC
names(gene_list) <- rownames(healthy_vs_asthma_eos_downsample)
gene_list

m_df <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

m_df$symbol <- getSYMBOL(as.character(m_df$entrez_gene), data = 'org.Hs.eg')

H_t2g <- m_df %>% dplyr::select(gs_name, symbol)

em2 <- GSEA(gene_list, TERM2GENE = H_t2g, eps = 1e-300, pvalueCutoff = 0.05)
df <- as.data.frame(em2@result)

df$direction <- ifelse(df$NES > 0, yes = "healthy", no = "asthma")
df$Description <- factor(df$Description, levels = df$Description[order(df$NES)])
ggplot(df, aes(x = Description, y = NES, fill = direction)) + geom_bar(stat="identity") + coord_flip() + theme_classic() + scale_fill_manual("legend", values = c("healthy" = "#008ECE", "asthma" = "red")) + theme(text = element_text(size=20))

ifng <- enrichplot::gseaplot2(em2, geneSetID = 1, title = "HALLMARK INTERFERON GAMMA RESPONSE")
ifna <- enrichplot::gseaplot2(em2, geneSetID = 2, title = "HALLMARK INTERFERON ALPHA RESPONSE")
inflam <- enrichplot::gseaplot2(em2, geneSetID = 4, title = "HALLMARK INFLAMMATORY RESPONSE")

ggsave(ifng, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/IFNG_line.svg")
ggsave(ifna, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/IFNA_line.svg")
ggsave(inflam, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/inflam_line.svg")

###################
#### FIGURE 4B ###
##################

IFNA_response <- df[df$Description == "HALLMARK_INTERFERON_ALPHA_RESPONSE",]$core_enrichment
IFNA_response_genes <- unlist(strsplit(IFNA_response, split = "/"))
cat(IFNA_response_genes, file = "IFNA_response_genes.txt", sep = ",")

IFNG_response <- df[df$Description == "HALLMARK_INTERFERON_GAMMA_RESPONSE",]$core_enrichment
IFNG_response_genes <- unlist(strsplit(IFNG_response, split = "/"))
cat(IFNG_response_genes, file = "IFNG_response_genes.txt", sep = ",")

InflammatoryResponse <- df[df$Description == "HALLMARK_INFLAMMATORY_RESPONSE",]$core_enrichment
InflammatoryResponseGenes <- unlist(strsplit(InflammatoryResponse, split = "/"))
cat(InflammatoryResponseGenes, file="InflammatoryResponseGenes.txt", sep = ",")

library(GOplot)
circ <- readxl::read_xlsx("Inflam_forCord.xlsx", sheet = 1)


healthy_vs_asthma_eos_downsample <- healthy_vs_asthma_eos_downsample[healthy_vs_asthma_eos_downsample$p_val_adj < 0.05,]

genes <- healthy_vs_asthma_eos_downsample[rownames(healthy_vs_asthma_eos_downsample) %in% IFNA_response_genes | rownames(healthy_vs_asthma_eos_downsample) %in% IFNG_response_genes | rownames(healthy_vs_asthma_eos_downsample) %in% InflammatoryResponseGenes,]
genes$gene <- rownames(genes)
rownames(genes) <- NULL
genes <- genes %>% dplyr::select(gene, avg_log2FC)
colnames(genes) <- c("ID", "logFC")
genes <- na.omit(genes)

circ <- circle_dat(circ, genes)

genes <- data.frame(circ$genes, circ$logFC) %>% unique()

process <- circ$term %>% unique()

chord <- chord_dat(circ, genes, process)
chord <- na.omit(chord)

a <- GOChord(chord, gene.size = 4, ribbon.col = c("#80CDC1", "#B2ABD2","#FDB863"))
ggsave(a, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/HealthyVAsthma_Cord.svg")


#######################
# mild vs severe only #
#######################

mild_vs_severe <- subset(obj, severity == "healthy", invert = T)

mild_vs_severe <- NormalizeData(mild_vs_severe)

mild_vs_severe <- FindNeighbors(mild_vs_severe, reduction = "integrated.cca", dims = 1:20)
mild_vs_severe <- FindClusters(mild_vs_severe, resolution = 0.2)
mild_vs_severe <- RunUMAP(mild_vs_severe, dims = 1:20, reduction = "integrated.cca")

######################
### FIGURE 5A AND 5B #
######################
umap1 <- DimPlot(mild_vs_severe, pt.size=1) + ggtitle("") + theme_void() + theme(legend.position="none")
ggsave(g, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/mildVsevere_umap1.png")

umap2 <- DimPlot(mild_vs_severe, group.by = "severity", cols = c("#6B00A0", "darkorange"), shuffle=T, pt.size=1) + ggtitle("") + theme_void()
ggsave(umap2, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/mildVsevere_umap2.png")

################
### FIGURE 5C #
################
i <- dittoSeq::dittoBarPlot(mild_vs_severe, var = "severity", group.by = "RNA_snn_res.0.2", color.panel = c("#6B00A0", "darkorange"), xlab="condition", ylab="Fraction of Cells", legend.title="Cluster", main = NULL)
ggsave(i, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/mildVsevere_clusterBar.svg")


################
### FIGURE 5D #
################
mild_vs_severe_deg <- FindMarkers(mild_vs_severe, ident.1="mild", ident.2="severe", group.by = "severity", min.pct = 0.1, logfc.threshold = 0.25)
mild_vs_severe_deg <- mild_vs_severe_deg[mild_vs_severe_deg$p_val_adj < 0.05,]
write.table(mild_vs_severe_deg, file = "MildvsSevere_deg.tsv", sep = "\t", quote=F)

keyvals <- ifelse(
  mild_vs_severe_deg$avg_log2FC < 0, 'darkorange',
  ifelse(mild_vs_severe_deg$avg_log2FC > 0, '#6B00A0',
         'black'))

names(keyvals)[keyvals == 'darkorange'] <- 'Severe'
names(keyvals)[keyvals == '#6B00A0'] <- 'Mild'

j <- EnhancedVolcano::EnhancedVolcano(mild_vs_severe_deg, pointSize=5, lab = rownames(mild_vs_severe_deg), x='avg_log2FC',y = 'p_val_adj',FCcutoff = 0.25,pCutoff = 0.05, title = 'Mild vs Severe', subtitle = "", legendPosition = 'bottom', colCustom = keyvals, gridlines.major = F, gridlines.minor = F)
ggsave(j, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/MildvSevere_Volcano.svg")

####################
#### FIGURE 6A ####
####################

# GSEA
mild_vs_severe_deg <- FindMarkers(mild_vs_severe, ident.1="mild", ident.2="severe", group.by = "severity", min.pct = 0.1, logfc.threshold = 0.25)
mild_vs_severe_deg <- mild_vs_severe_deg[order(-mild_vs_severe_deg$avg_log2FC),]

gene_list <- mild_vs_severe_deg$avg_log2FC
names(gene_list) <- rownames(mild_vs_severe_deg)
gene_list

m_df <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

m_df$symbol <- getSYMBOL(as.character(m_df$entrez_gene), data = 'org.Hs.eg')

H_t2g <- m_df %>% dplyr::select(gs_name, symbol)

em2 <- GSEA(gene_list, TERM2GENE = H_t2g, eps = 1e-300, pvalueCutoff = 0.05, pAdjustMethod = "bonferroni")
df <- as.data.frame(em2@result)

df$direction <- ifelse(df$NES > 0, yes = "mild", no = "severe")
df$Description <- factor(df$Description, levels = df$Description[order(df$NES)])
bar <- ggplot(df, aes(x = Description, y = NES, fill = direction)) + geom_bar(stat="identity") + coord_flip() + theme_classic() + scale_fill_manual("legend", values = c("mild" = "#6B00A0", "severe" = "orange"))
ggsave(bar, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/MildVsevere_GSEA.svg")

IFNA_response <- df[df$Description == "HALLMARK_INTERFERON_ALPHA_RESPONSE",]$core_enrichment
IFNA_response_genes <- unlist(strsplit(IFNA_response, split = "/"))

IFNG_response <- df[df$Description == "HALLMARK_INTERFERON_GAMMA_RESPONSE",]$core_enrichment
IFNG_response_genes <- unlist(strsplit(IFNG_response, split = "/"))

mild_vs_severe <- AddModuleScore(mild_vs_severe, features = list(IFNA_response_genes), name = "IFNA_response")
mild_vs_severe <- AddModuleScore(mild_vs_severe, features = list(IFNG_response_genes), name = "IFNG_response")

####################
#### FIGURE 6B ####
####################
ifna <- enrichplot::gseaplot2(em2, geneSetID = 1, title = "Hallmark IFNA response")
ggsave(ifna, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/MildvSevere_IFNA.svg")

ifna_feat <- FeaturePlot(mild_vs_severe, features = "IFNA_response1", order=T, pt.size = 0.5) + scale_color_gradientn(colours = rev(brewer.pal(n=11, name="RdBu"))) + theme_void() + ggtitle("")
ggsave(ifna_feat, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/IFNA_feature.svg")

####################
#### FIGURE 6C ####
####################
ifng <- enrichplot::gseaplot2(em2, geneSetID = 2, title = "Hallmark IFNG response")
ggsave(ifng, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/MildvSevere_IFNG.svg")

ifng_feat <- FeaturePlot(mild_vs_severe, features = "IFNG_response1", order=T, pt.size = 0.5) + scale_color_gradientn(colours = rev(brewer.pal(n=11, name="RdBu"))) + theme_void() + ggtitle("")
ggsave(ifng_feat, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/IFNG_feature.svg")


####################
#### FIGURE 6D ####
####################
FeaturePlot(mild_vs_severe, features = "IFI6", order=T, pt.size = 0.5) + scale_color_viridis(option = "magma")
FeaturePlot(mild_vs_severe, features = "MX1", order=T, pt.size = 0.5) + scale_color_viridis(option = "magma")
FeaturePlot(mild_vs_severe, features = "ISG15", order=T, pt.size = 0.5) + scale_color_viridis(option = "magma")
FeaturePlot(mild_vs_severe, features = "XAF1", order=T, pt.size = 0.5) + scale_color_viridis(option = "magma")
FeaturePlot(mild_vs_severe, features = "RSAD2", order=T, pt.size = 0.5) + scale_color_viridis(option = "magma")
FeaturePlot(mild_vs_severe, features = "IFIT3", order=T, pt.size = 0.5) + scale_color_viridis(option = "magma")

isgs <- FeaturePlot(mild_vs_severe, features = c("IFI6", "MX1", "ISG15", "XAF1", "RSAD2", "IFIT3"), order=T, pt.size=0.7, ncol=3) & scale_color_viridis(option = "magma") & theme_void() & theme(legend.position="none") & ggtitle("")
ggsave(isgs, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/ISGs_feat.png")

##############################################################################################################

# HEALTHY MNC ANALYSIS #

eos <- readRDS("EOS_integrated_updated.rds")
healthy_eos <- subset(eos, severity == "healthy")

# EOS #
# for each file we will make sure to only include the same EOS we used in the downstream analysis
# load patient C01 data and create Seurat object
c01.data <- Read10X("C01_EOS/filtered_feature_bc_matrix/")
c01.eos <- CreateSeuratObject(counts=c01.data, project = 'C01_EOS', min.cells = 3, min.features = 0)
keep <- read.csv("C01_EOS/C01_EOS_keep.csv")
c01.eos <- subset(c01.eos, cells = keep$Barcode)
c01.eos[["percent.mt"]] <- PercentageFeatureSet(c01.eos, pattern = "^MT-")
VlnPlot(c01.eos, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
c01.eos <- subset(c01.eos, subset = nFeature_RNA > 100 & nFeature_RNA < 1200 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 2000)
c01.eos$annot <- "EOS"
colnames(c01.eos) <- paste("c01", colnames(c01.eos), sep = "_")
c01.eos <- subset(c01.eos, cells = colnames(healthy_eos))

# load patient C02 data and create Seurat object
c02.data <- Read10X("C02_EOS/filtered_feature_bc_matrix/")
c02.eos <- CreateSeuratObject(counts = c02.data, project = "C02_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("C02_EOS/C02_EOS_keep.csv")
c02.eos <- subset(c02.eos, cells = keep$Barcode)
c02.eos[["percent.mt"]] <- PercentageFeatureSet(c02.eos, pattern = "^MT-")
VlnPlot(c02.eos, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
c02.eos <- subset(c02.eos, subset = nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 3000)
c02.eos$annot <- "EOS"
colnames(c02.eos) <- paste("c02", colnames(c02.eos), sep = "_")
c02.eos <- subset(c02.eos, cells = colnames(healthy_eos))

# load patient C03 data and create Seurat object
c03.data <- Read10X("C03_EOS/filtered_feature_bc_matrix/")
c03.eos <- CreateSeuratObject(counts = c03.data, project = "C03_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("C03_EOS/C03_EOS_keep.csv")
c03.eos <- subset(c03.eos, cells = keep$Barcode)
c03.eos[["percent.mt"]] <- PercentageFeatureSet(c03.eos, pattern = "^MT-")
VlnPlot(c03.eos, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
c03.eos <- subset(c03.eos, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 2000)
c03.eos$annot <- "EOS"
colnames(c03.eos) <- paste("c03", colnames(c03.eos), sep = "_")
c03.eos <- subset(c03.eos, cells = colnames(healthy_eos))

# load patient C04 data and create Seurat object
c04.data <- Read10X("C04_EOS/filtered_feature_bc_matrix/")
c04.eos <- CreateSeuratObject(counts = c04.data, project = "C04_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("C04_EOS/C04_EOS_keep.csv")
c04.eos <- subset(c04.eos, cells = keep$Barcode)
c04.eos[["percent.mt"]] <- PercentageFeatureSet(c04.eos, pattern = "^MT-")
VlnPlot(c04.eos, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
c04.eos <- subset(c04.eos, subset = nFeature_RNA > 100 & nFeature_RNA < 1500 & percent.mt < 10 & nCount_RNA > 200 & nCount_RNA < 2000)
c04.eos$annot <- "EOS"
colnames(c04.eos) <- paste("c04", colnames(c04.eos), sep = "_")
c04.eos <- subset(c04.eos, cells = colnames(healthy_eos))

# MNCs #
# get only healthy patient cells and grab cell annotations
mnc <- readRDS("MNC_integrated.rds")
mnc <- subset(mnc, severity == "healthy")

meta <- mnc[[]]
annot <- meta %>% dplyr::select(updated_clusters)

annot[c('patient', 'cellID')] <- stringr::str_split_fixed(rownames(annot), '_', 2)
annot$group <- "mnc"

annot$patient.group <- paste(annot$patient, annot$group, sep = ".")
annot$patient.group.cellID <- paste(annot$patient.group, annot$cellID, sep = "_")

# load patient C02 data
c02.data <- Read10X("C02_MNC/filtered_feature_bc_matrix/")
# create seurat object
c02.mnc <- CreateSeuratObject(counts = c02.data, project = "C02_MNC", min.cells = 3, min.features = 200)
# quality control
c02.mnc[["percent.mt"]] <- PercentageFeatureSet(c02.mnc, pattern = "^MT-")
VlnPlot(c02.mnc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
c02.mnc <- subset(c02.mnc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15 & nCount_RNA > 100 & nCount_RNA < 20000)

c02.annot <- annot[annot$patient == "c02",]
c02.annot <- c02.annot %>% remove_rownames() %>% column_to_rownames(var="cellID")
c02.annot <- c02.annot %>% dplyr::select(updated_clusters)
c02.mnc <- AddMetaData(c02.mnc, metadata = c02.annot, col.name = "annot")

# load patient C03 data
c03.data <- Read10X("C03_MNC/filtered_feature_bc_matrix/")
# create Seurat object
c03.mnc <- CreateSeuratObject(counts = c03.data, project = "C03_MNC", min.cells = 3, min.features = 200)
# quality control
c03.mnc[["percent.mt"]] <- PercentageFeatureSet(c03.mnc, pattern = "^MT-")
VlnPlot(c03.mnc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
c03.mnc <- subset(c03.mnc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15 & nCount_RNA > 100 & nCount_RNA < 25000)

c03.annot <- annot[annot$patient == "c03",]
c03.annot <- c03.annot %>% remove_rownames() %>% column_to_rownames(var="cellID")
c03.annot <- c03.annot %>% dplyr::select(updated_clusters)
c03.mnc <- AddMetaData(c03.mnc, metadata = c03.annot, col.name = "annot")

# load patient C04 data
c04.data <- Read10X("C04_MNC/filtered_feature_bc_matrix/")
# create Seurat object
c04.mnc <- CreateSeuratObject(counts = c04.data, project = "C04_MNC", min.cells = 3, min.features = 200)
# quality control
c04.mnc[["percent.mt"]] <- PercentageFeatureSet(c04.mnc, pattern = "^MT-")
VlnPlot(c04.mnc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
c04.mnc <- subset(c04.mnc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15 & nCount_RNA > 100 & nCount_RNA < 10000)

c04.annot <- annot[annot$patient == "c04",]
c04.annot <- c04.annot %>% remove_rownames() %>% column_to_rownames(var="cellID")
c04.annot <- c04.annot %>% dplyr::select(updated_clusters)
c04.mnc <- AddMetaData(c04.mnc, metadata = c04.annot, col.name = "annot")

c01.eos$patient <- "C01"

c02.eos$patient <- "C02"
c02.mnc$patient <- "C02"

c03.eos$patient <- "C03"
c03.mnc$patient <- "C03"

c04.eos$patient <- "C04"
c04.mnc$patient <- "C04"

# Neutrophils #
# load sample 1 data
neut1.data <- Read10X("../AsthmaSeq-EOS-paper/neutrophils/sample1/filtered_feature_bc_matrix/")
# create seurat object
neut1 <- CreateSeuratObject(counts = neut1.data, project = "neut1", min.cells = 3, min.features = 0)
# quality control
neut1[["percent.mt"]] <- PercentageFeatureSet(neut1, pattern = "^MT-")
VlnPlot(neut1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
neut1 <- subset(neut1, subset = nFeature_RNA > 200 & nFeature_RNA < 1200 & percent.mt < 10 & nCount_RNA > 100 & nCount_RNA < 3000)

# load sample 2 data
neut2.data <- Read10X("../AsthmaSeq-EOS-paper/neutrophils/sample2/filtered_feature_bc_matrix/")
# create seurat object
neut2 <- CreateSeuratObject(counts = neut2.data, project = "neut2", min.cells = 3, min.features = 0)
# quality control
neut2[["percent.mt"]] <- PercentageFeatureSet(neut2, pattern = "^MT-")
VlnPlot(neut2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
neut2 <- subset(neut2, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 & percent.mt < 10 & nCount_RNA > 100 & nCount_RNA < 3000)

neut1$patient <- "N1"
neut2$patient <- "N2"

neut1$annot <- "neutrophil"
neut2$annot <- "neutrophil"

# analyze combined datasets
obj <- merge(c01.eos, y = c(c02.eos, c02.mnc, c03.eos, c03.mnc, c04.eos, c04.mnc, neut1, neut2), add.cell.ids = c("c01.eos", "c02.eos", "c02.mnc", "c03.eos", "c03.mnc", "c04.eos", "c04.mnc", "n1", "n2"))

obj <- subset(obj, annot == "RBC" | annot == "platelet" | annot == "endothelial" | annot == "granulocyte", invert = T)

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$patient)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

ElbowPlot(obj, ndims = 30) #ndims 25

obj <- FindNeighbors(obj, dims = 1:25)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:25)

DimPlot(obj)
DimPlot(obj, group.by = "patient")
DimPlot(obj, group.by = "annot")

obj <- IntegrateLayers(object = obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca")

obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:25)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:25, reduction = "integrated.cca")

#############
# Figure 1B #
#############
a <-  DimPlot(obj, group.by = "annot", label=F, pt.size=0.5, cols = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#E68613", "#0CB702", "#00B8E7", "#ED68ED", "#CD9600", "#00BE67", "#00A9FF"), shuffle = T) + ggtitle("") + theme_void()
ggsave(a, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthyEOS_MNC_umap.png")

obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

#############
# Figure 1C #
#############
Idents(obj) <- obj$annot
my_levels <- c("CD14+ monocyte", "CD16+ monocyte", "ILC2", "CD4", "CD8", "NK", "EOS", "neutrophil", "B", "DC", "pDC")
obj@active.ident <- factor(x = obj@active.ident, levels = my_levels)

markers <- FindMarkers(obj, ident.1 = "EOS", only.pos = T, logfc.threshold = 0.25)
markers <- markers[markers$p_val_adj < 0.05,]
write.table(markers, file = "healthyEOS_vs_healthyMNC_markers.tsv", quote = F, sep = "\t")

top500 <- markers %>% top_n(500, avg_log2FC)

avg <- AverageExpression(obj, features = rownames(top500), return.seurat = T)
heat <- DoHeatmap(avg, features = rownames(top500), draw.lines = F) + theme(axis.text.y = element_text(size = 0)) + theme(legend.position="none")
ggsave(heat, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/EOS_pbmc_top500heat.svg")

#############
# Figure 1D #
#############
ccr3_violin <- VlnPlot(obj, features = "CCR3", cols = c("#7CAE00", "#00BFC4","#ED68ED", "#E68613", "#C77CFF","#0CB702", "#00A9FF", "#CD9600", "#F8766D", "#00BE67",  "#00B8E7")) + theme(legend.position = "none")
ggsave(ccr3_violin, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthyEOS_CCR3_violin.svg")

il5ra_violin <- VlnPlot(obj, features = "IL5RA", cols = c("#7CAE00", "#00BFC4","#ED68ED", "#E68613", "#C77CFF","#0CB702", "#00A9FF", "#CD9600", "#F8766D", "#00BE67",  "#00B8E7")) + theme(legend.position = "none")
ggsave(il5ra_violin, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthyEOS_IL5RA_violin.svg")

siglec8_violin <- VlnPlot(obj, features = "SIGLEC8", cols = c("#7CAE00", "#00BFC4","#ED68ED", "#E68613", "#C77CFF","#0CB702", "#00A9FF", "#CD9600", "#F8766D", "#00BE67",  "#00B8E7")) + theme(legend.position = "none")
ggsave(siglec8_violin, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthyEOS_SIGLEC8_violin.svg")

adgre1_violin <- VlnPlot(obj, features = "ADGRE1", cols = c("#7CAE00", "#00BFC4","#ED68ED", "#E68613", "#C77CFF","#0CB702", "#00A9FF", "#CD9600", "#F8766D", "#00BE67",  "#00B8E7")) + theme(legend.position = "none")
ggsave(adgre1_violin, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthyEOS_ADGRE1_violin.svg")

itgam_violin <- VlnPlot(obj, features = "ITGAM", cols = c("#7CAE00", "#00BFC4","#ED68ED", "#E68613", "#C77CFF","#0CB702", "#00A9FF", "#CD9600", "#F8766D", "#00BE67",  "#00B8E7")) + theme(legend.position = "none")
ggsave(itgam_violin, filename = "../AsthmaSeq-EOS-paper/figures/healthyEOS_ITGAM_violin.svg")

ceacam8_violin <- VlnPlot(obj, features = "CEACAM8", cols = c("#7CAE00", "#00BFC4","#ED68ED", "#E68613", "#C77CFF","#0CB702", "#00A9FF", "#CD9600", "#F8766D", "#00BE67",  "#00B8E7")) + theme(legend.position = "none")
ggsave(ceacam8_violin, filename = "../AsthmaSeq-EOS-paper/figures/healthyEOS_CEACAM8_violin.svg")

#############
# Figure 1E #
#############
markers_entrez <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = rownames(markers), columns = c("ENTREZID"), keytype = "GENENAME")

library(clusterProfiler)
ggo <- enrichGO(gene = markers_entrez$ENTREZID,
                OrgDb = org.Hs.eg.db,
                pvalueCutoff = 0.05,
                ont = "BP", 
                readable = T,
                pAdjustMethod = "bonferroni")
df.go <- ggo@result
df.go <- df.go[df.go$p.adjust < 0.05,]

library(ReactomePA)
reac <- enrichPathway(gene=markers_entrez$ENTREZID, pvalueCutoff = 0.05, readable=T, pAdjustMethod = "bonferroni")
df.reac <- reac@result
df.reac <- df.reac[df.reac$p.adjust < 0.05,]

write.table(df.reac, file = "EOS_Reactome_pathways.tsv", sep = "\t", quote = F)

degran.reac <- df.reac[grepl("degranulation", df.reac$Description),]
degran.genes <- unique(unlist(strsplit(as.character(degran.reac$geneID), "/")))

markers1 <- markers[rownames(markers) %in% degran.genes,]
cat(rownames(markers1), file = "degran_genes.txt", sep = "\n")

avg <- AverageExpression(obj, features = rownames(markers1), return.seurat=T)
heat <- DoHeatmap(avg, features = rownames(markers1), draw.lines = F) + theme(axis.text.y = element_text(size = 0)) + theme(legend.position="none")

ggsave(heat, filename = '../AsthmaSeq-EOS-paper/REVISION/updated_figures/degran_heat.svg')

neut_mark <- FindMarkers(obj, ident.1 = "neutrophil", only.pos = T, logfc.threshold = 0.25, group.by = "annot")
neut_mark <- neut_mark[neut_mark$p_val_adj < 0.05,]
markers_entrez <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = rownames(neut_mark), columns = c("ENTREZID"), keytype = "GENENAME")
reac <- enrichPathway(gene=markers_entrez$ENTREZID, pvalueCutoff = 0.05, readable=T)
reac_bar <- barplot(reac, showCategory = 10) + ggtitle("Top 10 Reactome pathways") + theme_classic()
reac.df <- reac@result
degran.reac2 <- reac.df[grepl("degranulation", reac.df$Description),]
degran.genes2 <- unique(unlist(strsplit(as.character(degran.reac2$geneID), "/")))
intersect(degran.genes2, degran.genes)

markers[rownames(markers) %in% degran.genes,]
markers2 <- neut_mark[rownames(neut_mark) %in% degran.genes2,]


##########################################
########### healthy eos only #############
##########################################

# read in EOS object
obj <- readRDS("EOS_integrated_updated.rds")

healthy_eos <- subset(obj, severity == "healthy")

healthy_eos[["RNA"]] <- split(healthy_eos[["RNA"]], f = healthy_eos$patient)

healthy_eos <- NormalizeData(healthy_eos)
healthy_eos <- FindVariableFeatures(healthy_eos)
healthy_eos <- ScaleData(healthy_eos)
healthy_eos <- RunPCA(healthy_eos)

ElbowPlot(healthy_eos, ndims = 25) #ndims 15

healthy_eos <- FindNeighbors(healthy_eos, dims = 1:15)
healthy_eos <- FindClusters(healthy_eos, resolution = 0.3)
healthy_eos <- RunUMAP(healthy_eos, dims = 1:15)

DimPlot(healthy_eos)
DimPlot(healthy_eos, group.by = "patient")

healthy_eos <- IntegrateLayers(object = healthy_eos, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca")
healthy_eos[["RNA"]] <- JoinLayers(healthy_eos[["RNA"]])

healthy_eos <- FindNeighbors(healthy_eos, reduction = "integrated.cca", dims = 1:15)
healthy_eos <- FindClusters(healthy_eos, resolution = 0.3)
healthy_eos <- RunUMAP(healthy_eos, dims = 1:15, reduction = "integrated.cca")

DimPlot(healthy_eos)
DimPlot(healthy_eos, group.by = "patient")

FeaturePlot(healthy_eos, features = "percent.mt")
FeaturePlot(healthy_eos, features = "nFeature_RNA")
FeaturePlot(healthy_eos, features = "nCount_RNA")

# get cell IDs for cluster 3
healthy_inflammatory <- WhichCells(healthy_eos, idents = '3')

# umap for figure 3_
healthy_inflamHighlight <- DimPlot(healthy_eos, cells.highlight = healthy_inflammatory, pt.size=1) + theme_void() + theme(legend.position="none") 
ggsave(healthy_inflamHighlight, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/Healthy_InflamHighlight.png")

clusterMarkers <- FindAllMarkers(healthy_eos, logfc.threshold = 0.25, min.pct = 0.1, only.pos=T)
clusterMarkers <- clusterMarkers[clusterMarkers$p_val_adj < 0.05,]

cluster0 <- clusterMarkers[clusterMarkers$cluster == '0',] %>% top_n(10, wt = avg_log2FC)
cluster1 <- clusterMarkers[clusterMarkers$cluster == '1',] %>% top_n(10, wt = avg_log2FC)
cluster2 <- clusterMarkers[clusterMarkers$cluster == '2',] %>% top_n(10, wt = avg_log2FC)
cluster3 <- clusterMarkers[clusterMarkers$cluster == '3',] %>% top_n(10, wt = avg_log2FC)

####################
#### FIGURE 2A ####
###################
umap <- DimPlot(healthy_eos, pt.size=0.7) + ggtitle("") + theme_void() + theme(legend.position = "none")
ggsave(umap, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthy_UMAP.png")

umap2 <- DimPlot(healthy_eos, group.by = "patient", shuffle=T, cols = c("palegreen3", "pink1", "violet" , "lightblue"), pt.size=1) + theme_void() + theme(legend.position='none') + ggtitle("")
ggsave(umap2, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthy_UMAP2.png")

###################
#### FIGURE 2B ####
##################
clusters <- FindAllMarkers(healthy_eos, logfc.threshold = 0.25, min.pct = 0.1, only.pos = T)
clusters <- clusters[clusters$p_val_adj < 0.05,]
write.table(clusters, file = "healthyEOS_clusterMarkers.tsv", sep = "\t", quote=F)

healthy_eos <- ScaleData(healthy_eos, features = rownames(clusters), assay = "RNA")

top20 <- clusters %>%
  group_by(cluster) %>%
  top_n(20, avg_log2FC)

avg <- AverageExpression(healthy_eos, features = top20$gene, return.seurat=T, group.by = "RNA_snn_res.0.3")
heat <- DoHeatmap(avg, features = top20$gene, group.by = "RNA_snn_res.0.3", draw.lines = F, group.bar = F) + theme(legend.position="none")
ggsave(heat, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthyEOS_heat.png")

###################
#### FIGURE 2C ####
##################
# Pathway analysis
up_c0 <- clusters[clusters$cluster == "0",]
up_c1 <- clusters[clusters$cluster == "1",]
up_c2 <- clusters[clusters$cluster == "2",]
up_c3 <- clusters[clusters$cluster == "3",]

# Return the Entrez ID for the set of gene symbols
c0_entrez <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                   keys = rownames(up_c0),
                                   columns = "ENTREZID",
                                   keytype = "GENENAME")

c1_entrez <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                   keys = rownames(up_c1),
                                   columns = "ENTREZID",
                                   keytype = "GENENAME")

c2_entrez <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                   keys = rownames(up_c2),
                                   columns = "ENTREZID",
                                   keytype = "GENENAME")

c3_entrez <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                   keys = rownames(up_c3),
                                   columns = "ENTREZID",
                                   keytype = "GENENAME")


# Compare pathways
list <- list(c0_entrez$ENTREZID, c1_entrez$ENTREZID, c2_entrez$ENTREZID, c3_entrez$ENTREZID)
names(list) <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")

library(clusterProfiler)
compGO <- compareCluster(geneCluster = list,
                         fun = "enrichGO",
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP", 
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "bonferroni", readable = T)
goBp_dot <- dotplot(compGO, title = "GO Enrichment Analysis", by = "count")
ggsave(goBp_dot, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthyEOS_pathwayDot.svg")
goBp_DF <- compGO@compareClusterResult
goBp_DF <- goBp_DF[goBp_DF$p.adjust < 0.05,]
write.table(goBp_DF, file = "healthyEOS_Clusters_GOBP.tsv", sep = "\t", quote=F)

###################
#### FIGURE 2D ####
##################
vln <- VlnPlot(healthy_eos, features = c("RNASE2", "RNASE3", 'MBP', "CLC", "HLA-A", "HLA-B", "HLA-C"), flip=T, stack=T) + theme(legend.position="none")
ggsave(vln, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthy_effectorViolin.svg")

rnase2 <- VlnPlot(healthy_eos, features = "RNASE2") + theme(legend.position="none")
rstatix::kruskal_test(RNASE2~ident, data = rnase2$data)
ggsave(rnase2, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthy_rnase2.svg")

rnase3 <- VlnPlot(healthy_eos, features = "RNASE3") + theme(legend.position="none")
rstatix::kruskal_test(RNASE3~ident, data = rnase3$data)
ggsave(rnase3, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthy_rnase3.svg")

hlaa <- VlnPlot(healthy_eos, features = "HLA-A") + theme(legend.position="none")
rstatix::kruskal_test(`HLA-A`~ident, data = hlaa$data)
ggsave(hlaa, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthy_HLAA.svg")

hlab <- VlnPlot(healthy_eos, features = "HLA-B") + theme(legend.position="none")
rstatix::kruskal_test(`HLA-B`~ident, data = hlab$data)

hlac <- VlnPlot(healthy_eos, features = "HLA-C") + theme(legend.position="none")
rstatix::kruskal_test(`HLA-C`~ident, data = hlac$data)
ggsave(hlac, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthy_HLAC.svg")

mbp <- VlnPlot(healthy_eos, features = "MBP") + theme(legend.position="none")
rstatix::kruskal_test(MBP~ident, data = mbp$data)

clc <- VlnPlot(healthy_eos, features = "CLC") + theme(legend.position="none")
rstatix::kruskal_test(CLC~ident, data = clc$data)

###################
#### FIGURE 2E ###
##################
isg <- DotPlot(healthy_eos, features = c("IFI6", "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", 'MX2', "EPSTI1", "IFITM2"), dot.scale = 8, cols = c("grey", "red")) + coord_flip()
ggsave(isg, filename = "../AsthmaSeq-EOS-paper/REVISION/updated_figures/healthy_ISG.svg")