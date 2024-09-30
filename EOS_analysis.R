library(Seurat)
library(dplyr)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(DESeq2)
library(tibble)
library(msigdbr)
library(annotate)
library(clusterProfiler)

set.seed(123)
setwd("/Users/kyndalgoss/Desktop/MELANIE/")

# read patient A01 data
a01.data <- Read10X("A01_EOS/filtered_feature_bc_matrix/")
# create seurat object
a01 <- CreateSeuratObject(counts = a01.data, project = "A01_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("A01_EOS/keep.csv")
a01 <- subset(a01, cells = keep$Barcode)

# quality control
a01[["percent.mt"]] <- PercentageFeatureSet(a01, pattern = "^MT-", assay = "RNA")
VlnPlot(a01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(a01, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
a01 <- subset(a01, subset = nFeature_RNA > 100 & nFeature_RNA < 1500 & percent.mt < 20 & nCount_RNA > 100 & nCount_RNA < 2000)

# read patient A02 data
a02.data <- Read10X("A02_EOS/filtered_feature_bc_matrix/")
# create seurat object
a02 <- CreateSeuratObject(counts = a02.data, project = "A02_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("A02_EOS/keep.csv")
a02 <- subset(a02, cells = keep$Barcode)
# quality control
a02[["percent.mt"]] <- PercentageFeatureSet(a02, pattern = "^MT-")
VlnPlot(a02, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(a02, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
a02 <- subset(a02, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 15 & nCount_RNA > 100 & nCount_RNA < 4000)

# laod patient A03 data
a03.data <- Read10X("A03_EOS/filtered_feature_bc_matrix/")
# create seurat object
a03 <- CreateSeuratObject(counts = a03.data, project = "A03_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("A03_EOS/keep.csv")
a03 <- subset(a03, cells = keep$Barcode)
# quality control
a03[["percent.mt"]] <- PercentageFeatureSet(a03, pattern = "^MT-")
VlnPlot(a03, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(a03, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
a03 <- subset(a03, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 15 & nCount_RNA > 100 & nCount_RNA < 5000)

#  laod patient A04 data
a04.data <- Read10X("A04_EOS/filtered_feature_bc_matrix/")
# create seurat object
a04 <- CreateSeuratObject(counts = a04.data, project = "A04_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("A04_EOS/keep.csv")
a04 <- subset(a04, cells = keep$Barcode)
# quality control
a04[["percent.mt"]] <- PercentageFeatureSet(a04, pattern = "^MT-")
VlnPlot(a04, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(a04, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
a04 <- subset(a04, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 15 & nCount_RNA > 100 & nCount_RNA < 4000)

#  laod patient A05 data
a05.data <- Read10X("A05_EOS/filtered_feature_bc_matrix/")
# create seurat object
a05 <- CreateSeuratObject(counts = a05.data, project = "A05_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("A05_EOS/keep.csv")
a05 <- subset(a05, cells = keep$Barcode)
# quality control
a05[["percent.mt"]] <- PercentageFeatureSet(a05, pattern = "^MT-")
VlnPlot(a05, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(a05, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
a05 <- subset(a05, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 20 & nCount_RNA > 100 & nCount_RNA < 3000)

#  laod patient A06 data
a06.data <- Read10X("A06_EOS/filtered_feature_bc_matrix/")
# create seurat object
a06 <- CreateSeuratObject(counts = a06.data, project = "A06_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("A06_EOS/keep.csv")
a06 <- subset(a06, cells = keep$Barcode)
# quality control
a06[["percent.mt"]] <- PercentageFeatureSet(a06, pattern = "^MT-")
VlnPlot(a06, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(a06, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
a06 <- subset(a06, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 20 & nCount_RNA > 100 & nCount_RNA < 2000)

# read patient b01 data
b01.data <- Read10X("B01_EOS/filtered_feature_bc_matrix/")
# create seurat object
b01 <- CreateSeuratObject(counts = b01.data, project = "B01_EOS", min.cells = 3, min.features = 0)
# quality control
b01[["percent.mt"]] <- PercentageFeatureSet(b01, pattern = "^MT-")
VlnPlot(b01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
b01 <- subset(b01, subset = nFeature_RNA > 200 & percent.mt < 15 & nCount_RNA > 100 & nCount_RNA < 2000)

# read patient b02 data
b02.data <- Read10X("B02_EOS/filtered_feature_bc_matrix/")
# create seurat object
b02 <- CreateSeuratObject(counts = b02.data, project = "B02_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("B02_EOS/keep.csv")
b02 <- subset(b02, cells = keep$Barcode)
# quality control
b02[["percent.mt"]] <- PercentageFeatureSet(b02, pattern = "^MT-")
VlnPlot(b02, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
b02 <- subset(b02, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 15 & nCount_RNA > 100 & nCount_RNA < 2000)

# load patient B03 data
b03.data <- Read10X("B03_EOS/filtered_feature_bc_matrix/")
# create seurat object
b03 <- CreateSeuratObject(counts = b03.data, project = "B03_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("B03_EOS/keep.csv")
b03 <- subset(b03, cells = keep$Barcode)
# quality control
b03[["percent.mt"]] <- PercentageFeatureSet(b03, pattern = "^MT-")
VlnPlot(b03, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
b03 <- subset(b03, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 15 & nCount_RNA > 100 & nCount_RNA < 4000)

# load patient B05 data and create Seurat object
b05.data <- Read10X("B05_EOS/filtered_feature_bc_matrix/")
b05 <- CreateSeuratObject(counts = b05.data, project = "B05_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("B05_EOS/B05_EOS_keep.csv")
b05 <- subset(b05, cells = keep$Barcode)
# quality control
b05[["percent.mt"]] <- PercentageFeatureSet(b05, pattern = "^MT-")
VlnPlot(b05, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
b05 <- subset(b05, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA > 100 & nCount_RNA < 3000)

# load patient C01 data and create Seurat object
c01.data <- Read10X("C01_EOS/filtered_feature_bc_matrix/")
c01 <- CreateSeuratObject(counts=c01.data, project = 'C01_EOS', min.cells = 3, min.features = 0)
keep <- read.csv("C01_EOS/keep.csv")
c01 <- subset(c01, cells = keep$Barcode)
c01[["percent.mt"]] <- PercentageFeatureSet(c01, pattern = "^MT-")
VlnPlot(c01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
c01 <- subset(c01, subset = nFeature_RNA > 100 & nFeature_RNA < 1200 & percent.mt < 15 & nCount_RNA > 100 & nCount_RNA < 2000)

# load patient C02 data and create Seurat object
c02.data <- Read10X("C02_EOS/filtered_feature_bc_matrix/")
c02 <- CreateSeuratObject(counts = c02.data, project = "C02_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("C02_EOS/keep.csv")
c02 <- subset(c02, cells = keep$Barcode)
c02[["percent.mt"]] <- PercentageFeatureSet(c02, pattern = "^MT-")
VlnPlot(c02, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
c02 <- subset(c02, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 15 & nCount_RNA > 100 & nCount_RNA < 15000)

# load patient C03 data and create Seurat object
c03.data <- Read10X("C03_EOS/filtered_feature_bc_matrix/")
c03 <- CreateSeuratObject(counts = c03.data, project = "C03_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("C03_EOS/keep.csv")
c03 <- subset(c03, cells = keep$Barcode)
c03[["percent.mt"]] <- PercentageFeatureSet(c03, pattern = "^MT-")
VlnPlot(c03, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
c03 <- subset(c03, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 15 & nCount_RNA > 100 & nCount_RNA < 2000)

# load patient C04 data and create Seurat object
c04.data <- Read10X("C04_EOS/filtered_feature_bc_matrix/")
c04 <- CreateSeuratObject(counts = c04.data, project = "C04_EOS", min.cells = 3, min.features = 0)
keep <- read.csv("C04_EOS/keep.csv")
c04 <- subset(c04, cells = keep$Barcode)
c04[["percent.mt"]] <- PercentageFeatureSet(c04, pattern = "^MT-")
VlnPlot(c04, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
c04 <- subset(c04, subset = nFeature_RNA > 100 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA > 100 & nCount_RNA < 5000)


#############################################################################################################
##################### COMBINED ANALYSIS ###########################
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

#a01$batch <- "3"
#a02$batch <- "2"
#a03$batch <- "3"
#a04$batch <- "2"
#a05$batch <- "2"
#a06$batch <- "2"
#b01$batch <- "1"
#b02$batch <- "3"
#b03$batch <- "3"
#b05$batch <- "4"

# analyze combined datasets
obj <- merge(a01, y = c(a02, a03, a04, a05, a06, b01, b02, b03, b05, c01, c02, c03, c04), add.cell.ids = c("a01", "a02", "a03", "a04", "a05", "a06", "b01", "b02", "b03", "b05", "c01", "c02", "c03", "c04"), project = "EOS")

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$patient)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

ElbowPlot(obj, ndims = 25) #ndims 20

obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.01)
obj <- RunUMAP(obj, dims = 1:20)

DimPlot(obj)
DimPlot(obj, group.by = "patient")
DimPlot(obj, group.by = "severity")
#DimPlot(obj, group.by = "batch")

obj <- IntegrateLayers(object = obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca")

obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:20)
obj <- FindClusters(obj, resolution = 0.1)
obj <- RunUMAP(obj, dims = 1:20, reduction = "integrated.cca")

DimPlot(obj)
DimPlot(obj, group.by = "patient")
DimPlot(obj, group.by = "severity")
DimPlot(obj, group.by = "batch")
DimPlot(obj, split.by = "batch")
DimPlot(obj, split.by = "severity")
dittoSeq::dittoBarPlot(obj, var = "severity", group.by = "seurat_clusters")

# re-join layers after integration
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
#saveRDS(obj, file = "EOS_integrated.rds")

# read in integrated EOS object
obj <- readRDS("EOS_integrated.rds")

obj$severity1 <- ifelse(obj$severity == "healthy", yes = "healthy", no = "asthma")


ISGs <- list(c("ADAR", "APOBEC3", "BST2", "CD74","DDIT4", "DDX58", "DDX60", "EIF2AK2", "GPB1", "GPB2", "HPSE", "IFI44L", "IFI6", "IFIH1", "IFIT1", "ITIT2", "IFIT3", "IFIT5", "IFITM1", "IFITM2", "IFITM3", "IRF1", "IRF7", "ISG15", "ISG20", "MAP3K14", "MOV10", "MS4A4A", "MOV10", "MX1", "MX2", "NAMPT", "NT5C3", "OAS1", "OAS2", "OAS3", "OASL", "P2RY6", "PHF15", "PML", "RSAD2", "RTP4", "SLC15A3", "SLC25A28", "SSBP3", "TRIM5", "SUN2", "ZC3HAV1"))
# ISG list obtained from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3274382/

obj <- AddModuleScore(obj, features = ISGs, name = "ISG")
a <- VlnPlot(obj, features = "ISG1", group.by = "severity1") + ggtitle("Module score - ISGs")
b <- ks.test(ISG1~ident, data = a$data)

c <- VlnPlot(obj, features = "ISG1", group.by = "severity") + ggtitle("Module Score - ISGs")
d <- kruskal_test(ISG1~ident, data = c$data)
FSA::dunnTest(ISG1~ident, data = c$data, method="bonferroni")


umap1 <- DimPlot(obj) + ggtitle("") + theme_void() + theme(legend.position = "none")
ggsave(umap1, filename = "../AsthmaSeq-EOS-paper/figures/healthyVSasthma_umap1.svg")

umap2 <- DimPlot(obj, group.by = "severity1", cols = c("#008ECE", "red")) + ggtitle("") + theme_void() + theme(legend.position = "none")
ggsave(umap2, filename = "../AsthmaSeq-EOS-paper/figures/healthyVSasthma_umap2.svg")

DimPlot(obj, split.by = "severity", group.by = "RNA_snn_res.0.1")
DimPlot(obj, group.by = "severity1", split.by = "RNA_snn_res.0.1")
z <- dittoSeq::dittoBarPlot(obj, var = "RNA_snn_res.0.1", group.by = "severity1", color.panel = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"), xlab="condition", ylab="Fraction of Cells", legend.title="Cluster", main = NULL)
ggsave(z, filename = "../AsthmaSeq-EOS-paper/figures/healthyVasthma_clusters_bar.svg")

# look at canonical EOS markers
FeaturePlot(obj, features = "CCR3", order=T) 
FeaturePlot(obj, features = "SIGLEC8", order=T) 
FeaturePlot(obj, features = "IL5RA", order=T) 
FeaturePlot(obj, features = "EPX", order=T) 
FeaturePlot(obj, features = "PRG2", order=T)

# platelet marker
FeaturePlot(obj, features = "PPBP", order=T) 
FeaturePlot(obj, features = "PF4", order=T) 
FeaturePlot(obj, features = "TUBB1", order=T)

# Cluster markers
Idents(obj) <- obj$RNA_snn_res.0.1
markers <- FindAllMarkers(obj, only.pos = T, min.pct = 0.1)
markers <- markers[markers$p_val_adj < 0.05,]

cluster2_vs_cluster3 <- FindMarkers(obj, ident.1="2", ident.2="3", logfc.threshold = 0.25)
cluster2_vs_cluster3 <- cluster2_vs_cluster3[cluster2_vs_cluster3$p_val_adj < 0.05,]
write.table(cluster2_vs_cluster3, file="Cluster2_vs_Cluster3.tsv", quote = F, sep = "\t", row.names = T, col.names = T)


keyvals <- ifelse(
  cluster2_vs_cluster3$avg_log2FC > 0, '#00BFC4',
  ifelse(cluster2_vs_cluster3$avg_log2FC < 0, '#C77CFF',
         'black'))

names(keyvals)[keyvals == '#00BFC4'] <- 'Cluster 2 (Asthma)'
names(keyvals)[keyvals == '#C77CFF'] <- 'Cluster 3 (Healthy)'

volcano <- EnhancedVolcano::EnhancedVolcano(cluster2_vs_cluster3, lab = rownames(cluster2_vs_cluster3), labSize = 0, x='avg_log2FC',y = 'p_val_adj',FCcutoff = 0.25,pCutoff = 0.05, title = 'Cluster 2 vs Cluster 3', subtitle = "", legendPosition = 'bottom', colCustom = keyvals, gridlines.major = F, gridlines.minor = F)
ggsave(volcano, filename = "../AsthmaSeq-EOS-paper/figures/cluster2_vs_cluster3_volcano.svg")

cluster2_up <- cluster2_vs_cluster3[cluster2_vs_cluster3$avg_log2FC > 0,]
cluster3_up <- cluster2_vs_cluster3[cluster2_vs_cluster3$avg_log2FC < 0,]

top100_c2 <- cluster2_up %>%
  top_n(100, avg_log2FC)

top100_c3 <- cluster3_up %>%
  top_n(100, avg_log2FC)

obj <- ScaleData(obj, features = rownames(obj))

# pathway analysis cluster 2 vs cluster 3
cluster2_entrez <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                        keys = rownames(cluster2_up),
                                        columns = "ENTREZID",
                                        keytype = "GENENAME")

cluster3_entrez <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                       keys = rownames(cluster3_up),
                                       columns = "ENTREZID",
                                       keytype = "GENENAME")


# GSEA on cluster 2 vs cluster 3
cluster2_vs_cluster3 <- FindMarkers(obj, ident.1="2", ident.2="3", logfc.threshold = 0.25)

cluster2_vs_cluster3 <- cluster2_vs_cluster3[order(-cluster2_vs_cluster3$avg_log2FC),]

gene_list <- cluster2_vs_cluster3$avg_log2FC
names(gene_list) <- rownames(cluster2_vs_cluster3)
gene_list

m_df <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

m_df$symbol <- getSYMBOL(as.character(m_df$entrez_gene), data = 'org.Hs.eg')

H_t2g <- m_df %>% dplyr::select(gs_name, symbol)

c2 <- msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::select(gs_name, gene_symbol)

em2 <- GSEA(gene_list, TERM2GENE = H_t2g, eps = 1e-300)
df <- as.data.frame(em2@result)

df$direction <- ifelse(df$NES > 0, yes = "cluster 2 (asthma)", no = "cluster 3 (healthy)")
df$Description <- factor(df$Description, levels = df$Description[order(df$NES)])
gsea_bar <- ggplot(df, aes(x = Description, y = NES, fill = direction)) + geom_bar(stat="identity") + coord_flip() + theme_classic() + scale_fill_manual("legend", values = c("cluster 2 (asthma)" = '#00BFC4', "cluster 3 (healthy)" = '#C77CFF'))
ggsave(gsea_bar, filename = "../AsthmaSeq-EOS-paper/figures/c2_v_c3_GSEAbar.svg")

enrichplot::gseaplot2(em2, geneSetID = 1, title = "Hallmark Interferon Alpha response")





############################################
# Find DEGs between healthy vs asthma EOS #
Idents(obj) <- obj$severity1

#a <- FindAllMarkers(obj, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25)

b <- FindMarkers(obj, ident.1 = "healthy", ident.2 = "asthma", group.by = "severity1", logfc.threshold = 0.25, min.pct = 0.1)

# GSEA
b <- b[order(-b$avg_log2FC),]

gene_list <- b$avg_log2FC
names(gene_list) <- rownames(b)
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

d <- enrichplot::gseaplot2(em2, geneSetID = 1, title = "HALLMARK INTERFERON GAMMA RESPONSE")
ggsave(d, filename = "../AsthmaSeq-EOS-paper/figures/IFNG_response_line.svg")

e <- enrichplot::gseaplot2(em2, geneSetID = 2, title = "HALLMARK INTERFERON ALPHA RESPONSE")
ggsave(e, filename = "../AsthmaSeq-EOS-paper/figures/IFNA_response_line.svg")

f <- enrichplot::gseaplot2(em2, geneSetID = 3, title = "HALLMARK INFLAMMATORY RESPONSE")
ggsave(f, filename = "../AsthmaSeq-EOS-paper/figures/InfResponse_line.svg")

enrichplot::gseaplot2(em2, geneSetID = c(1,2,3))

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

genes <- b[rownames(b) %in% IFNA_response_genes | rownames(b) %in% IFNG_response_genes | rownames(b) %in% InflammatoryResponseGenes,]
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
ggsave(a, filename = "../AsthmaSeq-EOS-paper/figures/healthyVasthma_Cord.svg")


b <- b[b$p_val_adj < 0.05,]
write.table(b, file="Healthy_vs_Asthma_EOS_DEG.tsv", quote = F, sep = "\t", row.names = T, col.names = T)

keyvals <- ifelse(
  b$avg_log2FC < 0, 'red',
  ifelse(b$avg_log2FC > 0, '#008ECE',
         'black'))

names(keyvals)[keyvals == 'red'] <- 'Asthma'
names(keyvals)[keyvals == '#008ECE'] <- 'Healthy'

volcano <- EnhancedVolcano::EnhancedVolcano(b, lab = rownames(b), x='avg_log2FC',y = 'p_val_adj',FCcutoff = 0.25,pCutoff = 0.05, title = 'Healthy vs Asthma', subtitle = "", legendPosition = 'bottom', colCustom = keyvals, gridlines.major = F, gridlines.minor = F)
ggsave(volcano, filename = "../AsthmaSeq-EOS-paper/figures/healthyVasthma_volcano.svg")

up_healthy <- b[b$avg_log2FC > 0,]
up_asthma <- b[b$avg_log2FC < 0,]


# Return the Entrez ID for the set of gene symbols
Healthy_entrez <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                        keys = rownames(up_healthy),
                                        columns = "ENTREZID",
                                        keytype = "GENENAME")

Asthma_entrez <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                     keys = rownames(up_asthma),
                                     columns = "ENTREZID",
                                     keytype = "GENENAME")


# Compare pathways
list <- list(Healthy_entrez$ENTREZID, Asthma_entrez$ENTREZID)
names(list) <- c("Healthy", "Asthma")

library(clusterProfiler)
compGO <- compareCluster(geneCluster = list,
                         fun = "enrichGO",
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP", 
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
dotplot(compGO, title = "GO Enrichment Analysis", by = "count", showCategory = 20)
#write.table(compGO@compareClusterResult, file = "EOS_compareGO_pathways.tsv", sep = "\t")

compReactome <- compareCluster(geneCluster = list,
                               fun = "enrichPathway", 
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH")
df <- compReactome@compareClusterResult
dotplot(compReactome, title = "Reactome Pathway Enrichment Analysis", by = "count")
#write.table(compReactome@compareClusterResult, file = "EOS_compareReactome_pathways.tsv", sep = "\t")



clc <- VlnPlot(obj, "CLC", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(CLC~ident, data = clc$data)
ggsave(clc, filename = "healthyVasthma_CLC.svg")

rnase2 <- VlnPlot(obj, "RNASE2", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(RNASE2~ident, data = rnase2$data)
ggsave(rnase2, filename = "healthyVasthma_RNASE2.svg")

rnase3 <- VlnPlot(obj, "RNASE3", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(RNASE3~ident, data = rnase3$data)
ggsave(rnase3, filename = "healthyVasthma_RNASE3.svg")

ccr3 <- VlnPlot(obj, "CCR3", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(CCR3~ident, data = ccr3$data)
ggsave(ccr3, filename = "healthyVasthma_CCR3.svg")

b2m <- VlnPlot(obj, "B2M", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(B2M~ident, data = b2m$data)
ggsave(b2m, filename = "healthyVasthma_B2M.svg")

tapbp <- VlnPlot(obj, "TAPBP", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(TAPBP~ident, data = tapbp$data)
ggsave(tapbp, filename = "healthyVasthma_TAPBP.svg")

cd74 <- VlnPlot(obj, "CD74", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(CD74~ident, data = cd74$data)
ggsave(cd74, filename = "healthyVasthma_CD74.svg")

hla_a <- VlnPlot(obj, "HLA-A", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(`HLA-A`~ident, data = hla_a$data)
ggsave(hla_a, filename = "healthyVasthma_HLAA.svg")

hla_b <- VlnPlot(obj, "HLA-B", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(`HLA-B`~ident, data = hla_b$data)
ggsave(hla_b, filename = "healthyVasthma_HLAB.svg")

hla_c <- VlnPlot(obj, "HLA-C", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(`HLA-C`~ident, data = hla_c$data)
ggsave(hla_c, filename = "healthyVasthma_HLAC.svg")

hla_drb1 <- VlnPlot(obj, "HLA-DRA", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(`HLA-DRA`~ident, data = hla_drb1$data)

hla_dqb1 <- VlnPlot(obj, "HLA-DQB1", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(`HLA-DQB1`~ident, data = hla_dqb1$data)

hla_dpb1 <- VlnPlot(obj, "HLA-DPB1", cols = c("red", "#008ECE"), group.by = "severity1", pt.size = 0.01) + theme(legend.position = "none")
ks.test(`HLA-DPB1`~ident, data = hla_dpb1$data)

#########################################################
# mild vs severe only
sub <- subset(obj, severity == "healthy", invert = T)

sub <- NormalizeData(sub)

sub <- FindNeighbors(sub, reduction = "integrated.cca", dims = 1:20)
sub <- FindClusters(sub, resolution = 0.1)
sub <- RunUMAP(sub, dims = 1:20, reduction = "integrated.cca")

g <- DimPlot(sub) + ggtitle("") + theme_void()
ggsave(g, filename = "../AsthmaSeq-EOS-paper/figures/mildVsevere_umap1.svg")

h <- DimPlot(sub, group.by = "severity", cols = c("#6B00A0", "orange")) + ggtitle("") + theme_void()
ggsave(h, filename = "../AsthmaSeq-EOS-paper/figures/mildVsevere_umap2.svg")

i <- dittoSeq::dittoBarPlot(sub, var = "RNA_snn_res.0.1", group.by = "severity", color.panel = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"), xlab="condition", ylab="Fraction of Cells", legend.title="Cluster", main = NULL)
ggsave(i, filename = "../AsthmaSeq-EOS-paper/figures/mildVsevere_clusterBar.svg")

j <- dittoSeq::dittoBarPlot(sub, group.by = "RNA_snn_res.0.1", var = "severity")
i$data

c <- FindMarkers(sub, ident.1="mild", ident.2="severe", group.by = "severity", min.pct = 0.1, logfc.threshold = 0.25)

# GSEA
c <- c[order(-c$avg_log2FC),]

gene_list <- c$avg_log2FC
names(gene_list) <- rownames(c)
gene_list

m_df <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

m_df$symbol <- getSYMBOL(as.character(m_df$entrez_gene), data = 'org.Hs.eg')

H_t2g <- m_df %>% dplyr::select(gs_name, symbol)

em2 <- GSEA(gene_list, TERM2GENE = H_t2g, eps = 1e-300, pvalueCutoff = 0.05, pAdjustMethod = "bonferroni")
df <- as.data.frame(em2@result)

ifna <- enrichplot::gseaplot2(em2, geneSetID = 1, title = "Hallmark IFNA response")
ggsave(ifna, filename = "../AsthmaSeq-EOS-paper/figures/mildVsevere_IFNA_line.svg")

ifng <- enrichplot::gseaplot2(em2, geneSetID = 2, title = "Hallmark IFNG response")
ggsave(ifng, filename = "../AsthmaSeq-EOS-paper/figures/mildVsevere_IFNG_line.svg")

df$direction <- ifelse(df$NES > 0, yes = "mild", no = "severe")
df$Description <- factor(df$Description, levels = df$Description[order(df$NES)])
bar <- ggplot(df, aes(x = Description, y = NES, fill = direction)) + geom_bar(stat="identity") + coord_flip() + theme_classic() + scale_fill_manual("legend", values = c("mild" = "#6B00A0", "severe" = "orange"))
ggsave(bar, filename = "../AsthmaSeq-EOS-paper/figures/mildVsevere_GSEAbar.svg")

c <- c[c$p_val_adj < 0.05,]
write.table(c, file="Mild_vs_Severe_EOS_DEG.tsv", quote = F, sep = "\t", row.names = T, col.names = T)


keyvals <- ifelse(
  c$avg_log2FC < 0, 'orange',
  ifelse(c$avg_log2FC > 0, '#6B00A0',
         'black'))

names(keyvals)[keyvals == 'orange'] <- 'Severe'
names(keyvals)[keyvals == '#6B00A0'] <- 'Mild'

j <- EnhancedVolcano::EnhancedVolcano(c, lab = rownames(c), x='avg_log2FC',y = 'p_val_adj',FCcutoff = 0.25,pCutoff = 0.05, title = 'Mild vs Severe', subtitle = "", legendPosition = 'bottom', colCustom = keyvals, gridlines.major = F, gridlines.minor = F)
ggsave(j, filename = "../AsthmaSeq-EOS-paper/figures/mildVsevere_volcano.svg")

EnhancedVolcano::EnhancedVolcano(c, lab = rownames(c), labSize = 0.5, x='avg_log2FC',y = 'p_val_adj',FCcutoff = 0.25,pCutoff = 0.05, title = 'Mild vs Severe', subtitle = "", legendPosition = 'bottom')

mild <- c[c$avg_log2FC > 0,]
severe <- c[c$avg_log2FC < 0,]

IFNA_response <- df[df$Description == "HALLMARK_INTERFERON_ALPHA_RESPONSE",]$core_enrichment
IFNA_response_genes <- unlist(strsplit(IFNA_response, split = "/"))

IFNG_response <- df[df$Description == "HALLMARK_INTERFERON_GAMMA_RESPONSE",]$core_enrichment
IFNG_response_genes <- unlist(strsplit(IFNG_response, split = "/"))

sub <- AddModuleScore(sub, features = list(IFNA_response_genes), name = "IFNA_response")
sub <- AddModuleScore(sub, features = list(IFNG_response_genes), name = "IFNG_response")

ifna_feat <- FeaturePlot(sub, features = "IFNA_response1", order=T, pt.size = 0.5) + scale_color_gradientn(colours = rev(brewer.pal(n=11, name="RdBu"))) + theme_void() + ggtitle("")
ggsave(ifna_feat, filename = "../AsthmaSeq-EOS-paper/figures/IFNA_feature.svg")

ifng_feat <- FeaturePlot(sub, features = "IFNG_response1", order=T, pt.size = 0.5) + scale_color_gradientn(colours = rev(brewer.pal(n=11, name="RdBu"))) + theme_void() + ggtitle("")
ggsave(ifng_feat, filename = "../AsthmaSeq-EOS-paper/figures/IFNG_feature.svg")


IFI6 <- FeaturePlot(sub, features = "IFI6", order=T, pt.size = 0.5) + scale_color_viridis(option = "magma") + ggtitle("") + theme_void()
ggsave(IFI6, filename = "../AsthmaSeq-EOS-paper/figures/IFI6_feature.svg")

MX1 <- FeaturePlot(sub, features = "MX1", order=T, pt.size = 0.5) + scale_color_viridis(option = "magma") + ggtitle("") + theme_void()
ggsave(MX1, filename = "../AsthmaSeq-EOS-paper/figures/MX1_feature.svg")

ISG15 <- FeaturePlot(sub, features = "ISG15", order=T, pt.size = 0.5) + scale_color_viridis(option = "magma") + ggtitle("") + theme_void()
ggsave(ISG15, filename = "../AsthmaSeq-EOS-paper/figures/ISG15_feature.svg")

XAF1 <- FeaturePlot(sub, features = "XAF1", order=T, pt.size = 0.5) + scale_color_viridis(option = "magma") + ggtitle("") + theme_void()
ggsave(XAF1, filename = "../AsthmaSeq-EOS-paper/figures/XAF1_feature.svg")

RSAD2 <- FeaturePlot(sub, features = "RSAD2", order=T, pt.size = 0.5) + scale_color_viridis(option = "magma") + ggtitle("") + theme_void()
ggsave(XAF1, filename = "../AsthmaSeq-EOS-paper/figures/RSAD2_feature.svg")

IFIT3 <- FeaturePlot(sub, features = "IFIT3", order=T, pt.size = 0.5) + scale_color_viridis(option = "magma") + ggtitle("") + theme_void()
ggsave(IFIT3, filename = "../AsthmaSeq-EOS-paper/figures/IFIT3_feature.svg")




d <- FindMarkers(obj, ident.1 = "mild", ident.2 = "healthy", min.pct = 0.1, logfc.threshold = 0.25, group.by = "severity")

# GSEA
d <- d[order(-d$avg_log2FC),]

gene_list <- d$avg_log2FC
names(gene_list) <- rownames(d)
gene_list

m_df <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

m_df$symbol <- getSYMBOL(as.character(m_df$entrez_gene), data = 'org.Hs.eg')

H_t2g <- m_df %>% dplyr::select(gs_name, symbol)

em2 <- GSEA(gene_list, TERM2GENE = H_t2g, eps = 1e-300, pvalueCutoff = 0.05, pAdjustMethod = "bonferroni")
df <- as.data.frame(em2@result)



