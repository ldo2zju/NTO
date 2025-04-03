####pre-processing single-cell RNA-seq data
###R version=3.6
library(Seurat)
library(monocle)
###construction of seurat object
NTO_count <- Read10X(data.dir = "10x_matrix_fold_path")
NTO <- CreateSeuratObject(counts = NTO_count, project = "NTO", min.cells = 3, min.features = 200)
NTO[["percent.mt"]] <- PercentageFeatureSet(NTO, pattern = "^MT")
VlnPlot(NTO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(NTO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(NTO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
NTO <- subset(NTO, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5)
NTO <- NormalizeData(NTO)
NTO <- FindVariableFeatures(NTO, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(NTO)
all.genes <- rownames(NTO)
NTO <- ScaleData(NTO, features = all.genes) ##select all.genes can be visualized by dotplot 
NTO <- RunPCA(NTO, features = VariableFeatures(object = NTO))
print(NTO[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(NTO, dims = 1:2, reduction = "pca")
DimPlot(NTO, reduction = "pca")
NTO <- JackStraw(NTO, num.replicate = 100)
NTO <- ScoreJackStraw(NTO, dims = 1:20)
JackStrawPlot(NTO, dims = 1:20)
ElbowPlot(NTO)
NTO <- FindNeighbors(NTO, dims = 1:15)
NTO <- FindClusters(NTO, resolution = 1)
NTO <- RunUMAP(NTO, dims = 1:15,return.model = T)
NTO <- RunTSNE(NTO, dims = 1:15)
DimPlot(NTO, reduction = "umap",label = T)
NTO.markers <- FindAllMarkers(NTO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(NTO.markers,"path/NTO.markers.csv")
saveRDS(NTO,"path/NTO.rds")

###annotation of cell clusters
NTO_markers <- c("BMP4","KRT8","STMN2","ONECUT2","TFAP2B","SOX10","SOX2","SOX3")

NTO_cell_types <- c("Progenitor","Progenitor","Progenitor","Progenitor","Progenitor",
                    "Progenitor","Epidermis","Progenitor","Neuron","Progenitor","Neural Crest")
NTO_color <- c("#D5C6CF","#D52918","#87C798","#A8B7DE")
names(NTO_cell_types) <- levels(NTO)
NTO <- RenameIdents(NTO,NTO_cell_types)
NTO@meta.data$cell_types <- as.character(Idents(NTO))
DimPlot(NTO,cols = NTO_color,pt.size = 1)   

##neural progenitor
NTO_progenitor <- subset(NTO, cells = rownames(NTO@meta.data[which(NTO@meta.data$cell_types=="Progenitor")]))
NTO_progenitor <- ScaleData(NTO_progenitor, vars.to.regress = "percent.mt")
NTO_progenitor <- RunPCA(NTO_progenitor, features = human_spinalcord_markers_pca)
print(NTO_progenitor[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(NTO_progenitor, dims = 1:2, reduction = "pca")
DimPlot(NTO_progenitor, reduction = "pca")
NTO_progenitor <- JackStraw(NTO_progenitor, num.replicate = 100)
NTO_progenitor <- ScoreJackStraw(NTO_progenitor, dims = 1:20)
JackStrawPlot(NTO_progenitor, dims = 1:20)
ElbowPlot(NTO_progenitor)
NTO_progenitor <- FindNeighbors(NTO_progenitor, dims = 1:20)
NTO_progenitor <- FindClusters(NTO_progenitor, resolution = 1)
NTO_progenitor <- RunUMAP(NTO_progenitor, dims = 1:20,return.model = T)
NTO_progenitor <- RunTSNE(NTO_progenitor, dims = 1:20)
DimPlot(NTO_progenitor, reduction = "umap",label = T)
##annotation neural progenitor identity based on DV_markers
DV_markers <- c("SOX2","LMX1A","MSX1","MSX2","PAX3","WNT1","IRX3","IRX5","OLIG3","PAX6","ASCL1","GBX2","GSX2",
                "PAX7","GSX1","DBX2","DBX1","SP8","NKX6-2","SHH")

###cellchat analysis
library("CellChat")
library(Seurat)
library(liana)

cellchat <- createCellChat(object = NTO_progenitor, meta = NTO_progenitor@meta.data, group.by = "cell_types")

# Setup database
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data 
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

# Compute the commnunication probability
cellchat <- computeCommunProb(cellchat)
cell.levels <- c("RP","dp1","dp2","dp3","dp4","dp5","dp6",
                 "p0","p1","p2",
                 "pMN","p3","FP")
cellchat <- updateClusterLabels(cellchat, new.order = cell.levels)

# Filter out by cell numbers
cellchat <- filterCommunication(cellchat, min.cells = 5)



cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

NTO_DV_cell_communication.net <- subsetCommunication(cellchat)

mat <- cellchat@net$weight
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


pathways.show <- c("BMP") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",thresh = 0.1,color.use = c("#FF4500","#E762D7FF","#AD002AFF","#E89242FF","#7E6148FF","#FDAF91FF","#ADB6B6FF",
                                                                                                      "#917C5DFF", "#91D1C2FF","#82491EFF",
                                                                                                      "#3C5488FF","#0099B4FF","#42B540FF"))

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

pathways.show <- c("HH") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",thresh = 0.1)

pathways.show <- c("WNT") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",thresh = 0.1,color.use = c("#FF4500","#E762D7FF","#AD002AFF","#E89242FF","#7E6148FF","#FDAF91FF","#ADB6B6FF",
                                                                                                      "#917C5DFF", "#91D1C2FF","#82491EFF",
                                                                                                      "#3C5488FF","#0099B4FF","#42B540FF"))

netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",signaling = c("WNT"),color.use = c("#FF4500","#E762D7FF","#AD002AFF","#E89242FF","#7E6148FF","#FDAF91FF","#ADB6B6FF",
                                                                                                    "#917C5DFF", "#91D1C2FF","#82491EFF",
                                                                                                    "#3C5488FF","#0099B4FF","#42B540FF"))

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",color.use = c("#FF4500","#E762D7FF","#AD002AFF","#E89242FF","#7E6148FF","#FDAF91FF","#ADB6B6FF",
                                                                                      "#917C5DFF", "#91D1C2FF","#82491EFF",
                                                                                      "#3C5488FF","#0099B4FF","#42B540FF"))
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",color.use = c("#FF4500","#E762D7FF","#AD002AFF","#E89242FF","#7E6148FF","#FDAF91FF","#ADB6B6FF",
                                                                                      "#917C5DFF", "#91D1C2FF","#82491EFF",
                                                                                      "#3C5488FF","#0099B4FF","#42B540FF"))
ht1 + ht2

netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",signaling = c("BMP","HH"))

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)

###TCTN2 vKO and dKO analysis
object.list <- c(NTO,NTO_vKO,NTO_dKO)
object.list <- lapply(X = object.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = object.list)
combined.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = features)
object.combined <- IntegrateData(anchorset = combined.anchors)
DefaultAssay(object.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
object.combined <- ScaleData(object.combined, verbose = FALSE)
object.combined <- RunPCA(object.combined, npcs = 50, verbose = FALSE)
object.combined <- RunUMAP(object.combined, reduction = "pca", dims = 1:50)
object.combined <- FindNeighbors(object.combined, reduction = "pca", dims = 1:50)
object.combined <- FindClusters(object.combined, resolution = 1)

# Visualization
DimPlot(object.combined, reduction = "umap", group.by = "orig.ident")
DimPlot(object.combined, reduction = "umap", label = TRUE, repel = TRUE,split.by = "orig.ident")


DefaultAssay(object.combined) <- "RNA"
object.combined.markers <- FindAllMarkers(object.combined,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(object.combined.markers,"object.combined.markers.csv")
object.combined.markers_top20 <- object.combined.markers %>% group_by(cluster) %>% top_n(20,avg_log2FC)
write.csv(object.combined.markers_top20,"object.combined.markers_top20.csv")
DimPlot(object.combined,group.by = "orig.ident")
DimPlot(object.combined,split.by = "orig.ident",label = T)
object.combined@meta.data$orig.ident <- factor(object.combined@meta.data$orig.ident,levels = c("NTO_human","NTO_vKO","NTO_dKO"))

NTO_combined_cell_types <- c("Progenitor","Progenitor","Progenitor","Progenitor","Progenitor_SHH_high",
                             "Progenitor","Progenitor","Neuron","Epidermis","Progenitor","Progenitor",
                             "Progenitor","Progenitor","Progenitor","Progenitor","Progenitor",
                             "Progenitor_proliferation","Progenitor_SHH_high","Neuron","Neural Crest")

names(NTO_combined_cell_types) <- levels(object.combined)
object.combined <- RenameIdents(object.combined,NTO_combined_cell_types)
object.combined@meta.data$cell_types_combined <- as.character(Idents(object.combined))
levels(object.combined) <- c("Progenitor","Epidermis","Neuron","Neural Crest","Progenitor_proliferation","Progenitor_SHH_high")
NTO_combined_color <- c("#F39800","#D52918","#87C798","#A8B7DE","#BFCB92","#6FAC04")
DimPlot(object.combined,label = F,cols = NTO_combined_color,split.by = "orig.ident")

###reanalyze progenitor cells
Progenitor_cells <- object.combined@meta.data %>% filter(cell_types_combined =="Progenitor") %>% rownames()
object.combined_progenitor <- subset(object.combined,cells = Progenitor_cells)
human_NTO_progenitor.list <- SplitObject(object.combined_progenitor,split.by = "orig.ident")

human_NTO_progenitor.list <- lapply(X = human_NTO_progenitor.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = human_NTO_progenitor.list, anchor.features = features)
# this command creates an 'integrated' data assay
Human_NTO_progenitor.combined <- IntegrateData(anchorset = anchors)

DefaultAssay(Human_NTO_progenitor.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
Human_NTO_progenitor.combined <- ScaleData(Human_NTO_progenitor.combined, verbose = FALSE)
Human_NTO_progenitor.combined <- RunPCA(Human_NTO_progenitor.combined, npcs = 50, verbose = FALSE)
Human_NTO_progenitor.combined <- RunUMAP(Human_NTO_progenitor.combined, reduction = "pca", dims = 1:50)
Human_NTO_progenitor.combined <- FindNeighbors(Human_NTO_progenitor.combined, reduction = "pca", dims = 1:50)
Human_NTO_progenitor.combined <- FindClusters(Human_NTO_progenitor.combined, resolution = 2)
DimPlot(Human_NTO_progenitor.combined, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(Human_NTO_progenitor.combined, reduction = "umap", label = TRUE, repel = TRUE,split.by = "orig.ident")
DefaultAssay(Human_NTO_progenitor.combined) <- "RNA"
FeaturePlot(Human_NTO_progenitor.combined,features = c("FOXA2","NKX2-2","OLIG2","PAX6","OLIG3","MSX1","PAX3","DBX1","DBX2"),pt.size = 0.2)

DefaultAssay(Human_NTO_progenitor.combined) <- "RNA"
Human_NTO_progenitor.combined.markers <- FindAllMarkers(Human_NTO_progenitor.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Human_NTO_progenitor.combined.markers_top20 <- Human_NTO_progenitor.combined.markers %>% group_by(cluster) %>% top_n(20,avg_log2FC)
write.csv(Human_NTO_progenitor.combined.markers,"Human_NTO_progenitor.combined.markers.csv")
write.csv(Human_NTO_progenitor.combined.markers_top20,"Human_NTO_progenitor.combined.markers_top20.csv")

Human_NTO_combined_celltypes <- c("Progenitor_SC","Progenitor_SC","Progenitor_brain","Progenitor_brain",
                                  "Progenitor_SC","Progenitor_SC","Progenitor_brain","Progenitor_SC",
                                  "Progenitor_SC","Progenitor_SC","Progenitor_SC","Progenitor_SC",
                                  "Progenitor_SC","Progenitor_SC","Progenitor_SC","Progenitor_SC",
                                  "Progenitor_SC","Progenitor_SC","Progenitor_SC","Progenitor_brain",
                                  "Progenitor_SC","Progenitor_SC","Progenitor_SC","Progenitor_SC")

names(Human_NTO_combined_celltypes) <- levels(Human_NTO_progenitor.combined)
Human_NTO_progenitor.combined <- RenameIdents(Human_NTO_progenitor.combined,Human_NTO_combined_celltypes)
Human_NTO_progenitor.combined@meta.data$cell_types_combined_2 <- as.character(Idents(Human_NTO_progenitor.combined))
DimPlot(Human_NTO_progenitor.combined,cols = c("#A65628","#FFFF33"))
FeaturePlot(Human_NTO_progenitor.combined,features = c("HOXB4","GDF15","SLC2A3","CDKN1A"),cols = c("grey","red"))
DotPlot(Human_NTO_progenitor.combined,features = c("HOXB4","GDF15","SLC2A3","CDKN1A"),cols = c("grey","red"))
saveRDS(Human_NTO_progenitor.combined,"path/Human_NTO_progenitor.combined.rds")
Human_NTO_DV_progenitor.combined_cells <- Human_NTO_progenitor.combined@meta.data %>% filter(cell_types_combined_2=="Progenitor_SC") %>% rownames()
Human_NTO_DV_progenitor.combined <- subset(Human_NTO_progenitor.combined,cells = Human_NTO_DV_progenitor.combined_cells)
###analyzing DV patterned progenitors
features <- unique(human_spinalcord_DV_markers$gene)[unique(human_spinalcord_DV_markers$gene) %in% Human_NTO_DV_progenitor.combined@assays$RNA@counts@Dimnames[[1]]]
human_NTO_DV_progenitor.list <- SplitObject(Human_NTO_DV_progenitor.combined,split.by = "orig.ident")
human_NTO_DV_progenitor.list <- lapply(X = human_NTO_DV_progenitor.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = human_NTO_DV_progenitor.list, anchor.features = features)
# this command creates an 'integrated' data assay
Human_NTO_DV_progenitor.combined <- IntegrateData(anchorset = anchors)

DefaultAssay(Human_NTO_DV_progenitor.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
Human_NTO_DV_progenitor.combined <- ScaleData(Human_NTO_DV_progenitor.combined, verbose = FALSE)
Human_NTO_DV_progenitor.combined <- RunPCA(Human_NTO_DV_progenitor.combined, npcs = 50, verbose = FALSE)
Human_NTO_DV_progenitor.combined <- RunUMAP(Human_NTO_DV_progenitor.combined, reduction = "pca", dims = 1:50)
Human_NTO_DV_progenitor.combined <- FindNeighbors(Human_NTO_DV_progenitor.combined, reduction = "pca", dims = 1:50)
Human_NTO_DV_progenitor.combined <- FindClusters(Human_NTO_DV_progenitor.combined, resolution = 2)
DimPlot(Human_NTO_DV_progenitor.combined, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(Human_NTO_DV_progenitor.combined) <- "RNA"
FeaturePlot(Human_NTO_DV_progenitor.combined,features = c("FOXA2","NKX2-2","OLIG2","PAX6","OLIG3","MSX1","PAX3","DBX1","DBX2"),pt.size = 0.2)
Human_NTO_DV_progenitor.combined.markers <- FindAllMarkers(Human_NTO_DV_progenitor.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Human_NTO_DV_progenitor.combined.markers_top20 <- Human_NTO_DV_progenitor.combined.markers %>% group_by(cluster) %>% top_n(20,avg_log2FC)
write.csv(Human_NTO_DV_progenitor.combined.markers,"Human_NTO_DV_progenitor.combined.markers_res1.csv")
write.csv(Human_NTO_DV_progenitor.combined.markers_top20,"Human_NTO_DV_progenitor.combined.markers_top20_res1.csv")
###change resolution 
DefaultAssay(Human_NTO_DV_progenitor.combined) <- "integrated"
Human_NTO_DV_progenitor.combined <- FindClusters(Human_NTO_DV_progenitor.combined, resolution = 4)
DimPlot(Human_NTO_DV_progenitor.combined, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(Human_NTO_DV_progenitor.combined) <- "RNA"
FeaturePlot(Human_NTO_DV_progenitor.combined,features = c("FOXA2","NKX2-2","OLIG2","PAX6","OLIG3","MSX1","PAX3","DBX1","DBX2"),pt.size = 0.2)
Human_NTO_DV_progenitor.combined.markers_res4 <- FindAllMarkers(Human_NTO_DV_progenitor.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Human_NTO_DV_progenitor.combined.markers_top20_res4 <- Human_NTO_DV_progenitor.combined.markers_res4 %>% group_by(cluster) %>% top_n(20,avg_log2FC)
write.csv(Human_NTO_DV_progenitor.combined.markers_res4,"Human_NTO_DV_progenitor.combined.markers_res4.csv")
write.csv(Human_NTO_DV_progenitor.combined.markers_top20_res4,"Human_NTO_DV_progenitor.combined.markers_top20_res4.csv")

##FP
DotPlot(Human_NTO_DV_progenitor.combined,features = c("SHH","FOXA2"))
##P3
DotPlot(Human_NTO_DV_progenitor.combined,features = c("NKX2-2","NKX2-8","NKX6-2"))
##pMN
DotPlot(Human_NTO_DV_progenitor.combined,features = c("OLIG2","NKX6-1"))
##p2
DotPlot(Human_NTO_DV_progenitor.combined,features = c("FOXN4","NKX6-1","MOXD1"))


cell_types_annotated <- c("dp5","FP","p1","p1","p0","pMN","dp4","p3","RP","p0","p2",
                          "pMN","dp2","dp3","dp1","dp6","RP","dp6","dp2","dp4","dp6","p0",
                          "RP","dp1","pMN","dp4","p2","p3","dp6","p2","FP","RP","dp3",
                          "dp2","dp4")

names(cell_types_annotated) <- levels(Human_NTO_DV_progenitor.combined)
Human_NTO_DV_progenitor.combined <- RenameIdents(Human_NTO_DV_progenitor.combined,cell_types_annotated)
Human_NTO_DV_progenitor.combined@meta.data$cell_types_annotated <- as.character(Idents(Human_NTO_DV_progenitor.combined))


NT_DV_color <- c("#FF4500","#E762D7FF","#AD002AFF","#E89242FF","#7E6148FF","#FDAF91FF","#ADB6B6FF",
                 "#917C5DFF", "#91D1C2FF","#82491EFF",
                 "#3C5488FF","#0099B4FF","#42B540FF")

levels(Human_NTO_DV_progenitor.combined) <- c("RP","dp1","dp2","dp3","dp4","dp5","dp6",
                                              "p0","p1","p2",
                                              "pMN","p3","FP")
Human_NTO_DV_progenitor.combined@meta.data$orig.ident <- factor(Human_NTO_DV_progenitor.combined@meta.data$orig.ident,levels = c("NTO_human","NTO_WT_day11","NTO_TKO_day11"))
FeaturePlot(Human_NTO_DV_progenitor.combined,features = c("FOXA2","NKX2-2","OLIG2","PAX6","OLIG3","MSX1","PAX3","DBX1","DBX2"))

DimPlot(Human_NTO_DV_progenitor.combined,cols = NT_DV_color,split.by = "orig.ident",pt.size =0.5,label = T)




















































































































