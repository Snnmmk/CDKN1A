library(Seurat)
library(tidyverse)
library(CellChat)
library(progeny)
library(viridis)
library(monocle3)
library(ComplexHeatmap)
library(CellChat)
library(ggalluvial)
library(patchwork)
library(igraph)
library(svglite)
library(ggsignif)
library(TCellSI)
#All data could be downloaded in NCBI GEO: GSE205013, GSE155698, GSE149813, GSE176078; Code Oceanic capsule: 10.24433/CO.0121060.v1; Beijing Institute of Genomics (BIG) : HRA001130 

####Extra Step 01#####

meta_data <- sce@meta.data
meta_data$CDKN1A_expression <- FetchData(sce.all.int, vars = "CDKN1A")[,1]
threshold <- median(meta_data$CDKN1A_expression, na.rm = TRUE)
meta_data$CDKN1A_expression_level <- ifelse(meta_data$CDKN1A_expression > threshold, "High", "Low")
sce$CDKN1A_expression_level<-meta_data$CDKN1A_expression_level

####Extra Step 02#####
malignant_indices <- which(sce$Malignant == "Malignant Epithelial")
sce$SubLabel[malignant_indices] <- ifelse(
  sce$CDKN1A_expression_level[malignant_indices] == "Low",
  "Malignant Epithelial LOW",
  "Malignant Epithelial HIGH"
)

####Figure 4A 4B  ####
DimPlot(sce, label = F, group.by = "CellType")  &
  theme_pubr(base_size = 8) &
  #NoLegend() &
  theme(plot.title = element_blank(), 
        legend.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        legend.position = "right", legend.key.size = unit(2, "pt"),
        legend.text = element_text(size = 12) ,
        legend.justification = c(0,0.5)) &
  scale_color_manual(values = custom_colors)#you can change the color as you wish

###Figure 4C 4D  ####
custom_colors <- c('darkred','lightblue')
DimPlot(sce, label = F, group.by = "CDKN1A_expression_level")  &
  theme_pubr(base_size = 8) &
  #NoLegend() &
  theme(plot.title = element_blank(), 
        legend.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        legend.position = "right", legend.key.size = unit(2, "pt"),
        legend.text = element_text(size = 12) ,
        legend.justification = c(0,0.5)) &
  scale_color_manual(values = custom_colors)

###Figure 4E & 5F####
sce@meta.data %>%
  ggplot(aes(x = CDKN1A_expression_level, fill = CellType)) +
  theme_pubr(base_size = 8) +
  theme(
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.key.size = unit(3, "pt"),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),  
    legend.position = "right"
  ) +
  facet_grid(~Group, scales = "free_x", space = "free_x") +
  geom_bar(position = "fill", alpha = 0.9) +
  scale_fill_manual(values = custom_colors_v2) +
  ylab("'Cancer Type Fraction'")

####Extra Step 03#####
sce <- progeny(sce, scale=FALSE, organism="Human", top=500, perm=1, 
                return_assay = TRUE)
sce <- Seurat::ScaleData(sce, assay = "progeny") 
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(sce, layer = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 
CellsClusters <- data.frame(Cell = names(Idents(sce)), 
                            CellType = as.character(sce$SubLabel),
                            stringsAsFactors = FALSE)
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, SubLabel) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

####Figure 5G 5H ####
progeny_hmap = pheatmap(t(summarized_progeny_scores_df),
                        fontsize=12, 
                        fontsize_row = 10, color=myColor,
                        breaks = progenyBreaks, main = "Activating Pathways",
                        angle_col = '45' ,treeheight_col = 4.5, border_color = NA)

####Extra Step 04#####
data <- GetAssayData(sce,assay = "RNA",layer = "counts")
cell_metadata <- sce@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,cell_metadata = cell_metadata,gene_metadata = gene_annotation)
cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds,cluster_method = "louvain")
cds <- learn_graph(cds)
cds <- order_cells(cds, root_cells = root_cells)

###Figure 5A 5C ####
plot_cells(cds, color_cells_by = "SubLabel",label_roots = FALSE, label_cell_groups=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE, graph_label_size=1.5)

#####Figure 5B 5D####
plot_cells(cds, color_cells_by = "partition",label_groups_by_cluster=FALSE, label_leaves=FALSE,label_branch_points=FALSE)

#####Figure 6G ####
plot_cells(cds, genes=c("CDKN1A","MYC","EPCAM","KRAS","PMS2","MSH2","MSH6","MLH1"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)



####Supplementary Figure6A####
plot_genes_in_pseudotime(cds[c("CDKN1A","MYC","EPCAM","KRAS","PMS2","MSH2","MSH6","MLH1"),],
                         color_cells_by="Malignant_Epithelial_CDKN1A",
                         min_expr=0.5)

####Extra Step 03#####
cellchat <- createCellChat(object =GetAssayData(sce, assay = "RNA", layer = "data"), meta = sce@meta.data, group.by = "SubLabel")

cellchat <- addMeta(cellchat, meta = sce@meta.data)
cellchat <- setIdent(cellchat, ident.use = "SubLabel")

CellChatDB.use <- CellChatDB.human
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)


#####Figure 7####
selected_pathways <- c("SPP1", "FN1","COLLAGEN")

netAnalysis_signalingRole_network(cellchat, signaling = selected_pathways, width = 14, height = 6, font.size = 15)

netVisual_aggregate(
  object = cellchat,
  signaling = selected_pathways,
  top = 1,
  layout = "circle"
)
