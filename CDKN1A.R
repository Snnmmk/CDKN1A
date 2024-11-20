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

#All data has been uploaded in BaiDu Netdisk

sce_pancreas<-readRDS('./pancreas.rds')
sce_lung<-readRDS('./lung.rds')

####Figure 5A & 5C####
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

###Figure 5B & 5D####
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

###Figure 5E & 5F####
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
  xlab("CDKN1A Expression Level") + 
  ylab("Fraction")

####Figure 5G & 5H####
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
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
progeny_hmap = pheatmap(t(summarized_progeny_scores_df),
                        fontsize=12, 
                        fontsize_row = 10, color=myColor,
                        breaks = progenyBreaks, main = "Activating Pathways",
                        angle_col = '45' ,treeheight_col = 4.5, border_color = NA)

###Figure 8A & 8C####
data <- GetAssayData(sce,assay = "RNA",layer = "counts")
cell_metadata <- sce@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,cell_metadata = cell_metadata,gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 7) #when processes pancreas, num_dim=5
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds,cluster_method = "louvain")
cds <- learn_graph(cds)
cds <- order_cells(cds, root_cells = root_cells)

plot_cells(cds, color_cells_by = "CellType",label_roots = FALSE,label_groups_by_cluster=FALSE, label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cds, color_cells_by = "SubLabel",label_roots = FALSE, label_cell_groups=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE, graph_label_size=1.5)

#####Figure 8B & 8D####
plot_cells(cds, color_cells_by = "partition",label_groups_by_cluster=FALSE, label_leaves=FALSE,label_branch_points=FALSE)

#####Figure 8E & 8F####
plot_cells(cds, genes=c("CDKN1A","MYC","EPCAM","KRAS","PMS2","MSH2","MSH6","MLH1"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

####Supplementary Figure6A####
plot_genes_in_pseudotime(cds[c("CDKN1A","MYC","EPCAM","KRAS","JAK1","JAK2","STAT1","STAT2","STAT3"),],
                         color_cells_by="Malignant_Epithelial_TFAP4",
                         min_expr=0.5)

####Supplementary Figure7####

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



selected_pathways <- c("SPP1", "CDH1", "FN1","COLLAGEN")

netAnalysis_signalingRole_network(cellchat, signaling = selected_pathways, width = 14, height = 6, font.size = 15)

netVisual_aggregate(
  object = cellchat,
  signaling = selected_pathways,
  top = 1,
  layout = "circle"
)




