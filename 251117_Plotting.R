
library(qs)
library(dplyr)
library(tibble)
library(Seurat)
library(ggpubr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(magrittr)

# 1.1 -----
library(Seurat)
library(ggplot2)
library(patchwork)
seuratobj <- qread("D:/2025 PainHY/data/snRNAseq_Mm_HYiPain_XiaoYu/processed/250712.qs")

# Define color palettes
color_brain <- c(
  "LH" = "#66C2A5",
  "PVN" = "#8DA0CB"
  # Add more regions as needed
)

color_group <- c(
  "Young" = "#4DAF4A",
  "Aging" = "#984EA3"
  # Adjust based on your actual group names
)

# Plot BrainRegion
p1 <- DimPlot(
  seuratobj,
  group.by = "BrainRegion",
  label = FALSE,
  pt.size = 3,
  cols = color_brain
) +
  ggtitle("Brain Region")

# Plot Group
p2 <- DimPlot(
  seuratobj,
  group.by = "Group",
  label = FALSE,
  pt.size = 3,
  cols = color_group
) +
  ggtitle("Group")

# Combine plots
p1 | p2

# 1.2 -----
FeaturePlot(seuratobj,features = c(
  "Slc17a6"),pt.size = 2,keep.scale = "all",order = TRUE) +
  scale_colour_gradient(low = "lightgrey", high = "red")

FeaturePlot(seuratobj,features = c(
  "Slc32a1"),pt.size = 2,keep.scale = "all",order = TRUE) +
  scale_colour_gradient(low = "lightgrey", high = "red")

# 1.3 -----
FeaturePlot(
  seuratobj,
  order = TRUE,
  ncol = 4,
  features = c("Aqp4","Mobp","Mog","Pdgfra"),
  cols = c("lightgrey","red"),
  pt.size = 2,
  keep.scale = "all")

# 2.1 -----
FeaturePlot(seuratobj,
  order = TRUE,
  ncol = 4,
  features = c("Sim1","Hcrt","Crh","Rfx4"),
  cols = c("lightgrey","red"),
  pt.size = 2,
  keep.scale = "all")

# 2.2 -----
color_ct <- c(
  "LH_Hcrt+"   = "#66C2A5",  # Original LH color
  "LH_Htr1b+"  = "#308e73ff",  # Darker shade of LH color
  "PVN_Npy2r+" = "#8DA0CB",  # Original PVN color
  "PVN_Avp+"   = "#263994ff"   # Darker shade of PVN color
)

DimPlot(
  seuratobj,
  group.by = "celltype_coarse",
  label = FALSE,
  pt.size = 3,
  cols = color_ct
)

# 2.3 -----
DotPlot(seuratobj,features = c(
  "Hcrt","Rfx4","Htr1b","Npy2r","Avp"),group.by = "celltype_coarse") +
  scale_colour_gradient(low = "lightgrey", high = "red")

# 3.1 -----
library(Seurat)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)

# Extract data for plotting
plot_data <- FetchData(
  seuratobj,
  vars = c("celltype_coarse", "Group", "Sen_Mayo1")
) %>%
  mutate(
    Group = factor(Group, levels = c("Young", "Aging")),
    celltype_coarse = factor(celltype_coarse)
  )

# Calculate statistics for each cell type
stat_test <- plot_data %>%
  group_by(celltype_coarse) %>%
  t_test(Sen_Mayo1 ~ Group) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "celltype_coarse", dodge = 0.8)

# Create plot
p <- ggplot(plot_data, aes(x = celltype_coarse, y = Sen_Mayo1, fill = Group)) +
  # Violin plot
  geom_violin(
    position = position_dodge(width = 0.8),
    trim = TRUE,
    alpha = 0.6,
    scale = "width"
  ) +
  # Boxplot overlay
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.2,
    alpha = 0.8,
    outlier.shape = NA
  ) +
  # Scatter points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.8
    ),
    size = 1.5,
    alpha = 0.3
  ) +
  # Add statistical brackets and p-values
  stat_pvalue_manual(
    stat_test,
    label = "p.adj.signif",
    tip.length = 0.01,
    bracket.size = 0.5,
    step.increase = 0.05
  ) +
  # Colors
  scale_fill_manual(values = color_group) +
  # Labels
  labs(
    title = "Senescence Score (SenMayo) by Cell Type",
    x = "Cell Type",
    y = "Sen_Mayo Score",
    fill = "Group"
  ) +
  # Theme
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_blank()
  )

p

# 3.2 -----
### Augur
# library(Augur)
# augur_res = calculate_auc(
#   seuratobj,
#   cell_type_col = "celltype_coarse",label_col = "Group",
#   n_subsamples = 50,subsample_size = 20,folds = 3,
#   rf_params = list(mtry = 2,trees = 50,importance = "accuracy"),
#   n_threads = 1)
augur_res_df <- data.frame(
  cell_type = c("PVN_Npy2r+", "LH_Htr1b+", "LH_Hcrt+", "PVN_Avp+"),
  auc = c(0.6233333, 0.5600000, 0.5541667, 0.5366667),
  stringsAsFactors = FALSE
)

### scDist
library(scDist)
scDist_res <- scDist(
  normalized_counts = GetAssayData(seuratobj,assay = "RNA",layer = "data"),meta.data = seuratobj@meta.data,
  fixed.effects = "Group",clusters="celltype_coarse")
# DistPlot(scDist_res)
scdist_df <- scDist_res$results %>% 
  rownames_to_column("cell_type")

### Plot
combined_df <- augur_res_df %>%
  left_join(scdist_df, by = "cell_type") %>%
  mutate(cell_type = factor(cell_type, levels = c("LH_Hcrt+", "LH_Htr1b+", 
                                                    "PVN_Npy2r+", "PVN_Avp+")))

cor_result <- cor.test(combined_df$auc, combined_df$Dist., method = "pearson")
cor_value <- cor_result$estimate
p_value <- cor_result$p.value

# Create correlation plot
ggplot(combined_df, aes(x = auc, y = Dist.)) +
  # Points with cell type colors
  geom_point(aes(fill = cell_type), 
             size = 7, 
             shape = 21, 
             color = "black",
             stroke = 1) +
  # Add cell type labels
  ggrepel::geom_text_repel(
    aes(label = cell_type),
    size = 4,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  # Add correlation annotation
  annotate("text",
           x = min(combined_df$auc) + 0.01,
           y = max(combined_df$scDist) * 0.95,
           label = sprintf("Pearson r = %.3f\np = %.3f", cor_value, p_value),
           hjust = 0,
           vjust = 1,
           size = 4.5,
           fontface = "bold") +
  # Colors
  scale_fill_manual(values = color_ct) +
  # Labels
  labs(
    title = "Correlation between Augur AUC and scDist",
    subtitle = "Cell Type Perturbation Metrics",
    x = "Augur AUC",
    y = "scDist (D̂)",
    fill = "Cell Type"
  ) +
  # Theme
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# 4 -----
seuratobj <- qread("D:/2025 PainHY/data/snRNAseq_Mm_HYiPain_XiaoYu/processed/250712.qs")

con <- DBI::dbConnect(RSQLite::SQLite(), dbname = 'D:/2025 PainHY/bin/msigdb_v2025.1.Mm.db/msigdb_v2025.1.Mm.db')
DBI::dbListTables(con)

# define tables we want to combine
geneset_db <- dplyr::tbl(con, 'gene_set')                                              # standard_name, collection_name
details_db <- dplyr::tbl(con, 'gene_set_details')                                      # description_brief, description_full
geneset_genesymbol_db <- dplyr::tbl(con, 'gene_set_gene_symbol')                       # meat and potatoes
genesymbol_db <- dplyr::tbl(con, 'gene_symbol')                                        # mapping from ids to gene symbols
collection_db <- dplyr::tbl(con, 'collection') %>% dplyr::select(collection_name, full_name)  # collection metadata

# join tables
msigdb <- geneset_db %>%
  dplyr::left_join(details_db, by = c('id' = 'gene_set_id')) %>%
  dplyr::left_join(collection_db, by = 'collection_name') %>%
  dplyr::left_join(geneset_genesymbol_db, by = c('id' = 'gene_set_id')) %>%
  dplyr::left_join(genesymbol_db, by = c('gene_symbol_id' = 'id')) %>%
  dplyr::select(collection = collection_name, subcollection = full_name, geneset = standard_name, description = description_brief, symbol) %>%
  dplyr::as_tibble() 

# clean up
DBI::dbDisconnect(con)
unique(msigdb$collection)
unique(msigdb$subcollection)

# convert to list[Hallmark] required by irGSEA package
msigdb.h <- msigdb %>% 
  #dplyr::filter(geneset %in% pathway.list) %>% 
  dplyr::filter(collection=="M5:GO:BP") %>% 
  dplyr::select(c("geneset", "symbol"))
msigdb.h$geneset <- factor(msigdb.h$geneset)
msigdb.h <- msigdb.h %>% 
  dplyr::group_split(geneset, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(symbol) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.h$geneset))

# filter gluN
library(irGSEA)
irGSEA_res <- irGSEA.score(object = seuratobj,assay = "RNA",slot = "data",
                            seeds = 777,ncores = 1,
                            custom = T, geneset = msigdb.h, 
                            #method = c("AUCell"),
                            #method = c("AUCell","UCell","singscore","ssgsea","JASMINE", "viper"),
                            method = c("singscore","ssgsea"),
                            kcdf = 'Gaussian')

result.dge <- irGSEA.integrate(object = irGSEA_res,
                               group.by = "celltype_coarse",
                               #method = c("AUCell")
                               #method = c("AUCell","UCell","singscore","ssgsea","JASMINE","viper")
                               method = c("singscore","ssgsea")
)

genesets <- result.dge$RRA %>%
  group_by(cluster) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  pull(Name)

color_ct <- c(
  "LH_Hcrt+"   = "#66C2A5",  # Original LH color
  "LH_Htr1b+"  = "#308e73ff",  # Darker shade of LH color
  "PVN_Npy2r+" = "#8DA0CB",  # Original PVN color
  "PVN_Avp+"   = "#263994ff"   # Darker shade of PVN color
)

irGSEA.heatmap(object = result.dge,
               method = "RRA",
               cluster.levels = c("LH_Hcrt+","LH_Htr1b+","PVN_Npy2r+","PVN_Avp+"),
               cluster.color = color_ct,
               direction.color = c("#D73027","#4575B4"),
               significance.color = c("#DDDDDD","#F0E442"),
               top = 40,heatmap.width = 20,heatmap.heigh = 15,
               show.geneset = genesets)

# 5.1 -----
library(fgsea)
fgseaRes_list <- list()
pathways <- gmtPathways("Z:/2023 EpiAdo/bin/MSigDB_Mm/msigdb_v2023.2.Mm_files_to_download_locally/msigdb_v2023.2.Mm_files_to_download_locally/msigdb_v2023.2.Mm_GMTs/msigdb.v2023.2.Mm.symbols.gmt")
seuratobj <- qread("D:/2025 PainHY/data/snRNAseq_Mm_HYiPain_XiaoYu/processed/250712.qs")
options(mc.cores = 1)
for(celltype in unique(seuratobj$celltype_coarse)) {
  cat(celltype,"\n")
  markers <- read.csv(glue::glue("D:/2025 PainHY/src/{celltype}.csv"),row.names = 1)
 
  gene_rank <- markers %>%
    tibble::rownames_to_column("gene") %>%
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene, avg_log2FC)
  rank <- setNames(gene_rank$avg_log2FC, gene_rank$gene)
  
  fgseaRes <- fgsea(pathways,rank,minSize=1,maxSize=500, nproc=1)
  fgseaRes_list[[celltype]] <- fgseaRes
}

res.tp <- fgseaRes_list$`LH_Hcrt+` %>%
  dplyr::filter(str_detect(pathway, "GOBP_")) %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Description = str_remove(pathway, "GOBP_") %>% str_to_title(),
    pathwayGenes = paste(unlist(leadingEdge), collapse = "/"),
    colorBy = NES,
    nodeSize = size
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(Description, pathwayGenes, colorBy, nodeSize) %>%
  as.data.frame()
aPEAR::enrichmentNetwork(
  res.tp,
  repelLabels = TRUE,fontSize = 5,minClusterSize = 3,
  colorBy = "colorBy",nodeSize = "nodeSize",drawEllipses = TRUE)

# 5.2 -----
fgseaRes.GO <- fgseaRes_list$`LH_Hcrt+` %>%
  dplyr::filter(str_detect(pathway, "GOBP_"))
topPathwaysUp <- fgseaRes.GO[NES > 0][head(order(NES,decreasing = TRUE), n=10), pathway]
topPathwaysDown <- fgseaRes.GO[ES < 0][head(order(NES), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], rank, fgseaRes.GO, 
              gseaParam=0.5,
              pathwayLabelStyle = list(size = 8),colwidths = c(5,3,1,1,1))

# 6.1 -----
res.tp <- fgseaRes_list$`LH_Htr1b+` %>%
  dplyr::filter(str_detect(pathway, "GOBP_")) %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Description = str_remove(pathway, "GOBP_") %>% str_to_title(),
    pathwayGenes = paste(unlist(leadingEdge), collapse = "/"),
    colorBy = NES,
    nodeSize = size
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(Description, pathwayGenes, colorBy, nodeSize) %>%
  as.data.frame()
aPEAR::enrichmentNetwork(
  res.tp,
  repelLabels = TRUE,fontSize = 5,minClusterSize = 3,
  colorBy = "colorBy",nodeSize = "nodeSize",drawEllipses = TRUE)

# 6.2 -----
fgseaRes.GO <- fgseaRes_list$`LH_Htr1b+` %>%
  dplyr::filter(str_detect(pathway, "GOBP_"))
topPathwaysUp <- fgseaRes.GO[NES > 0][head(order(NES,decreasing = TRUE), n=10), pathway]
topPathwaysDown <- fgseaRes.GO[ES < 0][head(order(NES), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], rank, fgseaRes.GO, 
              gseaParam=0.5,
              pathwayLabelStyle = list(size = 8),colwidths = c(5,3,1,1,1))

# 7.1 -----
res.tp <- fgseaRes_list$`PVN_Npy2r+` %>%
  dplyr::filter(str_detect(pathway, "GOBP_")) %>%
  dplyr::filter(pval < 0.01) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Description = str_remove(pathway, "GOBP_") %>% str_to_title(),
    pathwayGenes = paste(unlist(leadingEdge), collapse = "/"),
    colorBy = NES,
    nodeSize = size
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(Description, pathwayGenes, colorBy, nodeSize) %>%
  as.data.frame()
aPEAR::enrichmentNetwork(
  res.tp,
  repelLabels = TRUE,fontSize = 5,minClusterSize = 3,
  colorBy = "colorBy",nodeSize = "nodeSize",drawEllipses = TRUE)

# 7.2 -----
fgseaRes.GO <- fgseaRes_list$`PVN_Npy2r+` %>%
  dplyr::filter(str_detect(pathway, "GOBP_"))
topPathwaysUp <- fgseaRes.GO[NES > 0][head(order(NES,decreasing = TRUE), n=10), pathway]
topPathwaysDown <- fgseaRes.GO[ES < 0][head(order(NES), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], rank, fgseaRes.GO, 
              gseaParam=0.5,
              pathwayLabelStyle = list(size = 8),colwidths = c(5,3,1,1,1))

# 8.1 -----
res.tp <- fgseaRes_list$`PVN_Avp+` %>%
  dplyr::filter(str_detect(pathway, "GOBP_")) %>%
  dplyr::filter(pval < 0.01) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Description = str_remove(pathway, "GOBP_") %>% str_to_title(),
    pathwayGenes = paste(unlist(leadingEdge), collapse = "/"),
    colorBy = NES,
    nodeSize = size
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(Description, pathwayGenes, colorBy, nodeSize) %>%
  as.data.frame()
aPEAR::enrichmentNetwork(
  res.tp,
  repelLabels = TRUE,fontSize = 5,minClusterSize = 3,
  colorBy = "colorBy",nodeSize = "nodeSize",drawEllipses = TRUE)

# 8.2 -----
fgseaRes.GO <- fgseaRes_list$`PVN_Avp+` %>%
  dplyr::filter(str_detect(pathway, "GOBP_"))
topPathwaysUp <- fgseaRes.GO[NES > 0][head(order(NES,decreasing = TRUE), n=10), pathway]
topPathwaysDown <- fgseaRes.GO[ES < 0][head(order(NES), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], rank, fgseaRes.GO, 
              gseaParam=0.5,
              pathwayLabelStyle = list(size = 8),colwidths = c(5,3,1,1,1))
