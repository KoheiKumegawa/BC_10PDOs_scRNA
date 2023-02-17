#----------------------------------------------------------------------------
# Rcode_for_BC_10PDOs_scRNA
#----------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(org.Hs.eg.db)
library(ggridges)
'%ni%' <- Negate('%in%')

seu <- readRDS("rds/seu.rds")
sampleName <- names(table(seu$Sample))

#Normalize, Find variable features, scale data, and Run PCA
seu_sample <- lapply(sampleName, function(x){
  out <- seu[, which(seu$Sample == x)]
  out <- NormalizeData(out, normalization.method = "LogNormalize", scale.factor = 10000) %>%
         FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>% 
         ScaleData(., features = rownames(.)) %>%
         RunPCA(.)
})
names(seu_sample) <- sampleName

#analyze PC
PCFeatures_ls <- lapply(names(seu_sample), function(x){
  seu_tmp <- seu_sample[[x]]
  out <- lapply(c(1:50), function(i) TopFeatures(seu_tmp, dim = i, nfeatures = 40, balanced = T) %>% unlist()) %>% do.call(cbind, .)
  colnames(out) <- paste0("PC", c(1:50))
  write.csv(out, paste0("output/Tables/GEx_PC_topFeatures_", x,".csv"), quote = F)
  return(NULL)
})

#https://www.genenames.org/data/genegroup/#!/group/1054
rbproteins <- data.table::fread("ref/ribosomeproteins.txt") %>% data.frame() %>% .[,2]
usePC <- lapply(names(seu_sample), function(x){
  seu_tmp <- seu_sample[[x]]
  out <- lapply(c(1:50), function(i){
    f <- TopFeatures(seu_tmp, dim = i, nfeatures = 40, balanced = T)
    n1 <- length(c(which(f$positive %in% rbproteins), grep("^MT.", f$positive)))
    n2 <- length(c(which(f$negative %in% rbproteins), grep("^MT.", f$negative)))
    if(n1 >= 10 | n2 >= 10){return(FALSE)} else {return(TRUE)}
    }) %>% unlist
  return(out)
})
names(usePC) <- sampleName

#removing ribosomal protein weighted PCs, and clustering
seu_sample <- lapply(names(seu_sample), function(x){
  seu_tmp <- seu_sample[[x]]
  useDimention <- c(1:50)[usePC[[x]]] %>% .[c(1:20)]
  out <- FindNeighbors(seu_tmp, dims = useDimention) %>% 
         FindClusters(., resolution = 0.4) %>% 
         RunUMAP(., dims = useDimention, return.model = T)
  out$Clusters <- paste0("C", as.numeric(out$seurat_clusters))
  return(out)
})
names(seu_sample) <- sampleName

saveRDS(seu_sample, "rds/seu_sample_v2.rds")

#UMAP plot
p1 <- lapply(sampleName, function(x){
  tmp <- seu_sample[[x]]
  cluster_cols <- ArchR::ArchRPalettes$kelly[seq_along(names(table(tmp$Clusters)))] %>% 
    `names<-`(., sort(factor(names(table(tmp$Clusters)), levels = paste0("C", c(1:16)))))
  p <- DimPlot(tmp, reduction = "umap", label = F, group.by = "Clusters", cols = cluster_cols) + ggtitle(x) + ArchR::theme_ArchR()
  return(p)
})
pdf("output/Plots/GEx_UMAP_sample_v2.pdf", width = 4, height = 5)
p1
dev.off()

#identify cluster-specific GEx signatures
DEG_ls <- parallel::mclapply(sampleName, function(x){
  tmp <- seu_sample[[x]]
  Idents(tmp) <- tmp$Clusters
  DEG <- FindAllMarkers(tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)
  return(DEG)
}, mc.cores = 10)
names(DEG_ls) <- sampleName

clustGS_ls <- lapply(DEG_ls, function(x){
  out <- x[which(x$p_val_adj < 0.05),]
  out <- out[grep("^MT.", out$gene, invert = T),]
  out <- out[which(out$gene %ni% rbproteins),]
  return(out)
}) #significant genes (adjP < 0.05) and removing ribosome and MT genes

clustGS <- lapply(sampleName, function(x){
  tmp <- clustGS_ls[[x]]
  idx <- paste0("C", seq_along(names(table(tmp$cluster))))
  top50 <- lapply(idx, function(y){tmp[tmp$cluster == y, "gene"][c(1:25)]})
  names(top50) <- paste0(x, "#", idx)
  out <- do.call(cbind, top50)
  return(out)
}) %>% do.call(cbind, .)

clustGS_robust <- clustGS[, apply(clustGS, 2, function(x) length(which(is.na(x) == TRUE))) == 0]
write.csv(clustGS_robust, "output/Tables/GEx_ClustGS.csv", quote = F)

#visualization
cluster_cols <- ArchR::ArchRPalettes$kelly[c(1:15)] %>% `names<-`(., paste0("C", c(1:15)))
col_fun1 <- colorRamp2(c(-2,0,2), c("#3361A5", "white", "#A31D1D"))

p2 <- lapply(sampleName, function(x){
  idx <- clustGS_robust[, grep(x, colnames(clustGS_robust))] %>% as.vector()
  mtx <- seu_sample[[x]]@assays$RNA@scale.data[idx,]
  ha1 <- HeatmapAnnotation(Cluster = seu_sample[[x]]$Clusters, col = list(Cluster = cluster_cols))
  ht1 <- Heatmap(mtx, name = "Relative expression", col = col_fun1, top_annotation = ha1,
                 cluster_rows = F, cluster_columns = F, show_column_names = F,
                 row_names_gp = gpar(fontsize = 4),
                 row_title_gp = gpar(fontsize = 8),
                 column_split = seu_sample[[x]]$Clusters, 
                 row_split = lapply(grep(x, colnames(clustGS_robust), value = T), function(y) rep(y, 25)) %>% unlist(.),
                 column_title = x, use_raster = T)
  out <- draw(ht1)
  return(out)
})
pdf("output/Plots//GEx_ClustGS_Heatmap_Sample.pdf", width = 8, height = 6)
p2
dev.off()

#-------------- Similarity -------------- #
#jaccard index
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

jaI <- NULL
for(i in 1:ncol(clustGS_robust)){
  a <- clustGS_robust[,i]
  tmp2 <- NULL
  for(n in 1:ncol(clustGS_robust)){
    tmp <- jaccard(a, clustGS_robust[, n])
    tmp2 <- c(tmp2, tmp)
  }
  jaI <- rbind(jaI, tmp2)
}
rownames(jaI) <- colnames(clustGS_robust)
colnames(jaI) <- colnames(clustGS_robust)

col_fun2 <- colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8), c("white", rev(viridis(4, option = "G"))))
fh <- function(x) hclust(dist(x), method = "ward.D2")
subtype <- c("ob154" = "ER+/HER2+",
             "ob155" = "ER-/HER2+",
             "ob165" = "ER+/HER2-",
             "ob166" = "ER-/HER2+",
             "ob180" = "ER+/HER2-",
             "ob195" = "ER+/HER2-",
             "ob202" = "ER+/HER2-",
             "ob203" = "ER+/HER2-",
             "ob207P" = "ER+/HER2-",
             "ob210" = "ER+/HER2-")
subtype_cols <- c("ER+/HER2-" = "blue", "ER-/HER2+" = "orange", "ER+/HER2+" = "yellow")
sample_cols  <- ArchR::ArchRPalettes$stallion2[c(1:10)] %>% `names<-`(., names(table(seu$Sample)))

ra1 <- rowAnnotation(sample = stringr::str_split(colnames(jaI), pattern = "#", simplify = T)[,1], 
                     subtype = subtype[stringr::str_split(colnames(jaI), pattern = "#", simplify = T)[,1]],
                     col = list(sample = sample_cols, subtype = subtype_cols))
ht1 <- Heatmap(jaI, name = "Jaccard index", col = col_fun2, right_annotation = ra1,
               cluster_rows = fh, cluster_columns = fh,
               show_column_names = F, show_column_dend = F, show_row_names = T, show_row_dend = T, row_names_gp = gpar(fontsize = 8), 
               column_title = "DEG programs", column_title_side = "bottom")
p3 <- draw(ht1)
pdf("output/Plots//GEx_ClustGS_jaccardHeatmap.pdf", width = 7, height = 6)
p3
dev.off()
write.csv(clustGS_robust[, column_order(p3)], "output/Tables//GEx_ClustGS_genes.csv")
write.csv(round(jaI[row_order(p3), column_order(p3)], digits = 2), "output/Tables//GEx_ClustGS_jaccard.csv")

#-------------- Meta-clusters -------------- #
MC <- list(c(1:9), c(10:12), c(13:17), c(18:19), c(20:23), c(24:26), c(29:31)) %>% `names<-`(., paste0("MC", c(1:7)))
clustGS_robust_ordered <- clustGS_robust[, column_order(p3)]
MC_commonGEx <- lapply(seq_along(MC), function(x){
  tmp <- table(as.vector(as.matrix(clustGS_robust_ordered[, MC[[x]]])))
  out <- sort(names(tmp)[tmp >= 2])
}) %>% `names<-`(., names(MC))
lapply(names(MC_commonGEx), function(x) write.table(data.frame(gene = MC_commonGEx[[x]]), paste0("output/Tables/GEx_MC_gene_", x, ".txt"), quote = F, row.names = F))

#UMAPoverlay
p4 <- lapply(names(MC_commonGEx), function(x){
  gene <- MC_commonGEx[[x]]
  out <- lapply(names(seu_sample), function(i){
    seu_tmp <- seu_sample[[i]]
    seu_tmp$signature <- colMeans(GetAssayData(seu_tmp, slot = "data")[gene,])
    #seu_tmp$signature[which(seu_tmp$signature > 2)] <- 2
    p <- FeaturePlot(seu_tmp, features = "signature", order = T) + ggtitle(paste0(i, ": ", x)) +
      scale_colour_gradientn(colours = viridis(256, option = "D")) + ArchR::theme_ArchR()
    return(p)
  })
  return(out)
})
pdf("output/Plots/GEx_UMAPoverlay_MCsig.pdf", width = 5, height = 6)
p4
dev.off()

#violinplot
p5 <- lapply(names(MC_commonGEx), function(x){
  df <- data.frame(colMeans(GetAssayData(seu, slot = "data")[MC_commonGEx[[x]],]), seu$Sample)
  colnames(df) <- c(x, "sample")
  df <- reshape2::melt(df)
  p <- ggplot(df, aes(x = sample, y = value, fill = variable)) + geom_violin() + ArchR::theme_ArchR() +
    labs(x = "PDO", y = "Expression") + scale_fill_manual(values = "#481769FF") + ggtitle(x)
})
pdf("output/Plots/GEx_VlnPlot_MCsig.pdf", width = 6, height = 3)
p5
dev.off()

#-------------- Meta-clusters annotation --------------#
MC_commonGEx_entrez <- lapply(MC_commonGEx, function(x) AnnotationDbi::select(org.Hs.eg.db, keys = x, keytype = "SYMBOL", columns = c("ENTREZID")) %>% .$ENTREZID)

install.packages("msigdbr")
library(msigdbr)
msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
MC_commonGEx_msigDB <- lapply(MC_commonGEx_entrez, function(x) clusterProfiler::enricher(x, TERM2GENE = msig_hallmark, pAdjustMethod = "fdr", pvalueCutoff = 1, qvalueCutoff = 0.3))
p6 <- lapply(MC_commonGEx_msigDB, function(x) mutate(x, qscore = -log(p.adjust, base=10)) %>% barplot(x="qscore") + ArchR::theme_ArchR())
pdf("output/Plots//GEx_MC_enrich_MSigDB.pdf", width = 8, height = 4)
p6
dev.off()

#only PDO165
PDO165_clustGS <- clustGS_robust[, grep("ob165", colnames(clustGS_robust))] %>% as.data.frame %>% as.list
PDO165_clustGS_entrez <- lapply(PDO165_clustGS, function(x) AnnotationDbi::select(org.Hs.eg.db, keys = x, keytype = "SYMBOL", columns = c("ENTREZID")) %>% .$ENTREZID)
PDO165_clustGS_msigDB <- lapply(PDO165_clustGS_entrez, function(x) clusterProfiler::enricher(x, TERM2GENE = msig_hallmark, pAdjustMethod = "fdr", pvalueCutoff = 1, qvalueCutoff = 0.3))
p7 <- lapply(PDO165_clustGS_msigDB, function(x) mutate(x, qscore = -log(p.adjust, base=10)) %>% barplot(x="qscore") + ArchR::theme_ArchR())
pdf("output/Plots//GEx_PDO165_MSigDB.pdf", width = 8, height = 4)
p7
dev.off()

#PDO specific genes
PDO_clustGS <- clustGS_robust_ordered[,c(24:38)]
PDO_clustGS <- PDO_clustGS[,order(colnames(PDO_clustGS))]
PDO_clustGS <- as.data.frame(PDO_clustGS) %>% as.list
PDO_clustGS_entrez <- lapply(PDO_clustGS, function(x) AnnotationDbi::select(org.Hs.eg.db, keys = x, keytype = "SYMBOL", columns = c("ENTREZID")) %>% .$ENTREZID)
PDO_clustGS_msigDB <- lapply(PDO_clustGS_entrez, function(x) clusterProfiler::enricher(x, TERM2GENE = msig_hallmark, pAdjustMethod = "fdr", pvalueCutoff = 1, qvalueCutoff = 0.3))
p8 <- lapply(PDO_clustGS_msigDB, function(x) mutate(x, qscore = -log(p.adjust, base=10)) %>% barplot(x="qscore") + ArchR::theme_ArchR())
pdf("output/Plots//GEx_PDOspecific_MSigDB.pdf", width = 8, height = 4)
p8
dev.off()

#-------------- Meta-clusters annotation --------------#
p9 <- lapply(names(PDO_clustGS), function(x){
  gene <- PDO_clustGS[[x]]
  df <- data.frame(colMeans(GetAssayData(seu, slot = "data")[gene,]), seu$Sample)
  colnames(df) <- c("value", "sample")
  p <- ggplot(df, aes(x = value, y = sample, fill = sample)) + geom_density_ridges(alpha = 0.5) +
    scale_fill_manual(values = sample_cols) + ArchR::theme_ArchR() + ggtitle(x)
  return(p)
})
pdf("output/Plots/GEx_RidgePlots_MeansPDOclustGS.pdf", width = 4, height = 5)
p9
dev.off()

markers <- c("CHI3L1", "CST3", "VIM", "CYP4Z1")
p10 <- lapply(markers, function(x){
  out <- lapply(names(seu_sample), function(i){
    seu_tmp <- seu_sample[[i]]
    p <- FeaturePlot(seu_tmp, features = x, order = F, pt.size = 0.5) + ggtitle(paste0(i, ": ", x)) +
      scale_colour_gradientn(colours = viridis(256, option = "D")) + ArchR::theme_ArchR()
    return(p)
  })
  return(out)
})
pdf("output/Plots/GEx_UMAPoverlay_markerGenes_v1.pdf", width = 5, height = 6)
p10
dev.off()

p10_2 <- lapply(names(seu_sample), function(i){
    seu_tmp <- seu_sample[[i]]
    p <- FeaturePlot(seu_tmp, features = "CCND1", order = F, pt.size = 0.5) + ggtitle(paste0(i, ": CCND1")) +
      scale_colour_gradientn(colours = viridis(256, option = "D")) + ArchR::theme_ArchR()
    return(p)
  })
pdf("output/Plots/GEx_UMAPoverlay_CCND1.pdf", width = 5, height = 6)
p10_2
dev.off()

df <- data.frame(GetAssayData(seu, slot = "data")["CCND1",], seu$Sample)
colnames(df) <- c("value", "sample")
p10_3 <- ggplot(df, aes(x = value, y = sample, fill = sample)) + geom_density_ridges(alpha = 0.5) +
  scale_fill_manual(values = sample_cols) + ArchR::theme_ArchR() + ggtitle("CCND1")
pdf("output/Plots/GEx_RidgePlots_CCND1.pdf", width = 4, height = 5)
p10_3
dev.off()

MC_mtx <- lapply(seq_along(MC), function(i){
  res <- colnames(clustGS_robust_ordered)[MC[[i]]]
  res <- stringr::str_split(res, "#", simplify = T)[,1] %>% unique %>% sort
  out <- data.frame(PDO = res, Meta = i)
  return(out)
}) %>% do.call(rbind, .) %>% reshape2::dcast(., formula = PDO ~ Meta)
rownames(MC_mtx) <- MC_mtx[,1]
MC_mtx <- MC_mtx[, -1]
MC_mtx[is.na(MC_mtx)] <- 0
MC_mtx[MC_mtx > 1] <- 1

subtype_cols <- c("ER+/HER2-" = "blue", "ER-/HER2+" = "orange", "ER+/HER2+" = "yellow")
sample_cols  <- ArchR::ArchRPalettes$stallion2[c(1:10)] %>% `names<-`(., names(table(seu$Sample)))
col_fun3 <- colorRamp2(c(0,1), c("white", "black"))
ra2 <- rowAnnotation(sample = rownames(MC_mtx), subtype = subtype[rownames(MC_mtx)], col = list(sample = sample_cols, subtype = subtype_cols))
ht2 <- Heatmap(MC_mtx, name = "MetaCluster", col = col_fun3, right_annotation = ra2,
               cluster_rows = F, cluster_columns = F, row_names_gp = gpar(fontsize = 8), rect_gp = gpar(col = "gray", lwd = 2))
p12 <- draw(ht2)
pdf("output/Plots//GEx_MetaCluster_Overlap.pdf", width = 5, height = 5)
p12
dev.off()

#-------------- QC --------------#
Idents(seu) <- "Sample"
p11 <- VlnPlot(seu, features = c("nFeature_RNA", "percent.mt"), ncol = 2, pt.size = 0, cols = sample_cols)
mean(seu$nFeature_RNA)
mean(seu$nCount_RNA)
pdf("output/Plots/GEx_QC_vlnPlots.pdf", width = 6, height = 3)
p11
dev.off()

sample_cols  <- ArchR::ArchRPalettes$stallion2[c(1:10)] %>% `names<-`(., names(table(seu$Sample)))
p12 <- DimPlot(seu, reduction = "umap", label = F, group.by = "Sample", cols = sample_cols) + ArchR::theme_ArchR()
p13 <- FeaturePlot(seu, features = "EPCAM", order = T) + scale_colour_gradientn(colours = viridis(256, option = "A")) + ArchR::theme_ArchR()
p14 <- FeaturePlot(seu, features = "KRT8", order = T) + scale_colour_gradientn(colours = viridis(256, option = "A")) + ArchR::theme_ArchR()
p15 <- FeaturePlot(seu, features = "KRT5", order = T) + scale_colour_gradientn(colours = viridis(256, option = "A")) + ArchR::theme_ArchR()

df <- data.frame(GetAssayData(seu_sample$ob155, slot = "data")["CST3",], paste0("GEx_", seu_sample$ob155$Clusters))
colnames(df) <- c("value", "clusters")
p16 <- ggplot(df, aes(x = value, y = clusters, fill = clusters)) + geom_density_ridges(alpha = 0.5) +
  scale_fill_manual(values = cluster_cols) + ArchR::theme_ArchR() + ggtitle("CST3")

pdf("output/Plots/GEx_UMAP_all.pdf", width = 5, height = 6)
p12
p13
p14
p15
dev.off()
