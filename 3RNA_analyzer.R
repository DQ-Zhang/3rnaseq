#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


# Check arguments ----------------------------------------------------------
if(length(args) < 3) {
  stop("
Required arguments:
1) Count matrix file (TSV)
2) Metadata file (CSV)
3) Species (h/m/r/custom) or orgDB name
4) [Optional] Output directory (default: current)", call.=FALSE)
}

# Initialize environment ---------------------------------------------------
counts_file <- args[1]
metadata_file <- args[2]
species <- tolower(args[3])
output_dir <- ifelse(length(args)>3, args[4], getwd())

dir.create(file.path(output_dir, "results"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(output_dir, "plots"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(output_dir, "rdata"), recursive=TRUE, showWarnings=FALSE)

# Setup logging
log_file <- file.path(output_dir, "analysis.log")
sink(log_file, split=TRUE)
on.exit(sink())

# Species configuration ----------------------------------------------------
orgdb <- NULL
valid_species <- c("h", "m", "r")

if(species %in% valid_species) {
  orgdb <- switch(species,
                 "h" = "org.Hs.eg.db",
                 "m" = "org.Mm.eg.db",
                 "r" = "org.Rn.eg.db")
} else {
  orgdb <- species
}

if(!require(orgdb, character.only=TRUE, quietly=TRUE)) {
  if(species %in% valid_species) {
    BiocManager::install(orgdb)
  } else {
    stop("Custom organism package ", orgdb, " not found!")
  }
}
suppressPackageStartupMessages(library(orgdb, character.only=TRUE))

#Load libraries
required_packages <- c("dplyr","data.table","DESeq2","limma","corrplot",
                       "ggplot2","pheatmap","EnhancedVolcano","clusterProfiler",
                       "enrichplot","ggsci","tidyr","stringr","RColorBrewer","gridExtra","grid","svglite",
                       "GSVA","ggridges")
suppressPackageStartupMessages({
  for(pkg in required_packages) {
    library(pkg, character.only = TRUE)
  }
})

# Load local databases
tryCatch({
  MsigDb <- readRDS(file.path("localdata/Msigdb", paste0(species, ".m_df.RDS"))) #2026.01.02
  collectri <- readRDS(file.path("localdata/TF", 
                                paste0(c("human","mouse","rat")[match(species, c("h","m","r"))], 
                                "_collectri.RDS")))[,c(3,4)]
}, error=function(e) warning("Local database loading failed: ", e$message))

# Data loading -------------------------------------------------------------
message("\n[1/8] Loading data...")
counts <- as.matrix(read.delim(counts_file, row.names=1, check.names=FALSE))
metadata <- read.csv(metadata_file, row.names=1)

if(!all(rownames(metadata) %in% colnames(counts)))
  stop("Metadata samples don't match count matrix columns")
counts <- counts[,rownames(metadata)]

#DESeq2 analysis -----------------------------------------------------------------
message("\n[2/8] Generating DESeq2 object...")
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata, 
  design=if("biorep" %in% colnames(metadata)) ~biorep + group else ~group
)%>% DESeq()

normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file.path(output_dir, "results/normalized_counts.csv"))

# QC Plots -----------------------------------------------------------------
message("\n[3/8] Generating QC plots...")
try({
  vsd <- varianceStabilizingTransformation(dds)
  
  pdf(file.path(output_dir, "plots/QC_plots.pdf"), width=11, height=8.5)
  
  # Sample correlation heatmap
  sampleDists <- dist(t(assay(vsd)))
  pheatmap(as.matrix(sampleDists),
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           main = "Sample Distance Matrix",
           color = colorRampPalette(brewer.pal(9, "Blues"))(255))
  
  # Count distribution
  plot_data <- log10(counts[rowMeans(counts)>1,])
  print(ggplot(reshape2::melt(plot_data), aes(x=value)) + 
          geom_density(fill="steelblue", alpha=0.7) +
          labs(x="log10(counts+1)", y="Density", 
               title="Gene Expression Distribution"))
  dev.off()

  print("creating PCA plot")
  dds_norm <- vst(dds)
  pdf(file.path(output_dir, "plots/PCA_plots_1.pdf"), width=5, height=4)
  if("biorep" %in% colnames(metadata)){
   PCA <- plotPCA(dds_norm,intgroup=c("group","biorep"),returnData=T)
    PCAplot <- ggplot(PCA,aes(x=PC1,y=PC2,color = group, shape = biorep))+geom_point()
    print(PCAplot)
    PCAplot2 <- ggplot(PCA,aes(x=PC1,y=PC2,color = group, shape = biorep))+geom_point()+geom_label_repel(aes(label = name))
    print(PCAplot2)
    assay(dds_norm) <- limma::removeBatchEffect(assay(dds_norm), dds_norm$biorep)
    PCA <- plotPCA(dds_norm,intgroup=c("group","biorep"),returnData=T)
    PCAplot <- ggplot(PCA,aes(x=PC1,y=PC2,color = group.1,shape = biorep))+geom_point()
    print(PCAplot)
    PCAplot2 <- ggplot(PCA,aes(x=PC1,y=PC2,color = group, shape = biorep))+geom_point()+geom_label_repel(aes(label = name))
    print(PCAplot2)
    print(plotPCA(dds_norm,intgroup=c("group","biorep")))
  }else{
    PCA <- plotPCA(dds_norm,intgroup=c("group"))
    print(PCA)
    PCA2 <- PCA + geom_label(aes(label = name))
    print(PCA2)
  }
  dev.off()
  print("PCA plot created")

  # PCA with biorep handling
  pdf(file.path(output_dir, "plots/PCA_plots_2.pdf"), width=5, height=4)
  pca_data <- plotPCA(vsd, intgroup=if("biorep" %in% colnames(metadata)) c("group","biorep") else "group", 
                     returnData=TRUE)
  
  if ("biorep" %in% colnames(metadata)) {
    pca_plot1 <- ggplot(pca_data, aes(PC1, PC2, color = group, shape = biorep)) +
      labs(shape = "secondary factor", color ="Group") +
      geom_point(size = 3) +
      theme_minimal(base_size = 14) +
      scale_color_brewer(palette = "Set1") +
      ggtitle("Principal Component Analysis")
  } else {
    pca_plot1 <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
      labs(color ="Group") +
      geom_point(size = 3) +
      theme_minimal(base_size = 14) +
      scale_color_brewer(palette = "Set1") +
      ggtitle("Principal Component Analysis")
  }

  pca_plot2 <- pca_plot1 +
    geom_label_repel(aes(label = name))

  
  if("biorep" %in% colnames(metadata)) {
    assay(vsd) <- removeBatchEffect(assay(vsd), vsd$biorep)
    pca_corrected <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
    pca_plot3 <- 
      ggplot(pca_corrected, aes(PC1, PC2, color = group)) +
      geom_point(size = 3) +
      theme_minimal() +
      ggtitle("Batch Corrected PCA")
    pca_plot4 <- pca_plot3 +
      geom_label_repel(aes(label = name))
  }
  print(pca_plot1)
  print(pca_plot2)
  if(exists("pca_plot3")){
    print(pca_plot3)
    print(pca_plot4)
    }
  dev.off()
})



# Function to run enrichment for a gene set
run_enrichment <- function(genes, prefix, out_dir) {
  if (length(genes) < 10) {
    message(paste0("Gene set ", prefix, " has too few genes for enrichment analysis"))
    return(NULL)
  }
  
  # Main function to run and save enrichment results
  run_and_save_enrichment <- function(enrich_result, analysis_name, title) {
    if (is.null(enrich_result) || nrow(enrich_result) == 0) {
      message(paste("No results for", analysis_name, "- skipping plots"))
      return(FALSE)
    }
    
    # Create specific directory for this analysis
    analysis_dir <- file.path(out_dir, paste0(prefix, "-", analysis_name))
    dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save CSV results
    write.csv(enrich_result, file.path(analysis_dir, paste0(prefix, "_", analysis_name, ".csv")))
    
    # Create dotplot
    p <- dotplot(enrich_result, showCategory = min(nrow(enrich_result), 15)) + 
      ggtitle(title)
    
    # Save in multiple formats
    try({
      # PDF
      pdf(file.path(analysis_dir, paste0(prefix, "_", analysis_name, ".pdf")), 
          width = 10, height = 8)
      print(p)
      dev.off()
      
      # TIFF
      tiff(file.path(analysis_dir, paste0(prefix, "_", analysis_name, ".tiff")), 
           width = 2000, height = 1600, res = 300)
      print(p)
      dev.off()
      
      # PNG
      png(file.path(analysis_dir, paste0(prefix, "_", analysis_name, ".png")), 
          width = 2000, height = 1600, res = 300)
      print(p)
      dev.off()
      
      #SVG
      ggsave(file.path(analysis_dir, paste0(prefix, "_", analysis_name, ".svg")), plot = p, 
             width = 10, height = 8, unit = "in")
    }, silent = TRUE)
    
    return(TRUE)
  }
  
  tryCatch({
    # GO:BP
    message(paste("Running GOBP enrichment analysis for", prefix))
    go_bp <- enrichGO(genes, orgdb, keyType = "SYMBOL", ont = "BP")
    run_and_save_enrichment(go_bp, "GOBP", 
                            paste(prefix, "ORA: Gene ontology - Biological process"))
    
    # GO:MF
    message(paste("Running GOMF enrichment analysis for", prefix))
    go_mf <- enrichGO(genes, orgdb, keyType = "SYMBOL", ont = "MF")
    run_and_save_enrichment(go_mf, "GOMF", 
                            paste(prefix, "ORA: Gene ontology - Molecular function"))
    
    # GO:CC
    message(paste("Running GOCC enrichment analysis for", prefix))
    go_cc <- enrichGO(genes, orgdb, keyType = "SYMBOL", ont = "CC")
    run_and_save_enrichment(go_cc, "GOCC", 
                            paste(prefix, "ORA: Gene ontology - Cellular component"))
    
    # GO:ALL
    message(paste("Running GOALL enrichment analysis for", prefix))
    go_all <- enrichGO(genes, orgdb, keyType = "SYMBOL", ont = "ALL")
    run_and_save_enrichment(go_all, "GOALL", 
                            paste(prefix, "ORA: Gene ontology - All"))
    
    # MsigDB-hallmark
    if (exists("MsigDb")) {
      message(paste("Running MsigDB hallmark enrichment analysis for", prefix))
      T2G_hallmark <- dplyr::filter(MsigDb, gs_collection == "H") %>% 
        dplyr::distinct(gs_name, gene_symbol)
      msig_hallmark <- enricher(genes, TERM2GENE = T2G_hallmark)
      run_and_save_enrichment(msig_hallmark, "MsigDB-hallmark", 
                              paste(prefix, "ORA: hallmark genesets"))
      
      # MsigDB-CGP
      message(paste("Running MsigDB CGP enrichment analysis for", prefix))
      T2G_cgp <- dplyr::filter(MsigDb, gs_subcollection == "CGP") %>% 
        dplyr::distinct(gs_name, gene_symbol)
      msig_cgp <- enricher(genes, TERM2GENE = T2G_cgp)
      run_and_save_enrichment(msig_cgp, "MsigDB-CGP", 
                              paste(prefix, "ORA: chemical and genetic perturbation genesets"))
      
      # MsigDB-KEGG_MEDICUS
      message(paste("Running MsigDB KEGG medicus enrichment analysis for", prefix))
      T2G_kegg <- dplyr::filter(MsigDb, gs_subcollection == "CP:KEGG_MEDICUS") %>% 
        dplyr::distinct(gs_name, gene_symbol)
      msig_kegg <- enricher(genes, TERM2GENE = T2G_kegg)
      run_and_save_enrichment(msig_kegg, "MsigDB-KEGGmedicus", 
                              paste(prefix, "ORA: KEGG medicus pathway genesets"))
      
      # MsigDB-REACTOME
      message(paste("Running MsigDB REACTOME enrichment analysis for", prefix))
      T2G_reactome <- dplyr::filter(MsigDb, gs_subcollection == "CP:REACTOME") %>% 
        dplyr::distinct(gs_name, gene_symbol)
      msig_reactome <- enricher(genes, TERM2GENE = T2G_reactome)
      run_and_save_enrichment(msig_reactome, "MsigDB-REACTOME", 
                              paste(prefix, "ORA: REACTOME pathway genesets"))
      
      # MsigDB-miRNA
      message(paste("Running MsigDB miRNA enrichment analysis for", prefix))
      T2G_mirna <- dplyr::filter(MsigDb, gs_subcollection == "MIR:MIRDB") %>% 
        dplyr::distinct(gs_name, gene_symbol)
      msig_mirna <- enricher(genes, TERM2GENE = T2G_mirna)
      run_and_save_enrichment(msig_mirna, "MsigDB-miRNA", 
                              paste(prefix, "ORA: microRNA target genesets"))
      
      # MsigDB-ONCOsig
      message(paste("Running MsigDB ONCOsig enrichment analysis for", prefix))
      T2G_onco <- dplyr::filter(MsigDb, gs_collection == "C6") %>% 
        dplyr::distinct(gs_name, gene_symbol)
      msig_onco <- enricher(genes, TERM2GENE = T2G_onco)
      run_and_save_enrichment(msig_onco, "MsigDB-ONCOsig", 
                              paste(prefix, "ORA: oncogenic signature genesets"))
      
      # MsigDB-IMMUNESIGDB
      message(paste("Running MsigDB IMMUNESIGDB enrichment analysis for", prefix))
      T2G_immune <- dplyr::filter(MsigDb, gs_subcollection == "IMMUNESIGDB") %>% 
        dplyr::distinct(gs_name, gene_symbol)
      msig_immune <- enricher(genes, TERM2GENE = T2G_immune)
      run_and_save_enrichment(msig_immune, "MsigDB-IMMUNESIGDB", 
                              paste(prefix, "ORA: ImmunesigDB genesets"))
      
      # MsigDB-ALL
      message(paste("Running MsigDB ALL enrichment analysis for", prefix))
      T2G_all <- dplyr::distinct(MsigDb, gs_name, gene_symbol)
      msig_all <- enricher(genes, TERM2GENE = T2G_all)
      run_and_save_enrichment(msig_all, "MsigDB-ALL", 
                              paste(prefix, "ORA: All MsigDB genesets"))
    }
    
    # TF
    if (exists("collectri")) {
      message(paste("Running TF enrichment analysis for", prefix))
      tf <- enricher(genes, TERM2GENE = collectri)
      run_and_save_enrichment(tf, "TF", 
                              paste(prefix, "ORA: transcription factor target genesets"))
    }
    
    message(paste("Enrichment analysis completed for", prefix))
    
  }, error = function(e) {
    message(paste("Enrichment analysis failed for", prefix, ":", e$message))
  })
}

gsea_plots <- function(gsea_result, prefix, out_dir) {
  
  if (is.null(gsea_result) || nrow(gsea_result) == 0) {
    message(paste("No GSEA results for", prefix, "- skipping plots"))
    return(NULL)
  }
  
  message(paste("Creating GSEA plots for", prefix, "with", nrow(gsea_result), "gene sets"))
  
  # Create directory for this prefix's GSEA results
  gsea_dir <- file.path(out_dir, paste0(prefix, "-GSEA"))
  dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(gsea_result, file.path(gsea_dir, paste0(prefix, ".csv")))
  
  dot_p <- try(dotplot(gsea_result, showCategory=min(nrow(gsea_result),15)))
  cnet_p <- try(cnetplot(gsea_result, categorySize="pvalue"))
  gsea_p <- try(gseaplot2(gsea_result, geneSetID=1:min(nrow(gsea_result),10)))
  
  pdf(file.path(gsea_dir, paste0(prefix, "_plots.pdf")), width=11, height=8.5)
  try(print(dot_p))
  try(print(cnet_p))
  try(print(gsea_p))
  for (j in 1:min(nrow(gsea_result),100)){
    term <- gsea_result$Description[j]
    print(paste("plotting GSEA plot for: ",term, sep=""))
    print(gseaplot(gsea_result, by = "all", title = term, geneSetID = j))
  }
  dev.off()
  
  #save dotplot in multiple formats
  tiff(file.path(gsea_dir, paste0(prefix, "_dotplot.tiff")), 
       width = 11 * 300, height = 8.5 * 300, res = 300)
  print(dot_p)
  dev.off()
  
  png(file.path(gsea_dir, paste0(prefix, "_dotplot.png")), 
      width = 11 * 150, height = 8.5 * 150, res = 150)
  print(dot_p)
  dev.off()
  
  ggsave(file.path(gsea_dir, paste0(prefix, "_dotplot.svg")), plot = dot_p, 
                   width = 11, height = 8.5, unit = "in")
  
  #save network plot in multiple formats
  tiff(file.path(gsea_dir, paste0(prefix, "_cnetplot.tiff")), 
       width = 11 * 600, height = 8.5 * 600, res = 600)
  print(cnet_p)
  dev.off()
  
  png(file.path(gsea_dir, paste0(prefix, "_cnetplot.png")), 
      width = 11 * 600, height = 8.5 * 600, res = 600)
  print(cnet_p)
  dev.off()
  
  ggsave(file.path(gsea_dir, paste0(prefix, "_cnetplot.svg")), plot = cnet_p, 
                   width = 11, height = 8.5, unit = "in")
  
  #save gsea plot in multiple formats
  tiff(file.path(gsea_dir, paste0(prefix, "_gseaplot2.tiff")), 
       width = 11 * 300, height = 8.5 * 300, res = 300)
  print(gsea_p)
  dev.off()
  
  png(file.path(gsea_dir, paste0(prefix, "_gseaplot2.png")), 
      width = 11 * 300, height = 8.5 * 300, res = 300)
  print(gsea_p)
  dev.off()
  
  ggsave(file.path(gsea_dir, paste0(prefix, "_gseaplot2.svg")), plot = gsea_p, 
                   width = 11, height = 8.5, unit = "in")
}

# Global DEG Analysis & Heatmaps -------------------------------------------------
# Reinstated 2026/01/02
message("\n[4/8] Performing global analysis...")

# Global DEG analysis with batch correction if needed
if("biorep" %in% colnames(metadata) && length(unique(metadata[,"biorep"])) > 1 ) {
  corrected.dds <- dds
  assay(corrected.dds) <- limma::removeBatchEffect(assay(corrected.dds), corrected.dds$biorep)
  GlobalDEG <- results(corrected.dds)
} else {
  GlobalDEG <- results(dds)
}

message("Summary of global differential gene analysis:")
summary(results(dds))

GlobalDEG <- as.data.frame(GlobalDEG)

# Save global DEG results
write.csv(GlobalDEG, file.path(output_dir, "results/global_DEG_results.csv"))

GlobalDEG <- GlobalDEG[GlobalDEG$padj < 0.05 & !is.na(GlobalDEG$padj),]
GlobalDEG <- GlobalDEG[order(GlobalDEG$padj, decreasing = FALSE),]
GlobalDEG <- na.omit(GlobalDEG)

# Save global significant DEG results
write.csv(GlobalDEG, file.path(output_dir, "results/global_DEG_results_significant.csv"))

# Function to determine optimal number of clusters based on data complexity
determine_clusters <- function(n_genes, n_groups, max_clusters = 8) {
  # Base clusters on number of groups
  base_clusters <- min(max_clusters, 2^(n_groups - 1))
  
  # Adjust based on gene count complexity
  if (n_genes < 100) {
    clusters <- min(3, base_clusters)
  } else if (n_genes < 500) {
    clusters <- min(4, base_clusters)
  } else if (n_genes < 1000) {
    clusters <- min(6, base_clusters)
  } else {
    clusters <- base_clusters
  }
  
  # Ensure at least 2 clusters
  return(max(2, clusters))
}

# Create heatmaps for top DEGs at different cutoffs
top_n_list <- c(20, 50, 100)
n_groups <- length(unique(metadata$group))

# Heatmap plotting function
plot_heatmap <- function(data_matrix, title, show_rownames = FALSE, cluster_rows = TRUE, 
                         annotation_col = metadata, annotation_row = NULL) {
  # Create annotation colors
  ann_colors = list()
  
  # Group colors
  group_colors <- brewer.pal(min(8, n_groups), "Set1")[1:n_groups]
  names(group_colors) <- unique(metadata$group)
  ann_colors[["group"]] <- group_colors
  
  # Biorep colors if present
  if("biorep" %in% colnames(metadata)) {
    biorep_colors <- brewer.pal(min(8, length(unique(metadata$biorep))), "Set2")[1:length(unique(metadata$biorep))]
    names(biorep_colors) <- unique(metadata$biorep)
    ann_colors[["biorep"]] <- biorep_colors
  }
  
  pheatmap(
    data_matrix,
    cluster_rows = cluster_rows,
    cluster_cols = TRUE,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = ann_colors,
    show_rownames = show_rownames,
    show_colnames = TRUE,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    fontsize_row = ifelse(show_rownames, 6, 0),
    fontsize_col = 8,
    main = title,
    border_color = NA,
    silent = FALSE
  )
}

# Plot heatmaps for top N DEGs
for (top_n in top_n_list) {
  if (nrow(GlobalDEG) >= top_n) {
    top_genes <- rownames(GlobalDEG)[1:top_n]
    hmaptable <- normalized_counts[top_genes, ]
    
    # Determine if we should show gene names
    show_names <- top_n <= 20
    
    pdf(file.path(output_dir, paste0("plots/heatmap_top", top_n, "_DEGs.pdf")), 
        width = 10, height = max(6, 0.1 * top_n))
    
    hm <- plot_heatmap(
      data_matrix = hmaptable,
      title = paste("Top", top_n, "Differentially Expressed Genes"),
      show_rownames = show_names,
      cluster_rows = TRUE
    )
    
    dev.off()
  }
}

# Global heatmap with clustering (all significant DEGs)
if (nrow(GlobalDEG) >= 10) {
  hmaptable_all <- normalized_counts[rownames(GlobalDEG), ]
  
  # Determine optimal number of clusters
  n_clusts <- determine_clusters(nrow(GlobalDEG), n_groups)
  message(paste("Number of clusters for heatmap:", n_clusts))
  
  # Heatmap without clusters
  pdf(file.path(output_dir, "plots/heatmap_all_DEGs_noclust.pdf"), 
      width = 10, height = max(8, 0.05 * nrow(GlobalDEG)))
  hm_all <- plot_heatmap(
    data_matrix = hmaptable_all,
    title = "All Significant DEGs",
    show_rownames = FALSE,
    cluster_rows = TRUE
  )
  dev.off()
  
  # Get clusters from hierarchical clustering
  row_dend <- as.dendrogram(hm_all$tree_row)
  clusts <- cutree(hm_all$tree_row, k = n_clusts)
  clusts_df <- data.frame(cluster = as.character(clusts))
  rownames(clusts_df) <- names(clusts)
  
  # Heatmap with cluster annotation
  pdf(file.path(output_dir, "plots/heatmap_all_DEGs_clustered.pdf"), 
      width = 11, height = max(8, 0.05 * nrow(GlobalDEG)))
  hm_clustered <- plot_heatmap(
    data_matrix = hmaptable_all,
    title = "All Significant DEGs (Clustered)",
    show_rownames = FALSE,
    cluster_rows = TRUE,
    annotation_row = clusts_df
  )
  dev.off()
  
  # Save cluster assignments
  write.csv(clusts_df, file.path(output_dir, "results/gene_cluster_assignments.csv"))
  
  # Enrichment analysis for each cluster
  message("Performing enrichment analysis for gene clusters...")
  clust_dir <- file.path(output_dir, "results/gene_clusters_enrichment")
  dir.create(clust_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  
  # Process each cluster
  for (clust_num in 1:n_clusts) {
    clust_name <- paste("cluster", clust_num, sep = "_")
    cluster_genes <- names(clusts[clusts == clust_num])
    
    message(paste("Processing", clust_name, "with", length(cluster_genes), "genes"))
    
    # Save gene list
    write.table(cluster_genes, 
                file = file.path(clust_dir, paste0(clust_name, "_genes.txt")),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # Run enrichment
    clust_n_dir <- file.path(clust_dir, clust_name)
    dir.create(clust_n_dir)
    run_enrichment(cluster_genes, clust_name, clust_n_dir)
  }
}

# Pairwise Comparisons -----------------------------------------------------
message("\n[5/8] Performing pairwise comparisons...")
groups <- unique(metadata$group)
comps <- combn(groups, 2, simplify=FALSE)

#Revised 2026/01/02
plot_top_genes <- function(genes, title_suffix) {
  # Ensure we don't have more than 12 genes for 4x3 grid
  genes <- head(genes, 12)
  
  plots <- lapply(genes, function(g){
    plot_data <- plotCounts(dds, gene=g, intgroup="group", returnData=TRUE)
    # Convert group to factor with levels in the order they appear in metadata
    plot_data$group <- factor(plot_data$group, levels = group_levels)
    
    ggplot(plot_data, aes(x=group, y=count, fill=group)) +
      geom_boxplot() +
      geom_jitter(width=0.2) +
      ggtitle(g) +
      theme_minimal() +
      scale_x_discrete(limits = group_levels) +  # Ensure x-axis order matches metadata
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Use grid.arrange instead of marrangeGrob
  do.call("grid.arrange", c(plots, nrow=4, ncol=3, top=title_suffix))
}

for(comp in comps) {
  tryCatch({
    contrast <- c("group", comp[2], comp[1])
    res_name <- paste0(contrast[3], "_vs_", contrast[2])
    message(paste0("Performing pairwise comparisons for:", res_name))
    res_dir <- file.path(output_dir, "results", res_name)
    dir.create(file.path(res_dir, "Overrepresentation analysis"), recursive=TRUE, showWarnings=FALSE)
    dir.create(file.path(res_dir, "GSEA"), recursive=TRUE, showWarnings=FALSE)
    dir.create(file.path(res_dir, "plots"), recursive=TRUE, showWarnings=FALSE)
    
    # Differential expression
    message(paste0("Performing differential expression analysis for:", res_name))
    res <- results(dds, contrast=contrast) %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("gene") %>%
      filter(baseMean >= 5) %>% #not verified 260208
      arrange(padj)
    write.csv(res, file.path(res_dir, "DEG_results.csv"))
    
    # Volcano plot
    message(paste0("Volcano plot for:", res_name))
    volcano <- tryCatch(
      EnhancedVolcano(res,
                      lab=res$gene,
                      x="log2FoldChange",
                      y="padj",
                      title=res_name,
                      pCutoff=0.05,
                      FCcutoff=1)
    )
    ggsave(file.path(res_dir, "plots/volcano.pdf"), volcano, width=8, height=10)
    
    # Top DEG expression plots
    sig_genes <- res %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)
    up_genes <- sig_genes %>% filter(log2FoldChange > 0)
    top_up_genes <- head(up_genes$gene, min(nrow(up_genes), 12))  # Limit to 12 for 4x3 grid
    down_genes <- sig_genes %>% filter(log2FoldChange < 0)
    top_down_genes <- head(down_genes$gene, min(nrow(down_genes), 12))  # Limit to 12 for 4x3 grid
    write.csv(sig_genes, file.path(res_dir, "sig_DEG_results.csv"))
    write.csv(up_genes, file.path(res_dir, "sig_up_DEG_results.csv"))
    write.csv(down_genes, file.path(res_dir, "sig_down_DEG_results.csv"))
    
    group_levels <- unique(metadata$group)
    #DEG barplots
    pdf(file.path(res_dir, "plots/top_DEGs_barplots.pdf"), width=11, height=8.5)
    if(length(top_up_genes) > 0) {
      tryCatch({
        print(plot_top_genes(top_up_genes, "Top Upregulated Genes"))
      }, error = function(e) {
        message("Error in upregulated genes plot: ", e$message)
      })
    }
    if(length(top_down_genes) > 0) {
      tryCatch({
        print(plot_top_genes(top_down_genes, "Top Downregulated Genes"))
      }, error = function(e) {
        message("Error in downregulated genes plot: ", e$message)
      })
    }
    dev.off()
    
    # Overrepresentation analysis
    
    if(nrow(sig_genes)>=10) {
      ORA_out_dir <- file.path(res_dir, "Overrepresentation analysis")
      run_enrichment(sig_genes$gene, "all_DEGs", ORA_out_dir)
      run_enrichment(up_genes$gene, "upregulated", ORA_out_dir)
      run_enrichment(down_genes$gene, "downregulated", ORA_out_dir)
    }
    
    # GSEA analysis
    if(nrow(res)>=50) {
      gene_list <- res$log2FoldChange
      names(gene_list) <- res$gene
      gene_list <- na.omit(gene_list) 
      gene_list <- sort(gene_list,decreasing = TRUE)
      
      gsea_out_dir <- file.path(res_dir, "GSEA")
      dir.create(gsea_out_dir)
      
      try({
        gsea_go <- gseGO(gene_list, ont = "BP", OrgDb=orgdb, keyType="SYMBOL")
        gsea_plots(gsea_go, "GOBP", gsea_out_dir)
      })
      
      try({
        gsea_go <- gseGO(gene_list, ont = "MF", OrgDb=orgdb, keyType="SYMBOL")
        gsea_plots(gsea_go, "GOMF", gsea_out_dir)
      })
      
      try({
        gsea_go <- gseGO(gene_list, ont = "CC", OrgDb=orgdb, keyType="SYMBOL")
        gsea_plots(gsea_go, "GOCC", gsea_out_dir)
      })
      
      try({
        gsea_go <- gseGO(gene_list, ont = "ALL", OrgDb=orgdb, keyType="SYMBOL")
        gsea_plots(gsea_go, "GOALL", gsea_out_dir)
      })
      
      try({
        T2G <- dplyr::filter(MsigDb, gs_collection == "H") %>% dplyr::distinct(gs_name, gene_symbol)
        gsea_msig <- GSEA(gene_list, TERM2GENE=T2G)
        gsea_plots(gsea_msig, "MsigDB-hallmark", gsea_out_dir)
      })
      
      try({
        T2G <- dplyr::filter(MsigDb, gs_subcollection == "CGP") %>% dplyr::distinct(gs_name, gene_symbol)
        gsea_msig <- GSEA(gene_list, TERM2GENE=T2G)
        gsea_plots(gsea_msig, "MsigDB-CGP", gsea_out_dir)
      })
      
      try({
        T2G <- dplyr::filter(MsigDb, gs_subcollection == "CP:KEGG_MEDICUS") %>% dplyr::distinct(gs_name, gene_symbol)
        gsea_msig <- GSEA(gene_list, TERM2GENE=T2G)
        gsea_plots(gsea_msig, "MsigDB-KEGGmedicus", gsea_out_dir)
      })
      
      try({
        T2G <- dplyr::filter(MsigDb, gs_subcollection == "CP:REACTOME") %>% dplyr::distinct(gs_name, gene_symbol)
        gsea_msig <- GSEA(gene_list, TERM2GENE=T2G)
        gsea_plots(gsea_msig, "MsigDB-REACTOME", gsea_out_dir)
      })
      
      try({
        T2G <- dplyr::filter(MsigDb, gs_subcollection == "MIR:MIRDB") %>% dplyr::distinct(gs_name, gene_symbol)
        gsea_msig <- GSEA(gene_list, TERM2GENE=T2G, gsea_out_dir)
        gsea_plots(gsea_msig, "MsigDB-miRNA")
      })
      
      try({
        T2G <- dplyr::filter(MsigDb, gs_collection == "C6") %>% dplyr::distinct(gs_name, gene_symbol)
        gsea_msig <- GSEA(gene_list, TERM2GENE=T2G)
        gsea_plots(gsea_msig, "MsigDB-ONCOsig", gsea_out_dir)
      })
      
      try({
        T2G <- dplyr::filter(MsigDb, gs_subcollection == "IMMUNESIGDB") %>% dplyr::distinct(gs_name, gene_symbol)
        gsea_msig <- GSEA(gene_list, TERM2GENE=T2G)
        gsea_plots(gsea_msig, "MsigDB-IMMUNESIGDB", gsea_out_dir)
      })
      
      try({
        T2G <- dplyr::distinct(MsigDb, gs_name, gene_symbol)
        gsea_msig <- GSEA(gene_list, TERM2GENE=T2G)
        gsea_plots(gsea_msig, "MsigDB-ALL", gsea_out_dir)
      })
        
      try({
        gsea_tf <- GSEA(gene_list, TERM2GENE=collectri)
        gsea_plots(gsea_tf, "TF", gsea_out_dir)
      })
    }
  }, error=function(e) message("Comparison ", res_name, " failed: ", e$message))
}




# Gene Set Variation Analysis (GSVA) -------------------------------------------
message("\n[6/8] Performing Gene Set Variation Analysis (GSVA)...")

# 1. Prepare Input Data
# GSVA requires a normalized, non-negative expression matrix.
# The variance-stabilized transformation (VST) from DESeq2 is suitable[citation:6].
message("  Preparing normalized expression matrix...")
try({
  # Use the variance-stabilized transformation, created earlier in the QC plots section
  # Extract the expression matrix: genes as rows, samples as columns
  expr_matrix <- assay(vst(dds))
  
  # Ensure no missing or infinite values
  expr_matrix <- expr_matrix[is.finite(rowSums(expr_matrix)), ]
}, silent = TRUE)

# 2. Load & Prepare Gene Set Collections from Local Data
message("  Loading and preparing gene set collections from localdata...")
gene_sets <- list() # This list will store all our gene set collections

tryCatch({
  # --- GO Gene Sets from clusterProfiler (using the specified org.db) ---
  message("    Building GO gene sets from clusterProfiler...")
  
  # Alternative: More efficient helper function
  get_GO_genesets <- function(ontology, orgdb_var) {
    # Get the loaded organism database
    OrgDb <- get(orgdb_var)
    
    # Get GO annotations directly from organism database
    go_annot <- AnnotationDbi::select(
      OrgDb,
      keys = keys(OrgDb, keytype = "SYMBOL"),
      columns = c("GO", "ONTOLOGY"),
      keytype = "SYMBOL"
    )
    
    # Remove NAs
    go_annot <- go_annot[!is.na(go_annot$GO) & !is.na(go_annot$ONTOLOGY), ]
    
    # Filter by ontology if needed
    if (ontology != "ALL") {
      go_annot <- go_annot[go_annot$ONTOLOGY == ontology, ]
    }
    
    # Get term names from GO.db
    if (nrow(go_annot) > 0) {
      # Get unique GO IDs
      unique_go_ids <- unique(go_annot$GO)
      
      # Get term names for these GO IDs
      go_term_names <- AnnotationDbi::select(GO.db::GO.db,
                                             keys = unique_go_ids,
                                             columns = "TERM",
                                             keytype = "GOID")
      
      # Merge term names with annotations
      go_annot <- merge(go_annot, go_term_names, by.x = "GO", by.y = "GOID")
      
      # Create descriptive names
      descriptive_names <- paste0(go_annot$GO, ": ", go_annot$TERM)
      
      # Convert to named list
      go_list <- split(go_annot$SYMBOL, descriptive_names)
      
      message(paste("        Mapped", length(unique(go_annot$SYMBOL)), 
                    "genes to", length(go_list), "GO", ontology, "terms"))
      
      return(go_list)
    } else {
      return(list())
    }
  }
  
  # GO: Biological Process (BP)
  message("      Extracting GO Biological Process terms...")
  go_bp_sets <- get_GO_genesets("BP", orgdb)
  if (length(go_bp_sets) > 0) {
    gene_sets[["GO_BP"]] <- go_bp_sets
    message(paste("        Loaded", length(go_bp_sets), "GO BP gene sets."))
  }
  
  # GO: Molecular Function (MF)
  message("      Extracting GO Molecular Function terms...")
  go_mf_sets <- get_GO_genesets("MF", orgdb)
  if (length(go_mf_sets) > 0) {
    gene_sets[["GO_MF"]] <- go_mf_sets
    message(paste("        Loaded", length(go_mf_sets), "GO MF gene sets."))
  }
  
  # GO: Cellular Component (CC)
  message("      Extracting GO Cellular Component terms...")
  go_cc_sets <- get_GO_genesets("CC", orgdb)
  if (length(go_cc_sets) > 0) {
    gene_sets[["GO_CC"]] <- go_cc_sets
    message(paste("        Loaded", length(go_cc_sets), "GO CC gene sets."))
  }
  
  # GO: ALL - Combine all three ontologies
  message("      Combining GO ALL terms...")
  go_all_sets <- c(go_bp_sets, go_mf_sets, go_cc_sets)
  if (length(go_all_sets) > 0) {
    # Remove any potential duplicates if a GO term appears in multiple ontologies
    gene_sets[["GO_ALL"]] <- go_all_sets[!duplicated(names(go_all_sets))]
    message(paste("        Combined GO ALL gene sets:", 
                  length(gene_sets[["GO_ALL"]]), "unique terms"))
  }
  
  # MSigDB Hallmark (H)
  hallmark_df <- dplyr::filter(MsigDb, gs_collection == "H")
  hallmark_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
  gene_sets[["MSigDB_Hallmark"]] <- hallmark_sets
  
  # MSigDB CGP (Chemical & Genetic Perturbations)
  cgp_df <- dplyr::filter(MsigDb, gs_subcollection == "CGP")
  cgp_sets <- split(cgp_df$gene_symbol, cgp_df$gs_name)
  gene_sets[["MSigDB_CGP"]] <- cgp_sets
  
  # MSigDB KEGG_MEDICUS
  keggm_df <- dplyr::filter(MsigDb, gs_subcollection == "CP:KEGG_MEDICUS")
  keggm_sets <- split(keggm_df$gene_symbol, keggm_df$gs_name)
  gene_sets[["MSigDB_KEGG_MEDICUS"]] <- keggm_sets
  
  # MSigDB REACTOME
  reactome_df <- dplyr::filter(MsigDb, gs_subcollection == "CP:REACTOME")
  reactome_sets <- split(reactome_df$gene_symbol, reactome_df$gs_name)
  gene_sets[["MSigDB_REACTOME"]] <- reactome_sets
  
  # MSigDB miRNA Targets (from miRDB)
  mirna_df <- dplyr::filter(MsigDb, gs_subcollection == "MIR:MIRDB")
  mirna_sets <- split(mirna_df$gene_symbol, mirna_df$gs_name)
  gene_sets[["MSigDB_miRNA"]] <- mirna_sets
  
  # MSigDB Oncogenic Signatures (C6)
  onco_df <- dplyr::filter(MsigDb, gs_collection == "C6")
  onco_sets <- split(onco_df$gene_symbol, onco_df$gs_name)
  gene_sets[["MSigDB_ONCOsig"]] <- onco_sets
  
  # MSigDB IMMUNESIGDB
  immune_df <- dplyr::filter(MsigDb, gs_subcollection == "IMMUNESIGDB")
  immune_sets <- split(immune_df$gene_symbol, immune_df$gs_name)
  gene_sets[["MSigDB_IMMUNESIGDB"]] <- immune_sets
  
  # MSigDB ALL (All available sets - use with caution, can be very large)
  # gene_sets[["MSigDB_ALL"]] <- split(MsigDb$gene_symbol, MsigDb$gs_name)
  
  # Transcription Factor Targets from Collectri
  tf_sets <- split(collectri$source_genesymbol, collectri$target_genesymbol)
  gene_sets[["TF_Targets"]] <- tf_sets
  
  # Report loaded collections
  for (gs_name in names(gene_sets)) {
    message(paste("    Loaded", length(gene_sets[[gs_name]]), gs_name, "gene sets."))
  }
}, error = function(e) {
  warning("Failed to load or process local gene set databases: ", e$message)
})

# 3. Run GSVA
message("  Running GSVA enrichment scoring...")
gsva_results <- list()
gsva_scores_all <- NULL

dir.create(file.path(output_dir, "GSVA"))

if (length(gene_sets) > 0 && exists("expr_matrix")) {
  for (gs_name in names(gene_sets)) {
    message(paste("    Analyzing", gs_name, "..."))
    tryCatch({
      # 1. Build the parameter object
      gsvapar <- gsvaParam(expr = expr_matrix,
                           geneSets = gene_sets[[gs_name]])
      # 2. Run GSVA using the parameter object
      gsva_res <- gsva(gsvapar)
      
      gsva_results[[gs_name]] <- gsva_res
      write.csv(gsva_res, file.path(output_dir, "results", paste0("GSVA_scores_", gs_name, ".csv")))
      
      # Combine scores (optional)
      if (is.null(gsva_scores_all)) {
        gsva_scores_all <- gsva_res
      } else {
        gsva_scores_all <- rbind(gsva_scores_all, gsva_res)
      }
      
      message(paste("      Completed GSVA for", gs_name))
      
    }, error = function(e) {
      message(paste("      GSVA failed for", gs_name, ":", e$message))
    })
  }
  
  if (!is.null(gsva_scores_all)) {
    write.csv(gsva_scores_all, file.path(output_dir, "GSVA", "GSVA_scores_COMBINED.csv"))
  }
} else {
  message("  GSVA skipped.")
}


# 4. Identify Differentially Active Pathways
de_pathways_results <- list()
if (length(gsva_results) > 0 && length(unique(metadata$group)) >= 2) {
  message("  Identifying differentially active pathways between groups...")
  
  for (gs_name in names(gsva_results)) {
    gsva_mat <- gsva_results[[gs_name]]
    if (nrow(gsva_mat) < 3) next # Need enough pathways for stats
    
    # ---- BUILD DESIGN MATRIX ----
    # Match the DESeq2 design structure
    # Get sample metadata in same order as GSVA matrix columns
    sample_meta <- metadata[colnames(gsva_mat), , drop = FALSE]
    
    # Check if 'biorep' exists and has multiple levels
    if ("biorep" %in% colnames(sample_meta) && 
        length(unique(sample_meta$biorep)) > 1) {
      # Design with batch correction: ~biorep + group
      design_formula <- ~0 + group + biorep
      message(paste("    Using design: ~group + biorep for", gs_name))
    } else {
      # Simple design without batch: ~group
      design_formula <- ~0 + group
      message(paste("    Using design: ~group for", gs_name))
    }
    
    # Create design matrix
    design <- model.matrix(design_formula, data = sample_meta)
    
    # ---- FIT LINEAR MODEL AND TEST CONTRASTS ----
    fit <- lmFit(gsva_mat, design)
    
    for (comp in comps) {
      contrast_name <- paste(comp[2], "vs", comp[1], sep = "_")
      
      contrast_formula <- paste0("group", comp[2], " - group", comp[1])
      
      # If a group is the reference (no column), adjust
      if (!paste0("group", comp[1]) %in% colnames(design)) {
        contrast_formula <- paste0("group", comp[2])  # comp[1] is reference
      }
      if (!paste0("group", comp[2]) %in% colnames(design)) {
        contrast_formula <- paste0("-group", comp[1])  # comp[2] is reference
      }
      
      dir.create(file.path(output_dir, "GSVA", contrast_name))
      
      contrast_vec <- makeContrasts(contrasts = contrast_formula, levels = design)
      
      fit2 <- contrasts.fit(fit, contrast_vec)
      fit2 <- eBayes(fit2)
      
      # Extract results
      res_path <- topTable(fit2, number = Inf, adjust.method = "BH")
      res_path$Pathway <- rownames(res_path)
      res_path$GeneSet_Collection <- gs_name
      
      # Store and save results
      result_key <- paste0(gs_name, "_", contrast_name)
      de_pathways_results[[result_key]] <- res_path
      write.csv(res_path,
                file.path(output_dir, "GSVA", contrast_name, paste0("GSVA_DE_", result_key, ".csv")),
                row.names = FALSE)
      
      # Export heatmap for differential pathways
      sig_pathways <- res_path[res_path$adj.P.Val < 0.05, ]
      
      if (nrow(sig_pathways) >= 5) {
        top_pathways <- head(sig_pathways$Pathway, min(20, nrow(sig_pathways)))
        plot_data <- gsva_mat[top_pathways, , drop = FALSE]
        
        # Create annotation for samples
        annotation_col <- data.frame(Group = metadata[colnames(plot_data), "group"])
        rownames(annotation_col) <- colnames(plot_data)
        
        gsva_diffpath_hm <- plot_heatmap(
          data_matrix = plot_data,
          title =  paste0("GSVA_Heatmap_", gs_name, "_", contrast_name),
          show_rownames = TRUE,
          cluster_rows = TRUE
        )
        
        pdf(file.path(output_dir,  "GSVA", contrast_name, 
                      paste0("GSVA_Heatmap_", gs_name, "_", contrast_name, ".pdf")), 
            width = 10, height = 8)
        print(gsva_diffpath_hm)
        dev.off()
        
        tiff(file.path(output_dir,  "GSVA", contrast_name, 
                       paste0("GSVA_Heatmap_", gs_name, "_", contrast_name, ".tiff")), 
             width = 10 * 300, height = 8 * 300, res = 300)
        print(gsva_diffpath_hm)
        dev.off()
        
        png(file.path(output_dir,  "GSVA", contrast_name, 
                      paste0("GSVA_Heatmap_", gs_name, "_", contrast_name, ".png")), 
            width = 10 * 300, height = 8 * 300, res = 300)
        print(gsva_diffpath_hm)
        dev.off()
        
        ggsave(file.path(output_dir,  "GSVA", contrast_name, 
                         paste0("GSVA_Heatmap_", gs_name, "_", contrast_name, ".png")), plot = gsva_diffpath_hm, 
               width = 10, height = 8, unit = "in")
      }
     
    # Using EnhancedVolcano 
    if (nrow(res_path) > 5) {
      title_suffix = paste(gs_name, "-", contrast_name)
      gsva_diffpath_volcano <- EnhancedVolcano(res_path,
                           lab = res_path$Pathway,
                           x = 'logFC',
                           y = 'adj.P.Val',
                           title = paste('GSVA Pathway Analysis:', title_suffix),
                           subtitle = paste('Significant: adj.P.Val < 0.05', 
                                            '& |logFC| > 0.25'),
                           xlab = bquote(~Log[2]~ 'fold change (pathway activity)'),
                           ylab = bquote(~-Log[10]~ 'adjusted p-value'),
                           
                           # Visual parameters optimized for pathway names
                           pointSize = 4.0,  # Larger points for pathways
                           labSize = 3.5,    # Slightly larger text for pathway names
                           boxedLabels = TRUE, # Box around labels for better readability
                           
                           # Cutoff parameters
                           pCutoff = 0.05,
                           FCcutoff = 0.25,
                           
                           # Coloring and theming
                           colAlpha = 0.8,
                           col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                           cutoffLineType = 'dashed',
                           cutoffLineWidth = 0.5,
                           
                           # Label selection
                           selectLab = c(),  # Can specify particular pathways to label
                           drawConnectors = TRUE,  # Connectors help with overlapping labels
                           widthConnectors = 0.5,
                           max.overlaps = 20,  # Increased for pathway names
                           
                           # Axis adjustments
                           xlim = c(min(res_path$logFC, na.rm = TRUE) - 0.5,
                                    max(res_path$logFC, na.rm = TRUE) + 0.5),
                           ylim = c(0, max(-log10(res_path$adj.P.Val), na.rm = TRUE) + 1),
                           
                           # Grid and theme
                           gridlines.major = TRUE,
                           gridlines.minor = FALSE,
                           caption = paste('Total pathways:', nrow(res_path),
                                           '| Significant:', 
                                           sum(res_path$adj.P.Val < 0.05 & 
                                                 abs(res_path$logFC) > 0.25)),
                           
                           legendPosition = 'right',
                           legendLabSize = 12,
                           legendIconSize = 4.0) +
        theme(plot.title = element_text(size = 16, face = 'bold'),
              plot.subtitle = element_text(size = 12),
              axis.title = element_text(size = 14),
              axis.text = element_text(size = 12))
      
      ggsave(file.path(output_dir,  "GSVA", contrast_name, 
                       paste0("GSVA_Volcano_", gs_name, "_", contrast_name, ".png")),
             gsva_diffpath_volcano, width = 10, height = 8)
      
      ggsave(file.path(output_dir,  "GSVA", contrast_name, 
                       paste0("GSVA_Volcano_", gs_name, "_", contrast_name, ".pdf")),
             gsva_diffpath_volcano, width = 10, height = 8)
      
      ggsave(file.path(output_dir,  "GSVA", contrast_name, 
                       paste0("GSVA_Volcano_", gs_name, "_", contrast_name, ".tiff")),
             gsva_diffpath_volcano, width = 10, height = 8)
      
      ggsave(file.path(output_dir,  "GSVA", contrast_name, 
                       paste0("GSVA_Volcano_", gs_name, "_", contrast_name, ".svg")),
             gsva_diffpath_volcano, width = 10, height = 8, unit = "in")
      
    } 
    
    # Ridge/Density Plot for a few selected pathways
    if (nrow(sig_pathways) > 0) {
      selected_paths <- head(sig_pathways$Pathway, min(6, nrow(sig_pathways)))
      plot_data <- as.data.frame(t(gsva_mat[selected_paths, ]))
      plot_data$Group <- metadata[rownames(plot_data), "group"]
        
      plot_data_long <- tidyr::pivot_longer(plot_data, -Group, names_to = "Pathway")
        
      gsva_diffpath_ridge <- ggplot(plot_data_long, aes(x = value, y = Group, fill = Group)) +
          ggridges::geom_density_ridges(alpha = 0.7) +
          facet_wrap(~Pathway, scales = "free_x", ncol = 2) +
          theme_minimal()
        
      ggsave(file.path(output_dir,  "GSVA", contrast_name, 
                       paste0("GSVA_Ridge_", gs_name, "_", contrast_name, ".pdf")),
             gsva_diffpath_ridge, width = 10, height = 8)
      
      ggsave(file.path(output_dir,  "GSVA", contrast_name, 
                       paste0("GSVA_Ridge_", gs_name, "_", contrast_name, ".png")),
             gsva_diffpath_ridge, width = 10, height = 8)
      
      ggsave(file.path(output_dir,  "GSVA", contrast_name, 
                       paste0("GSVA_Ridge_", gs_name, "_", contrast_name, ".tiff")),
             gsva_diffpath_ridge, width = 10, height = 8)
      
      ggsave(file.path(output_dir,  "GSVA", contrast_name, 
                       paste0("GSVA_Ridge_", gs_name, "_", contrast_name, ".svg")),
             gsva_diffpath_ridge, width = 10, height = 8, unit = "in")
    } 
    }
  }
}

# 5. Generate Publication-Quality Global Plots
message("  Generating publication-quality global plots...")

if (length(gsva_results) > 0) {
  GSVA_global_out <- file.path(output_dir, "GSVA", "GSVA_global")
  dir.create(GSVA_global_out)
  
  for (gs_name in names(gsva_results)) {
    message(paste("    Creating global plots for:", gs_name))
    
    gs_GSVA_global_out <- file.path(GSVA_global_out, gs_name)
    dir.create(gs_GSVA_global_out)
    
    gsva_mat <- gsva_results[[gs_name]]
    
    # A. HEATMAP: Global activity patterns across all samples
    try({
      # Select top variable pathways for clarity (adjust 'top_n' as needed)
      pathway_variance <- apply(gsva_mat, 1, var)
      top_n <- min(50, nrow(gsva_mat))
      top_var_pathways <- names(sort(pathway_variance, decreasing = TRUE))[1:top_n]
      plot_mat <- gsva_mat[top_var_pathways, ]
      
      # Prepare sample annotation
      sample_anno <- metadata[colnames(plot_mat), "group", drop = FALSE]
      if ("biorep" %in% colnames(metadata)) {
        sample_anno$Biorep <- metadata[colnames(plot_mat), "biorep"]
      }
      
      # Color palette (RdBu is standard for expression/activity divergence)
      heatmap_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100)
      
      # Draw and save heatmap
      dir.create(file.path(gs_GSVA_global_out, "GSVA_global_heatmaps"))
      heatmap_filename <- file.path(gs_GSVA_global_out, "GSVA_global_heatmaps",
                                    paste0("GSVA_Global_Heatmap_", gs_name, ".pdf"))
      
      
      gsva_global_hm <- pheatmap(plot_mat,
                                 scale = "row",
                                 color = heatmap_palette,
                                 border_color = NA,
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 annotation_col = sample_anno,
                                 show_rownames = TRUE,
                                 show_colnames = ncol(plot_mat) <= 30,
                                 fontsize_row = 7,
                                 main = paste("GSVA:", gs_name, "\nTop", top_n, "Variable Pathways"),
                                 width = 14, height = max(7, 0.25 * top_n))
      
      ggsave(file.path(gs_GSVA_global_out, "GSVA_global_heatmaps",
                       paste0("GSVA_Global_Heatmap_", gs_name, ".pdf")),
             gsva_global_hm, width = 10, height = 8)
      
      ggsave(file.path(gs_GSVA_global_out, "GSVA_global_heatmaps",
                       paste0("GSVA_Global_Heatmap_", gs_name, ".png")),
             gsva_global_hm, width = 10, height = 8, dpi = 600)
      
      ggsave(file.path(gs_GSVA_global_out, "GSVA_global_heatmaps",
                       paste0("GSVA_Global_Heatmap_", gs_name, ".tiff")),
             gsva_global_hm, width = 10, height = 8, dpi = 600)
      
      ggsave(file.path(gs_GSVA_global_out, "GSVA_global_heatmaps",
                       paste0("GSVA_Global_Heatmap_", gs_name, ".svg")),
             gsva_global_hm, width = 10, height = 8, unit = "in")
      
      
      message(paste("      Saved heatmap:", heatmap_filename))
      
    }, silent = TRUE)
  }
  
  message("  Global plotting complete.")
  
} else {
  message("  No GSVA results available for plotting.")
}

# Generate report -------------------------------------------------------------
message("\n[7/8] Generating final report...")
#Not complete

message("\nAnalysis complete! Results saved to: ", normalizePath(output_dir))

# Save session -------------------------------------------------------------
message("\n[8/8] Saving session info...")
save.image(file.path(output_dir, "rdata/analysis_session.RData"))
message("\nAnalysis complete! Results saved to: ", normalizePath(output_dir))


