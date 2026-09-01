#' @title Baseline Dimension Reduction of RNA-seq Data
#'
#' @name baseline_dimreduction
#'
#' @description
#' Takes a formatted \code{DESeqDataSet} object and performs several
#' exploratory analyses to evaluate sample similarity and major sources
#' of variation within the dataset.
#' Performs baseline exploratory dimension reduction and sample-level
#' visualization of a \code{DESeqDataSet}. Variance-stabilized expression
#' values are used to calculate sample distances and principal components.
#' The most variable genes are used for PCA and downstream visualization.
#'
#' The analysis:
#' \enumerate{
#' \item applies a variance-stabilizing transformation to the expression
#' matrix,
#' \item identifies the top variable genes,
#' \item calculates Euclidean distances between samples,
#' \item performs principal component analysis (PCA) using the top variable
#' genes,
#' \item generates sample embedding plots,
#' \item generates an elbow plot showing the variance explained by the
#' leading principal components, and
#' \item performs UMAP for nonlinear dimensionality reduction.
#' }
#'
#' @param bulk_dataset A formatted \code{DESeqDataSet} (\code{DDS}) object
#' containing the RNA-seq expression data and associated sample metadata.
#'
#' @param top_var Integer specifying the number of genes with the highest
#' variance to retain for PCA and downstream dimension reduction.
#' For example, \code{top_var = 2000} uses the 2,000 most variable genes.
#'
#' @param grouping Character string specifying the sample metadata variable
#' used to annotate or group samples in the distance and dimension
#' reduction visualizations. The variable should correspond to a column
#' in the \code{colData} of \code{bulk_dataset}.
#'
#'
#'@returns
#'\itemize{
#'  \item{ \code{distance_plot}: euclidean distance plot,}
#'  \item{\code{elbow_plot}: shows which PC has the majority of the variance within the data }
#'  \item{\code{PCA plot}: scatterplot of the left unitary vectors}
#'  \item{\code{UMAP_plot}: A UMAP visualization of the sample embeddings.}
#'
#'  results$distance_plot
#' results$PCA_plot
#' results$elbow_plot
#' results$umap_plot
#'}
#'
#'@examples
#'\dontrun{ baseline_dimreduction(dds_object, top_var = 2000, grouping = NA)}
#'@export
baseline_dimreduction <- function(bulk_dataset,
                                  top_var = 2000,
                                  grouping = NA){

  library(ggplot2)
  require(umap)
  require(prcomp)

  if (!is.na(grouping) && is.character(grouping)) {
    if (!grouping %in% colnames(colData(bulk_dataset))) {
      stop(paste("The grouping variable", grouping, "was not found in your colData."))
    }
    anno_variable <- c("condition", grouping)
    color_variable <- grouping
  } else {
    anno_variable <- c("condition")
    color_variable <- "condition"
  }




  #prioritize variance stable assay within the bulk_dataset
  assays(bulk_dataset) <- assays(bulk_dataset)[c("var_stable", "counts", "log_counts")]
  # takes the euclidan distance of a variance stabilized dataset
  # row by row euclidean distance against each other
  # plots the visual distance between each sample onto heatmap
  distance_plot <- plot_sample_clustering(bulk_dataset, anno_vars = anno_variable, distance = "euclidean")


  #____Extracting Top Variable Gene per sample -----#
  # Setting Variance Stable Counts as Priority Assay
  vst_mat <- assay(bulk_dataset, "var_stable")
  #Select top variable gene
  top_var_genes <- head(order(rowVars(vst_mat), decreasing = TRUE), top_var)
  vst_top <- vst_mat[top_var_genes, ]
  # Run PCA via R
  # Remember PCA is taking the right unitary martrix (VT) (feature space) of SVD
  # VT tells us the principal directions of the data set
  # multiply the the original matrix by the right unitary which is equivalent to the
  # U time signma (left unitary) * Singular values (Size of the PCs)
  pca_res <- prcomp(t(vst_top))


  # Plotting first two PCAs
  pca_plot <- plot_pca(bulk_dataset, PC_x = 1, PC_y = 2, color_by = color_variable)



  # Calculating Percentage of Variance Per PC for the Scree Plot

  # Calculate percentage variance explained
  var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100

  # Create a data frame for plotting (first 10 PCs)
  scree_data <- data.frame(
    PC = paste0("PC", 1:10),
    Variance = var_explained[1:10]
  )

  # Ensure PC is treated as a factor in order
  scree_data$PC <- factor(scree_data$PC, levels = scree_data$PC)

  # Scree Plot ensures that we are plotting PCs that capture the most amount of variability within our data
  scree_plot <- ggplot(scree_data, aes(x = PC, y = Variance, group = 1)) +
    geom_line(color = "red", size = 1) +
    geom_point(color = "red", size = 2) +
    labs(title = "Scree Plot: Variance Explained per PC",
         x = "Principal Component",
         y = "Percentage of Variance (%)") +
    theme_minimal()


  #Input the Left Unitary * Singular Values
  # 1. Run PCA
  pca_res <- prcomp(t(assay(bulk_dataset, "var_stable")), rank. = 50)

  # 2. Use the 'x' slot (the PC scores) as input for UMAP
  # We only take the number of PCs identified by your scree plot (e.g., 1:15)
  umap_out <- umap(pca_res$x[, 1:15])


  umap_df <- data.frame(
    Sample_IDS = rownames(pca_res$x),
    UMAP_1 = umap_out$layout[, 1],
    UMAP_2 = umap_out$layout[, 2],
    ColorGroup = as.factor(colData(bulk_dataset)[[color_variable]]))

  umap_plot <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = ColorGroup)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "UMAP of COVID-19 Samples", subtitle = "Input: Top 15 PCs from VST Data")



  dimreduction_visuals <- c(distance_plot, scree_plot, pca_plot$plot, umap_plot)
  return(dimreduction_visuals)
}
