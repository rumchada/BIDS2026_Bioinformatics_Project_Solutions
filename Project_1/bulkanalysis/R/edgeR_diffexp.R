#' edgeR Differential Expression Analysis
#'
#' Performs differential expression analysis on a formatted DESeq2 object.
#'
#' @name edgeR_diffexp
#' @aliases edgeR_diffexp
#'
#' @param dds_object A formatted DESeq2 object where each sample and its
#'   counts are aligned with the respective metadata in `@colData`.
#'
#' @param condition The name of the condition column in the metadata
#'   table in `colData`.
#'
#' @param reference_group_name The name of the control group to compare
#'   against for each comparison. For example, `"healthy"` or `"control"`.
#'
#' @param initial_pvalue_cutoff The initial p-value cutoff used for
#'   filtering.
#'
#'  @param adjust.method choose a method of pval correct Family Wise Error Rate reccomended
#'
#' @return A list of data frames containing every possible combination
#'   of binary comparisons (control vs. experimental group) within your
#'   metadata.
#'
#' @description
#' Performs differential expression analysis on a formatted DESeq2 object.
#'
#' @details
#' INPUTS:
#'
#' A DESeq2 object with metadata columns corresponding to the samples
#' supplied.
#'
#' `condition`: A string giving the name of the condition column within
#' the DESeq2 object's `@colData` metadata.
#'
#' `reference_group_name`: A string specifying the control group for
#' the experiment.
#'
#' OUTPUTS:
#'
#' A list of data frames containing every possible combination of binary
#' comparisons (control vs. experimental group) within the supplied
#' metadata.
#'
#' The resulting data frames contain the following columns:
#'
#' * `log2FC`: The log2-scaled fold change of a given gene across the
#'   specified samples.
#'
#' * `P-value`: The probability of obtaining data as extreme as the
#'   observed values assuming that the null hypothesis is true.
#'
#' * `Q-value (FDR)`: The proportion of false positives after accounting
#'   for multiple testing.


edgeR_diffexp <- function(dds_object,
                          condition_col,
                          ref_group_name,
                          init_pval_cutof = 1,
                          adjust.method = "bonfferonni") {
  require(edgeR)
  #For each condition combination, runs an indepedent differential expression analysis
  #0-Intercept

  #Converting to a DGEList
  raw_counts <- dds_object@assays@data$counts
  metadata <- dds_object@colData
  y <- DGEList(counts = raw_counts , group = metadata$condition)

  #Filter by low counts
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  #Performs TMM Normalization, Corrects for systemic biases between samples.
  y <- calcNormFactors(y)

  # reaffirm healthy as level reference
  y$samples$condition <- relevel(
    factor(y$samples$group),
    ref = ref_group_name
  )

  # ~0+ not to include an intercept column and instead to include a column for each group
  # Always use a model intecept term because then the data will no fit to the glm line as closely
  # A model without an intercept term would only be recommended in cases where there is a strong biological reason why a zero covariate should be associated with a zero expression value (Law, et.al., F10000Research, 2020, PMID: 33604029)

  design <- model.matrix(~0+condition, data = y$samples)
  conds <- gsub("^condition", "", colnames(design))

  contrast_strings <- combn(conds, 2, FUN = function(x) paste0(x[1], "-", x[2]))


  contrast_strings
  colnames(design) <- gsub(glue("^{condition_col}"), "", colnames(design))

  my.contrasts <- makeContrasts(
    contrasts = contrast_strings,
    levels = colnames(design))

  #Test Using a GLM
  results <- list()

  # Run Dispersion on the Data
  y <- estimateDisp(y, design)

  #fit model once
  fit <- glmQLFit(y, design)

  # run the test for each contrast
  for(i in seq_len(ncol(my.contrasts))) {

    qlf <- glmQLFTest( fit, contrast = my.contrasts[, i])

    results[[colnames(my.contrasts)[i]]] <-
      topTags(qlf, n = nrow(y), p.value = init_pval_cutoff, adjust.method = adjust.method)$table %>%
      as.data.frame() %>%
      tibble::rownames_to_column("geneid") %>%
      dplyr::rename_with( ~ "log2foldchange", .cols = "logFC") %>%
      dplyr::rename_with( ~ "p.val_adj", .cols = dplyr::any_of(c("FWER", "FDR")) )
  }

  return(results)
}
