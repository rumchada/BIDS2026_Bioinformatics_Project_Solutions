# Project 1 Part 1 and Part 3 Functions
#QC functions and dim_reduction functions.

pc_var_association <- function(bulk_ds){
  library(pheatmap)
  vst_mat <- assay(bulk_ds, "var_stable")
  #Select top variable gene
  top_var_genes <- head(order(rowVars(vst_mat), decreasing = TRUE), 2000)
  vst_top <- vst_mat[top_var_genes, ]
  # Run PCA via R
  # Remember PCA is taking the right unitary martrix (VT) (feature space) of SVD
  # VT tells us the principal directions of the data set
  # multiply the the original matrix by the right unitary which is equivalent to the
  # U time signma (left unitary) * Singular values (Size of the PCs)
  pca_res <- prcomp(t(vst_top))

  top_pcs <- pca_res$x[, 1:10]

  #creating a list of possible covariate to compare with PCs
  meta_data <- as.data.frame(bulk_ds@colData)

  meta_data[] <- lapply(meta_data, factor)

  covars <- split(meta_data,
                  rep(1:ncol(meta_data),
                      each = nrow(meta_data)))

  covars <- lapply(seq_len(ncol(meta_data)),
                   function(x) meta_data[ , x])

  names(covars) <- colnames(meta_data)


  get_p_value <- function(pc_vector, var) {


    valid <- !is.na(pc_vector) & !is.na(var)

    pc_vector <- pc_vector[valid]
    var_vector <- var[valid]


    if (is.numeric(var_vector)) {
      # Continuous variable: use linear regression/correlation test
      return(summary(lm(pc_vector ~ var_vector))$coefficients[2, 4])
    } else {
      # Categorical variable: use ANOVA
      fit <- aov(pc_vector ~ as.factor(var_vector))

      anova_sum <- summary(fit)[[1]]
      p_val <- anova_sum[1, "Pr(>F)"]

      if (length(p_val) != 1 || is.na(p_val) || is.nan(p_val) || is.infinite(p_val)) {
        return(NA)
      }
      return(p_val)

    }
  }

  p_matrix <- matrix(NA,
                     nrow = length(covars),
                     ncol = ncol(top_pcs),
                     dimnames = list(names(covars), colnames(top_pcs)))


  for(vars in names(covars)){

    current_var <- covars[[vars]]

    if (is.factor(current_var)) {
      n_levels <- nlevels(current_var)
    } else {
      n_levels <- length(unique(current_var))
    }
    if (n_levels < 2) next

    for (pcs in colnames(top_pcs)) {

      current_pc <- top_pcs[, pcs]

      p_matrix[vars, pcs] <- get_p_value(current_pc, current_var)
    }
  }

  log_p_matrix <- -log10(p_matrix)

  max_finite <- max(log_p_matrix[is.finite(log_p_matrix)])

  log_p_matrix[is.infinite(log_p_matrix)] <- max_finite + 1

  log_p_matrix <- log_p_matrix[!apply(is.na(p_matrix), 1 , all), ]

  # 5. Plot the heatmap
  p_val_map <- pheatmap(log_p_matrix,
                        legend_labels = "log10(pvalue)",
                        cluster_rows=F,
                        cluster_cols=F,
                        main = "PC-Metadata Associations (-log10(p-value))")
  return(p_val_map)
}

baseline_dimreduction <- function(bulk_dataset,
                                  grouping = NA,
                                  object_type = "dds"){

  library(umap)
  library(ggplot2)

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
  top_var_genes <- head(order(rowVars(vst_mat), decreasing = TRUE), 2000)
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



# Part 2 Intra-Dataset Differntial Expression Functions
edgeR_diffexp <- function(dds_object,
                          condition_col,
                          ref_group_name,
                          init_pval_cutoff = 1,
                          adjust.method = "bonferroni") {
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
      topTags(qlf, n = nrow(y), p.value = init_pval_cutoff, adjust.method = 'bonferroni')$table %>%
      as.data.frame() %>%
      tibble::rownames_to_column("geneid") %>%
      dplyr::rename_with( ~ "log2foldchange", .cols = "logFC") %>%
      dplyr::rename_with( ~ "p.val_adj", .cols = dplyr::any_of(c("FDR", "FWER")))
  }

  return(results)
}


# Volcano Plot Visualizations
volcano_plot <- function(diffexp_df,
                         #data-alterations
                         log2fc_thresh = 0,
                         p.val_adj_thresh = 0.05,
                         #visual parts
                         lower_xlim = -5,
                         upper_xlim = 5,
                         step = 1) {

  library(ggrepel)
  library(ggplot2)
  library(dplyr)

volcano_df <- diffexp_df %>%
    mutate(
      color = case_when(
        p.val_adj < p.val_adj_thresh & log2foldchange >  log2fc_thresh  ~ "Upregulated",
        p.val_adj < p.val_adj_thresh & log2foldchange < -log2fc_thresh  ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    )

    deg_df <- volcano_df %>%
      dplyr::filter(color != "Not Significant")


  volplot <- ggplot(volcano_df) +

    aes(x = log2foldchange,y = -log10(p.val_adj), color = color, label = geneid) +

    geom_point(size = 2) +

    geom_vline(xintercept = c(-0.6, 0.6), col = "gray",linetype = "dashed"
    ) +

    geom_hline(yintercept = -log10(0.05), col = "gray",linetype = "dashed") +

    theme_classic(base_size = 12) +

    theme(axis.title.y = element_text(margin = margin(0, 20, 0, 0),size = rel(1.1),color = "black"
    ),

    axis.title.x = element_text(hjust = 0.5, margin = margin(20, 0, 0, 0),size = rel(1.1),
                                color = "black"
    ),

    plot.title = element_text(hjust = 0.5)
    ) +

    scale_color_manual(values = c("Upregulated" = "red","Downregulated" = "blue","Not Significant" = "gray"
    )
    ) +

    geom_text(nudge_y = 0.5, check_overlap = TRUE) +

    coord_cartesian(xlim = c(lower_xlim, upper_xlim)) +

    scale_x_continuous(
      breaks = seq(lower_xlim, upper_xlim, step)
    )

  return(list(deg_df, volplot))
}


#unpacks enrichGO results, converts them to the user's desired gene IDs, then aligns them to their log2fc results

enrichgo_unpack_ver2 <- function(intial_table = filtered_results$`healthy-covid19`,
                                 key_table = enrichgo_results$`healthy-covid19`,
                                 chosen_desc,
                                 #assume user does not want to convert geneids
                                 #arguments below are for the conversion helperfunction
                                 convert = TRUE,
                                 from_type = "ensembl",
                                 to_type = "symbol",
                                 ensembl_dataset = "hsapiens_gene_ensembl"){



  # Loop will be placed here.
  cluster <- tibble(key@result) %>%
    dplyr::filter(str_detect(Description, chosen_desc)) %>%
    dplyr::mutate(Description = as.character(Description))
  cluster


  ### helper functionf or converting geneidA to geneidB
  ### user
  gene_id_converter <- function(vector, from_type, to_type, ensembl_dataset){
    library(biomaRt)
    # geneid pulling for correct attribute
    id_map <- c(
      ensembl = "ensembl_gene_id",
      entrez  = "entrezgene_id",
      symbol  = "external_gene_name"
    )
    # p
    from_attr <- id_map[[from_type]]
    to_attr <-id_map[[to_type]]

    ensembl <- useMart(biomart = "ensembl", dataset = ensembl_dataset)

    return_df <- getBM( attributes = c(from_attr, to_attr), filters = from_attr, values = vector, mart = ensembl
    )
    return(return_df)

  }


  ###creates a helper function to unlist and strsplits the enriched Gene IDs per term
  unlist_converter <- function(str) {
    if (grepl("^ENSG|\\d", str)){
      bruh <- as.character(unlist(strsplit(str, "/")))
    }else{
      bruh <- as.integer(unlist(strsplit(str, "/")))
    }
    return(bruh)
  }

  #apply this function to the geneID column
  extract_result <- lapply(cluster$geneID, unlist_converter)
  #Next check the length of the new vector to see if it matches the Count column in the GSEA DF.
  length(extract_result[[1]])

  #create an empty vector
  path_vector <- c()
  #for each index in the first column of numbers
  for (i in seq_along(cluster$Description)) {
    #replicate the name of the Description for the length of the GO term in column 1
    path_vector <- c(path_vector, rep(cluster$Description[i], length(extract_result[[i]])))
  }

  # make new tibble that pairs the name of the Pathway from the GOTerm and enriched Entrez ID from the "post-processed" GeneID column
  pathway_tibble <- tibble(gene_desc = path_vector, geneid = unlist(extract_result))
  pathway_tibble

  if(convert == FALSE){
    pathway_genenames <- pathway_tibble
  }else{
    pathway_genenames <- gene_id_converter(pathway_tibble$geneid,
                                           from_type = from_type,
                                           to_type = to_type,
                                           ensembl_dataset = ensembl_dataset)

    new_colname <- paste0(to_type,"_", "geneid")
    pathway_genenames <- pathway_genenames %>%
      rename(geneid = colnames(pathway_genenames)[1]) %>%
      rename(new_colname = colnames(pathway_genenames)[2])
  }


  final_path_tibble <- inner_join(pathway_tibble, pathway_genenames, by = "geneid")
  final_path_tibble <- final_path_tibble %>% mutate(geneid = as.character(geneid))
  final_path_tibble


  final_path_fc_tibble <- inner_join(final_path_tibble, initial_table, by = "geneid")

  fpfc_tibble <- final_path_fc_tibble %>%
    relocate(new_colname, .after = geneid)

  return(fpfc_tibble)
}

### helper functionf or converting geneidA to geneidB
### user
gene_id_converter <- function(vector, from_type, to_type, ensembl_dataset){
  library(biomaRt)
  # geneid pulling for correct attribute
  id_map <- c(
    ensembl = "ensembl_gene_id",
    entrez  = "entrezgene_id",
    symbol  = "external_gene_name"
  )
  # p
  from_attr <- id_map[[from_type]]
  to_attr <-id_map[[to_type]]

  ensembl <- useMart(biomart = "ensembl", dataset = ensembl_dataset)

  return_df <- getBM(
    attributes = c(from_attr, to_attr),
    filters = from_attr,
    values = vector,
    mart = ensembl
  )
  return(return_df)

}

gene_id_converter_ver2 <- function(vector, from_type, to_type, ensembl_datset){
  require(biomaRt)
  # geneid pulling for correct attribute
  id_map <- c(
    ensembl = "ensembl_gene_id",
    entrez  = "entrezgene_id",
    symbol  = "external_gene_name"
  )
  # Pass the the ID mapping index over to the new attribute string
  from_attr <- id_map[[from_type]]
  to_attr <-id_map[[to_type]]

  #initialize the return and parameters
  #baseline NULL intialization
  return_df <- NULL
  max_attempts <- 2
  attempt <- 1
  org_db <- org.Hs.eg.db


  #Fix 1: There is a huge API crash almost 50% of the time that I use getBM() function to connect to the API
  #Implementing a while loop/tryCatch function to loop API calls just in case that the getBM work
  #Solution:After 10 attempts,

  #while the current attemp is less than or equal to the max attempts and while the return DF is null
  while(attempt <= max_attempts && is.null(return_df)){
    return_df <- tryCatch({

      message(sprintf("biomaRt attempt %d of %d...", attempt, max_attempts))

      #call the getBM ensemble API
      ensembl <- useMart(biomart = "ensembl", dataset = ensembl_dataset)
      getBM( attributes = c(from_attr, to_attr), filters = from_attr, values = vector, mart = ensembl)

      #If there is any kind of error
    }, error=function(e){
      #message for that error
      message(sprintf("  -> biomaRt attempt %d failed: %s", attempt, trimws(e$message)))
      #Return a NULL DF
      return(NULL)
    })

    # During that current attempt, If the return_df is still null
    if (is.null(return_df)) {
      # add on the attempt counter
      attempt <- attempt + 1
      # If the attempt counter is still less than the max attempts
      if (attempt <= max_attempts) {
        Sys.sleep(2) # Brief pause to prevent hammering the server
      }
    }# end of second if statement
  }#end of while loop

  if(is.null(return_df)){
    dbi_map <- c(
      ensembl = "ENSEMBL",
      entrez  = "ENTREZID",
      symbol  = "SYMBOL"
    )

    dbi_from <- dbi_map[[from_type]]
    dbi_to   <- dbi_map[[to_type]]

    return_df <- AnnotationDbi::mapIds(org_db,
                                       keys = vector,
                                       keytype = dbi_from,
                                       column =   dbi_to,
                                       multiVals = "first")

    return_df <- data.frame(
      from_attr = names(return_df),
      to_attr = unname(return_df))

    colnames(return_df) <- c(from_attr, to_attr)
  }
  return(return_df)
}#end of whole function



heatmap_function <- function(initial_table, table_list) {

  require(colorRamp2)
  require(ComplexHeatmap)
  require(circlize)

  # Initialize the list to store heatmap objects
  heatmap_list <- list()

  # Extract the raw count matrix once outside the loop for efficiency
  # Assumes a SingleCellExperiment or SummarizedExperiment-like structure
  expr_raw <- as.matrix(initial_table@assays@data$log_counts)

  #Early stopping points just in case
  if (any(is.na(colnames(expr_raw)))) stop("NA column names in expr_raw")
  if (any(duplicated(colnames(expr_raw)))) stop("Duplicate column names in expr_raw")

  # Loop through the list of filtered differential expression results
  for (i in seq_along(table_list)) {

    # Grab the ith table
    filtered_diff_exp_results <- table_list[[i]]

    if (!"geneid" %in% colnames(filtered_diff_exp_results)) {
      stop(glue::glue("Missing geneid column in {names(table_list)[i]}"))
      print(colnames(filtered_diff_exp_results))
    }

    # Split the names of the ctrl-treatment comparison
    comparison_name <- names(table_list)[i]
    ctrl_treatment <- strsplit(comparison_name, "-")[[1]]
    #itemizing the names of the
    ctrl_name      <- ctrl_treatment[1]
    treatment_name <- ctrl_treatment[2]



    # Message tracking
    message(glue::glue("Processing: {ctrl_name} vs {treatment_name}"))



    # Pull the differentially expressed gene IDs
    diffexp_geneids <- unique(as.character(filtered_diff_exp_results$geneid))

    print(length(diffexp_geneids))

    genes_found <- diffexp_geneids %in% rownames(expr_raw)

    validated_geneids <- diffexp_geneids[genes_found]


    if (sum(genes_found) == 0) {
      stop(glue::glue(
        "No DE genes found in expression matrix for {comparison_name}"
      ))
    }

    if (sum(genes_found) < length(diffexp_geneids) * 0.5) {
      warning(glue::glue(
        "Less than 50% of DE genes found in matrix for {comparison_name}"
      ))
    }

    # Filter the matrix for the target genes and matching sample columns
    # Note: Added anchors to ensure exact suffix matching
    col_pattern <- glue::glue("_{ctrl_name}$|_{treatment_name}$")
    matching_cols <- grepl(col_pattern, colnames(expr_raw))

    expr <- expr_raw[validated_geneids, matching_cols, drop = FALSE]

    message(glue("Comparison:{ctrl_name} vs {treatment_name} has {nrow(expr)} DEGs"))

    # Scaling and z-score normalizing the counts matrix
    expr_z <- t(scale(t(expr)))
    expr_z[is.na(expr_z)] <- 0
    expr_z[!is.finite(expr_z)] <- 0
    rownames(expr_z) <- rownames(expr)
    colnames(expr_z) <- colnames(expr)

    # Gene ID conversion
    #Quality check of the number of DEGs matched from raw_expr table that were z-scaled
    message(glue("Comparison:{ctrl_name} vs {treatment_name} has {nrow(expr_z)} DEGs"))


    conversion_table <- gene_id_converter_ver2(
      rownames(expr_z),
      "ensembl",
      "symbol",
      "hsapiens_gene_ensembl"
    )

    gene_map <- setNames(
      conversion_table$external_gene_name,
      conversion_table$ensembl_gene_id
    )

    new_names <- gene_map[rownames(expr_z)]

    # Fallback to Ensembl ID if no symbol is found
    rownames(expr_z) <- ifelse(is.na(new_names), rownames(expr_z), new_names)


    # Clustering data (Optional: ComplexHeatmap does this natively if cluster_rows = TRUE,
    # but kept if you intend to pass hclust objects explicitly)


    # THIS IS A RISKY MOVE BUT IT WILL MAKE THE HEATMAP A LOT CLEANER
    # PROS: IT WILL MAKE THE HEATMAP CLEANER
    # CONS: IT MAY REMOVE BIOLOGICAL SIGNAL FOR A GENE: COULD BE EXPLAINED BY VARIANTS OF THE SAME GENE
    expr_z <- expr_z[!duplicated(rownames(expr_z)), ]


    gene_dist <- dist(expr_z)
    gene_hclust <- hclust(gene_dist)

    #____________________________________________________#
    #Sample Color Annotation Block
    sample_condition <- c(sub(".*_", "", colnames(expr_z)))


    condition_levels <- unique(sample_condition)


    #Setting colors for the chosen disease states
    condition_cols <- setNames(
      scales::hue_pal()(length(condition_levels)),
      condition_levels
    )

    #Using HeatmapAnnotation
    HA <- HeatmapAnnotation(
      Condition = sample_condition,
      col = list(Condition = condition_cols)
    )

    #____________________________________________________#

    # Define color palette
    col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

    # THIS IS A RISKY MOVE BUT IT WILL MAKE THE HEATMAP A LOT CLEANER

    # Gene ID conversion
    #Quality check of the number of DEGs matched from raw_expr table that were z-scaled
    message(glue("Comparison:{ctrl_name} vs {treatment_name} has {nrow(expr_z)} UNIQUE DEGs"))

    # Generate Heatmap object
    p <- Heatmap(
      expr_z,
      name            = "Z-score",
      col             = col_fun,
      cluster_columns = TRUE,
      cluster_rows    = TRUE, # To use your manual cluster, change to: cluster_rows = gene_hclust
      row_names_gp    = gpar(fontsize = 12),
      top_annotation  = HA,
      show_row_dend   =  FALSE,
      show_column_names = FALSE,
      row_gap         = unit(5, "mm"),
      column_title    = glue::glue("{ctrl_name} - {treatment_name} Diff Exp Genes")
    )

    # Store the plot object into the initialized list
    heatmap_list[[comparison_name]] <- p
  }

  return(heatmap_list)
}

ora_dotplot <- function(ora_results, gspval_cutoff, org_by, num_disp, graph_title, file_name){
  if(org_by %in% 'Count'){
    #inside if
    dotplot <- ora_results%>%
      dplyr::arrange(desc(Count))%>%
      utils::head(num_disp)%>%
      dplyr::mutate(as.factor(Description))%>%
      ggplot() +
      aes(y = Count, x = fct_reorder(Description, Count)) +
      scale_colour_gradient(low = "red", high = "blue") +
      geom_point(aes(size = GeneRatio, color = pvalue)) +
      labs(title = graph_title, x = "Enriched Term", y = "Gene Count") +
      coord_flip() +
      theme(text = element_text(face = "bold") , plot.title = element_text(hjust = 1))

  }
  if(org_by %in% 'pvalue'){
    #inside if
    dotplot <- ora_results %>%
      dplyr::filter(pvalue < gspval_cutoff) %>%
      utils::head(num_disp) %>%
      dplyr::mutate(as.factor(Description)) %>%
      dplyr::mutate(pvalue = -1 * log(pvalue)) %>%
      ggplot() +
      aes(x = pvalue, y = fct_reorder(Description, pvalue))+
      scale_colour_gradient(low = "red", high = "blue") +
      geom_point(aes(size = GeneRatio, color = Count)) +
      labs(title = graph_title, x = "-log10(P)", y = "Enriched Term") +
      theme(text = element_text(face = "bold"), plot.title = element_text(hjust = 1))


  }
  return(list(dotplot, ggsave(file_name, device = "png", width = 8, height = 6, units = "in")))
}


ora_enrichgo <- function(filtered_results,
                         direction = c("up", "down"),
                         log2fc_thresh = 1,
                         pval_thresh = 0.05,
                         top_terms = 15,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENSEMBL",
                         ont = "BP"){

  direction <- match.arg(direction)

  # store enrichment results
  enrichgo_results <- vector("list", length(filtered_results))
  names(enrichgo_results) <- names(filtered_results)

  # run enrichGO
  for(i in seq_along(filtered_results)){

    if(direction == "up"){

      genes <- filtered_results[[i]]$geneid[
        filtered_results[[i]]$log2foldchange > log2fc_thresh
      ]

    } else if(direction == "down"){

      genes <- filtered_results[[i]]$geneid[
        filtered_results[[i]]$log2foldchange < -log2fc_thresh
      ]

    }

    enrichgo_results[[i]] <- clusterProfiler::enrichGO(
      gene = genes,
      OrgDb = OrgDb,
      keyType = keyType,
      ont = ont,
      pAdjustMethod = "BH"
    )
  }

  # store visualization outputs
  enrichgo_visuals <- vector("list", length(enrichgo_results))
  names(enrichgo_visuals) <- names(enrichgo_results)

  # generate plots
  for(i in seq_along(enrichgo_results)){

    ora_dotplot_count <- ora_dotplot(
      enrichgo_results[[i]],
      pval_thresh,
      "Count",
      top_terms,
      glue::glue(
        "ORA Results Ranked by Count: {names(enrichgo_results)[i]}"
      ),
      glue::glue(
        "ORA Results Ranked by Count: {names(enrichgo_results)[i]}.png"
      )
    )

    ora_dotplot_pvalue <- ora_dotplot(
      enrichgo_results[[i]],
      pval_thresh,
      "pvalue",
      top_terms,
      glue::glue(
        "ORA Results Ranked by PValue: {names(enrichgo_results)[i]}"
      ),
      glue::glue(
        "ORA Results Ranked by PValue: {names(enrichgo_results)[i]}.png"
      )
    )

    enrichgo_visuals[[i]] <- list(
      ora_count_ranked = ora_dotplot_count,
      ora_pvalue_ranked = ora_dotplot_pvalue
    )
  }

  return(list(
    enrichgo_results = enrichgo_results,
    enrichgo_visuals = enrichgo_visuals
  ))
}



enrichgo_unpack_ver3 <- function(
    # Pass the whole list so we can grab comparisons
  initial_table_list = filtered_results,
  key = ora_output_down$enrichgo_results,
  pval_adj_threshold = 0.001,

  # Assume user does not want to convert geneids
  convert = FALSE,
  from_type = "ensembl",
  to_type = "symbol",
  ensembl_dataset = "hsapiens_gene_ensembl"){
  require(tibble)
  require(dplyr)

  ### helper functionf or converting geneidA to geneidB
  ### user
  gene_id_converter <- function(vector, from_type, to_type, ensembl_dataset){
    require(biomaRt)
    # geneid pulling for correct attribute
    id_map <- c(
      ensembl = "ensembl_gene_id",
      entrez  = "entrezgene_id",
      symbol  = "external_gene_name"
    )
    # p
    from_attr <- id_map[[from_type]]
    to_attr <-id_map[[to_type]]

    ensembl <- useMart(biomart = "ensembl", dataset = ensembl_dataset)

    return_df <- getBM( attributes = c(from_attr, to_attr), filters = from_attr, values = vector, mart = ensembl
    )
    return(return_df)

  }



  ###creates a helper function to unlist and strsplits the enriched Gene IDs per term
  unlist_converter <- function(str) {
    if (grepl("^ENSG|\\d", str)){
      bruh <- as.character(unlist(strsplit(str, "/")))
    }else{
      bruh <- as.integer(unlist(strsplit(str, "/")))
    }
    return(bruh)
  }


  # FIX 2: Corrected the spelling typo here
  unpacking_results <- list()

  # FIX 1: Loop over names(key) so 'comparison' is a valid character string identifier
  for (comparison in names(key)) {

    # Safeguard check: If the comparison doesn't exist in the cluster results, skip
    if (is.null(key[[comparison]])) next

    cluster <- tibble::as_tibble(key[[comparison]]@result) %>%
      dplyr::mutate(Description = as.character(Description))

    if (!is.null(pval_adj_threshold)) {
      significant_clusters <- cluster %>%
        dplyr::filter(p.adjust < pval_adj_threshold)
    } else {
      significant_clusters <- cluster
    }

    labeled_gene_lists <- list()

    # FIX 3: Pull the specific comparison dataframe from your filtered results list
    current_initial_table <- initial_table_list[[comparison]] %>%
      dplyr::mutate(geneid = as.character(geneid))

    # If there are no significant clusters, jump to saving the empty list or next comparison
    if (nrow(significant_clusters) == 0) {
      unpacking_results[[comparison]] <- labeled_gene_lists
      next
    }

    for (i in seq_len(nrow(significant_clusters))) {

      go_term <- significant_clusters$Description[i]

      # Extract genes for current GO term
      genes <- unlist_converter(significant_clusters$geneID[i])

      pathway_tibble <- tibble::tibble(
        gene_desc = rep(go_term, length(genes)),
        geneid = as.character(genes)
      )

      if (convert == FALSE) {
        pathway_genenames <- pathway_tibble
      } else {
        pathway_genenames <- gene_id_converter(
          pathway_tibble$geneid,
          from_type = from_type,
          to_type = to_type,
          ensembl_dataset = ensembl_dataset
        )

        new_colname <- paste0(to_type, "_geneid")

        # FIX 4: Cleaner, bulletproof column renaming by position
        colnames(pathway_genenames)[1] <- "geneid"
        colnames(pathway_genenames)[2] <- new_colname
      }

      # Join converted names
      final_path_tibble <- dplyr::inner_join(pathway_tibble, pathway_genenames, by = "geneid") %>%
        dplyr::mutate(geneid = as.character(geneid))

      # Join fold changes / statistics using our dynamically updated initial table
      final_path_fc_tibble <- dplyr::inner_join(final_path_tibble, current_initial_table, by = "geneid")

      # Relocate converted gene column
      if (convert == TRUE) {
        fpfc_tibble <- final_path_fc_tibble %>%
          dplyr::relocate(dplyr::all_of(new_colname), .after = geneid)
      } else {
        fpfc_tibble <- final_path_fc_tibble
      }

      # Store result under the GO Term name
      labeled_gene_lists[[go_term]] <- fpfc_tibble
    }

    # Store the list of terms into the master list under the Comparison name
    unpacking_results[[comparison]] <- labeled_gene_lists
  }
  return(unpacking_results)
}





#Rendering Patchwork Images bypassing R Windows
#Helper function used to render the spatial feature plots and bypass RStudip window contrainsts
gg_patchwork <- function(plot, filename, width = 8, height = 6, dpi = 300, ...) {
  if (!grepl("\\.png$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".png")
  }
  # Open a device and print the plot explicitly to bypass RStudio window constraints
  grDevices::png(filename, width = width, height = height, units = "in", res = dpi)
  print(plot)  # works for ggplot OR patchwork
  dev.off()
  message("Saved: ", normalizePath(filename))
}




run_de_pipeline <- function(dds_object,
                            edgeR_results,
                            #controls the intial filtering of DEGs directly after the edgeR differential expression

                            deg_log2fc_thresh = 1,
                            #controls the intial p.val thresh of DEGs directly after the edgeR differential expression

                            deg_pval_adj_thresh = 0.05,
                            # controls the volcano plot's lower x-axis limit on the image
                            vol_plot_lower_xlim = -10,
                            #controls the volcano plot's upper x-axis limit on the image
                            vol_plot_upper_xlim = 10,
                            # ticker step of the plot on both x and y axis
                            cartesian_step = 5,
                            #over-representation test pvalue threshold
                            ora_pval_adj_threshold = 0.05) {

  require(glue)
  require(clusterProfiler)
  require(ggplot2)
  require(AnnotationDbi)
  require(DESeq2)
  #save vector the length of the edgeR results with comparisons
  diffexp_results <- vector("list", length(edgeR_results))
  names(diffexp_results) <- names(edgeR_results)

  # for each comparison
  for (i in seq_along(edgeR_results)) {

    volcano_res <- volcano_plot(
      edgeR_results[[i]],
      log2fc_thresh = deg_log2fc_thresh,
      p.val_adj_thresh = deg_pval_adj_thresh,
      lower_xlim = vol_plot_lower_xlim,
      upper_xlim = vol_plot_upper_xlim,
      step = cartesian_step
    )
    #for each comparison within the edgeR diff exp analysis
    # this function will save the filtered DEG table

    #save directly to the initialized vector
    diffexp_results[[i]] <- list(
      #saving the first subcript content of the volcano_res function
      data = volcano_res[[1]],

      vol_plot = volcano_res[[2]] +
        labs(title = glue::glue(
          "Volcano Plot: {names(edgeR_results)[i]}"
        ))

    )
  }

    # save the filtered Differential Expression Results
    filtered_results <- lapply(diffexp_results , function(x) x$data)

    #report for each comparison the number of up and down regulated DEGS
    for(i in seq_along(filtered_results)){

      num_sig_up_degs <- unique(filtered_results[[i]]$geneid[filtered_results[[i]]$log2foldchange > deg_log2fc_thresh & filtered_results[[i]]$p.val_adj < deg_pval_adj_thresh])

      num_sig_down_degs <-unique(filtered_results[[i]]$geneid[filtered_results[[i]]$log2foldchange < -deg_log2fc_thresh & filtered_results[[i]]$p.val_adj < deg_pval_adj_thresh])

      total <- length(num_sig_down_degs) + length(num_sig_up_degs)

      message(
        glue::glue(
          "Within Dataset {unique(colData(dds_object)$ds_origin)}, ",
          "Comparison: {names(filtered_results)[i]} has ",
          "{length(num_sig_up_degs)} significantly up-regulated DEGs and ",
          "{length(num_sig_down_degs)} down-regulated DEGs",
          "Making a total of {total} DEGS"
        )
      )
    }


  #Plot the Heatmap of the Differential Expression Anaalysis
  heatmap_visuals <- heatmap_function(dds_object, filtered_results)


  #saving visuals for Volcano plot
  results_vol_plots <- lapply(diffexp_results, function(x) x$vol_plot)


  #perform an EnrichGO ORA
  #ORA only takes in significant DEGs Answers the question:
  #Which pathways are overrepresented in my list of significant genes?

  ora_output_up <- ora_enrichgo(
    # Output of the differential expression results.
    filtered_results,
    # Which Direction of the DEGs (Upregulated("up")/Down("down"))
    direction = "up",
    #How strict do you want your DEGs entered into the ORA?
    log2fc_thresh = deg_log2fc_thresh,
    # what is the pval_thresh of the ORA
    pval_thresh = ora_pval_adj_threshold,
    #How many terms do you wantv visualized in your dotplots?
    top_terms = 15
  )

  ora_output_down <- ora_enrichgo(
    # Output of the differential expression results.
    filtered_results,
    # Which Direction of the DEGs (Upregulated("up")/Down("down"))
    direction = "down",
    #How strict do you want your DEGs entered into the ORA?
    log2fc_thresh = deg_log2fc_thresh,
    # what is the pval_thresh of the ORA
    pval_thresh = ora_pval_adj_threshold,
    #How many terms do you wantv visualized in your dotplots?
    top_terms = 15
  )

  #Unpacking ORA from enrichGO

  #Unpacking ORA from GO and align to your differential expression results for the input condition

  #ORA only takes in significant DEGs Answers the question:

  #Which pathways are overrepresented in my list of significant genes?


  # making a data table for each Description Term that is consider statistically significant within the context of the ORA

  significant_degs <- filtered_results

  key_down <- ora_output_down$enrichgo_results

  key_up <-ora_output_up$enrichgo_results


  unpacked_results_down <-  enrichgo_unpack_ver3(initial_table_list = significant_degs,
                                                 key = key_down,
                                                 pval_adj_threshold = ora_pval_adj_threshold,
                                                 convert = FALSE,
                                                 from_type = "ensembl",
                                                 to_type = "symbol",
                                                 ensembl_dataset = "hsapiens_gene_ensembl")

  unpacked_results_up <-  enrichgo_unpack_ver3(initial_table_list = significant_degs,
                                               key = key_up,
                                               pval_adj_threshold = ora_pval_adj_threshold,
                                               convert = FALSE,
                                               from_type = "ensembl",
                                               to_type = "symbol",
                                               ensembl_dataset = "hsapiens_gene_ensembl")
  return(
    list(
      filtered_results      = filtered_results,
      volcano_plots         = results_vol_plots,
      heatmap               = heatmap_visuals,
      ora_up_raw            = ora_output_up,
      ora_down_raw          = ora_output_down,
      unpacked_results_up   = unpacked_results_up,
      unpacked_results_down = unpacked_results_down
    )
  )
}


