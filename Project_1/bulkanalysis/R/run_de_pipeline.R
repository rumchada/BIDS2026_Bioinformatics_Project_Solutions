#' @title Run all Downstream analysis pipeline
#'
#' @description
#' Uses the edgeR differentiall expression results and the formatted Deseq2Dataset
#' Performs volcano_plot.R which filters for the chosen log2fc and pval_adj
#'
#'
#'
#' @param edgeR_results dataframe of log2fc results from edge
#'
#' @param deg_log2fc_thresh absolute value cutoff of log2fc
#'
#' @param vol_plot_lower_xlim Numeric value specifying the lower limit of the
#'   x-axis. Should be less than 0. For example, `-5`.
#'
#' @param vol_plot_upper_xlim Numeric value specifying the upper limit of the
#'   x-axis. Should be greater than 0. For example, `5`.
#'
#' @param cartesian_step Numeric value specifying the interval between x-axis ticks.
#'
#' @param ora_pval_adj_threshold pvalue cutoff of the over-representation test
#'
#'@return
#'A list of the following objects
#'\itemize{
#'    \item{filtered_results}: The initial filtering of DEGs log2fc and pval_adj_threshold
#'    \item{volcano_plots}: visualization of the filtered results on a volcano plot
#'    \item{ora_up_raw}: Result of the over-representation analysis on upregulated DEGs via clusterProfiler
#'    \item{ora_down_raw}: Result of the over-represenation analysis on downregulated DEGs via clusterProfiler
#'    \item{unpacked_results_up}: unpacked results table consisting of a nested list of ORA results on the upregulated DEGs
#'    \item{unpacked_results_down}: unpacked results table consisting of a nested list of ORA results on the down-regulated DEGs
#'}
#'@export
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
                            ora_convert_ids = FALSE,
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
                                                 convert = ora_convert_ids,
                                                 from_type = "ensembl",
                                                 to_type = "symbol",
                                                 ensembl_dataset = "hsapiens_gene_ensembl")

  unpacked_results_up <-  enrichgo_unpack_ver3(initial_table_list = significant_degs,
                                               key = key_up,
                                               pval_adj_threshold = ora_pval_adj_threshold,
                                               convert = ora_convert_ids,
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
