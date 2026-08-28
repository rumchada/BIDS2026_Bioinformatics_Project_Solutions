#' Over-Representation Analysis Using enrichGO
#'
#' Performs Gene Ontology (GO) over-representation analysis (ORA) using
#' \code{\link[clusterProfiler]{enrichGO}} on differential expression results.
#' Genes can be analyzed separately based on their direction of differential
#' expression (upregulated or downregulated).
#'
#'@name ora_enrichgo
#'@aliases ora_enrichgo
#'
#'
#'@param results list of tables or a table of differential expression results
#'
#'@param direction Which direction +log2fc -log2fc should be looked at ? e.g. \code{"up"}, \code{"down"}
#'
#'@param pval_thesh after running fisher exact test for ORA, what should the distribution of pvalues be?
#'
#'@param top_terms How many resulting terms should be display in the visuals.
#'
#'@param OrgDb rganism database. Default is human: `org.Hs.eg.db`.
#'   Other organism annotation databases are available from
#'   https://bioconductor.org/packages/3.23/data/annotation/.
#'
#'@param keyType which type of geneid are being entered into the cluster
#'
#'@param ont which ontological DB will you be using ? options \code{c("BP", "CC", "MF")}
#'
#' @return A list containing:
#' \itemize{
#' \item \code{enrichgo_results}: Results from the GO over-representation
#' analysis returned by \code{\link[clusterProfiler]{enrichGO}}.
#' \item \code{visuals}: A collection of visualizations summarizing the
#' enriched GO terms, including dotplots and barplots.
#' }
#'
#' @details
#' Genes are first separated according to their direction of differential
#' expression and filtered using the supplied log2 fold-change and p-value
#' thresholds. Over-representation analysis is then performed using
#' \code{\link[clusterProfiler]{enrichGO}}.
#'
#'@examples
#'\dontrun{
#'     enrich_ora_results <- ora_enrichgo(filtered_results,
#'                                  direction = c("up", "down"),
#'                                      log2fc_thresh = 1,
#'                                       pval_thresh = 0.05,
#'                                       top_terms = 15,
#'                                       OrgDb = org.Hs.eg.db,
#'                                      keyType = "ENSEMBL",
#'                                        ont = "BP")
#'}


ora_enrichgo <- function(filtered_results,
                         direction = c("up", "down"),
                         log2fc_thresh = 1,
                         pval_thresh = 0.05,
                         top_terms = 15,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENSEMBL",
                         ont = "BP"){

  require(clusterProfiler)
  require(BiocManager)


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

