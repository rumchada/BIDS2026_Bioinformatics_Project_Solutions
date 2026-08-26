#' Volcano Plot Visualization
#'
#' Uses `ggplot2`, `ggrepel`, and `tidyverse` to construct a volcano
#' plot from differential expression results.
#'
#' @name volcano_plot
#' @alias volcano_plot
#'
#' @param diffexp_df A tibble or data frame containing differential
#'   expression results. Must contain, at minimum, `p.val_adj` and
#'   `log2FC` columns.
#'
#' @param log2fc_thresh Numeric threshold for the absolute log2 fold-change
#'   used to identify differentially expressed genes for display.
#'
#' @param p.val_adj_thresh Numeric threshold for the adjusted p-value used
#'   to identify statistically significant differentially expressed genes.
#'
#' @param lower_xlim Numeric value specifying the lower limit of the
#'   x-axis. Should be less than 0. For example, `-5`.
#'
#' @param upper_xlim Numeric value specifying the upper limit of the
#'   x-axis. Should be greater than 0. For example, `5`.
#'
#' @param step Numeric value specifying the interval between x-axis ticks.
#'
#' @param dot_size Numeric value specifying the size of points in the
#'   volcano plot.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `plot`: A ggplot object containing the volcano plot.
#'     \item `data`: A data frame containing the differential expression
#'       results used to construct the plot.
#'   }


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
