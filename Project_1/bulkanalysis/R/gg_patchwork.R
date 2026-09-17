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
