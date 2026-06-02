# --- Step 0b: Load all packages ---------------------------------------------
suppressPackageStartupMessages({
  library(vegan)
  library(devtools)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(ggnewscale)
  library(dplyr)
  library(tidyr)
  library(terra)
  library(sf)
  library(corrplot)
  library(RColorBrewer)
  library(viridis)
  library(Hmsc)
  library(coda)
  library(car)
  library(vegetarian)
})
# --- Step 0c: Helper function -----------------------------------------------
savefig <- function(p, name, w = 9, h = 7, dpi = 300) {
  ggsave(file.path(output_dir, name), plot = p,
         width = w, height = h, dpi = dpi, bg = "white")
  message("Saved -> ", name)
}

cat("\n=== Setup complete. Working dir:", getwd(), "===\n")