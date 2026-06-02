# =============================================================================
# SECTION 0: SETUP — paths, packages, output folder
# =============================================================================

data_dir   <- "E:/manchester/study/semister2/71922se/a2/Assessment2_Data_GEOG71922/Beetles"
output_dir <- file.path(data_dir, "output")
dir.create(output_dir, showWarnings = FALSE)
setwd(data_dir)

# --- Step 0a: Install packages -----
# repos set explicitly to avoid pop-up asking which mirror to use
options(repos = c(CRAN = "https://cloud.r-project.org"))

pkg_install <- function(pkg) {
  if (!pkg %in% rownames(installed.packages())) {
    message("Installing: ", pkg)
    install.packages(pkg, dependencies = TRUE,
                     quiet = TRUE, verbose = FALSE)
  }
}

cran_pkgs <- c(
  "vegan", "devtools",
  "ggplot2", "ggrepel", "patchwork", "ggnewscale",
  "dplyr", "tidyr",
  "terra", "sf",
  "corrplot", "RColorBrewer", "viridis",
  "Hmsc", "coda",
  "car",
  "remotes"           # needed by devtools for vegetarian
)

invisible(lapply(cran_pkgs, pkg_install))

if (!"vegetarian" %in% rownames(installed.packages())) {
  message("Installing: vegetarian (via remotes)")
  remotes::install_version("vegetarian", version = "1.2",
                           quiet = TRUE, upgrade = "never")
}
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
# =============================================================================
# SECTION 1: DATA IMPORT AND INSPECTION
# =============================================================================

comm_raw <- read.csv("scot_beetle_community.csv", row.names = 1,
                     check.names = FALSE)
env_raw  <- read.csv("scot_beetle_env.csv",       row.names = 1,
                     check.names = FALSE)

# Remove any non-species columns from community matrix
# The "Sites" column sometimes appears in community CSVs as a label column
non_species <- c("Sites", "sites", "site", "Site")
comm_raw <- comm_raw[, !names(comm_raw) %in% non_species]
cat("Community matrix columns after cleaning:", ncol(comm_raw), "\n")
cat("Column names:", paste(names(comm_raw), collapse = ", "), "\n")

# Verify site order matches
stopifnot("Site order mismatch between community and env data!" =
            all(rownames(comm_raw) == rownames(env_raw)))

cat("Community matrix:", nrow(comm_raw), "sites x", ncol(comm_raw), "species\n")
cat("Environment table:", nrow(env_raw), "sites x", ncol(env_raw), "variables\n")
cat("Matrix sparsity:",
    round(sum(comm_raw == 0) / prod(dim(comm_raw)) * 100, 1), "%\n")
# Separate coordinates and management from predictors
coords     <- env_raw[, c("X", "Y")]
management <- env_raw[["Management"]]
env_cont   <- env_raw[, !names(env_raw) %in% c("X", "Y", "Sites", "Management")]

# --- Fig 0: Species abundance barplot (data check) --------------------------
sp_tot <- sort(colSums(comm_raw), decreasing = TRUE)
p_sp <- ggplot(data.frame(sp = names(sp_tot), ab = sp_tot),
               aes(x = reorder(sp, ab), y = ab)) +
  geom_col(fill = "#1D9E75", width = 0.7) +
  coord_flip() +
  labs(title = "Total abundance per species (84 sites)",
       x = NULL, y = "Total abundance") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "italic"))
savefig(p_sp, "fig00_species_abundance.png", w = 8, h = 6)
# =============================================================================
# SECTION 2: STAGE 1 — BETA DIVERSITY + LCBD  (W8)
# =============================================================================

cat("\n=== STAGE 1: Beta diversity + LCBD ===\n")

# Hellinger transformation (W8 + W9)
comm_hel <- decostand(comm_raw, method = "hellinger")

# --- 2a: Alpha diversity using vegetarian (W8 method) -----------------------
alpha_q0 <- d(comm_raw, lev = "alpha", q = 0)   # species richness
alpha_q1 <- d(comm_raw, lev = "alpha", q = 1)   # Shannon
alpha_q2 <- d(comm_raw, lev = "alpha", q = 2)   # Simpson

cat("\nAlpha diversity (mean across sites):\n")
cat("  Species richness (q=0):", round(alpha_q0, 3), "\n")
cat("  Shannon index   (q=1):", round(alpha_q1, 3), "\n")
cat("  Simpson index   (q=2):", round(alpha_q2, 3), "\n")
# --- 2b: Beta diversity decomposition using vegetarian (W8 method) ----------
beta_q0 <- d(comm_raw, lev = "beta", q = 0)
beta_q1 <- d(comm_raw, lev = "beta", q = 1)
gamma_q0 <- d(comm_raw, lev = "gamma", q = 0)

cat("\nBeta diversity:\n")
cat("  Multiplicative beta (q=0):", round(beta_q0, 3), "\n")
cat("  Multiplicative beta (q=1):", round(beta_q1, 3), "\n")
cat("  Gamma richness     (q=0):", round(gamma_q0, 3), "\n")

# Beta diversity across Hill numbers (q=0 to 5) — W8 plot
qN  <- 0:5
qDat <- sapply(qN, function(q) {
  out <- d(comm_raw, lev = "beta", q = q, boot = TRUE)
  c(beta = out$D.Value, se = out$StdErr)
})
beta_df <- data.frame(q = qN, beta = qDat["beta",], se = qDat["se",])

p_beta_q <- ggplot(beta_df, aes(x = q, y = beta)) +
  geom_ribbon(aes(ymin = beta - se, ymax = beta + se),
              fill = "#1D9E75", alpha = 0.25) +
  geom_line(colour = "#1D9E75", linewidth = 1) +
  geom_point(colour = "#085041", size = 3) +
  labs(title = "Beta diversity across Hill number orders",
       subtitle = "Scottish carabid beetle communities (n = 84 sites)",
       x = "Order of diversity measure (q)",
       y = "Multiplicative beta diversity") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
savefig(p_beta_q, "fig01a_beta_hill_numbers.png", w = 7, h = 5)
# Save diversity summary
div_summary <- data.frame(
  Metric  = c("Alpha richness (q=0)", "Alpha Shannon (q=1)",
              "Alpha Simpson (q=2)", "Beta (q=0)", "Beta (q=1)",
              "Gamma richness"),
  Value   = round(c(alpha_q0, alpha_q1, alpha_q2,
                    beta_q0, beta_q1, gamma_q0), 4)
)
write.csv(div_summary,
          file.path(output_dir, "results_diversity_summary.csv"),
          row.names = FALSE)
# --- 2c: LCBD using teacher's W8 hand-written function ----------------------
# Direct replication of the function_lcbd() from the W8 practical
function_lcbd <- function(x) {
  # Hellinger transformation
  spe1 <- decostand(x, method = "hellinger")
  # Empty matrix to store squared deviations
  ss_mat <- spe1
  ss_mat[] <- 0
  # For each species: squared deviation from mean abundance across sites
  for (i in 1:ncol(spe1)) {
    sp.i       <- spe1[, i]
    col_mean   <- mean(sp.i)
    beta.i     <- sapply(sp.i, function(val) (val - col_mean)^2)
    ss_mat[, i] <- beta.i
  }
  ss_total   <- sum(ss_mat)
  site_LCBD  <- rowSums(ss_mat) / ss_total
  return(site_LCBD)
}
lcbd_vals <- function_lcbd(comm_raw)

cat("\nLCBD values computed for", length(lcbd_vals), "sites\n")
cat("Min:", round(min(lcbd_vals), 5),
    "  Max:", round(max(lcbd_vals), 5),
    "  Sum:", round(sum(lcbd_vals), 4), "\n")

# Identify sites with high LCBD (above mean + 1SD — same threshold logic as W8)
lcbd_threshold <- mean(lcbd_vals) + sd(lcbd_vals)
high_lcbd      <- lcbd_vals > lcbd_threshold
cat("Sites with high LCBD (>mean+1SD):", sum(high_lcbd), "\n")
# Save LCBD table
lcbd_df <- data.frame(
  site       = rownames(comm_raw),
  X          = coords$X,
  Y          = coords$Y,
  LCBD       = lcbd_vals,
  high_LCBD  = high_lcbd,
  Management = management
)
write.csv(lcbd_df,
          file.path(output_dir, "results_LCBD_table.csv"),
          row.names = FALSE)
# --- Export: LCBD as spatial vector points (shapefile) ----------------------
# Create sf point object in British National Grid (EPSG:27700)
lcbd_sf <- st_as_sf(lcbd_df,
                    coords = c("X", "Y"),
                    crs    = 27700)

# Split into high and low LCBD layers for easy GIS use
lcbd_high <- lcbd_sf[lcbd_sf$high_LCBD == TRUE,  ]
lcbd_low  <- lcbd_sf[lcbd_sf$high_LCBD == FALSE, ]

# Save as shapefiles
st_write(lcbd_sf,
         file.path(output_dir, "lcbd_all_points.shp"),
         delete_layer = TRUE, quiet = TRUE)
st_write(lcbd_high,
         file.path(output_dir, "lcbd_high_points.shp"),
         delete_layer = TRUE, quiet = TRUE)
st_write(lcbd_low,
         file.path(output_dir, "lcbd_low_points.shp"),
         delete_layer = TRUE, quiet = TRUE)

message("Saved -> lcbd_all_points.shp")
message("Saved -> lcbd_high_points.shp  (", nrow(lcbd_high), " sites)")
message("Saved -> lcbd_low_points.shp   (", nrow(lcbd_low),  " sites)")
# --- Export: LCBD as raster (TIFF) ------------------------------------------
# Rasterise LCBD point values to grid cells (nearest neighbour)
lcbd_vect <- vect(lcbd_sf)

if (file.exists(lcm_path)) {
  lcm_ref   <- rast(lcm_path)
  lcbd_rast <- rast(ext(lcm_ref),
                    resolution = res(lcm_ref)[1] * 10,
                    crs        = crs(lcm_ref))
} else {
  lcbd_rast <- rast(ext(lcbd_vect) + 10000,
                    resolution = 5000,
                    crs        = "EPSG:27700")
}

lcbd_tif <- rasterize(lcbd_vect,
                      lcbd_rast,
                      field = "LCBD",
                      fun   = "mean")

writeRaster(lcbd_tif,
            file.path(output_dir, "lcbd_raster.tif"),
            overwrite = TRUE)
message("Saved -> lcbd_raster.tif")
# --- Fig 1b: LCBD bubble map --------------------
png(file.path(output_dir, "fig01b_LCBD_map_baseR.png"),
    width = 1800, height = 1800, res = 250)
plot(coords,
     cex.axis = 0.8,
     pch      = 21,
     col      = "black",
     bg       = ifelse(high_lcbd, "#D85A30", "#1D9E75"),
     cex      = lcbd_vals * 120,
     main     = "Local Contribution to Beta Diversity (LCBD)",
     xlab     = "Easting (BNG)",
     ylab     = "Northing (BNG)")
legend("topright",
       legend = c("High LCBD (>mean+1SD)", "Low LCBD"),
       pt.bg  = c("#D85A30", "#1D9E75"),
       pch    = 21, pt.cex = 1.5, bty = "n", cex = 0.8)
dev.off()
message("Saved → fig01b_LCBD_map_baseR.png")

# --- Fig 1c: LCBD bubble map with land cover --------------------------------
lcm_path <- file.path(data_dir, "LCMUK_2000.tif")
lcm      <- rast(lcm_path)
lcm_agg  <- aggregate(lcm, fact = 10, fun = "modal")
lcm_df   <- as.data.frame(lcm_agg, xy = TRUE)
names(lcm_df)[3] <- "landcover"
lcm_df$landcover <- as.factor(lcm_df$landcover)
rm(lcm, lcm_agg); gc()

lcm_cols <- c("0" = "#FFFFFF", "1" = "#27500A", "2" = "#085041",
              "3" = "#FAC775", "4" = "#9FE1CB", "5" = "#B5D4F4",
              "6" = "#888780", "7" = "#5DCAA5",
              "8" = "#378ADD", "9" = "#185FA5")

p_lcbd_lcm <- ggplot() +
  geom_raster(data = na.omit(lcm_df),
              aes(x = x, y = y, fill = landcover),
              show.legend = FALSE) +
  scale_fill_manual(values = lcm_cols, na.value = "white") +
  ggnewscale::new_scale_fill() +
  geom_point(data = lcbd_df,
             aes(x = X, y = Y,
                 size   = LCBD,
                 colour = LCBD,
                 shape  = high_LCBD),
             alpha = 0.85) +
  scale_colour_gradient(low  = "#9FE1CB", high = "#085041",
                        name = "LCBD") +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18),
                     labels = c("Low LCBD", "High LCBD"),
                     name   = "") +
  coord_equal() +
  labs(title    = "LCBD overlaid on UK Land Cover Map 2000",
       subtitle = paste0("Diamonds: sites with LCBD > mean + 1 SD (n = ",
                         sum(high_lcbd), ")"),
       x = "Easting (BNG)", y = "Northing (BNG)",
       caption = "Coordinates: British National Grid EPSG:27700") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

savefig(p_lcbd_lcm, "fig01c_LCBD_landcover.png", w = 9, h = 8)
# =============================================================================
# SECTION 3: STAGE 2 — CCA CONSTRAINED ORDINATION  (W9)
# =============================================================================

cat("\n=== STAGE 2: CCA constrained ordination ===\n")

# --- 3a: DCA gradient length test — determines CCA vs RDA ------------------
dca_res   <- decorana(comm_raw)
print(dca_res)

# Gradient length of axis 1
dca_ax1 <- diff(range(scores(dca_res, display = "sites")[, 1]))
cat("\nDCA Axis 1 gradient length:", round(dca_ax1, 3), "SD units\n")

if (dca_ax1 > 3) {
  cat("→ > 3 SD: UNIMODAL response → CCA selected\n")
  use_cca <- TRUE
} else {
  cat("→ < 3 SD: LINEAR response → RDA selected\n")
  use_cca <- FALSE
}

# Save DCA result
dca_df <- data.frame(
  Axis       = paste0("DCA", 1:4),
  Eigenvalue = round(dca_res$evals, 4),
  GradLength = round(apply(scores(dca_res, display = "sites"), 2,
                           function(x) diff(range(x))), 3)
)
write.csv(dca_df,
          file.path(output_dir, "results_DCA_gradient.csv"),
          row.names = FALSE)