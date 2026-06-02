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
# --- 3b: VIF screening using vif.cca (W9 compatible) -----------------------
# Build a full CCA first, then check VIF per variable
# Iteratively remove highest VIF variable until all < 10

cat("\nVIF screening (threshold = 10):\n")
env_work <- env_cont
vif_log  <- data.frame()

repeat {
  cca_vif  <- cca(comm_raw ~ ., data = env_work)
  vif_vals <- vif.cca(cca_vif)
  max_vif  <- max(vif_vals)
  if (max_vif < 10) break
  drop_var <- names(which.max(vif_vals))
  cat("  Removing:", drop_var, "  VIF =", round(max_vif, 2), "\n")
  vif_log  <- rbind(vif_log,
                    data.frame(variable = drop_var,
                               VIF      = round(max_vif, 2)))
  env_work <- env_work[, names(env_work) != drop_var, drop = FALSE]
}

cat("Retained variables (n =", ncol(env_work), "):",
    paste(names(env_work), collapse = ", "), "\n")

final_vif <- round(vif.cca(cca(comm_raw ~ ., data = env_work)), 3)
cat("Final VIF values:\n"); print(final_vif)

# Save VIF results
write.csv(vif_log,
          file.path(output_dir, "results_VIF_dropped.csv"),
          row.names = FALSE)
write.csv(data.frame(variable = names(final_vif), VIF = final_vif),
          file.path(output_dir, "results_VIF_retained.csv"),
          row.names = FALSE)

# Assign final env dataset
env_sel <- env_work
# --- 3c: Build CCA model-----------
if (use_cca) {
  cca_mod <- cca(comm_raw ~ ., data = env_sel)
} else {
  cca_mod <- rda(comm_hel ~ ., data = env_sel)
}

method_label <- if (use_cca) "CCA" else "RDA"
cat("\n---", method_label, "model summary ---\n")
print(summary(cca_mod))

# Constrained variance
constr_pct <- cca_mod$CCA$tot.chi / cca_mod$tot.chi * 100
r2_adj     <- RsquareAdj(cca_mod)$adj.r.squared
cat("Constrained variance:", round(constr_pct, 2), "%\n")
cat("Adjusted R²:         ", round(r2_adj, 4), "\n")

# Permutation tests
set.seed(42)
anova_overall <- anova.cca(cca_mod, permutations = 999)
anova_axes    <- anova.cca(cca_mod, by = "axis",  permutations = 999)
anova_terms   <- anova.cca(cca_mod, by = "terms", permutations = 999)

cat("\nOverall model test:\n"); print(anova_overall)
cat("\nBy axis:\n");            print(anova_axes)
cat("\nBy term:\n");            print(anova_terms)

# envfit passive fitting
set.seed(42)
ef <- envfit(cca_mod, env_sel, perm = 9999)
print(ef)

# Save stats
cca_stats <- data.frame(
  Metric = c("Method", "n_sites", "n_species", "n_predictors",
             "Total_inertia", "Constrained_inertia_pct",
             "Adjusted_R2", "Model_F", "Model_p"),
  Value  = c(method_label, nrow(comm_raw), ncol(comm_raw),
             ncol(env_sel),
             round(cca_mod$tot.chi, 4),
             round(constr_pct, 2),
             round(r2_adj, 4),
             round(anova_overall$F[1], 3),
             round(anova_overall$`Pr(>F)`[1], 4))
)
write.csv(cca_stats,
          file.path(output_dir, "results_CCA_statistics.csv"),
          row.names = FALSE)
write.csv(as.data.frame(anova_terms),
          file.path(output_dir, "results_CCA_terms_anova.csv"))
write.csv(as.data.frame(anova_axes),
          file.path(output_dir, "results_CCA_axes_anova.csv"))
# --- Fig 2: CCA triplot ------------------------------------------
ax1 <- paste0(method_label, "1")
ax2 <- paste0(method_label, "2")

site_sc <- as.data.frame(scores(cca_mod, display = "sites",   scaling = 2))
sp_sc   <- as.data.frame(scores(cca_mod, display = "species", scaling = 2))
bp_sc   <- as.data.frame(scores(cca_mod, display = "bp",      scaling = 2))

site_sc$Management <- management
site_sc$LCBD       <- lcbd_vals
site_sc$high_lcbd  <- high_lcbd

ax1_pct <- round(summary(cca_mod)$cont$importance[2, 1] * 100, 1)
ax2_pct <- round(summary(cca_mod)$cont$importance[2, 2] * 100, 1)

p_triplot <- ggplot() +
  # Site points coloured by management, sized by LCBD
  geom_point(data = site_sc,
             aes(x      = .data[[ax1]],
                 y      = .data[[ax2]],
                 colour = Management,
                 size   = LCBD),
             alpha = 0.8) +
  # Ring around high-LCBD sites
  geom_point(data = subset(site_sc, high_lcbd),
             aes(x = .data[[ax1]], y = .data[[ax2]]),
             shape = 1, size = 6, colour = "#D85A30", stroke = 1.2) +
  # Species labels (italic, purple — W9 style)
  ggrepel::geom_text_repel(data = sp_sc,
                           aes(x     = .data[[ax1]] * 0.8,
                               y     = .data[[ax2]] * 0.8,
                               label = rownames(sp_sc)),
                           size = 2.8, colour = "#534AB7",
                           fontface = "italic", max.overlaps = 25) +
  # Environmental biplot arrows (W9 style)
  geom_segment(data = bp_sc,
               aes(x    = 0, y    = 0,
                   xend = .data[[ax1]],
                   yend = .data[[ax2]]),
               arrow    = arrow(length = unit(0.22, "cm"), type = "closed"),
               colour   = "black", linewidth = 0.8) +
  ggrepel::geom_text_repel(data = bp_sc,
                           aes(x     = .data[[ax1]] * 1.15,
                               y     = .data[[ax2]] * 1.15,
                               label = rownames(bp_sc)),
                           size = 3.5, colour = "black",
                           fontface = "bold", max.overlaps = 30) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey70", linewidth = 0.3) +
  scale_colour_viridis_c(name = "Management\nintensity", option = "D") +
  scale_size_continuous(name  = "LCBD", range = c(2, 6)) +
  labs(
    title    = paste0(method_label,
                      " triplot — carabid community vs environment"),
    subtitle = paste0("Constrained variance: ", round(constr_pct, 1),
                      "%  |  adj.R² = ", round(r2_adj, 3),
                      "  |  p = ",
                      round(anova_overall$`Pr(>F)`[1], 3)),
    x        = paste0(ax1, " (", ax1_pct, "%)"),
    y        = paste0(ax2, " (", ax2_pct, "%)"),
    caption  = paste0("Arrows = environmental predictors  |  ",
                      "Italic purple = species  |  ",
                      "Circle = high-LCBD site")
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40"))
savefig(p_triplot, "fig02_CCA_triplot.png", w = 11, h = 8)
# --- Fig 3: NMDS exploratory plot --------------------------------
set.seed(42)
nmds <- metaMDS(comm_raw, distance = "bray", k = 2,
                trymax = 100, trace = FALSE)
cat("\nNMDS stress:", round(nmds$stress, 4),
    if (nmds$stress < 0.1) "(excellent)"
    else if (nmds$stress < 0.2) "(acceptable)" else "(poor)", "\n")

nmds_sc  <- as.data.frame(scores(nmds)$sites)
nmds_sp  <- as.data.frame(scores(nmds)$species)
nmds_sc$Type <- cut(management, breaks = 3,
                    labels = c("Low", "Medium", "High"))

# Convex hull per management group 
hull_df <- nmds_sc |>
  group_by(Type) |>
  slice(chull(NMDS1, NMDS2))

# =============================================================================
# SECTION 4: STAGE 3 — PARTIAL CCA + VARIANCE PARTITIONING
# =============================================================================

cat("\n=== STAGE 3: Partial CCA + variance partitioning ===\n")

# --- 4a: Extract spatial vectors using pcnm ---------------------
dist_mat <- as.matrix(dist(coords[, c("X", "Y")]))
pcnm_res <- pcnm(dist_mat)

# Keep only positive eigenvalue PCNM vectors 
pos_pcnm   <- which(pcnm_res$values > 0)
pcnm_scores <- as.data.frame(scores(pcnm_res)[, pos_pcnm])
colnames(pcnm_scores) <- paste0("PCNM", pos_pcnm)

cat("Total PCNM vectors with positive eigenvalues:", ncol(pcnm_scores), "\n")

# Select significant PCNM vectors via forward selection
if (use_cca) {
  pcnm_null <- cca(comm_raw ~ 1, data = pcnm_scores)
  pcnm_full <- cca(comm_raw ~ ., data = pcnm_scores)
} else {
  pcnm_null <- rda(comm_hel ~ 1, data = pcnm_scores)
  pcnm_full <- rda(comm_hel ~ ., data = pcnm_scores)
}

set.seed(42)
pcnm_step <- ordistep(pcnm_null,
                      scope        = formula(pcnm_full),
                      direction    = "forward",
                      permutations = 999,
                      trace        = FALSE)

sel_pcnm_vars <- attr(terms(pcnm_step), "term.labels")
cat("Selected PCNM vectors:", paste(sel_pcnm_vars, collapse = ", "), "\n")

space_sel <- pcnm_scores[, sel_pcnm_vars, drop = FALSE]

# envfit vectors
set.seed(42)
ef_nmds    <- envfit(nmds, env_sel, perm = 9999)
ef_vec     <- as.data.frame(scores(ef_nmds, display = "vectors"))
ef_vec$var <- rownames(ef_vec)
ef_sig     <- ef_vec[ef_nmds$vectors$pvals < 0.05, ]

p_nmds <- ggplot(nmds_sc, aes(x = NMDS1, y = NMDS2)) +
  geom_polygon(data = hull_df,
               aes(fill = Type, group = Type),
               alpha = 0.2) +
  geom_point(aes(colour = Type, shape = Type), size = 3.5) +
  ggrepel::geom_text_repel(data = nmds_sp,
                           aes(x = NMDS1, y = NMDS2,
                               label = rownames(nmds_sp)),
                           size = 2.6, colour = "#534AB7",
                           fontface = "italic", max.overlaps = 20) +
  geom_segment(data = ef_sig,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow  = arrow(length = unit(0.2, "cm")),
               colour = "black", linewidth = 0.8) +
  ggrepel::geom_text_repel(data = ef_sig,
                           aes(x = NMDS1 * 1.1, y = NMDS2 * 1.1,
                               label = var),
                           colour = "black", size = 3.5,
                           fontface = "bold") +
  scale_colour_manual(values = c("#1D9E75", "#534AB7", "#D85A30"),
                      name = "Management") +
  scale_fill_manual(values   = c("#1D9E75", "#534AB7", "#D85A30"),
                    name = "Management") +
  scale_shape_manual(values  = c(16, 17, 15), name = "Management") +
  labs(title    = "NMDS ordination — carabid community composition",
       subtitle = paste0("Bray-Curtis dissimilarity  |  Stress = ",
                         round(nmds$stress, 3),
                         "  |  Arrows: sig. env vectors (p < 0.05)"),
       x = "NMDS1", y = "NMDS2") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.border = element_rect(colour = "black",
                                    fill = NA, linewidth = 0.8))
savefig(p_nmds, "fig03_NMDS.png", w = 10, h = 7)