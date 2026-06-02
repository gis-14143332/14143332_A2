# =============================================================================
# GEOG71922 Spatial Ecology — Assessment 2
# Scottish Ground Beetle Community Analysis
# Environmental and spatial drivers of carabid community composition
# Student ID: [YOUR_STUDENT_ID]
# =============================================================================


# =============================================================================
# SECTION 0: SETUP
# =============================================================================

data_dir   <- "E:/manchester/study/semister2/71922se/a2/Assessment2_Data_GEOG71922/Beetles"
output_dir <- file.path(data_dir, "output")
dir.create(output_dir, showWarnings = FALSE)
setwd(data_dir)

# Set CRAN mirror to avoid interactive pop-up during package installation
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install a package only if not already present
pkg_install <- function(pkg) {
  if (!pkg %in% rownames(installed.packages())) {
    message("Installing: ", pkg)
    install.packages(pkg, dependencies = TRUE, quiet = TRUE, verbose = FALSE)
  }
}

cran_pkgs <- c(
  "vegan", "devtools",
  "ggplot2", "ggrepel", "patchwork", "ggnewscale",
  "dplyr", "tidyr",
  "terra", "sf",
  "corrplot", "RColorBrewer", "viridis",
  "Hmsc", "coda",
  "car", "remotes"
)

invisible(lapply(cran_pkgs, pkg_install))

# vegetarian is not available on CRAN for R 4.x; install via remotes
if (!"vegetarian" %in% rownames(installed.packages())) {
  remotes::install_version("vegetarian", version = "1.2",
                           quiet = TRUE, upgrade = "never")
}

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

# Helper function to save ggplot figures to the output directory
savefig <- function(p, name, w = 9, h = 7, dpi = 300) {
  ggsave(file.path(output_dir, name), plot = p,
         width = w, height = h, dpi = dpi, bg = "white")
  message("Saved -> ", name)
}

# Path to land cover raster — defined here so it is available throughout
lcm_path <- file.path(data_dir, "LCMUK_2000.tif")


# =============================================================================
# SECTION 1: DATA IMPORT
# =============================================================================

comm_raw <- read.csv("scot_beetle_community.csv", row.names = 1,
                     check.names = FALSE)
env_raw  <- read.csv("scot_beetle_env.csv", row.names = 1,
                     check.names = FALSE)

# Remove any non-species label columns that may be present in the community CSV
comm_raw <- comm_raw[, !names(comm_raw) %in% c("Sites", "sites", "site", "Site")]

# Confirm that sites are in the same order in both matrices
stopifnot(all(rownames(comm_raw) == rownames(env_raw)))

# Separate spatial coordinates and management intensity from environmental predictors
coords     <- env_raw[, c("X", "Y")]
management <- env_raw[["Management"]]
env_cont   <- env_raw[, !names(env_raw) %in% c("X", "Y", "Sites", "Management")]

# Fig 00: species total abundance — used to inspect data sparsity
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
# SECTION 2: BETA DIVERSITY + LCBD
# =============================================================================

# Hellinger transformation standardises abundance by site totals,
# reducing the influence of dominant species
comm_hel <- decostand(comm_raw, method = "hellinger")

# Alpha diversity at three Hill number orders (q=0 richness, q=1 Shannon, q=2 Simpson)
alpha_q0 <- d(comm_raw, lev = "alpha", q = 0)
alpha_q1 <- d(comm_raw, lev = "alpha", q = 1)
alpha_q2 <- d(comm_raw, lev = "alpha", q = 2)

# Multiplicative beta diversity and gamma richness
beta_q0  <- d(comm_raw, lev = "beta",  q = 0)
beta_q1  <- d(comm_raw, lev = "beta",  q = 1)
gamma_q0 <- d(comm_raw, lev = "gamma", q = 0)

# Beta diversity profile across Hill numbers q = 0 to 5
qN   <- 0:5
qDat <- sapply(qN, function(q) {
  out <- d(comm_raw, lev = "beta", q = q, boot = TRUE)
  c(beta = out$D.Value, se = out$StdErr)
})
beta_df <- data.frame(q = qN, beta = qDat["beta", ], se = qDat["se", ])

p_beta_q <- ggplot(beta_df, aes(x = q, y = beta)) +
  geom_ribbon(aes(ymin = beta - se, ymax = beta + se),
              fill = "#1D9E75", alpha = 0.25) +
  geom_line(colour = "#1D9E75", linewidth = 1) +
  geom_point(colour = "#085041", size = 3) +
  labs(title    = "Beta diversity across Hill number orders",
       subtitle = "Scottish carabid beetle communities (n = 84 sites)",
       x = "Order of diversity measure (q)",
       y = "Multiplicative beta diversity") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
savefig(p_beta_q, "fig01a_beta_hill_numbers.png", w = 7, h = 5)

div_summary <- data.frame(
  Metric = c("Alpha richness (q=0)", "Alpha Shannon (q=1)",
             "Alpha Simpson (q=2)", "Beta (q=0)", "Beta (q=1)",
             "Gamma richness"),
  Value  = round(c(alpha_q0, alpha_q1, alpha_q2,
                   beta_q0, beta_q1, gamma_q0), 4)
)
write.csv(div_summary,
          file.path(output_dir, "results_diversity_summary.csv"),
          row.names = FALSE)

# LCBD: Local Contribution to Beta Diversity
# Computed as the site-level sum of squared deviations from species means,
# divided by the total sum of squares across the Hellinger-transformed matrix
function_lcbd <- function(x) {
  spe1   <- decostand(x, method = "hellinger")
  ss_mat <- spe1
  ss_mat[] <- 0
  for (i in 1:ncol(spe1)) {
    sp.i        <- spe1[, i]
    col_mean    <- mean(sp.i)
    beta.i      <- sapply(sp.i, function(val) (val - col_mean)^2)
    ss_mat[, i] <- beta.i
  }
  ss_total  <- sum(ss_mat)
  site_LCBD <- rowSums(ss_mat) / ss_total
  return(site_LCBD)
}

lcbd_vals      <- function_lcbd(comm_raw)
lcbd_threshold <- mean(lcbd_vals) + sd(lcbd_vals)
high_lcbd      <- lcbd_vals > lcbd_threshold

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

# Export LCBD as spatial vector layers (BNG EPSG:27700)
lcbd_sf   <- st_as_sf(lcbd_df, coords = c("X", "Y"), crs = 27700)
lcbd_high <- lcbd_sf[lcbd_sf$high_LCBD == TRUE,  ]
lcbd_low  <- lcbd_sf[lcbd_sf$high_LCBD == FALSE, ]

st_write(lcbd_sf,
         file.path(output_dir, "lcbd_all_points.shp"),
         delete_layer = TRUE, quiet = TRUE)
st_write(lcbd_high,
         file.path(output_dir, "lcbd_high_points.shp"),
         delete_layer = TRUE, quiet = TRUE)
st_write(lcbd_low,
         file.path(output_dir, "lcbd_low_points.shp"),
         delete_layer = TRUE, quiet = TRUE)

# Export LCBD values as a raster (nearest-neighbour rasterisation)
lcbd_vect <- vect(lcbd_sf)
lcm_ref   <- rast(lcm_path)
lcbd_rast <- rast(ext(lcm_ref), resolution = res(lcm_ref)[1] * 10,
                  crs = crs(lcm_ref))
lcbd_tif  <- rasterize(lcbd_vect, lcbd_rast, field = "LCBD", fun = "mean")
writeRaster(lcbd_tif,
            file.path(output_dir, "lcbd_raster.tif"),
            overwrite = TRUE)

# Fig 01b: base-R bubble map of LCBD values
png(file.path(output_dir, "fig01b_LCBD_map_baseR.png"),
    width = 1800, height = 1800, res = 250)
plot(coords,
     cex.axis = 0.8, pch = 21, col = "black",
     bg  = ifelse(high_lcbd, "#D85A30", "#1D9E75"),
     cex = lcbd_vals * 120,
     main = "Local Contribution to Beta Diversity (LCBD)",
     xlab = "Easting (BNG)", ylab = "Northing (BNG)")
legend("topright",
       legend = c("High LCBD (>mean+1SD)", "Low LCBD"),
       pt.bg  = c("#D85A30", "#1D9E75"),
       pch = 21, pt.cex = 1.5, bty = "n", cex = 0.8)
dev.off()

# Fig 01c: LCBD overlaid on land cover map
# Raster aggregated by factor 10 to reduce memory usage before conversion
lcm_agg <- aggregate(rast(lcm_path), fact = 10, fun = "modal")
lcm_df  <- as.data.frame(lcm_agg, xy = TRUE)
names(lcm_df)[3] <- "landcover"
lcm_df$landcover <- as.factor(lcm_df$landcover)
rm(lcm_agg); gc()

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
             aes(x = X, y = Y, size = LCBD,
                 colour = LCBD, shape = high_LCBD),
             alpha = 0.85) +
  scale_colour_gradient(low = "#9FE1CB", high = "#085041", name = "LCBD") +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18),
                     labels = c("Low LCBD", "High LCBD"), name = "") +
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
# SECTION 3: CCA CONSTRAINED ORDINATION
# =============================================================================

# DCA gradient length test: determines whether to use CCA (unimodal) or RDA (linear)
# Axis 1 length > 3 SD indicates unimodal species responses → CCA is appropriate
dca_res <- decorana(comm_raw)
dca_ax1 <- diff(range(scores(dca_res, display = "sites")[, 1]))

use_cca      <- dca_ax1 > 3
method_label <- if (use_cca) "CCA" else "RDA"

dca_df <- data.frame(
  Axis       = paste0("DCA", 1:4),
  Eigenvalue = round(dca_res$evals, 4),
  GradLength = round(apply(scores(dca_res, display = "sites"), 2,
                           function(x) diff(range(x))), 3)
)
write.csv(dca_df,
          file.path(output_dir, "results_DCA_gradient.csv"),
          row.names = FALSE)

# VIF screening: iteratively remove the variable with highest VIF until all < 10
# This eliminates severe multicollinearity before ordination
env_work <- env_cont
vif_log  <- data.frame()

repeat {
  cca_vif  <- cca(comm_raw ~ ., data = env_work)
  vif_vals <- vif.cca(cca_vif)
  max_vif  <- max(vif_vals)
  if (max_vif < 10) break
  drop_var <- names(which.max(vif_vals))
  vif_log  <- rbind(vif_log,
                    data.frame(variable = drop_var,
                               VIF      = round(max_vif, 2)))
  env_work <- env_work[, names(env_work) != drop_var, drop = FALSE]
}

final_vif <- round(vif.cca(cca(comm_raw ~ ., data = env_work)), 3)
write.csv(vif_log,
          file.path(output_dir, "results_VIF_dropped.csv"),
          row.names = FALSE)
write.csv(data.frame(variable = names(final_vif), VIF = final_vif),
          file.path(output_dir, "results_VIF_retained.csv"),
          row.names = FALSE)

env_sel <- env_work

# Fit constrained ordination model with all retained predictors
if (use_cca) {
  cca_mod <- cca(comm_raw ~ ., data = env_sel)
} else {
  cca_mod <- rda(comm_hel ~ ., data = env_sel)
}

constr_pct <- cca_mod$CCA$tot.chi / cca_mod$tot.chi * 100
r2_adj     <- RsquareAdj(cca_mod)$adj.r.squared

# Permutation tests for overall model significance, per axis, and per predictor
set.seed(42)
anova_overall <- anova.cca(cca_mod, permutations = 999)
anova_axes    <- anova.cca(cca_mod, by = "axis",  permutations = 999)
anova_terms   <- anova.cca(cca_mod, by = "terms", permutations = 999)

# Passive fitting of environmental vectors onto ordination space
set.seed(42)
ef <- envfit(cca_mod, env_sel, perm = 9999)

cca_stats <- data.frame(
  Metric = c("Method", "n_sites", "n_species", "n_predictors",
             "Total_inertia", "Constrained_inertia_pct",
             "Adjusted_R2", "Model_F", "Model_p"),
  Value  = c(method_label, nrow(comm_raw), ncol(comm_raw), ncol(env_sel),
             round(cca_mod$tot.chi, 4), round(constr_pct, 2),
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

# Fig 02: CCA triplot — sites coloured by management intensity, sized by LCBD
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
  geom_point(data = site_sc,
             aes(x = .data[[ax1]], y = .data[[ax2]],
                 colour = Management, size = LCBD),
             alpha = 0.8) +
  geom_point(data = subset(site_sc, high_lcbd),
             aes(x = .data[[ax1]], y = .data[[ax2]]),
             shape = 1, size = 6, colour = "#D85A30", stroke = 1.2) +
  ggrepel::geom_text_repel(data = sp_sc,
                           aes(x = .data[[ax1]] * 0.8,
                               y = .data[[ax2]] * 0.8,
                               label = rownames(sp_sc)),
                           size = 2.8, colour = "#534AB7",
                           fontface = "italic", max.overlaps = 25) +
  geom_segment(data = bp_sc,
               aes(x = 0, y = 0,
                   xend = .data[[ax1]], yend = .data[[ax2]]),
               arrow  = arrow(length = unit(0.22, "cm"), type = "closed"),
               colour = "black", linewidth = 0.8) +
  ggrepel::geom_text_repel(data = bp_sc,
                           aes(x = .data[[ax1]] * 1.15,
                               y = .data[[ax2]] * 1.15,
                               label = rownames(bp_sc)),
                           size = 3.5, colour = "black",
                           fontface = "bold", max.overlaps = 30) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey70", linewidth = 0.3) +
  scale_colour_viridis_c(name = "Management\nintensity", option = "D") +
  scale_size_continuous(name = "LCBD", range = c(2, 6)) +
  labs(title    = paste0(method_label, " triplot — carabid community vs environment"),
       subtitle = paste0("Constrained variance: ", round(constr_pct, 1),
                         "%  |  adj.R² = ", round(r2_adj, 3),
                         "  |  p = ", round(anova_overall$`Pr(>F)`[1], 3)),
       x       = paste0(ax1, " (", ax1_pct, "%)"),
       y       = paste0(ax2, " (", ax2_pct, "%)"),
       caption = "Arrows = environmental predictors  |  Italic purple = species  |  Circle = high-LCBD site") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40"))
savefig(p_triplot, "fig02_CCA_triplot.png", w = 11, h = 8)

# Fig 03: NMDS ordination — exploratory unconstrained ordination
# Sites grouped by management intensity using convex hulls
set.seed(42)
nmds <- metaMDS(comm_raw, distance = "bray", k = 2,
                trymax = 100, trace = FALSE)

nmds_sc  <- as.data.frame(scores(nmds)$sites)
nmds_sp  <- as.data.frame(scores(nmds)$species)
nmds_sc$Type <- cut(management, breaks = 3,
                    labels = c("Low", "Medium", "High"))

hull_df <- nmds_sc |>
  group_by(Type) |>
  slice(chull(NMDS1, NMDS2))

set.seed(42)
ef_nmds <- envfit(nmds, env_sel, perm = 9999)
ef_vec  <- as.data.frame(scores(ef_nmds, display = "vectors"))
ef_vec$var <- rownames(ef_vec)
ef_sig  <- ef_vec[ef_nmds$vectors$pvals < 0.05, ]

p_nmds <- ggplot(nmds_sc, aes(x = NMDS1, y = NMDS2)) +
  geom_polygon(data = hull_df,
               aes(fill = Type, group = Type), alpha = 0.2) +
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
                         round(nmds$stress, 3)),
       x = "NMDS1", y = "NMDS2") +
  theme_minimal(base_size = 11) +
  theme(plot.title   = element_text(face = "bold"),
        panel.border = element_rect(colour = "black",
                                    fill = NA, linewidth = 0.8))
savefig(p_nmds, "fig03_NMDS.png", w = 10, h = 7)


# =============================================================================
# SECTION 4: PARTIAL CCA + VARIANCE PARTITIONING
# =============================================================================

# Extract spatial structure using Principal Coordinates of Neighbourhood Matrices
# Positive-eigenvalue PCNM vectors capture broad-to-fine spatial gradients
dist_mat    <- as.matrix(dist(coords[, c("X", "Y")]))
pcnm_res    <- pcnm(dist_mat)
pos_pcnm    <- which(pcnm_res$values > 0)
pcnm_scores <- as.data.frame(scores(pcnm_res)[, pos_pcnm])
colnames(pcnm_scores) <- paste0("PCNM", pos_pcnm)

# Forward selection of significant PCNM vectors
if (use_cca) {
  pcnm_null <- cca(comm_raw ~ 1, data = pcnm_scores)
  pcnm_full <- cca(comm_raw ~ ., data = pcnm_scores)
} else {
  pcnm_null <- rda(comm_hel ~ 1, data = pcnm_scores)
  pcnm_full <- rda(comm_hel ~ ., data = pcnm_scores)
}

set.seed(42)
pcnm_step     <- ordistep(pcnm_null, scope = formula(pcnm_full),
                          direction = "forward", permutations = 999,
                          trace = FALSE)
sel_pcnm_vars <- attr(terms(pcnm_step), "term.labels")
space_sel     <- pcnm_scores[, sel_pcnm_vars, drop = FALSE]

# Variance partitioning: separates pure environment [a], pure space [b],
# shared [a∩b], and unexplained residual [c]
if (use_cca) {
  vp <- varpart(comm_raw, env_sel, space_sel)
} else {
  vp <- varpart(comm_hel, env_sel, space_sel)
}

# Permutation tests for the pure environment and pure space fractions
set.seed(42)
if (use_cca) {
  pCCA_env <- cca(comm_raw ~ . + Condition(as.matrix(space_sel)),
                  data = env_sel)
  pCCA_spa <- cca(comm_raw ~ . + Condition(as.matrix(env_sel)),
                  data = space_sel)
} else {
  pCCA_env <- rda(comm_hel ~ . + Condition(as.matrix(space_sel)),
                  data = env_sel)
  pCCA_spa <- rda(comm_hel ~ . + Condition(as.matrix(env_sel)),
                  data = space_sel)
}
anova_env <- anova.cca(pCCA_env, permutations = 999)
anova_spa <- anova.cca(pCCA_spa, permutations = 999)

vp_df <- data.frame(
  Fraction = c("[env]", "[space]", "[env+space]", "[residual]"),
  R2       = round(c(vp$part$fract$R.square[1],
                     vp$part$fract$R.square[2],
                     vp$part$indfract$R.square[3],
                     vp$part$indfract$R.square[4]), 4),
  Adj_R2   = round(c(vp$part$indfract$Adj.R.square[1],
                     vp$part$indfract$Adj.R.square[2],
                     vp$part$indfract$Adj.R.square[3],
                     vp$part$indfract$Adj.R.square[4]), 4)
)
write.csv(vp_df,
          file.path(output_dir, "results_varpart.csv"),
          row.names = FALSE)

vp_plot_df <- data.frame(
  fraction = c("Pure\nenvironment [a]", "Shared\n[a∩b]",
               "Pure\nspace [b]",       "Residual\n[c]"),
  adj_r2   = c(vp$part$indfract$Adj.R.square[1],
               vp$part$indfract$Adj.R.square[3],
               vp$part$indfract$Adj.R.square[2],
               vp$part$indfract$Adj.R.square[4]),
  fill_col = c("#1D9E75", "#9FE1CB", "#534AB7", "#D3D1C7")
)

p_vp <- ggplot(vp_plot_df,
               aes(x = reorder(fraction, adj_r2), y = adj_r2,
                   fill = fill_col)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(round(adj_r2 * 100, 1), "%")),
            hjust = -0.15, size = 3.5) +
  scale_fill_identity() +
  coord_flip() +
  ylim(min(vp_plot_df$adj_r2) - 0.05,
       max(vp_plot_df$adj_r2) + 0.08) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  labs(title    = "Variance partitioning: environment vs space",
       subtitle = paste0("Pure env p = ",
                         round(anova_env$`Pr(>F)`[1], 3),
                         "  |  Pure space p = ",
                         round(anova_spa$`Pr(>F)`[1], 3)),
       x = NULL, y = "Adjusted R²") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
savefig(p_vp, "fig04_variance_partitioning.png", w = 8, h = 5)


# =============================================================================
# SECTION 5: HMSC JSDM
# =============================================================================

# Remove rare species (present in fewer than 8 sites, i.e. < 10% of sites)
# to reduce model complexity and improve MCMC convergence
comm_hmsc <- comm_raw[, colSums(comm_raw > 0) >= 8]

# Use only CCA-significant predictors (p < 0.05) as fixed effects in Hmsc
# This reduces parameter space and directly links ordination and JSDM results
cca_terms_hmsc <- read.csv(
  file.path(output_dir, "results_CCA_terms_anova.csv"),
  row.names = 1
)

# Locate the p-value column robustly across R versions
# (Pr(>F) is sanitised differently depending on the system)
p_col    <- grep("Pr", names(cca_terms_hmsc), value = TRUE)[1]
sig_vars <- rownames(cca_terms_hmsc)[
  !is.na(cca_terms_hmsc[[p_col]]) &
    cca_terms_hmsc[[p_col]] < 0.05  &
    rownames(cca_terms_hmsc) != "Residual"
]

if (length(sig_vars) == 0) {
  warning("No variables significant at p<0.05 — relaxing threshold to p<0.10")
  sig_vars <- rownames(cca_terms_hmsc)[
    !is.na(cca_terms_hmsc[[p_col]]) &
      cca_terms_hmsc[[p_col]] < 0.10  &
      rownames(cca_terms_hmsc) != "Residual"
  ]
}

# Standardise predictors to improve MCMC sampling efficiency
env_hmsc <- env_sel[, sig_vars, drop = FALSE]
XData    <- data.frame(scale(env_hmsc))
XFormula <- as.formula(paste("~", paste(sig_vars, collapse = " + ")))

Y           <- as.matrix(comm_hmsc)
cn          <- as.factor(rownames(comm_hmsc))
studyDesign <- data.frame(sample = cn)
rL          <- HmscRandomLevel(units = studyDesign$sample)

# lognormal poisson distribution is appropriate for over-dispersed count data
# (abundance values range 0–2893, mean = 35.9, 60% zeros)
m <- Hmsc(Y           = Y,
          XData       = XData,
          XFormula    = XFormula,
          studyDesign = studyDesign,
          ranLevels   = list(sample = rL),
          distr       = "lognormal poisson")

# MCMC settings: thin=10 reduces autocorrelation; effective samples = 50000/10 = 5000 per chain
nChains <- 2
nIter   <- 50000
nBurn   <- 10000
thin    <- 10

set.seed(11)
fit <- sampleMcmc(m,
                  samples   = nIter,
                  transient = nBurn,
                  thin      = thin,
                  nChains   = nChains,
                  verbose   = 1000)

saveRDS(fit, file.path(output_dir, "hmsc_model_fitted.rds"))

# Convergence diagnostics: Gelman-Rubin PSRF should be < 1.1 for all parameters
mpost      <- convertToCodaObject(fit)
psrf_omega <- gelman.diag(mpost$Omega[[1]], multivariate = FALSE)$psrf
psrf_beta  <- gelman.diag(mpost$Beta,       multivariate = FALSE)$psrf

png(file.path(output_dir, "fig05a_PSRF_omega.png"),
    width = 1500, height = 900, res = 200)
hist(psrf_omega[, 1],
     main = "MCMC convergence: Gelman-Rubin PSRF (Omega)",
     xlab = "PSRF value", col = "#9FE1CB", border = "white", breaks = 20)
abline(v = 1.1, col = "#D85A30", lwd = 2, lty = 2)
legend("topright", legend = "Threshold = 1.1",
       col = "#D85A30", lty = 2, lwd = 2, bty = "n")
dev.off()

psrf_df <- data.frame(
  parameter = c("Beta", "Omega"),
  mean_PSRF = round(c(mean(psrf_beta[, 1]),  mean(psrf_omega[, 1])), 3),
  max_PSRF  = round(c(max(psrf_beta[, 1]),   max(psrf_omega[, 1])), 3)
)
write.csv(psrf_df,
          file.path(output_dir, "results_MCMC_convergence.csv"),
          row.names = FALSE)

# Fig 05b: species-specific responses to environmental predictors (posterior means)
postBeta <- getPostEstimate(fit, "Beta")
png(file.path(output_dir, "fig05b_plotBeta.png"),
    width = 2200, height = 1400, res = 200)
par(mar = c(6, 9, 2, 2))
plotBeta(fit, postBeta, param = "Mean",
         main = "Species responses to environmental predictors")
dev.off()

# Fig 05c: residual species co-occurrence matrix (Omega)
# Values represent pairwise associations after controlling for environmental effects
OmegaCor <- computeAssociations(fit)
toPlot   <- OmegaCor[[1]]$mean

png(file.path(output_dir, "fig05c_omega_cooccurrence.png"),
    width = 2000, height = 1900, res = 200)
corrplot(toPlot,
         method  = "color",
         col     = colorRampPalette(c("#185FA5", "white", "#D85A30"))(200),
         title   = "Residual species co-occurrence matrix (Omega)",
         type    = "lower", tl.col = "black", tl.cex = 0.65,
         mar     = c(0, 0, 3, 0))
dev.off()

write.csv(as.data.frame(toPlot),
          file.path(output_dir, "results_Omega_matrix.csv"))

# Fig 05d–05e: gradient predictions along the strongest environmental predictor
# sig_vars[1] = Org (soil organic matter), highest F-value in CCA (F=8.67, p=0.001)
focal_var <- sig_vars[1]
Gradient  <- constructGradient(fit, focalVariable = focal_var)
predY     <- predict(fit, Gradient = Gradient, expected = TRUE)

png(file.path(output_dir, "fig05d_gradient_richness.png"),
    width = 1600, height = 1000, res = 200)
plotGradient(fit, Gradient, pred = predY, measure = "S", index = 1,
             showData = TRUE,
             main = paste0("Predicted species richness along ", focal_var))
dev.off()

top_sp <- names(sort(colSums(comm_hmsc), decreasing = TRUE))[1]
sp_idx <- which(colnames(comm_hmsc) == top_sp)

png(file.path(output_dir, "fig05e_gradient_topspecies.png"),
    width = 1600, height = 1000, res = 200)
plotGradient(fit, Gradient, pred = predY, measure = "Y", index = sp_idx,
             showData = TRUE,
             main = paste0("Predicted abundance: ", top_sp,
                           " along ", focal_var))
dev.off()


# =============================================================================
# SECTION 6: MASTER SUMMARY
# =============================================================================

master <- data.frame(
  Metric = c(
    "n_sites", "n_species_CCA", "n_species_hmsc", "matrix_sparsity_pct",
    "alpha_richness_q0", "alpha_shannon_q1",
    "beta_q0", "beta_q1", "gamma_q0",
    "n_high_LCBD_sites",
    "DCA_axis1_SD", "ordination_method",
    "n_vars_after_VIF", "retained_variables",
    "constrained_variance_pct", "adjusted_R2",
    "model_F", "model_p",
    "NMDS_stress",
    "varpart_env_adjR2", "varpart_space_adjR2",
    "varpart_shared_adjR2", "varpart_residual_adjR2",
    "pure_env_p", "pure_space_p",
    "hmsc_sig_vars", "hmsc_n_fixed_effects",
    "MCMC_chains", "MCMC_samples", "MCMC_thin",
    "Beta_mean_PSRF", "Omega_mean_PSRF"
  ),
  Value = c(
    nrow(comm_raw), ncol(comm_raw), ncol(comm_hmsc),
    round(sum(comm_raw == 0) / prod(dim(comm_raw)) * 100, 1),
    round(alpha_q0, 3), round(alpha_q1, 3),
    round(beta_q0, 3), round(beta_q1, 3), round(gamma_q0, 3),
    sum(high_lcbd),
    round(dca_ax1, 3), method_label,
    ncol(env_sel), paste(names(env_sel), collapse = "; "),
    round(constr_pct, 2), round(r2_adj, 4),
    round(anova_overall$F[1], 3),
    round(anova_overall$`Pr(>F)`[1], 4),
    round(nmds$stress, 4),
    round(vp$part$indfract$Adj.R.square[1], 4),
    round(vp$part$indfract$Adj.R.square[2], 4),
    round(vp$part$indfract$Adj.R.square[3], 4),
    round(vp$part$indfract$Adj.R.square[4], 4),
    round(anova_env$`Pr(>F)`[1], 4),
    round(anova_spa$`Pr(>F)`[1], 4),
    paste(sig_vars, collapse = "; "), length(sig_vars),
    nChains, nIter, thin,
    round(mean(psrf_beta[, 1]),  3),
    round(mean(psrf_omega[, 1]), 3)
  )
)

write.csv(master,
          file.path(output_dir, "MASTER_results_summary.csv"),
          row.names = FALSE)

message("Analysis complete. All outputs saved to: ", output_dir)