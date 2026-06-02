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