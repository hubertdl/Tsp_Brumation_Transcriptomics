# =============================================================================
# maSigPro Time-Course Analysis - Testis
# =============================================================================
# Author:  David L. Hubert, Oregon State University (hubertd@oregonstate.edu)
# Project: Transcriptomic dynamics of brumation in Thamnophis sirtalis parietalis
# Journal: Journal of Experimental Biology
#
# Description:
#   Time-course differential expression analysis of testis RNA-seq data across
#   five brumation timepoints (weeks 0, 8, 16, 29, 31) using the maSigPro
#   package. Significant genes are clustered via hierarchical clustering, with
#   the number of clusters (k) selected at the elbow of the mean silhouette
#   width curve (k = 10 for testis).
#
# Input files (place in working directory):
#   raw_countsT.tab  - raw counts matrix (transcripts x samples), tab-delimited,
#                      first column is transcript ID
#   Tkey.tab         - maSigPro edesign key file with columns:
#                      Time, Replicate, <group columns as required by maSigPro>
#
# Output files:
#   checkpoints/fit_testis.rds     - saved p.vector() output (re-used on re-runs)
#   checkpoints/tstep_testis.rds   - saved T.fit() output  (re-used on re-runs)
#   silhouette_testis.png          - cluster number selection plot
#   k10/k10_cluster1.tab ...      - gene lists per cluster
#   Figure_3_testis_clusters_k10.tif / .png - cluster expression profiles (Fig. 3)
#   DEGOE_inputs/Testis_k10_Scores/ - binary cluster-membership files for DEGOE
#
# Packages required:
#   data.table, NOISeq, maSigPro, cluster
#
# Usage:
#   1. Set WORKING_DIR below to the folder containing your input files.
#   2. Run the script in full. The two slow steps (p.vector, T.fit) are
#      checkpointed: subsequent runs load saved results and skip them.
# =============================================================================

# ── 0. Configuration: edit these lines ───────────────────────────────────────

WORKING_DIR   <- ""   # <-- set this to your data directory
K_OVERRIDE    <- NA                     # NA = use the silhouette elbow (k = 10)
COVERAGE_THRD <- 5                      # min average reads per sample (filter)
MASIGPRO_Q    <- 0.05                   # FDR threshold for p.vector()
MASIGPRO_RSQ  <- 0.5                    # R² threshold for get.siggenes()
MIN_OBS       <- 7                      # min observations per gene for p.vector()
DEGREE        <- 4                      # polynomial degree (= n timepoints - 1)
K_EVAL_RANGE  <- 8:28                   # k values to evaluate for silhouette plot
TISSUE_COL    <- "#0072B2"              # testis colour (Okabe-Ito, colorblind safe)

# ── 1. Setup ──────────────────────────────────────────────────────────────────

setwd(WORKING_DIR)

suppressPackageStartupMessages({
  library(data.table)
  library(NOISeq)
  library(maSigPro)
  library(cluster)
})

dir.create("checkpoints", showWarnings = FALSE)

# ── 2. Load raw counts ────────────────────────────────────────────────────────

cat("Loading raw counts...\n")
rawcounts_dt           <- fread("raw_countsT.tab", data.table = FALSE)
rownames(rawcounts_dt) <- rawcounts_dt[, 1]
rawcounts              <- rawcounts_dt[, -1, drop = FALSE]

cat("  Raw transcripts:", nrow(rawcounts), "\n")

# ── 3. Filter low-expression transcripts ──────────────────────────────────────

fltrcounts <- rawcounts[rowSums(rawcounts) >= (ncol(rawcounts) * COVERAGE_THRD), ]
cat("  After coverage filter (min avg", COVERAGE_THRD, "reads):", nrow(fltrcounts), "\n")

# ── 4. TMM normalisation ──────────────────────────────────────────────────────

cat("Normalising (TMM)...\n")
data <- tmm(fltrcounts, long = 1000, lc = 0)
cat("  Transcripts after normalisation:", nrow(data), "\n")

# ── 5. Design matrix ──────────────────────────────────────────────────────────
# Timepoints are weeks from the start of brumation: 0, 8, 16, 29, 31

edesign <- read.table("Tkey.tab", header = TRUE)
design  <- make.design.matrix(edesign = edesign, degree = DEGREE)

# ── 6. p.vector: NB regression (slow; checkpointed) ──────────────────────────

if (file.exists("checkpoints/fit_testis.rds")) {
  cat("Loading saved p.vector output...\n")
  fit <- readRDS("checkpoints/fit_testis.rds")
} else {
  cat("Running p.vector(), this is the slowest step...\n")
  fit <- p.vector(data    = data,
                  design  = design,
                  Q       = MASIGPRO_Q,
                  counts  = TRUE,
                  min.obs = MIN_OBS)
  saveRDS(fit, "checkpoints/fit_testis.rds")
}
cat("  Significant transcripts:", fit$i, "\n")

# ── 7. T.fit: stepwise regression (slow; checkpointed) ───────────────────────

if (file.exists("checkpoints/tstep_testis.rds")) {
  cat("Loading saved T.fit output...\n")
  tstep <- readRDS("checkpoints/tstep_testis.rds")
} else {
  cat("Running T.fit()...\n")
  tstep <- T.fit(fit)
  saveRDS(tstep, "checkpoints/tstep_testis.rds")
}

# ── 8. Extract significant genes ──────────────────────────────────────────────

sigs.all   <- get.siggenes(tstep = tstep, rsq = MASIGPRO_RSQ, vars = "all")
sigs.group <- get.siggenes(tstep = tstep, rsq = MASIGPRO_RSQ, vars = "groups")
sigs.each  <- get.siggenes(tstep = tstep, rsq = MASIGPRO_RSQ, vars = "each")

cat("\n--- Significant gene summaries ---\n")
print(sigs.all$summary)
cat("  Transcripts with R² >=", MASIGPRO_RSQ, ":",
    nrow(sigs.all$sig.genes$sig.profiles), "\n")   # 5,822

# ── 9. Silhouette analysis: cluster number selection ─────────────────────────
# Evaluates k over K_EVAL_RANGE using mean silhouette width.
# Mirrors the internal clustering performed by see.genes() with cluster.data=1:
#   z-score scaling per gene → Euclidean distance → Ward D2 hierarchical clustering.
# k is selected at the initial inflection (elbow) of the silhouette curve,
# identified as the maximum of its second derivative.

cat("\nRunning silhouette analysis across k =",
    min(K_EVAL_RANGE), "to", max(K_EVAL_RANGE), "...\n")

expr_mat    <- sigs.all$sig.genes$sig.profiles
expr_scaled <- t(scale(t(expr_mat)))
expr_scaled <- expr_scaled[complete.cases(expr_scaled), ]

d  <- dist(expr_scaled, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

sil_scores <- sapply(K_EVAL_RANGE, function(k) {
  cuts <- cutree(hc, k = k)
  mean(silhouette(cuts, d)[, "sil_width"])
})

# Elbow detection via second derivative
d1      <- diff(sil_scores)
d2      <- diff(d1)
elbow_k <- K_EVAL_RANGE[which.max(d2) + 1]
K       <- if (is.na(K_OVERRIDE)) elbow_k else K_OVERRIDE

cat("  Elbow (second derivative) at k =", elbow_k, "\n")
cat("  Selected k =", K, "\n")

# Save silhouette plot
png("silhouette_testis.png", width = 1800, height = 1400, res = 200)
par(mfrow = c(2, 1), mar = c(4, 4, 2.5, 1))

plot(K_EVAL_RANGE, sil_scores, type = "b", pch = 19,
     xlab = "Number of clusters (k)", ylab = "Mean silhouette width",
     main = "Cluster number selection: Testis maSigPro analysis")
abline(v = K, lty = 2, col = "#0072B2")
legend("topright",
       legend = paste0("Selected (k = ", K, ")"),
       col = "#0072B2", lty = 2, bty = "n")

plot(K_EVAL_RANGE[2:(length(K_EVAL_RANGE) - 1)], d2,
     type = "b", pch = 19,
     xlab = "Number of clusters (k)", ylab = "Second derivative",
     main = "Rate of change in silhouette decline (elbow detection)")
abline(h = 0, lty = 3, col = "grey50")
abline(v = K, lty = 2, col = "#0072B2")
points(K, d2[K - min(K_EVAL_RANGE)], pch = 19, col = "#0072B2", cex = 1.5)

par(mfrow = c(1, 1))
dev.off()
cat("  Silhouette plot saved to silhouette_testis.png\n")

# ── 10. Cluster & visualise ───────────────────────────────────────────────────

cat("\nRunning see.genes() with k =", K, "...\n")
pdf(NULL)   # see.genes() plots are discarded; the figure is drawn below
k_result <- see.genes(sigs.all$sig.genes,
                      min.obs        = 10,
                      show.fit       = FALSE,
                      dis            = design$dis,
                      cluster.method = "hclust",
                      cluster.data   = 1,
                      k              = K)
dev.off()
cut <- k_result$cut

# x-axis labels: PlotGroups() centres every week label under its tick, so 29
# and 31 run together ("2931"). TIME_LABELS() is called right after each
# PlotGroups() plot to blank those labels and redraw them, splitting any pair
# that would crowd each other. Ticks, data and legend are not touched.
TIME_LABELS <- function(weeks, gap = 0.8, cex = 1) {
  weeks <- sort(unique(weeks))
  labs  <- as.character(weeks)
  usr   <- par("usr")
  m     <- par("mgp")[2]
  lh    <- par("mai")[1] / par("mar")[1]
  y_at  <- function(l) usr[3] - l * lh * diff(usr[3:4]) / par("pin")[2]
  op <- par(xpd = NA); on.exit(par(op))
  rect(usr[1] - 0.06 * diff(usr[1:2]), y_at(m + 1.25), usr[2], y_at(m - 0.05),
       col = "white", border = NA)
  w   <- strwidth(labs, cex = cex)
  g   <- gap * strwidth("0", cex = cex)
  at  <- weeks
  adj <- rep(0.5, length(weeks))
  for (i in seq_along(weeks)[-1]) {
    if (weeks[i] - weeks[i - 1] < (w[i] + w[i - 1]) / 2 + g) {
      mid <- (weeks[i] + weeks[i - 1]) / 2
      at[i - 1] <- mid - g / 2; adj[i - 1] <- 1
      at[i]     <- mid + g / 2; adj[i]     <- 0
    }
  }
  for (i in seq_along(weeks))
    mtext(labs[i], side = 1, at = at[i], adj = adj[i], line = m, cex = cex * par("cex"))
}

# Figure 3: cluster panels drawn as see.genes() draws them (median profiles),
# in one 3-column grid, saved as TIFF and PNG at 300 dpi
time    <- edesign[, 1]
repvect <- edesign[, 2]
groups  <- as.matrix(edesign[, c(3:ncol(edesign))])
colnames(groups) <- colnames(edesign)[3:ncol(edesign)]
prof <- as.data.frame(sigs.all$sig.genes$sig.profiles)
prof <- as.data.frame(prof[, (ncol(prof) - length(time) + 1):ncol(prof)])
prof <- prof[names(cut), ]

nrow_grid <- ceiling(K / 3)
palette(c("black", TISSUE_COL, palette.colors(8, "Okabe-Ito")[-1]))
for (fmt in c("tif", "png")) {
  f <- paste0("Figure_3_testis_clusters_k", K, ".", fmt)
  if (fmt == "tif") tiff(f, width = 250 / 25.4, height = 190 / 3 * nrow_grid / 25.4,
                         units = "in", res = 300, compression = "lzw")
  else              png(f,  width = 250 / 25.4, height = 190 / 3 * nrow_grid / 25.4,
                        units = "in", res = 300)
  par(mfrow = c(nrow_grid, 3), xpd = FALSE)
  for (j in seq_len(K)) {
    PlotGroups(data = prof[cut == j, ], show.fit = TRUE, dis = design$dis,
               step.method = "backward", min.obs = 10, alfa = 0.05,
               nvar.correction = FALSE, show.lines = TRUE, time = time,
               groups = groups, repvect = repvect, summary.mode = "median",
               xlab = "Time", main = paste("Cluster", j), ylim = NULL,
               legend = TRUE, groups.vector = design$groups.vector,
               item = "genes", cex.legend = 1, lty.legend = NULL)
    TIME_LABELS(time)
  }
  dev.off()
}
palette("default")
cat("  Cluster figure saved to Figure_3_testis_clusters_k", K, ".tif/.png\n", sep = "")

# ── 11. Write cluster gene lists ──────────────────────────────────────────────

out_dir <- paste0("k", K)
dir.create(out_dir, showWarnings = FALSE)
cat("Writing cluster gene lists to ", out_dir, "/...\n", sep = "")
for (i in seq_len(K)) {
  write.table(
    names(which(cut == i)),
    file      = paste0(out_dir, "/", out_dir, "_cluster", i, ".tab"),
    sep       = "\t",
    col.names = NA
  )
}

# ── 12. Write DEGOE input files ───────────────────────────────────────────────
# One file per cluster listing every clustered transcript: value = 1 if in the
# cluster, 0 if in another cluster (background for the Fisher's exact test).

score_dir <- file.path("DEGOE_inputs", paste0("Testis_k", K, "_Scores"))
dir.create(score_dir, showWarnings = FALSE, recursive = TRUE)
universe <- sort(names(cut))
for (i in seq_len(K)) {
  write.table(data.frame(gene = universe, value = as.integer(cut[universe] == i)),
              file.path(score_dir, paste0("testis_k", K, "_cluster", i, "_scores.tab")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("\nDone. Cluster files written to ", out_dir, "/\n", sep = "")
cat("Cluster sizes:\n")
print(table(cut))
