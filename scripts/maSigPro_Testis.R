# =============================================================================
# maSigPro Time-Course Analysis — Testis
# =============================================================================
# Author:  David L. Hubert, Oregon State University (hubertd@oregonstate.edu)
# Project: Transcriptomic dynamics of brumation in Thamnophis sirtalis parietalis
# Journal: Journal of Experimental Biology https://doi.org/10.1111%2F1365-2435.70357
#
# Description:
#   Time-course differential expression analysis of testis RNA-seq data across
#   five brumation timepoints using the maSigPro package. Significant genes are
#   clustered into k=18 groups via hierarchical clustering. 
#
# Input files (place in working directory):
#   raw_countsT.tab  — raw counts matrix (transcripts x samples), tab-delimited,
#                      first column is transcript ID
#   Tkey.tab         — maSigPro edesign key file with columns:
#                      Time, Replicate, <group columns as required by maSigPro>
#
# Output files:
#   checkpoints/fit_testis.rds    — saved p.vector() output (re-used on re-runs)
#   checkpoints/tstep_testis.rds  — saved T.fit() output  (re-used on re-runs)
#   k18/k18_cluster1.tab ...      — gene lists per cluster
#   silhouette_testis.png         — cluster number selection plot
#
# Packages required:
#   data.table, NOISeq, maSigPro, cluster
#
# Usage:
#   1. Set WORKING_DIR below to the folder containing your input files.
#   2. Run the script in full. The two slow steps (p.vector, T.fit) are
#      checkpointed: subsequent runs load saved results and skip them.
# =============================================================================

# ── 0. Configuration — edit these lines ──────────────────────────────────────

WORKING_DIR   <- ""   # <-- set this to your data directory
K             <- 18                     # number of clusters for see.genes()
COVERAGE_THRD <- 5                      # min average reads per sample (filter)
MASIGPRO_Q    <- 0.05                   # FDR threshold for p.vector()
MASIGPRO_RSQ  <- 0.5                    # R² threshold for get.siggenes()
MIN_OBS       <- 7                      # min observations per gene for p.vector()
DEGREE        <- 4                      # polynomial degree (= n timepoints - 1)
K_EVAL_RANGE  <- 8:28                   # k values to evaluate for silhouette plot

# ── 1. Setup ──────────────────────────────────────────────────────────────────

setwd(WORKING_DIR)

suppressPackageStartupMessages({
  library(data.table)
  library(NOISeq)
  library(maSigPro)
  library(cluster)
})

dir.create("checkpoints", showWarnings = FALSE)
dir.create("k18",         showWarnings = FALSE)

# ── 2. Load raw counts ────────────────────────────────────────────────────────

cat("Loading raw counts...\n")
rawcounts_dt      <- fread("raw_countsT.tab", data.table = FALSE)
rownames(rawcounts_dt) <- rawcounts_dt[, 1]
rawcounts         <- rawcounts_dt[, -1, drop = FALSE]

cat("  Raw transcripts:", nrow(rawcounts), "\n")

# ── 3. Filter low-expression transcripts ──────────────────────────────────────

fltrcounts <- rawcounts[rowSums(rawcounts) >= (ncol(rawcounts) * COVERAGE_THRD), ]
cat("  After coverage filter (min avg", COVERAGE_THRD, "reads):", nrow(fltrcounts), "\n")

# ── 4. TMM normalisation ──────────────────────────────────────────────────────

cat("Normalising (TMM)...\n")
data <- tmm(fltrcounts, long = 1000, lc = 0)
cat("  Transcripts after normalisation:", nrow(data), "\n")

# ── 5. Design matrix ──────────────────────────────────────────────────────────

edesign <- read.table("Tkey.tab", header = TRUE)
design  <- make.design.matrix(edesign = edesign, degree = DEGREE)

# ── 6. p.vector — NB regression (slow; checkpointed) ─────────────────────────

if (file.exists("checkpoints/fit_testis.rds")) {
  cat("Loading saved p.vector output...\n")
  fit <- readRDS("checkpoints/fit_testis.rds")
} else {
  cat("Running p.vector() — this is the slowest step...\n")
  fit <- p.vector(data    = data,
                  design  = design,
                  Q       = MASIGPRO_Q,
                  counts  = TRUE,
                  min.obs = MIN_OBS)
  saveRDS(fit, "checkpoints/fit_testis.rds")
}
cat("  Significant transcripts:", fit$i, "\n")

# ── 7. T.fit — stepwise regression (slow; checkpointed) ──────────────────────

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

# ── 9. Silhouette analysis — cluster number selection ─────────────────────────
# Evaluates k over K_EVAL_RANGE using mean silhouette width.
# Mirrors the internal clustering performed by see.genes() with cluster.data=1:
#   z-score scaling per gene → Euclidean distance → Ward D2 hierarchical clustering.

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
d1        <- diff(sil_scores)
d2        <- diff(d1)
elbow_k   <- K_EVAL_RANGE[which.max(d2) + 1]

cat("  Silhouette peak at k =", K_EVAL_RANGE[which.max(sil_scores)],
    "(score =", round(max(sil_scores), 4), ")\n")
cat("  Elbow (second derivative) at k =", elbow_k, "\n")
cat("  Silhouette at selected k =", K, ":",
    round(sil_scores[K_EVAL_RANGE == K], 4), "\n")

# Save silhouette plot
png("silhouette_testis.png", width = 1800, height = 1400, res = 200)
par(mfrow = c(2, 1), mar = c(4, 4, 2.5, 1))

plot(K_EVAL_RANGE, sil_scores, type = "b", pch = 19,
     xlab = "Number of clusters (k)", ylab = "Mean silhouette width",
     main = "Cluster number selection — Testis maSigPro analysis")
abline(v = elbow_k, lty = 2, col = "blue")
abline(v = K,       lty = 2, col = "red")
legend("topright",
       legend = c(paste0("Elbow (k = ", elbow_k, ")"),
                  paste0("Selected (k = ", K, ")")),
       col = c("blue", "red"), lty = 2, bty = "n")

plot(K_EVAL_RANGE[2:(length(K_EVAL_RANGE) - 1)], d2,
     type = "b", pch = 19,
     xlab = "Number of clusters (k)", ylab = "Second derivative",
     main = "Rate of change in silhouette decline (elbow detection)")
abline(h = 0, lty = 3, col = "grey50")
abline(v = elbow_k, lty = 2, col = "blue")
points(elbow_k, d2[which.max(d2)], pch = 19, col = "blue", cex = 1.5)

par(mfrow = c(1, 1))
dev.off()
cat("  Silhouette plot saved to silhouette_testis.png\n")

# ── 10. Cluster & visualise ───────────────────────────────────────────────────

cat("\nRunning see.genes() with k =", K, "...\n")
k_result <- see.genes(sigs.all$sig.genes,
                      min.obs        = 10,
                      show.fit       = TRUE,
                      dis            = design$dis,
                      cluster.method = "hclust",
                      cluster.data   = 1,
                      k              = K)

# ── 11. Write cluster gene lists ──────────────────────────────────────────────

cat("Writing cluster gene lists to k18/...\n")
for (i in seq_len(K)) {
  write.table(
    names(which(k_result$cut == i)),
    file      = paste0("k18/k18_cluster", i, ".tab"),
    sep       = "\t",
    col.names = NA
  )
}

cat("\nDone. Cluster files written to k18/\n")
cat("Cluster sizes:\n")
print(table(k_result$cut))
