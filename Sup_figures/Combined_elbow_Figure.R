library(maSigPro)
library(cluster)

#clear environment
rm(list=ls())+
  dev.off()

# ── Paths ─────────────────────────────────────────────────────────────────────
base <- ""

# add paths to the output from checkpoints directories made after running MaSigPro
paths <- list(
  Combined = file.path(base, ""),
  Liver    = file.path(base, ""),
  Testis   = file.path(base, "")
)

K_RANGE <- 8:28
K_SEL   <- 18

# ── Compute silhouette scores for each analysis ───────────────────────────────
compute_sil <- function(rds_path) {
  tstep    <- readRDS(rds_path)
  sigs.all <- get.siggenes(tstep = tstep, rsq = 0.5, vars = "all")
  expr     <- t(scale(t(sigs.all$sig.genes$sig.profiles)))
  expr     <- expr[complete.cases(expr), ]
  d        <- dist(expr, method = "euclidean")
  hc       <- hclust(d, method = "ward.D2")
  scores   <- sapply(K_RANGE, function(k) {
    mean(silhouette(cutree(hc, k), d)[, "sil_width"])
  })
  d1 <- diff(scores)
  d2 <- diff(d1)
  list(scores = scores, d2 = d2, elbow_k = K_RANGE[which.max(d2) + 1])
}

cat("Computing Combined...\n"); res_combined <- compute_sil(paths$Combined)
cat("Computing Liver...\n");    res_liver    <- compute_sil(paths$Liver)
cat("Computing Testis...\n");   res_testis   <- compute_sil(paths$Testis)

results <- list(Combined = res_combined, Liver = res_liver, Testis = res_testis)

# ── Plot ──────────────────────────────────────────────────────────────────────
out_path <- file.path(base, "silhouette_supplemental_figure.png")
png(out_path, width = 5400, height = 3200, res = 300)

par(mfrow = c(2, 3),
    mar   = c(4.5, 4.5, 3, 1),
    oma   = c(0, 0, 2, 0))

# Shared y-axis limits for cleaner comparison across panels
sil_ylim <- range(sapply(results, function(r) r$scores))
sil_ylim <- sil_ylim + c(-0.005, 0.005)
d2_ylim  <- range(sapply(results, function(r) r$d2))
d2_ylim  <- d2_ylim + c(-0.002, 0.002)

panel_labels <- c("A", "B", "C", "D", "E", "F")
label_idx    <- 1

for (name in c("Combined", "Liver", "Testis")) {
  r <- results[[name]]
  
  # Top row: silhouette curve
  plot(K_RANGE, r$scores, type = "b", pch = 19, cex = 0.8,
       ylim  = sil_ylim,
       xlab  = "Number of clusters (k)",
       ylab  = "Mean silhouette width",
       main  = paste(name, "analysis"))
  abline(v = r$elbow_k, lty = 2, col = "blue")
  abline(v = K_SEL,     lty = 2, col = "red")
  legend("topright",
         legend = c(paste0("Elbow (k = ", r$elbow_k, ")"),
                    paste0("Selected (k = ", K_SEL, ")")),
         col = c("blue", "red"), lty = 2, bty = "n", cex = 0.85)
  mtext(panel_labels[label_idx], side = 3, adj = 0, font = 2, cex = 1.1)
  label_idx <- label_idx + 1
}

for (name in c("Combined", "Liver", "Testis")) {
  r    <- results[[name]]
  k_d2 <- K_RANGE[2:(length(K_RANGE) - 1)]
  
  plot(k_d2, r$d2, type = "b", pch = 19, cex = 0.8,
       ylim  = d2_ylim,
       xlab  = "Number of clusters (k)",
       ylab  = "Second derivative",
       main  = "Elbow detection")
  abline(h = 0,         lty = 3, col = "grey50")
  abline(v = r$elbow_k, lty = 2, col = "blue")
  abline(v = K_SEL,     lty = 2, col = "red") 
  points(r$elbow_k, r$d2[which.max(r$d2)], pch = 19, col = "blue", cex = 1.5)
  mtext(panel_labels[label_idx], side = 3, adj = 0, font = 2, cex = 1.1)
  label_idx <- label_idx + 1
}

mtext("Cluster number selection — maSigPro hierarchical clustering",
      outer = TRUE, cex = 1.1, font = 2)

dev.off()
cat("Saved to:", out_path, "\n")

