# =============================================================================
# figure_generation.R — Brumation Manuscript Figure Generation
# =============================================================================
# Author:  David L. Hubert, Oregon State University (hubertd@oregonstate.edu)
# Project: Transcriptomic dynamics of brumation in Thamnophis sirtalis parietalis
# Journal: Journal of Experimental Biology https://doi.org/10.1111%2F1365-2435.70357
#
# Purpose:
#   Generate all manuscript figures at 300 dpi TIFF.
#   Uses PlotGroups() from maSigPro for individual gene expression panels
#   (Figures 4–7), and a custom cluster-mean visualization for Figure 2.
#
#   Run AFTER the main maSigPro analysis script has completed and checkpoint
#   files exist in checkpoints/.
#
# Input files (same directory as main analysis script):
#   raw_counts.tab                  — raw count matrix
#   key.tab                         — maSigPro edesign key
#   checkpoints/fit_combined.rds    — saved p.vector() output
#   checkpoints/tstep_combined.rds  — saved T.fit() output
#
# Output (figures_300dpi/ subdirectory):
#   Figure_2A_clusters1-9.tif
#   Figure_2B_clusters10-18.tif
#   Figure_4A_lipid_mobilization.tif
#   Figure_4B_lipid_transport.tif
#   Figure_4C_lipid_biosynthesis.tif
#   Figure_5Ai_glycogen.tif
#   Figure_5Aii_nitrogen_carbon.tif
#   Figure_5C_gluconeogenesis.tif
#   Figure_6A_HSPs.tif
#   Figure_6B_hypoxia.tif
#   Figure_6D_oxidative_stress.tif
#   Figure_7A_vitellogenin.tif
#   Figure_7B_HSDs.tif
#
#   Then run figure_assembly.py to combine sub-panels into final figures with
#   bold panel labels (A, B, C ...) and insert schematic panels for Figures 5 and 6.
# =============================================================================

# ── 0. Configuration — edit WORKING_DIR if needed ────────────────────────────
WORKING_DIR   <- "" # <- set working directory
OUTPUT_DIR    <- "figures_300dpi"
DPI           <- 300
K             <- 18
DEGREE        <- 4
COVERAGE_THRD <- 5

# JEB figure widths (mm)
W_FULL <- 174    # full page width
W_HALF <- 87     # single column

# ── 1. Setup ──────────────────────────────────────────────────────────────────
setwd(WORKING_DIR)
suppressPackageStartupMessages({
  library(data.table)
  library(NOISeq)
  library(maSigPro)
  library(cluster)
})
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ── 2. Load data ──────────────────────────────────────────────────────────────
cat("Loading raw counts...\n")
rawcounts_dt           <- fread("raw_counts.tab", data.table = FALSE)
rownames(rawcounts_dt) <- rawcounts_dt[, 1]
rawcounts              <- rawcounts_dt[, -1, drop = FALSE]
fltrcounts  <- rawcounts[rowSums(rawcounts) >= (ncol(rawcounts) * COVERAGE_THRD), ]
normcounts  <- tmm(fltrcounts, long = 1000, lc = 0)   # 'normcounts' = TMM-normalized data

edesign <- read.table("key.tab", header = TRUE)
design  <- make.design.matrix(edesign = edesign, degree = DEGREE)

cat("Loading maSigPro checkpoints...\n")
tstep    <- readRDS("checkpoints/tstep_combined.rds")
sigs.all <- get.siggenes(tstep = tstep, rsq = 0.5, vars = "all")

cat("Data loaded.\n")

# ── 3. Figure 2 — Cluster expression profiles ─────────────────────────────────
# Custom visualization: z-score profiles per gene (thin traces) +
# per-tissue group means (thick colored lines) for each of the k=18 clusters.
# Splits into two 3x3 panels (clusters 1-9 and 10-18).
cat("\n--- Figure 2: Cluster expression profiles ---\n")

# Compute cluster assignments (identical method to main analysis script)
expr_scaled <- t(scale(t(sigs.all$sig.genes$sig.profiles)))
expr_scaled <- expr_scaled[complete.cases(expr_scaled), ]
hc18 <- hclust(dist(expr_scaled, method = "euclidean"), method = "ward.D2")
k18  <- cutree(hc18, k = K)

time   <- edesign$Time
tissue <- ifelse(edesign$Liver == 1, "Liver", "Testis")
tps    <- sort(unique(time))

# Helper: plot one cluster panel
plot_cluster <- function(cl, show_legend = FALSE) {
  genes <- names(which(k18 == cl))
  genes <- genes[genes %in% rownames(sigs.all$sig.genes$sig.profiles)]
  if (length(genes) == 0) { plot.new(); return(invisible(NULL)) }
  
  mat   <- sigs.all$sig.genes$sig.profiles[genes, , drop = FALSE]
  mat_z <- t(scale(t(mat)))
  mat_z <- mat_z[complete.cases(mat_z), , drop = FALSE]
  
  lmean <- sapply(tps, function(t) mean(mat_z[, time == t & tissue == "Liver"],  na.rm = TRUE))
  tmean <- sapply(tps, function(t) mean(mat_z[, time == t & tissue == "Testis"], na.rm = TRUE))
  
  yr <- range(mat_z, na.rm = TRUE)
  yr <- yr + diff(yr) * c(-0.05, 0.05)   # small breathing room
  
  plot(NA, xlim = c(0.6, length(tps) + 0.4), ylim = yr,
       xaxt = "n", xlab = "", ylab = "z-score",
       main = paste0("Cluster ", cl, "  (n = ", nrow(mat_z), ")"),
       cex.main = 0.85, cex.axis = 0.75, cex.lab = 0.8)
  axis(1, at = seq_along(tps),
       labels = c("Pre", "Early", "Mid", "Late", "Post"), cex.axis = 0.75)
  abline(h = 0, lty = 3, col = "grey70", lwd = 0.7)
  
  # Individual gene traces (semi-transparent)
  for (g in seq_len(nrow(mat_z))) {
    gl <- sapply(tps, function(t) mean(mat_z[g, time == t & tissue == "Liver"],  na.rm = TRUE))
    gt <- sapply(tps, function(t) mean(mat_z[g, time == t & tissue == "Testis"], na.rm = TRUE))
    lines(seq_along(tps), gl, col = adjustcolor("#2166AC", 0.18), lwd = 0.4)
    lines(seq_along(tps), gt, col = adjustcolor("#D6604D", 0.18), lwd = 0.4)
  }
  
  # Group mean lines (thick)
  lines(seq_along(tps), lmean, col = "#2166AC", lwd = 2.2, type = "b", pch = 16, cex = 0.85)
  lines(seq_along(tps), tmean, col = "#D6604D", lwd = 2.2, type = "b", pch = 17, cex = 0.85)
  
  if (show_legend) {
    legend("topright", legend = c("Liver", "Testis"),
           col = c("#2166AC", "#D6604D"), lwd = 2, pch = c(16, 17),
           bty = "n", cex = 0.8)
  }
}

# Figure 2A: Clusters 1–9
cat("  Saving Figure_2A_clusters1-9.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_2A_clusters1-9.tif"),
     width = W_FULL, height = 190, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(3, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
for (cl in 1:9) plot_cluster(cl, show_legend = (cl == 1))
mtext("z-score", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

# Figure 2B: Clusters 10–18
cat("  Saving Figure_2B_clusters10-18.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_2B_clusters10-18.tif"),
     width = W_FULL, height = 190, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(3, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
for (cl in 10:18) plot_cluster(cl, show_legend = (cl == 10))
mtext("z-score", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

cat("  Figure 2 panels saved.\n")

# ── 4. Figure 4 — Lipid metabolism program ────────────────────────────────────
cat("\n--- Figure 4: Lipid metabolism ---\n")

# Panel A: Lipid mobilization and beta-oxidation (6 genes, 2x3)
cat("  Saving Figure_4A_lipid_mobilization.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_4A_lipid_mobilization.tif"),
     width = W_FULL, height = 130, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(2, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN30348_c1_g2_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           ylim = c(0, 800), main = "ATGL")
PlotGroups(normcounts[rownames(normcounts) == "DN27598_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           ylim = c(0, 60),  main = "FOXO1")
PlotGroups(normcounts[rownames(normcounts) == "DN31376_c0_g1_i2", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "PPARa")
PlotGroups(normcounts[rownames(normcounts) == "DN27065_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ACSL4")
PlotGroups(normcounts[rownames(normcounts) == "DN25931_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           ylim = c(0, 600), main = "CPT1A")
PlotGroups(normcounts[rownames(normcounts) == "DN3769_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "CPT2")
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

# Panel B: Lipid transport and acyl-CoA metabolism (4 genes, 2x2)
cat("  Saving Figure_4B_lipid_transport.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_4B_lipid_transport.tif"),
     width = W_FULL, height = 120, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(2, 2), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN27126_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ACSS3")
PlotGroups(normcounts[rownames(normcounts) == "DN27274_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ACADM")
PlotGroups(normcounts[rownames(normcounts) == "DN37813_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ApoB-100")
PlotGroups(normcounts[rownames(normcounts) == "DN128_c0_g1_i1", ],   edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "LDLR")
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

# Panel C: Lipid biosynthesis suppression (4 genes, 2x2)
cat("  Saving Figure_4C_lipid_biosynthesis.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_4C_lipid_biosynthesis.tif"),
     width = W_FULL, height = 120, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(2, 2), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN20955_c0_g1_i2", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           ylim = c(0, 60), main = "FASN")
PlotGroups(normcounts[rownames(normcounts) == "DN29654_c0_g1_i2", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ACAC")
PlotGroups(normcounts[rownames(normcounts) == "DN32362_c0_g1_i2", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ACAT2")
PlotGroups(normcounts[rownames(normcounts) == "DN33706_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "FITM2")
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

cat("  Figure 4 panels saved.\n")

# ── 5. Figure 5 — Carbohydrate and amino acid metabolism ──────────────────────
# Note: Panel B (Cahill cycle schematic) is not R-generated.
#       Insert it manually between panels A and C during final assembly.
cat("\n--- Figure 5: Carbohydrate and amino acid metabolism ---\n")

# Panel A-i: Glycogen metabolism (3 genes, 1x3)
cat("  Saving Figure_5Ai_glycogen.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_5Ai_glycogen.tif"),
     width = W_FULL, height = 80, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(1, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN31785_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HK1")
PlotGroups(normcounts[rownames(normcounts) == "DN27640_c0_g4_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "GBE1")
PlotGroups(normcounts[rownames(normcounts) == "DN27339_c0_g1_i2", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "GAA")
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

# Panel A-ii: Nitrogen and carbon cycling (4 genes, 2x2)
cat("  Saving Figure_5Aii_nitrogen_carbon.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_5Aii_nitrogen_carbon.tif"),
     width = W_FULL, height = 120, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(2, 2), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN30324_c1_g1_i3", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "GPT2L")
PlotGroups(normcounts[rownames(normcounts) == "DN1843_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "GDH1")
PlotGroups(normcounts[rownames(normcounts) == "DN32715_c1_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           ylim = c(0, 150), main = "SLC1A5")
PlotGroups(normcounts[rownames(normcounts) == "DN15192_c0_g3_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "GLUT2")
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

# Panel C: Gluconeogenesis (7 genes, 3x3 — last 2 panels blank)
cat("  Saving Figure_5C_gluconeogenesis.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_5C_gluconeogenesis.tif"),
     width = W_FULL, height = 190, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(3, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN20092_c0_g1_i2", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "CREBBP")
PlotGroups(normcounts[rownames(normcounts) == "DN11978_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ACLY")
PlotGroups(normcounts[rownames(normcounts) == "DN81280_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "PCK1")
PlotGroups(normcounts[rownames(normcounts) == "DN23586_c0_g1_i2", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "MDH1")
PlotGroups(normcounts[rownames(normcounts) == "DN31348_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           ylim = c(0, 75), main = "ME1")
PlotGroups(normcounts[rownames(normcounts) == "DN31143_c0_g1_i5", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ENO1")
PlotGroups(normcounts[rownames(normcounts) == "DN32454_c2_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "PKM")
plot.new()   # blank panel 8
plot.new()   # blank panel 9
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

cat("  Figure 5 panels saved.  (Panel B Cahill schematic inserted separately.)\n")

# ── 6. Figure 6 — Cellular stress response ────────────────────────────────────
# Note: Panel C (oxidative stress pathway schematic) is not R-generated.
#       Insert it manually between panels B and D during final assembly.
cat("\n--- Figure 6: Cellular stress response ---\n")

# Panel A: Heat shock proteins (9 genes, 3x3)
cat("  Saving Figure_6A_HSPs.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_6A_HSPs.tif"),
     width = W_FULL, height = 190, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(3, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN30203_c0_g5_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "HSP70")
PlotGroups(normcounts[rownames(normcounts) == "DN26946_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "HSC71")
PlotGroups(normcounts[rownames(normcounts) == "DN20278_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "HSP90A")
PlotGroups(normcounts[rownames(normcounts) == "DN28420_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "HSC90")
PlotGroups(normcounts[rownames(normcounts) == "DN27809_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "Hsp40A")
PlotGroups(normcounts[rownames(normcounts) == "DN29509_c0_g1_i2", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "HSP105")
PlotGroups(normcounts[rownames(normcounts) == "DN33801_c1_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "HSP30C")
PlotGroups(normcounts[rownames(normcounts) == "DN37739_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "HSF1")
PlotGroups(normcounts[rownames(normcounts) == "DN20226_c1_g2_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "Nst1")
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

# Panel B: Hypoxia signaling (3 genes, 1x3)
cat("  Saving Figure_6B_hypoxia.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_6B_hypoxia.tif"),
     width = W_FULL, height = 80, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(1, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN33682_c0_g1_i10", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "EP300")
PlotGroups(normcounts[rownames(normcounts) == "DN27032_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "HIF1A")
PlotGroups(normcounts[rownames(normcounts) == "DN30644_c1_g1_i3", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "HYOU1")
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

# Panel D: Oxidative stress response (9 genes, 3x3)
cat("  Saving Figure_6D_oxidative_stress.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_6D_oxidative_stress.tif"),
     width = W_FULL, height = 190, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(3, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN76450_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "SOD1")
PlotGroups(normcounts[rownames(normcounts) == "DN33384_c0_g1_i3", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "CAT")
PlotGroups(normcounts[rownames(normcounts) == "DN12366_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "XCT")
PlotGroups(normcounts[rownames(normcounts) == "DN33520_c0_g2_i3", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "GSH1")
PlotGroups(normcounts[rownames(normcounts) == "DN26062_c0_g1_i7", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "GSH2")
PlotGroups(normcounts[rownames(normcounts) == "DN24848_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "GPX3")
PlotGroups(normcounts[rownames(normcounts) == "DN16346_c0_g3_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "GLR1")
PlotGroups(normcounts[rownames(normcounts) == "DN19268_c3_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "PRDX1")
PlotGroups(normcounts[rownames(normcounts) == "DN32522_c0_g2_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector, main = "NRF2")
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

cat("  Figure 6 panels saved.  (Panel C oxidative stress schematic inserted separately.)\n")

# ── 7. Figure 7 — Reproductive biology ───────────────────────────────────────
cat("\n--- Figure 7: Reproductive biology ---\n")

# Panel A: Vitellogenins (3 genes, 1x3)
cat("  Saving Figure_7A_vitellogenin.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_7A_vitellogenin.tif"),
     width = W_FULL, height = 80, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(1, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN50748_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "VTG-A2")
PlotGroups(normcounts[rownames(normcounts) == "DN34958_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "VTG-2")
PlotGroups(normcounts[rownames(normcounts) == "DN77992_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "VTG (Fragment)")
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

# Panel B: Hydroxysteroid dehydrogenases (6 genes, 2x3)
cat("  Saving Figure_7B_HSDs.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_7B_HSDs.tif"),
     width = W_FULL, height = 130, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(2, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN45404_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSD3B1")
PlotGroups(normcounts[rownames(normcounts) == "DN11355_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSD17B11")
PlotGroups(normcounts[rownames(normcounts) == "DN33620_c3_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSD17B10")
PlotGroups(normcounts[rownames(normcounts) == "DN42745_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSD17B12")
PlotGroups(normcounts[rownames(normcounts) == "DN6372_c0_g2_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSD17B12 (iso2)")
PlotGroups(normcounts[rownames(normcounts) == "DN15294_c0_g1_i1", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSD11B2")
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

cat("  Figure 7 panels saved.\n")

# ── 8. Done ───────────────────────────────────────────────────────────────────
cat("\n=== All sub-panels saved to:", file.path(WORKING_DIR, OUTPUT_DIR), "===\n")
cat("\nNext step: run figure_assembly.py to combine sub-panels into final figures.\n")
cat("  Figure 2:  2A + 2B\n")
cat("  Figure 4:  4A + 4B + 4C\n")
cat("  Figure 5:  5Ai + 5Aii + [Cahill schematic PNG] + 5C\n")
cat("  Figure 6:  6A + 6B + [oxidative stress schematic PNG] + 6D\n")
cat("  Figure 7:  7A + 7B\n")