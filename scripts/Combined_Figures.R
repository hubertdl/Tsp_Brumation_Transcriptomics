# =============================================================================
# Combined_Figures.R - Brumation Manuscript Figure Generation
# =============================================================================
# Author:  David L. Hubert, Oregon State University (hubertd@oregonstate.edu)
# Project: Transcriptomic dynamics of brumation in Thamnophis sirtalis parietalis
# Journal: Journal of Experimental Biology
#
# Purpose:
#   Generate the gene-level expression panels for Figures 5-8 at 300 dpi TIFF.
#   Uses PlotGroups() from maSigPro for individual gene expression panels,
#   with liver (orange, #D55E00) and testis (blue, #0072B2) shown in
#   colorblind-safe Okabe-Ito colours.
#
#   The cluster figures (Figures 2, 3 and S2) are written by the maSigPro
#   scripts (maSigPro_liver.R, maSigPro_Testis.R, maSigPro_combined.R) and
#   Figure S1 by Combined_elbow_Figure.R.
#
# Input files (same directory as main analysis script):
#   raw_counts.tab                  - raw count matrix
#   key.tab                         - maSigPro edesign key (weeks 0, 8, 16, 29, 31)
#
# Output (figures_300dpi/ subdirectory):
#   Figure_5A_lipid_mobilization.tif
#   Figure_5B_lipid_transport.tif
#   Figure_5C_lipid_biosynthesis.tif
#   Figure_6A_carb_nitrogen.tif
#   Figure_6B_gluconeogenesis.tif
#   Figure_7A_HSPs.tif
#   Figure_7B_oxidative_stress.tif
#   Figure_7D_hypoxia.tif
#   Figure_8A_vitellogenin.tif
#   Figure_8B_HSDs.tif
#
#   Then run figure_assembly.py to combine sub-panels into final figures with
#   bold panel labels (A, B, C ...) and insert schematic panels for Figures 6
#   (Cahill cycle, panel C) and 7 (oxidative stress, panel C).
# =============================================================================

# ── 0. Configuration: edit WORKING_DIR if needed ─────────────────────────────
WORKING_DIR   <- "" # <- set working directory
OUTPUT_DIR    <- "figures_300dpi"
DPI           <- 300
DEGREE        <- 4
COVERAGE_THRD <- 5

# JEB figure widths (mm)
W_FULL <- 174    # full page width
W_HALF <- 85.7   # half width: (174 mm - 2.5 mm gap) / 2, for side-by-side panels

# ── 1. Setup ──────────────────────────────────────────────────────────────────
setwd(WORKING_DIR)
suppressPackageStartupMessages({
  library(data.table)
  library(NOISeq)
  library(maSigPro)
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
time    <- edesign[, 1]

# Palette index 2 = first group column in key.tab (testis), 3 = second (liver)
palette(c("black", "#0072B2", "#D55E00", palette.colors(8, "Okabe-Ito")[-1]))

cat("Data loaded.\n")

# ── 3. x-axis labels ──────────────────────────────────────────────────────────
# PlotGroups() centres every week label under its tick, so 29 and 31 (two weeks
# apart) run together ("2931"). TIME_LABELS() is called right after each
# PlotGroups() plot: it blanks the strip PlotGroups() labelled and redraws the
# labels, splitting any pair that would crowd each other so the left label ends
# and the right label starts a fixed gap either side of their midpoint.
# Ticks, data, the "Time" title and the legend are not touched.
#   gap = space between the two split labels, in widths of the character "0"
TIME_LABELS <- function(weeks, gap = 0.8, cex = 1) {
  weeks <- sort(unique(weeks))
  labs  <- as.character(weeks)
  usr   <- par("usr")
  m     <- par("mgp")[2]                                   # line PlotGroups() labels on
  lh    <- par("mai")[1] / par("mar")[1]                   # inches per margin line
  y_at  <- function(l) usr[3] - l * lh * diff(usr[3:4]) / par("pin")[2]
  op <- par(xpd = NA); on.exit(par(op))
  rect(usr[1] - 0.06 * diff(usr[1:2]), y_at(m + 1.25), usr[2], y_at(m - 0.05),
       col = "white", border = NA)                         # hide the original labels

  w   <- strwidth(labs, cex = cex)                         # label widths, user units
  g   <- gap * strwidth("0", cex = cex)
  at  <- weeks
  adj <- rep(0.5, length(weeks))
  for (i in seq_along(weeks)[-1]) {
    if (weeks[i] - weeks[i - 1] < (w[i] + w[i - 1]) / 2 + g) {   # would crowd
      mid <- (weeks[i] + weeks[i - 1]) / 2
      at[i - 1] <- mid - g / 2; adj[i - 1] <- 1
      at[i]     <- mid + g / 2; adj[i]     <- 0
    }
  }
  for (i in seq_along(weeks))
    mtext(labs[i], side = 1, at = at[i], adj = adj[i], line = m, cex = cex * par("cex"))
}

# ── 4. Figure 5: Lipid metabolism program ─────────────────────────────────────
# Laid out to fit one JEB page: panel A full width, panels B and C at half width
# placed side by side by figure_assembly.py. B and C are too narrow for the
# testis/liver legend, so it is left off there (panel A uses the same colours).
cat("\n--- Figure 5: Lipid metabolism ---\n")

# Panel A: Lipid mobilization and beta-oxidation (6 genes, 2x3)
cat("  Saving Figure_5A_lipid_mobilization.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_5A_lipid_mobilization.tif"),
     width = W_FULL, height = 118, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(2, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN30348_c1_g2_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           ylim = c(0, 800), main = "ATGL", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN27598_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           ylim = c(0, 60), main = "FOXO1", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN31376_c0_g1_i2", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "PPARa", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN27065_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ACSL4", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN25931_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           ylim = c(0, 600), main = "CPT1A", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN3769_c0_g1_i1", ],   edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "CPT2", cex.legend = 0.7); TIME_LABELS(time)
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

# Panel B: Substrate switching and lipoprotein transport (4 genes, 2x2, half width)
cat("  Saving Figure_5B_lipid_transport.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_5B_lipid_transport.tif"),
     width = W_HALF, height = 100, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(2, 2), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
par(cex = 0.66)   # same text size as the 3 x 3 figures
PlotGroups(normcounts[rownames(normcounts) == "DN27126_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ACSS3", legend = FALSE); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN27274_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ACADM", legend = FALSE); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN37813_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ApoB-100", legend = FALSE); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN128_c0_g1_i1", ],    edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "LDLR", legend = FALSE); TIME_LABELS(time)
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.7, line = 0.2)
dev.off()

# Panel C: Fatty acid biosynthesis and storage (4 genes, 2x2, half width)
cat("  Saving Figure_5C_lipid_biosynthesis.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_5C_lipid_biosynthesis.tif"),
     width = W_HALF, height = 100, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(2, 2), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
par(cex = 0.66)   # same text size as the 3 x 3 figures
PlotGroups(normcounts[rownames(normcounts) == "DN20955_c0_g1_i2", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           ylim = c(0, 60), main = "FASN", legend = FALSE); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN29654_c0_g1_i2", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ACAC", legend = FALSE); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN32362_c0_g1_i2", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ACAT2", legend = FALSE); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN33706_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "FITM2", legend = FALSE); TIME_LABELS(time)
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.7, line = 0.2)
dev.off()

cat("  Figure 5 panels saved.\n")

# ── 5. Figure 6: Carbohydrate and amino acid metabolism ───────────────────────
# Note: Panel C (Cahill cycle schematic) is not R-generated.
#       It is inserted below panels A and B by figure_assembly.py.
cat("\n--- Figure 6: Carbohydrate and amino acid metabolism ---\n")

# Panel A: Glycogen and glucose handling, alanine cycle (7 genes, 3x3; last 2 panels blank)
cat("  Saving Figure_6A_carb_nitrogen.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_6A_carb_nitrogen.tif"),
     width = W_FULL, height = 190, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(3, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN31785_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HK1", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN27640_c0_g4_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "GBE1", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN27339_c0_g1_i2", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "GAA", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN30324_c1_g1_i3", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "GPT2L", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN1843_c0_g1_i1", ],   edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "GDH1", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN32715_c1_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           ylim = c(0, 150), main = "SLC1A5", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN15192_c0_g3_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "GLUT2", cex.legend = 0.7); TIME_LABELS(time)
plot.new()   # blank panel 8
plot.new()   # blank panel 9
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

# Panel B: Gluconeogenesis (7 genes, 3x3; last 2 panels blank)
cat("  Saving Figure_6B_gluconeogenesis.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_6B_gluconeogenesis.tif"),
     width = W_FULL, height = 190, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(3, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN20092_c0_g1_i2", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "CREBBP", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN11978_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ACLY", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN81280_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "PCK1", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN23586_c0_g1_i2", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "MDH1", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN31348_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           ylim = c(0, 75), main = "ME1", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN31143_c0_g1_i5", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "ENO1", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN32454_c2_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "PKM", cex.legend = 0.7); TIME_LABELS(time)
plot.new()   # blank panel 8
plot.new()   # blank panel 9
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

cat("  Figure 6 panels saved.  (Panel C Cahill schematic inserted separately.)\n")

# ── 6. Figure 7: Cellular stress response ─────────────────────────────────────
# Note: Panel C (oxidative stress pathway schematic) is not R-generated.
#       It is inserted next to panel D by figure_assembly.py.
cat("\n--- Figure 7: Cellular stress response ---\n")

# Panel A: Heat shock proteins (9 genes, 3x3)
cat("  Saving Figure_7A_HSPs.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_7A_HSPs.tif"),
     width = W_FULL, height = 190, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(3, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN30203_c0_g5_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSP70", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN26946_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSC71", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN20278_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSP90A", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN28420_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSC90", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN27809_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "Hsp40A", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN29509_c0_g1_i2", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSP105", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN33801_c1_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSP30C", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN37739_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSF1", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN20226_c1_g2_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "Nst1", cex.legend = 0.7); TIME_LABELS(time)
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

# Panel B: Oxidative stress response (9 genes, 3x3)
cat("  Saving Figure_7B_oxidative_stress.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_7B_oxidative_stress.tif"),
     width = W_FULL, height = 190, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(3, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN76450_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "SOD1", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN33384_c0_g1_i3", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "CAT", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN12366_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "XCT", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN33520_c0_g2_i3", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "GSH1", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN26062_c0_g1_i7", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "GSH2", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN24848_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "GPX3", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN16346_c0_g3_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "GLR1", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN19268_c3_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "PRDX1", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN32522_c0_g2_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "NRF2", cex.legend = 0.7); TIME_LABELS(time)
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

# Panel D: Hypoxia signaling (3 genes, 1x3)
cat("  Saving Figure_7D_hypoxia.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_7D_hypoxia.tif"),
     width = W_FULL, height = 80, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(1, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN33682_c0_g1_i10", ], edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "EP300", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN27032_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HIF1A", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN30644_c1_g1_i3", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HYOU1", cex.legend = 0.7); TIME_LABELS(time)
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

cat("  Figure 7 panels saved.  (Panel C oxidative stress schematic inserted separately.)\n")

# ── 7. Figure 8: Reproductive biology ─────────────────────────────────────────
cat("\n--- Figure 8: Reproductive biology ---\n")

# Panel A: Vitellogenins (3 genes, 1x3)
cat("  Saving Figure_8A_vitellogenin.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_8A_vitellogenin.tif"),
     width = W_FULL, height = 80, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(1, 3), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN50748_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "VTG-A2", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN34958_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "VTG-2", cex.legend = 0.7); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN77992_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "VTG (Fragment)", cex.legend = 0.7); TIME_LABELS(time)
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

# Panel B: Hydroxysteroid dehydrogenases (6 genes, 3x2)
cat("  Saving Figure_8B_HSDs.tif...\n")
tiff(file.path(OUTPUT_DIR, "Figure_8B_HSDs.tif"),
     width = W_FULL, height = 190, units = "mm", res = DPI, compression = "lzw")
par(mfrow = c(3, 2), mar = c(3.8, 3.8, 2.2, 0.6),
    oma = c(0.5, 1.5, 0.5, 0.5), mgp = c(2.2, 0.7, 0))
PlotGroups(normcounts[rownames(normcounts) == "DN45404_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSD3B1"); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN11355_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSD17B11"); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN33620_c3_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSD17B10"); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN42745_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSD17B12"); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN6372_c0_g2_i1", ],   edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSD17B12 (iso2)"); TIME_LABELS(time)
PlotGroups(normcounts[rownames(normcounts) == "DN15294_c0_g1_i1", ],  edesign = edesign,
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector,
           main = "HSD11B2"); TIME_LABELS(time)
mtext("Normalized expression (TMM)", side = 2, outer = TRUE, cex = 0.85, line = 0.2)
dev.off()

cat("  Figure 8 panels saved.\n")
palette("default")

# ── 8. Done ───────────────────────────────────────────────────────────────────
cat("\n=== All sub-panels saved to:", file.path(WORKING_DIR, OUTPUT_DIR), "===\n")
cat("\nNext step: run figure_assembly.py to combine sub-panels into final figures.\n")
cat("  Figure 5:  5A / [5B | 5C]\n")
cat("  Figure 6:  [6A | 6B] / [Cahill schematic PNG]\n")
cat("  Figure 7:  [7A | 7B] / [oxidative stress schematic PNG | 7D]\n")
cat("  Figure 8:  8A / 8B\n")
