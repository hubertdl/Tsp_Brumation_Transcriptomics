#!/usr/bin/env python3
"""
figure_assembly.py - Brumation manuscript figure assembly
Hubert, Bentz and Mason, Journal of Experimental Biology

Combines the TIFF sub-panels written by Combined_Figures.R (figures_300dpi/)
with the two hand-drawn schematics into the final figures, laid out and
lettered as in the manuscript:

  Figure_5_lipid_metabolism.tif       A full width; B | C side by side below
  Figure_6_carbohydrate.tif           A | B side by side; C (Cahill schematic) full width below
  Figure_7_stress_response.tif        A | B side by side; C (oxidative-stress schematic) | D below
  Figure_8_vitellogenin_steroids.tif  A / B stacked, full width

Figures 2, 3 and S2 are written by the maSigPro scripts and Figure S1 by
Combined_elbow_Figure.R. Figures 1 and 4 are illustrations.

Usage: set WORKING_DIR below (the folder Combined_Figures.R was run in), then
    python figure_assembly.py

Requires: Pillow   (python -m pip install Pillow)
"""

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

# =============================================================================
# Configuration
# =============================================================================

WORKING_DIR = Path("")                       # <- set working directory
INPUT_DIR   = WORKING_DIR / "figures_300dpi"  # written by Combined_Figures.R
OUTPUT_DIR  = WORKING_DIR / "figures_final"   # assembled figures are written here

CAHILL_SCHEMATIC    = WORKING_DIR / "schematics" / "Cahill_schematic.png"             # Fig. 6C
OXIDATIVE_SCHEMATIC = WORKING_DIR / "schematics" / "oxidative_stress_schematic.png"   # Fig. 7C

DPI          = 300
FIG_WIDTH_MM = 174                                  # JEB full-page width
FIG_WIDTH_PX = round(FIG_WIDTH_MM / 25.4 * DPI)     # 2055 px at 300 dpi

COL_GAP_PX     = 30     # white space between side-by-side panels
ROW_GAP_PX     = 40     # white space between stacked panels
LABEL_SIZE_PX  = 64     # height of the bold panel letters (about 5.4 mm)
LABEL_OFFSET   = (6, 2) # (x, y) position of each letter within its panel

# =============================================================================
# Helpers
# =============================================================================

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def _font(size, bold=True):
    names = (["arialbd.ttf", "C:/Windows/Fonts/arialbd.ttf",
              "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
              "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
              "/System/Library/Fonts/Supplemental/Arial Bold.ttf"]
             if bold else
             ["arial.ttf", "C:/Windows/Fonts/arial.ttf",
              "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
              "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
              "/System/Library/Fonts/Supplemental/Arial.ttf"])
    for n in names:
        try:
            return ImageFont.truetype(n, size)
        except (OSError, IOError):
            continue
    print("  NOTE: TrueType font not found; using the default font.")
    return ImageFont.load_default()


LABEL_FONT = _font(LABEL_SIZE_PX)
NOTE_FONT  = _font(40, bold=False)


def load(filename):
    """Open a panel from INPUT_DIR, or return None with a warning."""
    return load_file(INPUT_DIR / filename)


def load_file(path):
    """Open an image from any path, or return None with a warning."""
    path = Path(path)
    if not path.exists():
        print(f"  WARNING: {path} not found.")
        return None
    return Image.open(path).convert("RGB")


def placeholder(width, height, text):
    """Grey-outlined white box marking a missing panel, so gaps are obvious."""
    img = Image.new("RGB", (width, height), "white")
    d = ImageDraw.Draw(img)
    d.rectangle([2, 2, width - 3, height - 3], outline=(170, 170, 170), width=4)
    d.text((30, height // 2 - 20), f"MISSING: {text}", fill=(150, 150, 150), font=NOTE_FONT)
    return img


def fit_width(img, width, missing_text="panel", missing_aspect=0.6):
    """Scale img to an exact width, keeping its aspect ratio."""
    if img is None:
        return placeholder(width, int(width * missing_aspect), missing_text)
    if img.width == width:
        return img
    h = round(img.height * width / img.width)
    return img.resize((width, h), Image.LANCZOS)


def label(img, letter):
    """Stamp a bold panel letter in the top-left corner."""
    if not letter:
        return img
    out = img.copy()
    ImageDraw.Draw(out).text(LABEL_OFFSET, letter, fill="black", font=LABEL_FONT)
    return out


def hstack(images, gap=COL_GAP_PX):
    """Place images side by side, top-aligned, on a white canvas."""
    w = sum(i.width for i in images) + gap * (len(images) - 1)
    h = max(i.height for i in images)
    canvas = Image.new("RGB", (w, h), "white")
    x = 0
    for i in images:
        canvas.paste(i, (x, 0))
        x += i.width + gap
    return canvas


def vstack(images, gap=ROW_GAP_PX):
    """Stack images vertically, left-aligned, on a white canvas."""
    images = [i for i in images if i is not None]
    w = max(i.width for i in images)
    h = sum(i.height for i in images) + gap * (len(images) - 1)
    canvas = Image.new("RGB", (w, h), "white")
    y = 0
    for i in images:
        canvas.paste(i, (0, y))
        y += i.height + gap
    return canvas


def row(blocks, total_width=FIG_WIDTH_PX, gap=COL_GAP_PX):
    """
    A row of equal-width labelled panels spanning total_width.
    blocks: list of (image_or_None, letter, name_for_warning).
    """
    n = len(blocks)
    w = (total_width - gap * (n - 1)) // n
    widths = [w] * (n - 1) + [total_width - gap * (n - 1) - w * (n - 1)]  # last absorbs rounding
    return hstack([label(fit_width(img, wi, name), letter)
                   for (img, letter, name), wi in zip(blocks, widths)], gap)


def full(img, letter, name):
    """One full-width labelled panel."""
    return label(fit_width(img, FIG_WIDTH_PX, name), letter)


def save(img, filename):
    if img is None:
        print(f"  Skipped {filename}: no panels found.")
        return
    out = OUTPUT_DIR / filename
    img.save(out, dpi=(DPI, DPI), compression="tiff_lzw")
    print(f"  Saved {out}  ({img.width} x {img.height} px, "
          f"{img.width / DPI * 25.4:.0f} x {img.height / DPI * 25.4:.0f} mm)")


# =============================================================================
# Figure 5: lipid metabolism
#   top:    A (lipid mobilization and beta-oxidation), full width
#   bottom: B (substrate switching, lipoprotein transport) | C (biosynthesis, storage)
# =============================================================================

print("\n=== Figure 5: lipid metabolism ===")
save(vstack([
        full(load("Figure_5A_lipid_mobilization.tif"), "A", "lipid mobilization"),
        row([(load("Figure_5B_lipid_transport.tif"),    "B", "lipid transport"),
             (load("Figure_5C_lipid_biosynthesis.tif"), "C", "lipid biosynthesis")]),
     ]),
     "Figure_5_lipid_metabolism.tif")

# =============================================================================
# Figure 6: carbohydrate metabolism
#   top:    A (glycogen, glucose handling, alanine cycle) | B (gluconeogenesis)
#   bottom: C Cahill schematic, full width
# =============================================================================

print("\n=== Figure 6: carbohydrate and amino acid metabolism ===")
save(vstack([
        row([(load("Figure_6A_carb_nitrogen.tif"),    "A", "carbohydrate panels"),
             (load("Figure_6B_gluconeogenesis.tif"),  "B", "gluconeogenesis")]),
        full(load_file(CAHILL_SCHEMATIC), "C", "Cahill schematic"),
     ]),
     "Figure_6_carbohydrate.tif")

# =============================================================================
# Figure 7: stress response
#   top:    A (heat shock proteins) | B (oxidative stress)
#   bottom: C (oxidative stress schematic) | D (hypoxia)
# =============================================================================

print("\n=== Figure 7: cellular stress response ===")
save(vstack([
        row([(load("Figure_7A_HSPs.tif"),             "A", "heat shock proteins"),
             (load("Figure_7B_oxidative_stress.tif"), "B", "oxidative stress")]),
        row([(load_file(OXIDATIVE_SCHEMATIC),         "C", "oxidative stress schematic"),
             (load("Figure_7D_hypoxia.tif"),          "D", "hypoxia")]),
     ]),
     "Figure_7_stress_response.tif")

# =============================================================================
# Figure 8: vitellogenin and steroid-modulating enzymes (A, B stacked)
# =============================================================================

print("\n=== Figure 8: vitellogenin and steroid metabolism ===")
save(vstack([
        full(load("Figure_8A_vitellogenin.tif"), "A", "vitellogenins"),
        full(load("Figure_8B_HSDs.tif"),         "B", "hydroxysteroid dehydrogenases"),
     ]),
     "Figure_8_vitellogenin_steroids.tif")

print(f"\nDone. Final figures are in: {OUTPUT_DIR.resolve()}")
