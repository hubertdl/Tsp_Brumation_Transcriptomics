# Tsp_Brumation_Transcriptomics

Code and processed data for:

**Cold, dark, and hungry: Dynamic transcriptional regulation across eight months of brumation**
David L. Hubert, Ehren J. Bentz, Robert T. Mason


Raw sequencing reads are deposited at NCBI SRA under BioProject PRJNA1473998.

---

## Repository structure

```
Tsp_Brumation_Transcriptomics/
├── data/
│   ├── raw_countsL.tab          # Raw count matrix, liver (transcripts × samples)
│   ├── raw_countsT.tab          # Raw count matrix, testis (transcripts × samples)
│   ├── raw_counts.tab           # Raw count matrix, liver + testis combined
│   ├── Lkey.tab                 # Sample IDs, timepoints (weeks 0, 8, 16, 29, 31), tissue for liver
│   ├── Tkey.tab                 # Sample IDs, timepoints (weeks 0, 8, 16, 29, 31), tissue for testis
│   └── key.tab                  # Sample IDs, timepoints (weeks 0, 8, 16, 29, 31), tissue for liver + testis combined
├── scripts/
│   ├── maSigPro_liver.R         # time-series differential expression and clustering, k = 11 (liver; Fig. 2)
│   ├── maSigPro_Testis.R        # time-series differential expression and clustering, k = 10 (testis; Fig. 3)
│   ├── maSigPro_combined.R      # time-series differential expression and clustering, k = 10 (liver + testis combined; Fig. S2)
│   ├── Combined_elbow_Figure.R  # Manuscript figure (Fig. S1, cluster number selection)
│   ├── Combined_Figures.R       # Manuscript figure panels (Figs. 5–8)
│   └── figure_assembly.py       # Assembles panels and schematics into final Figs. 5–8
└── README.md
```


## Contact

David L. Hubert - hubertd@oregonstate.edu
Department of Integrative Biology, Oregon State University
