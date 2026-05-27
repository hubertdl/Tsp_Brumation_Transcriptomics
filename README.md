# Tsp_Brumation_Transcriptomics

Code and processed data for:

**Active transcriptional remodeling during brumation in the red-sided garter snake (*Thamnophis sirtalis parietalis*)**
David L. Hubert, Ehren J. Bentz, Robert T. Mason


Raw sequencing reads and the reference transcriptome are deposited at NCBI SRA under BioProject [accession TBD].

---

## Repository structure

```
Tsp_Brumation_Transcriptomics/
├── data/
│   ├── raw_countsL.tab          # Raw count matrix, liver (transcripts × samples)
│   ├── raw_countsT.tab         # Raw count matrix, testis (transcripts × samples)
│   ├── raw_counts.tab       # Raw count matrix, liver + testis combined
│   ├── Lkey.tab       # Sample IDs, timepoints, tissue, collection dates for liver
│   ├── Tkey.tab       # Sample IDs, timepoints, tissue, collection dates for testis
│   └── key.tab       # Sample IDs, timepoints, tissue, collection dates for liver + testis combined
├── scripts/
│   ├── maSigPro_liver.R          # time-series differentail expression (liver)
│   ├── maSigPro_testis.R         # time-series differentail expression (testis)
│   ├── maSigPro_combined.R       # time-series differentail expression (liver + testis combined)
│   ├── Combined_Figures.R        # Manuscript figures (Figs. 2–8, S2)
│   └── Combined_elbow_Figure.R   # Manuscript figure (Figs. S1)           
└── README.md
```


## Contact

David L. Hubert — hubertd@oregonstate.edu
Department of Integrative Biology, Oregon State University
