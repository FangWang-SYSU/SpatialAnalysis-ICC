# SpatialAnalysis-ICC

Analysis code for integrated spatial transcriptomics, scRNA-seq, and TCR profiling in ICC.

## Project Description

This project accompanies the study **Mapping spatial transcriptional plasticity to clinical immune signatures in intrahepatic cholangiocarcinoma**. The main analyses include:

- Spatial data processing, cell-type annotation, and visualization
- Identification and scoring of spatial ecosystem programs (SEPs)
- CNV inference and spatial characterization of tumor clones
- Associations of SEPs with treatment response and TCR clonal expansion
- Periclonal microenvironment and cell-cell communication analyses

## System Requirements

Analyses use R and Python 3.11. CPU and memory requirements depend on dataset size, particularly for full-resolution Open-ST data.

## Software Dependencies

### R Packages

- **Seurat**, **UCell**, **scalop** — Single-cell analysis and signature scoring
- **NMF**, **infercnv**, **CytoTRACE** — Program discovery and clonal characterization
- **fgsea**, **clusterProfiler**, **mixtools** — Enrichment analysis and response-associated grouping
- **CytoSignal**, **CellChat** — Cell-cell communication analysis
- **ggplot2**, **ggsci**, **ComplexHeatmap** — Visualization

### Python Packages and Other Tools

- Scanpy, AnnData, NumPy, pandas, SciPy**, scikit-learn, Matplotlib
- Cell Ranger, Open-ST, Spacemake, FICTURE

Versions reported in the manuscript: Cell Ranger 6.1.2, Open-ST tools 0.2.3, Spacemake 0.8.0, FICTURE 0.0.5.0, Seurat 5.5.0, inferCNV 1.26.0, CytoTRACE 0.3.3, CytoSignal 0.2.2, and CellChat 2.2.0.

## Instructions for Use

1. **Prepare the environment:** Install the dependencies required by the selected analysis.
2. **Configure inputs:** Update sample identifiers, data paths, reference files, and output directories in the scripts or notebooks.
3. **Run the analyses:** Perform preprocessing and annotation before SEP and CNV analyses. Run downstream spatial and response-associated analyses after generating the required intermediate objects.
4. **Generate figures:** Run the corresponding R or Jupyter analysis sections with their processed inputs.

## Data Availability

Data identifiers reported in the manuscript:

- **Raw sequencing data:** GSA-Human accession [PRJCA080529](https://ngdc.cncb.ac.cn/)
- **Processed scRNA-seq and spatial data:** Zenodo DOI [10.5281/zenodo.24123314
        
        ](https://doi.org/10.5281/zenodo.24123314
        
        
        
        )

## Notes

- Scripts require study-specific data, metadata, and intermediate objects; update paths before running.
- Spatial R-like/NR-like labels are score-derived ecosystem states, distinct from clinically observed R/NR outcomes.
- Refer to the manuscript Methods for analysis parameters and statistical procedures.
