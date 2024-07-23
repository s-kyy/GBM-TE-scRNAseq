# Processing scRNA-seq datasets with Seurat

Run the following R scripts in order:

-   `1-createSeuratObj-GTE.R`
-   `2-qc-figs.R`
-   `3-normalize-integrate.R`

The original code was obtained from the following scripts:

-   [`seurat_brain-dataexpo.ipynb` from s-kyy/gete-gbm](https://app.reviewnb.com/s-kyy/gete-gbm/blob/scriptsonly/results%2FGBM_rxiv%2Fseurat_brain-dataexpo.ipynb)
    -   *Normalize and Find Variable Features*
-   [`neuroblastoma.ipynb` from s-kyy/gete-gbm](https://app.reviewnb.com/s-kyy/gete-gbm/blob/scriptsonly/results%2FGBMGSCTE%2Fneuroblastoma.ipynb)
    -   *Load Neuroblastoma Dataset* -\> 1
    -   *Compute Quality Metrics* -\> 1
    -   *Quality Control [Figures]* -\> 2
    -   ~~*Create Unique Identities for each Sample (Metadata)*~~ (todo)
    -   *Normalization, Find Variable Features (2021-06-15)* -\> 3
    -   *UMAP* -\> 3
