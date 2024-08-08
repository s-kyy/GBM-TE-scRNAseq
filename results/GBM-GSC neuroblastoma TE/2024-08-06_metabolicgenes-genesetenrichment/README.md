# GSEA Analysis for GBMSC PC20res0.3

## Prepare Expression dataset file (from 2021-08-30)

The `meanofratios2.R` script in `/bin`. 

The script uses 2 functions: `cellnamesbycluster()` and `meanofratios()`. 

The slurm file used to run the script was: `/data/2021-08-30_gbmscGSEA/slurm_meanofratios.sh`

The output files were  `GE_gsea-table.txt` and `GTE_gsea_table.txt` and associated output files to `ge_gbmsc.out` and `gte_gbmsc.out`

The following code was run on the command line to prepend the necessary information for the GSEA software:

```
$ cd $SCRATCH/gete-gbm/data/2021-08-30_gbmscGSEA/
$ sed -i '1s/^/#1.2\n31527\t11\n/' GE_gsea-table.txt
$ sed -i '1s/^/#1.2\n32501\t10\n/' GTE_gsea-table.txt
$ mv GE_gsea-table.txt GE_table.gct
$ mv GTE_gsea-table.txt GTE_table.gct
$ sed -i 's/"//g' GE_table.gct
$ sed -i 's/"//g' GTE_table.gct
$ sed -i 's/\tNA\t/\tna\t/g' GE_table.gct
$ sed -i 's/\tNA\t/\tna\t/g' GTE_table.gct
```

See `GE_gsea-table.gct` and `GTE_gsea-table.gct`.


## Phenotype Label


The phenotype label files were created manually. Two `.csl` files were written to make sure the number of clusters matched the corresponding datasets.

## Genesets

- H collection: Hallmark gene sets 2015
- C4 subcollection : Curated Cancer Cell Atlas (2023) with Weizmann [Brain | 3CA](https://www.weizmann.ac.il/sites/3CA/brain) 
	- Metamodules for: [[R- Neftel2019_IntegrativeModel|Neftel 2019]], [[R- Couturier2020_SinglecellRNAseq|Couturier 2020]], Tirosh 2016, etc
- C6 collection: oncogenic signature gene sets
- C7 collection: Immunological signature gene sets
- C8 cell type signature gene sets

Combined the genesets into one gmx file (where genesets are organized by columns) on personal PC and uploaded to the cluster as `genesets.gmx`.

```
c4.3ca.v2023.2.Hs.symbols.gmt
c6.all.v2023.2.Hs.symbols.gmt
c7.immunesigdb.v2023.2.Hs.symbols.gmt
c8.all.v2023.2.Hs.symbols.gmt
h.all.v2023.2.Hs.symbols.gmt
```

## Workflow

- copied previous GE_table.gct and GTE_table.gct files from previous run to current run. 
- Ran python script `apply_gsea4.py` with python3 and the same parameters as when I first performed GSEA analysis. Ensured the filepaths were written correctly.

GE : 
```
python3 ~/scratch/gsea/apply_gsea4.py /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls /home/samkyy/scratch/gsea/genesets/c4.3ca.v2023.2.Hs.symbols.gmt c4-3ca gene_set 500 5 160
python3 ~/scratch/gsea/apply_gsea4.py /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls /home/samkyy/scratch/gsea/genesets/c6.all.v2023.2.Hs.symbols.gmt c6-all gene_set 500 5 160
python3 ~/scratch/gsea/apply_gsea4.py /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls /home/samkyy/scratch/gsea/genesets/c7.immunesigdb.v2023.2.Hs.symbols.gmt c7-immune gene_set 500 5 160
python3 ~/scratch/gsea/apply_gsea4.py /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls /home/samkyy/scratch/gsea/genesets/c8.all.v2023.2.Hs.symbols.gmt c8-all gene_set 500 5 160
python3 ~/scratch/gsea/apply_gsea4.py /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls /home/samkyy/scratch/gsea/genesets/h.all.v2023.2.Hs.symbols.gmt c8-all gene_set 500 5 160
```

GTE : 
```
python3 ~/scratch/gsea/apply_gsea4.py /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls /home/samkyy/scratch/gsea/genesets/c4.3ca.v2023.2.Hs.symbols.gmt c4-3ca gene_set 500 5 160
python3 ~/scratch/gsea/apply_gsea4.py /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls /home/samkyy/scratch/gsea/genesets/c6.all.v2023.2.Hs.symbols.gmt c6-all gene_set 500 5 160
python3 ~/scratch/gsea/apply_gsea4.py /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls /home/samkyy/scratch/gsea/genesets/c7.immunesigdb.v2023.2.Hs.symbols.gmt c7-immune gene_set 500 5 160
python3 ~/scratch/gsea/apply_gsea4.py /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls /home/samkyy/scratch/gsea/genesets/c8.all.v2023.2.Hs.symbols.gmt c8-all gene_set 500 5 160
python3 ~/scratch/gsea/apply_gsea4.py /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls /home/samkyy/scratch/gsea/genesets/h.all.v2023.2.Hs.symbols.gmt c8-all gene_set 500 5 160
```

- Sbatch resulting scripts to create GSEA results. 
- Ran bash script `extGSEAReports.sh`
