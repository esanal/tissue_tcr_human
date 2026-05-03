# tissue_tcr_human

Analysis pipeline for human tissue TCR sequencing data.

The repository is organized in three main steps:

`0.mergeFASTQ`: FASTQ-level sample merging. Samples are combined based on sample IDs or GSIDs.

`1.align`: TCR alignment and clonotype assembly pipeline for merged or non-combined samples.

`2.postAnalysis`: Post-processing of clone tables, sample-level summary metrics, and figures for tissue-level diversity and clone sharing analyses.

## Repository structure

`data/`: Metadata tables, intermediate analysis objects, and exported summary tables.

`src/`: Scripts for sample merging, MiXCR/Snakemake alignment, and downstream R analysis.

`results/`: MiXCR outputs such as `.vdjca`, `.clns`, clone tables, reports, and QC plots.

`figures/`: Summary plots generated during post-analysis.

## Inputs

The pipeline expects paired-end FASTQ files together with metadata describing individuals, organs, subsets, and sequencing IDs.

Metadata used by the current workflow are stored under `data/metadata/`

- `metadata_gsid_cellNumber.csv`

## Step 0. Merge FASTQ files

FASTQ merging script are in `src/0.mergeFASTQ/`.

- `mergeFASTQ_GenomescanIDs.py`: exports FASTQ files per GSID.

Both scripts currently use hard-coded input and output paths so these should be adjusted before running on a new system.

## Step 1. Run MiXCR alignment

Alignment workflow is in `src/1.align/` and are organized as:

- `combinedGSIDs/`

This Snakemake workflow run MiXCR alignment, tag refinement, clonotype assembly, clone export, and QC export.

Current configs use the human preset:

- species: `hsa`
- preset: `takara-human-rna-tcr-umi-smartseq`
- UMI read threshold: `1`

Example:

```bash
cd src/1.align/nonCombined
snakemake -j <n>
```

Adjust `config.yaml` first, especially `results_dir`, memory settings and any data directory paths.

## Step 2. Post-analysis

Post-analysis scripts are in `src/2.postAnalysis/`.

`post_analysis.R` reads exported clone tables, splits TRA and TRB chains, summarizes clone statistics, joins metadata, and generates:

- sample-level metrics in `data/summary/`
- histogram panels in `figures/summary/`
- clones-per-cell plots in `figures/clones_per_cell/`
- clone distribution plots in `figures/clone_dist/`

The script currently reads clone tables from `./results/results_combinedGSIDs` and metadata from `./data/metadata/metadata_gsid_cellNumber.csv`.

## Software

Python dependencies are defined in `pyproject.toml`.

The workflow also depends on:

- `Snakemake`
- `MiXCR 4.7.0`
- R packages used by `src/2.postAnalysis/post_analysis.R`

## Notes

Several scripts still contain machine-specific absolute paths[TO BE FIXED]. Before re-running the full workflow, update these paths to match the local environment.
