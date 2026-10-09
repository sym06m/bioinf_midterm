# Project 18: Short Reads Versus Long Reads (E. coli K-12 MG1655)

Course project, Introduction to Bioinformatics, Astana IT University.
Authors: [Student 1], [Student 2].

**Question.** When does the extra per-base cost of long-read sequencing pay off for assembling a bacterial genome, and which genomic features does each technology fail on?

**Data (all public, retrieved 2026-10-08).**

| What | Accession | Platform |
|---|---|---|
| Reference genome | GCF_000005845.2 (NC_000913.3) | finished RefSeq assembly |
| Short reads | SRR1030394 | Illumina MiSeq, 2 x 251 bp, 285x |
| Long reads | ERR14686234 | Oxford Nanopore MinION, 36x |

Checksums and retrieval dates: `data/PROVENANCE.tsv` (written by `scripts/fetch_ena.py`).

## Repository layout

The Snakefile uses these paths, so keep this layout.

```
.
├── README.md
├── Snakefile                 workflow (rules mirror the Colab notebooks)
├── config.yaml               accessions, coverages, seeds, threads
├── environment.yml           conda environment (p18)
├── environment.lock.yml      exact versions used in Colab (from results/)
├── notebooks/
│   ├── 01_pipeline_part1.ipynb     download, strain check, QC, subsampling, assemblies, QUAST
│   └── 02_analysis_part2.ipynb     long_30x diagnosis, replicates, SPAdes hybrid, polishing, errors by context
├── scripts/
│   ├── fetch_ena.py          download an ENA run, verify md5, write PROVENANCE
│   ├── analysis.py           uncovered regions and indel/homopolymer analysis
│   ├── make_figures.py       figures 1-3
│   ├── cost_model.py         cost per assembly from published prices
│   └── make_test_data.sh     builds the 100 kb test dataset in tests/
├── results/                  small result tables (QC summaries, QUAST tables, CSVs)
├── figures/                  fig1_assembly_metrics.png, fig2_indel_homopolymers.png, fig3_missing_regions.png
├── tests/                    100 kb test dataset (created by make_test_data.sh)
└── report/                   report.pdf
```

Raw reads, subsamples and assemblies are not committed (see `.gitignore`).

## Reproduce

```bash
conda env create -f environment.yml
conda activate p18
snakemake --cores 2 -n      # dry run: check that all inputs resolve
snakemake --cores 2         # full run (about 2 to 3 hours on 2 cores)
```

[Before submission, fill in: has the full run been executed end to end? Date and machine.]

To check the pipeline quickly on a small input, build the test dataset once the raw data are present:

```bash
bash scripts/make_test_data.sh      # writes tests/ (100 kb of the reference and the reads that map to it)
```

[Before submission, fill in: has make_test_data.sh been run and the test data verified?]

## Pipeline

1. **Fetch** short and long reads from the ENA API (md5 checked) and the RefSeq reference.
2. **QC.** fastp for the short reads (`-q 20 -l 100`, 3' tail trimming); NanoPlot for the long reads (not filtered).
3. **Subsample** with rasusa to 10, 20 and 30x (seed 42; replicate seed 43 for long reads).
4. **Assemble.** SPAdes `--isolate` (short), Flye `--nano-raw` (long), SPAdes `--nanopore` (hybrid), Flye + Polypolish (polished).
5. **Evaluate** with QUAST against NC_000913.3, then `analysis.py` for uncovered regions and homopolymer context, `cost_model.py` for cost.

## Known issues

- Unicycler failed twice (SPAdes memory limit of 1 GB set by Unicycler), so the hybrid is SPAdes hybrid mode.
- The Flye 30x assembly with seed 42 is duplicated (9.46 Mb, duplication ratio 1.994); it is kept in the tables and excluded from curves and the cost model. Replicate seed 43 is normal.
- Peak memory was not recorded.

## AI assistance

Disclosed in Appendix A of the report.
