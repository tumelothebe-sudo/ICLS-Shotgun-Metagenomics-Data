# ICLS shotgun metagenomics: processed data and figure scripts

**Windows-friendly release:** short filenames; no extra top-level ZIP directory. Extract to `C:\ICLS`. The original-to-short filename mapping is in [`docs/filename_map.csv`](docs/filename_map.csv).

**Dataset:** *Shotgun metagenomics dataset of goat manure and maize rhizosphere from an integrated crop-livestock system*.  
**Associated manuscript:** *Data in Brief*, DIB-D-26-01861 (manuscript under revision).  
**Status:** Draft repository package; verify release permissions, license and public links before publishing.

## Dataset and intended reuse

This repository provides machine-readable **Micro_NR taxonomic**, **KEGG functional**,
**CARD ARO predicted resistance-associated annotation** matrices, sample-level
metadata, **27-sample sequencing QC base statistics**, and the Python scripts
that generated the final two non-normalized top-12 manuscript bar plots.

There are **27 sequenced sample units**: goat manure (**GM**, *n*=12), integrated
maize rhizosphere (**RS**, *n*=12), and non-integrated control rhizosphere (**CR**,
*n*=3). Each rhizosphere composite pools material from five maize plants.
The three control composites were sampled at three separate positions **within
one control field**; they do **not** constitute independent control fields.
The 12 integrated rhizosphere composites come from four quadrants with three
samples per quadrant. This is a single-farm observational dataset, not a
balanced treatment design.

Raw paired-end reads are deposited in NCBI SRA (not mirrored here):
[SRP684228](https://www.ncbi.nlm.nih.gov/sra/SRP684228),
[SRP684229](https://www.ncbi.nlm.nih.gov/sra/SRP684229), and
[SRP684230](https://www.ncbi.nlm.nih.gov/sra/SRP684230).

**Publication record and persistent code/data DOI:** `TO BE ADDED AFTER PUBLIC DEPOSITION`.

## Directory layout

```text
.
├── data/
│   ├── taxonomy/               # 4 full Micro_NR taxonomic matrices
│   ├── functional/             # 5 full KEGG abundance matrices
│   ├── resistance/             # 1 CARD ARO abundance matrix
│   ├── metadata/               # 27-sample S1 metadata (CSV)
│   ├── qc/                     # Supplementary Table S2 and original QC summaries
│   │   └── source/
│   └── figure_source/          # Exact 12-feature plotting and coverage tables
├── scripts/                   # Two final plotting scripts; reproduction & QC validation
├── figures/                   # Reference PNG files corresponding to the manuscript figures
├── docs/                      # Methodological details, limitations, provenance and release checklist
├── requirements.txt           # Figure-script Python dependencies
├── CITATION.cff               # Provisional author metadata; replace repository URL
└── LICENSE_PENDING.md         # Permission/licensing decisions before public release
```

All matrices are **sample-level** and contain a first column naming each taxon
or functional annotation, followed by 27 sample columns. Values are unitless
relative abundance **fractions**, not percentages (e.g. `0.01` = 1%). Sample
columns in the full matrices sum to approximately 1 within floating-point
rounding. `Others`/`NonARO` are categories already present in supplier/derived
source matrices, not newly inferred residual annotations.

See [data dictionary](docs/dictionary.md) and [provenance](docs/provenance.csv).

## Reproduce the main figures

Requires **Python 3.10+** with NumPy, pandas and Matplotlib.

```bash
python -m pip install -r requirements.txt
python scripts/validate.py
python scripts/reproduce.py
```

The two final plotting scripts are also individually runnable in Spyder:

- `scripts/species.py`
- `scripts/kegg_l3.py`

Figures will appear under `output/figures/` as PDF, SVG, and 600-dpi PNG. Their
cleaned top-12 tables and per-sample coverage CSVs will be written to the same
folder. Compare them with `data/figure_source/` or inspect `figures/`.

**Plotting rule:** top 12 displayed features are ranked by mean across all 27
samples; the supplied source-file `Others` remains eligible for selection. We
**do not** pool omitted lower-abundance features into `Others`, and we **do not**
rescale displayed features to 100%. The y-axis remains 0–100% so the plotted
subset is interpreted against the entire annotated profile. Some Micro_NR
species-level names are lineage-level or other non-binomial labels, not formal
species identifications.

## QC interpretation

`data/qc/qc_s2.csv` reports raw, cleaned and
non-host sequence **base totals (Gb)** and calculated base-retention metrics.
Its provider input totals are rounded to **0.01 Gb**; calculated differences
are approximate. A **0.00 Gb** difference does not establish that no reads
were removed. These files do **not** provide exact fastp-discarded read counts,
Bowtie2 host-aligned read counts or exact post-trimming read-pair counts.
Do not treat changes in bases as confirmed counts of host-mapped reads.

## Limits of computational reproducibility

The Python code here is the **downstream manuscript figure-generation code**,
not the vendor's complete pipeline for fastp/Bowtie2/MEGAHIT/gene prediction,
CD-HIT, Micro_NR/KEGG annotation or CARD/RGI profiling. The provider report
specified a workflow and parameters, but raw upstream commands/logs were not
supplied. See [provider parameters](docs/params.md).
The researchers confirmed that identical processing parameters were used for
**all 27 samples**; the provider-version differences from earlier manuscript
text have been recorded for editorial reconciliation.

CARD annotations are **sequence-based predictions** only. They do not demonstrate
phenotypic antimicrobial resistance, transmission, or validated surveillance.

## Citation, rights and release

The manuscript is under revision. Before depositing, confirm author order,
choose the code/data license, replace the placeholder URL in `CITATION.cff`,
publish GitHub release `v1.0.0`, connect Zenodo to obtain a version-specific
DOI (or archive in an approved data repository), and cite the DOI in the
manuscript. See [release checklist](docs/RELEASE_CHECKLIST.md).
