# ICLS shotgun metagenomics: processed data and figure scripts

**Windows-friendly release:** Extract to `C:\ICLS`. The original-to-short filename mapping is in [`docs/filename_map.csv`](docs/filename_map.csv).

**Dataset:** *Shotgun metagenomics dataset of goat manure and maize rhizosphere from an integrated crop-livestock system*.  
**Associated manuscript:** *Data in Brief*, DIB-D-26-01861 (manuscript under revision).  
**Status:** Repository package; open access license and public publishing.

## Dataset and intended reuse

This repository provides machine-readable **Micro_NR taxonomic**, **KEGG functional**,
**CARD ARO predicted resistance-associated annotation** matrices, sample-level
metadata, **27-sample sequencing QC base statistics**, and the Python scripts
that generated the final two top-12 manuscript bar plots.

There are **27 sequenced sample units**: goat manure (**GM**, *n*=12), integrated
maize rhizosphere (**RS**, *n*=12), and non-integrated control rhizosphere (**CR**,
*n*=3). Each rhizosphere composite pools material from five maize plants.
The three control composites were sampled at three separate positions **within
one control field**.
The 12 integrated rhizosphere composites come from four quadrants with three
samples per quadrant. This is a single-farm observational dataset.

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
Its provider input totals are rounded to **0.01 Gb**.

## Limits of computational reproducibility

The Python code here is the **downstream manuscript figure-generation code**,
not the complete pipeline for fastp/Bowtie2/MEGAHIT/gene prediction,
CD-HIT, Micro_NR/KEGG annotation or CARD/RGI profiling. The manuscript
specified a workflow and parameters. See [provider parameters](docs/params.md).
The identical processing parameters were used for
**all 27 samples**.

CARD annotations are **sequence-based predictions**.

## Citation, rights and release

The manuscript is under revision. See [release checklist](docs/RELEASE_CHECKLIST.md).

## Licensing

This repository uses separate licenses for code and data:

- **Python scripts (`scripts/`):** MIT License; see [`LICENSE`](LICENSE).
- **Researcher-produced processed data (`data/`) and figures (`figures/`):** Creative Commons Attribution 4.0 International (CC BY 4.0); see [`DATA_LICENSE.md`](DATA_LICENSE.md).
- **Original third-party/provider source files:** not automatically covered by either license. In particular, see the exception for `data/qc/source/` in `DATA_LICENSE.md`; redistribution and relicensing require applicable permission.
- **Raw SRA sequence reads:** maintained under the cited NCBI SRA accession records, rather than relicensed here.

Please cite the manuscript-associated archived release and follow the attribution instructions in `CITATION.cff` and `DATA_LICENSE.md`.

