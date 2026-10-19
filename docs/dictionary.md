# Data dictionary

## Sample labels

| Prefix or pattern | Interpretation | Samples |
|---|---|---:|
| `GM` | Goat manure; `GMMM`, `GMMF`, `GMYM`, `GMYF` encode mature/young and male/female | 12 |
| `RSCR` | Non-integrated control maize rhizosphere, three composites at separate sampling locations within one control field | 3 |
| `RSMN`, `RSME`, `RSMS`, `RSMW` | Integrated maize rhizosphere from N/E/S/W quadrants | 12 |

Each rhizosphere composite pooled five maize plants. Do not treat composite
sampling units as independent farm-level replicates or infer randomized treatment.

## Data matrices

The feature ID is in column 1. Remaining 27 columns are sample IDs. Values are
fractions, not percentage-point numbers, and per-sample feature sums should be
~1. Feature resolution and source-labelled residual categories vary:

| Directory | Files | Column-1 descriptor | Notes |
|---|---:|---|---|
| `data/taxonomy/` | 4 | Kingdom, Phylum, Genus, Species | Micro_NR taxonomy levels; species file includes partial/unresolved labels |
| `data/functional/` | 5 | KEGG Level 1, Level 2, Level 3, KO, Module | Level/category hierarchy and annotations as provided |
| `data/resistance/` | 1 | CARD ARO | Contains `NonARO` category; do not infer absolute ARG load |

## Sample metadata

`data/metadata/samples.csv` columns: `sample_id` (identifier),
`compartment` (GM/RS/CR), `material` (soil/manure), `location` (sampling subgroup or
field/pen label), `age`/`sex` (goat subgroups only; blank for soil), and `batch`
(provider QC source grouping, run1/run2). The `batch` label is not an independent
site or biological replicate.

## QC and sequencing data

`data/qc/qc_s2.csv` contains raw, cleaned and
non-host **bases** (Gb) per sample, plus reported clean Q20/Q30/GC and calculated
base retention or base reductions. Derived percentages use source totals rounded
to 0.01 Gb. `Effective (%)` is an original provider column, not a proxy for
overall retention calculated here. Text files in `data/qc/source/*.tsv` preserve
the original provider tables (original files carried `.xls` extensions but
contain tab-delimited text).

## Figure-source tables

`data/figure_source/` contains 12-row feature-by-sample abundance tables,
with feature label `Label` and per-sample proportions (not percent); 27-row
coverage tables have fields `Sample`, `Displayed_top12_sum_fraction`,
`Displayed_top12_sum_percent`, `Omitted_fraction`, `Omitted_percent`.
`Omitted` is a calculated coverage complement and is **not** added to plotted
`Others`. Figure scripts multiply fractions by 100 for the y-axis labels.
