# Documented upstream metagenomics workflow and limitations

**Source:** uploaded Novogene `X201SC25105229-Z01-F002_meta_Report.html`,
Appendix 5.1 (software table) and methods narrative. The research team has
confirmed that processing parameters were **identical across all 27 samples**,
including the seven additional rhizosphere samples. The report itself initially
listed 20 samples, with seven additional samples subsequently incorporated.

| Stage | Provider report software / version | Reported parameters (transcribed) |
|---|---|---|
| Quality filtering | fastp 0.23.1 | `-g -q 5 -u 50 -n 15 -l 150 --overlap_diff_limit 1 --overlap_diff_percent_limit 10` |
| Host screening / mapping | Bowtie2 2.5.4 | `--end-to-end --sensitive --no-hd --no-sq -I 200 -X 400 --threads 8` |
| Assembly | MEGAHIT 1.2.9 | `--presets meta-large -m 180000000000 -t 2` |
| Gene prediction | GeneMark.hmm 2.1 | `-a -d -f G -p 1` |
| ORF dereplication | CD-HIT 4.5.8 | `-T 6 -G 0 -aS 0.9 -g 1 -d 0 -c 0.95 -n 5 -M 8000` |
| Functional/taxonomic annotation | DIAMOND 2.1.9 | `blastp`, reported E-value ≤ 1e-5 |
| CARD profiling | RGI 6.0.2 / CARD 3.2.6 | Sequence-based predicted annotations |

**This table is a record of the provider's documented settings, not executable
command history.** Read paths, reference genome and Bowtie2 index identity,
fastp JSON logs and Bowtie2 counts, assembly intermediate files, database
release identifiers for some annotations and provider execution scripts were not
supplied. Do not claim that this repository reproduces raw-FASTQ-to-annotation
processing.

**Manuscript discrepancy to resolve before submission:** earlier draft states
fastp 0.23.4; MEGAHIT 1.1.2; MetaGeneMark 2.10; CD-HIT 4.6.8. The documented
Novogene report instead lists the versions above. The team confirms *identical
parameters across samples* but should approve the version numbers used in the
final Methods section. Full workflow execution logs would be definitive.
