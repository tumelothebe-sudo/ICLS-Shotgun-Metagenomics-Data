# Release checklist — complete before making this GitHub repository public

- [ ] Obtain co-author/institution approval to publicly distribute processed
      matrices and provider QC summaries (do not upload confidential vendor
      reports, raw FASTQ duplicates or reviewer correspondence).
- [ ] Verify `CITATION.cff` author order and manuscript details; replace the
      placeholder GitHub account/URL.
- [ ] Agree code and data licenses; add actual LICENSE and explicit data license
      (none has been granted by this draft package).
- [ ] Verify methods-version reconciliation with provider report and submitted
      manuscript; team confirms settings were identical for all 27 samples.
- [ ] From project root run `python scripts/validate.py`.
- [ ] Run `python scripts/reproduce.py` and inspect exported images.
- [ ] Confirm S1 and S2 supplementary tables and figure captions match the
      final submitted manuscript, including that the main plots are not
      normalized and lower-ranked features are not pooled into source `Others`.
- [ ] Upload **extracted files/folders**, not the ZIP as a single GitHub file.
- [ ] Create a GitHub release tag `v1.0.0` for the manuscript version.
- [ ] Archive the release in Zenodo to obtain a version-specific DOI (or upload
      processed datasets to an approved repository such as Figshare).
- [ ] Enter the permanent repository DOI(s), dataset availability dates, and
      GitHub URL into the manuscript and point-by-point reviewer response.
- [ ] Do not claim that full Novogene upstream code or exact Bowtie2-host read
      counts are deposited: only the documented settings and downstream plotting
      scripts are available in this package.
