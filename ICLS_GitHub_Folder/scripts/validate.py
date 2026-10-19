#!/usr/bin/env python3
"""Validate sample IDs, source matrices, QC, plot source data, and file hashes.

Uses only the Python standard library; can run before pip-installing requirements.
"""
from pathlib import Path
import csv
from collections import Counter
from hashlib import sha256
import math
import sys

ROOT = Path(__file__).resolve().parents[1]

def load(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def checks():
    sample = load(ROOT / "data/metadata/samples.csv")
    sample_ids = [r["sample_id"] for r in sample]
    assert len(sample_ids) == 27 and len(set(sample_ids)) == 27, "Expected 27 unique samples"
    assert Counter(r["compartment"] for r in sample) == Counter({"GM": 12, "RS": 12, "CR": 3})
    qc = load(ROOT / "data/qc/qc_s2.csv")
    assert len(qc) == 27 and {r["Sample ID"] for r in qc} == set(sample_ids)
    for r in qc:
        raw, clean, nonhost = [float(r[c]) for c in ("Raw bases (Gb)", "Clean bases (Gb)", "Non-host bases (Gb)")]
        assert raw > 0 and 0 <= nonhost <= clean <= raw
        assert abs(float(r["Overall base retention (%)"]) - 100*nonhost/raw) <= .025

    matrices = sorted((ROOT/"data/taxonomy").glob("*.csv")) + sorted((ROOT/"data/functional").glob("*.csv")) + sorted((ROOT/"data/resistance").glob("*.csv"))
    assert len(matrices) == 10, f"Expected 10 canonical full matrices, found {len(matrices)}"
    for path in matrices:
        with path.open(encoding="utf-8-sig", newline="") as f:
            r = csv.reader(f)
            header = next(r)
            assert len(header) == 28 and set(header[1:]) == set(sample_ids), path.name
            sums = [0.0]*27
            count = 0
            for row in r:
                if not row:
                    continue
                assert len(row) == len(header) and row[0].strip(), (path.name, count)
                for i, value in enumerate(row[1:]):
                    fval = float(value)
                    assert math.isfinite(fval) and 0 <= fval <= 1.00000001, (path.name, row[0], i, value)
                    sums[i] += fval
                count += 1
            assert count > 0 and all(abs(s - 1) < 0.0002 for s in sums), (path.name,min(sums),max(sums))
            print(f"PASS matrix: {path.name} ({count} rows; per-sample sums {min(sums):.6f}..{max(sums):.6f})")

    for prefix in ("f1", "f2"):
        p = ROOT/"data/figure_source"/f"{prefix}_plot.csv"
        c = ROOT/"data/figure_source"/f"{prefix}_coverage.csv"
        with p.open(encoding="utf-8-sig", newline="") as f:
            reader=csv.DictReader(f)
            h=reader.fieldnames
            assert h[0] == "Label" and set(h[1:]) == set(sample_ids)
            plot=list(reader)
        coverage = load(c)
        assert len(plot)==12 and len(coverage)==27
        coverage_map={r['Sample']:r for r in coverage}
        for s in h[1:]:
            actual=sum(float(row[s]) for row in plot)
            recorded=float(coverage_map[s]['Displayed_top12_sum_fraction'])
            assert abs(actual-recorded) < 1e-8, (prefix,s,actual,recorded)
        print(f"PASS plot source: {prefix} (12 features; 27 samples)")

    manifest = ROOT/"SHA256SUMS.txt"
    if manifest.exists():
        for line in manifest.read_text(encoding='utf8').splitlines():
            if not line.strip(): continue
            hexdigest, relative=line.split('  ',1)
            f=ROOT/relative
            assert f.is_file() and sha256(f.read_bytes()).hexdigest()==hexdigest, f"Checksum mismatch: {relative}"
        print('PASS SHA256SUMS.txt file integrity')
    print('PASS metadata, QC, matrices, figure sources and checksums')

if __name__ == '__main__':
    try:
        checks()
    except (AssertionError, ValueError, FileNotFoundError, KeyError) as exc:
        print('VALIDATION FAILED:',exc, file=sys.stderr)
        raise SystemExit(1)
