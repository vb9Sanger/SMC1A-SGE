#!/usr/bin/env python3
"""
excalibr_fit_quality.py

Computes ExCALIBR's own model-fit-quality metric for a finished run, because
`run_pipeline.py` evaluates it during fitting but does not serialise it into
the outputs it saves.

The metric is the normalised Yang et al. (2019) distance (p=2) between each
sample's empirical CDF and the CDF of the fitted skew-normal mixture for that
sample. Zeiberg, Stewart et al. 2026 treat a distance below 0.2 as a good fit
and report it per sample; across their 80 datasets, all 80 cleared 0.2 on the
gnomAD and synonymous samples (the reliable ones, being much larger) and 78
cleared it on every sample. The two that did not had fewer than five
pathogenic variants.

Each sample's fitted distribution is the same set of components with a
different mixture weight vector -- that shared-component structure is the
point of the multi-sample model -- so `component_params` is a list of
(a, loc, scale) skew-normal triples and `weights[sample]` is that sample's
mixing proportions over them.

Usage:
    python excalibr_fit_quality.py \\
        --fits DEE85_b20_bootstrap_fits.json.gz \\
        --dataset smc1a_DEE85_pathogenic.csv \\
        --sample_names Pathogenic Benign gnomAD Synonymous
"""
import argparse
import gzip
import json

import numpy as np
import pandas as pd
from scipy.stats import skewnorm


def yang_dist(x, y, p=2):
    """Normalised Yang distance, transcribed from the ExCALIBR implementation
    (`MulticomponentCalibrationModel.yang_dist`) so the number reported here is
    the same quantity the paper reports, not a substitute goodness-of-fit."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    gt = x >= y
    dP = ((x[gt] - y[gt]).sum() ** p + (y[~gt] - x[~gt]).sum() ** p) ** (1 / p)
    denom = sum(max(abs(xi), abs(yi), abs(xi - yi)) for xi, yi in zip(x, y))
    return dP / denom if denom else np.nan


def mixture_cdf(xs, params, w):
    return sum(wi * skewnorm.cdf(xs, a, loc, scale)
               for wi, (a, loc, scale) in zip(w, params))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fits", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--sample_names", nargs="+",
                    default=["Pathogenic", "Benign", "gnomAD", "Synonymous"])
    a = ap.parse_args()

    d = pd.read_csv(a.dataset)
    assign = d.sample_assignments.fillna("").astype(str)
    scores = {}
    for i, name in enumerate(a.sample_names):
        m = assign.apply(lambda s: str(i) in s.split(","))
        scores[i] = d.score[m].to_numpy()

    fits = json.load(gzip.open(a.fits))
    rows = []
    for _, per_component in fits.items():
        for comp, f in per_component.items():
            params = f["fit"]["component_params"]
            weights = f["fit"]["weights"]
            for i, name in enumerate(a.sample_names):
                s = np.sort(scores[i])
                if len(s) < 2 or i >= len(weights):
                    continue
                emp = np.arange(1, len(s) + 1) / len(s)
                rows.append(dict(component=comp, sample=name, n=len(s),
                                 dist=yang_dist(emp, mixture_cdf(s, params, weights[i]))))

    r = pd.DataFrame(rows)
    print(f"{len(fits)} bootstrap fits\n")
    g = r.groupby(["component", "sample"]).dist.agg(
        ["count", "median", "min", "max"]).reset_index()
    g["n_variants"] = g["sample"].map(
        {nm: len(scores[i]) for i, nm in enumerate(a.sample_names)})
    g["passes_0.2"] = g["median"].lt(0.2).map({True: "yes", False: "NO"})
    print(g.to_string(index=False))
    bad = r[r.dist >= 0.2]
    print(f"\nfits at or above the 0.2 threshold: {len(bad)}/{len(r)} "
          f"({len(bad) / len(r):.1%})")
    if len(bad):
        print(bad.groupby("sample").dist.agg(["count", "median"]).to_string())


if __name__ == "__main__":
    main()
