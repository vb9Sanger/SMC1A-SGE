#!/usr/bin/env python3
"""
excalibr_pathomechanism.py

Estimates P(M = 1 | Y = 1) -- the fraction of a disease arm's pathogenic
variants whose mechanism the assay actually measures -- by calling ExCALIBR's
own boundary estimator on bootstrap fits the pipeline has already produced.

**Why this is a separate script rather than a pipeline flag.** ExCALIBR
exposes this as `--pathomechanism-prior`, but that route cannot run on SMC1A:

  * Without `--manual-prior`, the population prior is estimated at ~2e-6 and
    discarded by the <= 1e-3 guard in `get_fit_prior`. Every downstream
    quantity is then NaN, including the pathomechanism prior
    P(Y=1, M=1) = P(Y=1) x P(M=1|Y=1), and the run aborts reporting
    "0/N fits had any assay-relevant excess/mask (degenerate on every fit)".
    That message is misleading here: the mechanism decomposition itself is
    fine, it is the P(Y=1) factor multiplying it that is NaN.
  * With `--manual-prior`, the pipeline warns that the pathomechanism path is
    ignored, and then raises `NameError: cannot access free variable
    'fit_priors_benign'` (visualize.py:676) -- the manual-prior branch never
    assigns that variable but the pathomechanism block reads it. So the two
    flags cannot currently be combined.

P(M=1|Y=1) is the useful half of that calculation and it does not depend on
the prior at all, so this recovers it directly. It calls
`compute_pathomechanism_pathogenic_density_boundary` unmodified, with the same
benign anchor the pipeline uses (`--benign-method avg`, the mean of the benign
and synonymous weight vectors), against the bootstrap fits already written by
`run_pipeline.py`. The number reported is therefore ExCALIBR's own estimator,
not a reimplementation -- only the plumbing around it differs.

**What it means.** ExCALIBR decomposes the pathogenic-labeled sample's score
density as

    f_P(x) = P(M=1|Y=1) . f_D(x)  +  (1 - P(M=1|Y=1)) . f_N(x)

with f_N fixed to the benign/synonymous anchor. A pathogenic variant whose
disease mechanism this assay cannot see will score like a benign variant, so
it contributes to the f_N term. P(M=1|Y=1) is thus the fraction of the arm
that is "on-mechanism" for this assay, estimated from the score distribution
alone without reference to any per-variant call.

For a gene causing two disorders by two mechanisms, run per arm, this is a
direct quantification of which disorder the assay is competent to report --
the question SMC1A poses and which the paper lists as unevaluated.

Usage:
    python excalibr_pathomechanism.py \\
        --excalibr_repo /path/to/exCALIBR \\
        --dataset smc1a_DEE85_pathogenic.csv \\
        --fits DEE85_b20_bootstrap_fits.json.gz --label DEE85
"""
import argparse
import gzip
import json
import sys

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--excalibr_repo", required=True,
                    help="checkout of github.com/rosstewart/exCALIBR, for its "
                         "own boundary estimator")
    ap.add_argument("--dataset", required=True, help="the BasicScoreset CSV")
    ap.add_argument("--fits", required=True,
                    help="*_bootstrap_fits.json.gz from run_pipeline.py")
    ap.add_argument("--component", default="3c")
    ap.add_argument("--benign_method", default="avg",
                    choices=["avg", "benign", "synonymous"],
                    help="must match the pipeline run that produced --fits "
                         "[%(default)s]")
    ap.add_argument("--label", default="")
    a = ap.parse_args()

    sys.path.insert(0, a.excalibr_repo)
    from src.assay_calibration.fit_utils.point_ranges import (
        compute_pathomechanism_pathogenic_density_boundary)

    d = pd.read_csv(a.dataset)
    assign = d.sample_assignments.fillna("").astype(str)

    def sample(i):
        return d.score[assign.apply(lambda s: str(i) in s.split(","))].to_numpy()

    labeled, population = sample(0), sample(2)
    if len(labeled) == 0:
        sys.exit("no pathogenic-labeled variants in this dataset")

    fits = json.load(gzip.open(a.fits))
    fracs = []
    for _, per_component in fits.items():
        f = per_component[a.component]["fit"]
        params, w = f["component_params"], f["weights"]
        if a.benign_method == "avg":
            w_N = (np.array(w[1]) + np.array(w[3])) / 2
        else:
            w_N = np.array(w[1 if a.benign_method == "benign" else 3])
        _, frac = compute_pathomechanism_pathogenic_density_boundary(
            labeled, population, params, w[0], w_N)
        fracs.append(frac)

    fr = np.array(fracs, dtype=float)
    ok = int(np.isfinite(fr).sum())
    print(f"{a.label or a.dataset}: pathogenic n={len(labeled)}, "
          f"population n={len(population)}, {len(fr)} bootstrap fits")
    print(f"  P(M=1|Y=1) estimable in {ok}/{len(fr)} fits")
    if ok:
        med = float(np.nanmedian(fr))
        print(f"  median  {med:.3f}   5th-95th  "
              f"{np.nanpercentile(fr, 5):.3f}-{np.nanpercentile(fr, 95):.3f}")
        print(f"  implies ~{med * len(labeled):.0f} of {len(labeled)} arm "
              f"variants are on-mechanism for this assay")


if __name__ == "__main__":
    main()
