#!/usr/bin/env python3
"""Test whether a D139-suppressed deposit's individual can be matched, one to
one, onto a specific individual of the publication it was deposited from.

Collapse is only safe where the match is unique. Where a publication's
individuals are all keyed to a single `:unlabelled` id, the deposit labels are
FINER than the publication's own, and collapsing would merge distinct patients
-- the opposite error.
"""
import re
import pandas as pd

ROOT = "/Users/vb9/Documents/Claude/SMC1A/thesis_workflow/curation_attempt2"
o = pd.read_csv(f"{ROOT}/data/processed/smc1a_observations.tsv",
                sep="\t", dtype=str, low_memory=False)


def token(s):
    """Reduce a patient label to its comparable core: drop a leading PMID,
    the PUB:<paper> prefix, and the Pat/Pt/P/PT patient-word, keep the rest."""
    s = str(s)
    s = re.sub(r"^PUB:[A-Za-z]+\d{4}", "", s)      # PUB:Deardorff2007...
    s = re.sub(r":unlabelled$", "", s)
    s = re.sub(r"^\d{7,8}[-.]?", "", s)            # leading PMID
    s = s.replace(":", "").replace(" ", "").replace(".", "").replace("-", "")
    s = re.sub(r"^(Pat|Pt|PT|P)(?=\d|Fam)", "", s)  # patient word before a number
    return s.lower()


sup = o[o.disease_attribution_suppressed.notna()].copy()
sup["parent"] = sup.citations.str.replace(" ", "_")

rows = []
for parent, grp in sup.groupby("parent"):
    par_obs = o[o.source == parent]
    par_keys = sorted(set(par_obs.individual_key.dropna()))
    par_tok = {}
    for k in par_keys:
        par_tok.setdefault(token(k), []).append(k)
    # a publication whose individuals are ALL one unlabelled key cannot be
    # matched against finer deposit labels
    degenerate = len(par_keys) == 1 and par_keys[0].endswith(":unlabelled") \
        and token(par_keys[0]) == ""
    for _, d in grp.iterrows():
        dt = token(d.individual_label_raw)
        cand = par_tok.get(dt, [])
        rows.append(dict(
            parent=parent, deposit_label=d.individual_label_raw,
            deposit_key=d.individual_key, deposit_token=dt,
            n_parent_individuals=len(par_keys),
            matched=(len(cand) == 1 and not degenerate),
            matched_to=cand[0] if len(cand) == 1 else "",
            degenerate_parent=degenerate))

m = pd.DataFrame(rows)
print(f"suppressed deposit observations: {len(m)}")
print()
summary = m.groupby("parent").agg(
    deposits=("deposit_label", "size"),
    parent_individuals=("n_parent_individuals", "max"),
    matched=("matched", "sum"),
    degenerate=("degenerate_parent", "max"))
summary["unmatched"] = summary.deposits - summary.matched
print(summary.to_string())
print()
print(f"TOTAL matched 1:1  : {m.matched.sum()}")
print(f"TOTAL unmatched    : {(~m.matched).sum()}")
print()
um = m[~m.matched]
if len(um):
    print("=== unmatched deposits (collapse would be unsafe) ===")
    print(um[["parent", "deposit_label", "deposit_token",
              "n_parent_individuals", "degenerate_parent"]].to_string(index=False))
m.to_csv("/private/tmp/claude-503/-Users-vb9-Documents-Claude-SMC1A-thesis-workflow/"
         "3868e5c0-820a-4d63-9a8f-9da30bba6c9e/scratchpad/deposit_individual_matching.tsv",
         sep="\t", index=False)
