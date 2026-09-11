"""
Kayla Brand / 2026 round
Derive the two benchmark-variant predictions tables from a validated table.

The eQTL benchmark is run in two variants this round:

  all_models      Only biosamples for which EVERY benchmarked model has predictions,
                  so each model is scored on an identical biosample-tissue set.
                  Membership is limited by the least-covered model (Open4Gene).

  all_biosamples  Only models that have predictions across the FULL biosample set,
                  so each model is scored on every pairing. Membership is limited by
                  the least-covered biosample.

Both are derived from the validated predictions table rather than hand-assembled,
because a cell pointing at the wrong file silently removes a method from the
benchmark (utils.smk selects the column by method name and drops NaN).

Re-run this after every re-validation; do not edit the outputs by hand.
"""

import argparse
import os
import sys

import pandas as pd


def summarize(df, label):
    tissues = sorted({t.strip() for g in df["GTExTissue"].dropna()
                      for t in str(g).split(",") if t.strip()})
    pairs = sum(len([t for t in str(g).split(",") if t.strip()])
                for g in df["GTExTissue"].dropna())
    print(f"  {label}: {len(df)} biosamples, {pairs} biosample-tissue pairs, "
          f"{len(tissues)} tissues {tissues}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--validated-table", required=True,
                        help="Validated predictions table (<name>_validated.tsv)")
    parser.add_argument("-o", "--out-dir", default="config",
                        help="Directory to write the two variant tables into")
    parser.add_argument("--all-models", nargs="+", required=True,
                        help="Models that must ALL be present for the all_models variant")
    parser.add_argument("--candidate-models", nargs="+", required=True,
                        help="Models considered for the all_biosamples variant")
    parser.add_argument("--exclude-biosample", action="append", default=[],
                        help="Biosample to drop from the all_biosamples variant "
                             "(repeatable). Use for biosamples whose missing models "
                             "would otherwise disqualify those models entirely.")
    parser.add_argument("--prefix", default="predictions_2026",
                        help="Output filename prefix")
    args = parser.parse_args()

    v = pd.read_csv(args.validated_table, sep="\t", index_col="biosample")

    # Only biosamples with at least one eQTL tissue pairing are benchmarkable; the
    # rest would consume full compute and never appear in a tissue map.
    matched = v[v["GTExTissue"].notna() & (v["GTExTissue"].astype(str).str.strip() != "")]
    print(f"Matched biosamples in {args.validated_table}: {len(matched)}")
    if matched.empty:
        sys.exit("ERROR: no biosample has a GTExTissue assignment.")

    missing = [m for m in set(args.all_models) | set(args.candidate_models)
               if m not in matched.columns]
    if missing:
        sys.exit(f"ERROR: models absent from the table: {sorted(missing)}")

    os.makedirs(args.out_dir, exist_ok=True)

    # --- all_models: intersect on biosamples -------------------------------------
    a = matched[matched[args.all_models].notna().all(axis=1)]
    a_out = os.path.join(args.out_dir, f"{args.prefix}_all_models.tsv")
    a[["GTExTissue"] + list(args.all_models)].to_csv(a_out, sep="\t")
    print("\nall_models variant")
    summarize(a, "kept")
    print(f"  models ({len(args.all_models)}): {list(args.all_models)}")
    print(f"  -> {a_out}")

    # --- all_biosamples: intersect on models -------------------------------------
    b = matched.drop(index=[i for i in args.exclude_biosample if i in matched.index])
    for excluded in args.exclude_biosample:
        if excluded not in matched.index:
            print(f"  WARNING: --exclude-biosample {excluded} is not in the table")

    full = [m for m in args.candidate_models if b[m].notna().all()]
    dropped = [m for m in args.candidate_models if m not in full]
    b_out = os.path.join(args.out_dir, f"{args.prefix}_all_biosamples.tsv")
    b[["GTExTissue"] + full].to_csv(b_out, sep="\t")
    print("\nall_biosamples variant")
    summarize(b, "kept")
    print(f"  models with full coverage ({len(full)}): {full}")
    if dropped:
        for m in dropped:
            blockers = b.index[b[m].isna()].tolist()
            print(f"  DROPPED {m}: missing for {blockers}")
    print(f"  -> {b_out}")


if __name__ == "__main__":
    main()
