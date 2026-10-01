# Checks run while the DAG is built, before any job starts, so config mistakes
# fail in seconds with a list of every problem instead of an hour later in a
# rule deep in the DAG.

import difflib
import gzip
import os
import re

import pandas as pd
from snakemake.exceptions import WorkflowError

REQUIRED_METHOD_COLUMNS = ["method", "boolean", "inverse_predictor", "threshold", "score_col", "geneUniverse"]
# nonneg_only is optional; a blank or absent cell means FALSE.
BOOLEAN_METHOD_COLUMNS = ["boolean", "inverse_predictor", "nonneg_only"]
# pandas turns true/false in any case into booleans. Any other text (yes, 1, T)
# reaches process_predictions.R as a string it reads as FALSE.
ACCEPTED_BOOLEAN_TEXT = {"TRUE", "FALSE", "True", "False"}
REQUIRED_PREDICTION_FILE_COLUMNS = ["ElementChr", "ElementStart", "ElementEnd", "GeneSymbol"]
MAX_FILES_LISTED = 5


def _did_you_mean(name, options):
    by_lower = {str(o).lower(): str(o) for o in options}
    close = difflib.get_close_matches(str(name).lower(), list(by_lower), n=1, cutoff=0.6)
    return f" Did you mean '{by_lower[close[0]]}'?" if close else ""


def _is_blank(value):
    return value is None or (isinstance(value, float) and pd.isna(value)) or str(value).strip() in ("", "None", "nan")


def _as_bool(value):
    if isinstance(value, bool):
        return value
    return str(value).strip() in ("TRUE", "True")


def _header_columns(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                return line.rstrip("\n").split("\t")
    return []


def _list_some(items):
    shown = ", ".join(items[:MAX_FILES_LISTED])
    extra = len(items) - MAX_FILES_LISTED
    return shown + (f" (and {extra} more)" if extra > 0 else "")


def _check_methods_table(methods_config, config, errors):
    mt = config["methodsTable"]
    methods = list(config.get("methods") or [])

    if not methods:
        errors.append("The main config's 'methods' list is empty. List the methods to benchmark.")
    repeated = sorted({m for m in methods if methods.count(m) > 1})
    if repeated:
        errors.append(f"The main config's 'methods' list repeats {repeated}. List each method once.")

    missing_cols = [c for c in REQUIRED_METHOD_COLUMNS if c not in methods_config.columns]
    if missing_cols:
        errors.append(f"Methods table {mt} is missing required column(s) {missing_cols}. "
                      f"Columns found: {list(methods_config.columns)}.")
        return False

    duplicated = sorted(set(methods_config.index[methods_config.index.duplicated()]))
    if duplicated:
        errors.append(f"Methods table {mt} has more than one row for {duplicated}. Keep one row per method.")

    for m in methods:
        if m not in methods_config.index:
            errors.append(f"Method '{m}' is listed in the main config but has no row in {mt}."
                          + _did_you_mean(m, methods_config.index))
            continue
        if m in duplicated:
            continue
        row = methods_config.loc[m]

        # A blank threshold is allowed: the method is then left out of every
        # thresholded output and appears only in the threshold-free ones.
        threshold = None
        if not _is_blank(row["threshold"]):
            try:
                threshold = float(row["threshold"])
            except (TypeError, ValueError):
                errors.append(f"{m}: threshold '{row['threshold']}' in {mt} is not a number. "
                              f"Leave it blank if this method has no threshold.")

        flags = {}
        for col in BOOLEAN_METHOD_COLUMNS:
            if col not in methods_config.columns:
                continue
            value = row[col]
            if _is_blank(value):
                if col == "nonneg_only":
                    flags[col] = False
                else:
                    errors.append(f"{m}: {col} is blank in {mt}. Write TRUE or FALSE.")
                continue
            if not isinstance(value, bool) and str(value).strip() not in ACCEPTED_BOOLEAN_TEXT:
                errors.append(f"{m}: {col} is '{value}' in {mt}. Write TRUE or FALSE; other values "
                              f"(e.g. 'yes', '1', 'T') are read as FALSE by the R scripts.")
                continue
            flags[col] = _as_bool(value)

        if threshold is not None and flags.get("boolean") and threshold != 1:
            errors.append(f"{m}: boolean is TRUE, so scores are 0/1 and the threshold should be 1, not {threshold:g}.")

        if _is_blank(row["score_col"]):
            errors.append(f"{m}: score_col is blank in {mt}. Name the prediction-file column holding the score.")

        gene_universe = row["geneUniverse"]
        if _is_blank(gene_universe):
            errors.append(f"{m}: geneUniverse is blank in {mt}.")
        elif not os.path.exists(gene_universe):
            errors.append(f"{m}: geneUniverse file '{gene_universe}' does not exist "
                          f"(paths are relative to {os.getcwd()}).")

        if "color" in methods_config.columns and not _is_blank(row["color"]) \
                and not re.fullmatch(r"#[0-9A-Fa-f]{6}", str(row["color"]).strip()):
            errors.append(f"{m}: color '{row['color']}' is not a hex code like #792374. Leave it blank to pick one automatically.")
    return True


def _check_predictions_table(methods_config, config, errors):
    pt = config["predictionsTable"]
    if not os.path.exists(pt):
        errors.append(f"Predictions table {pt} does not exist (paths are relative to {os.getcwd()}).")
        return
    key = pd.read_csv(pt, sep="\t", dtype=str)

    missing_cols = [c for c in ("biosample", "GTExTissue") if c not in key.columns]
    if missing_cols:
        errors.append(f"Predictions table {pt} is missing column(s) {missing_cols}.")
        return

    names = key["biosample"].dropna()
    padded = sorted(n for n in names if n != n.strip())
    if padded:
        errors.append(f"Predictions table {pt}: biosample name(s) {padded} have leading or trailing spaces.")
    repeated = sorted(names[names.duplicated()].unique())
    if repeated:
        errors.append(f"Predictions table {pt}: biosample name(s) {repeated} appear more than once. "
                      f"Names must be unique (e.g. prefix them with the dataset).")

    header_cache = {}
    for m in config.get("methods") or []:
        if m not in methods_config.index or isinstance(methods_config.loc[m], pd.DataFrame):
            continue
        if m not in key.columns:
            errors.append(f"Method '{m}' has no column in predictions table {pt}. The column name must "
                          f"match the method name exactly." + _did_you_mean(m, key.columns))
            continue

        with_file = key.dropna(subset=[m])
        if with_file.empty:
            errors.append(f"{m}: no biosample in {pt} has a prediction file for this method.")
            continue
        if not with_file["GTExTissue"].fillna("").str.strip().ne("").any():
            errors.append(f"{m}: none of the biosamples with a {m} file has a GTExTissue pairing in {pt}, "
                          f"so nothing would be benchmarked for this method.")

        score_col = str(methods_config.loc[m, "score_col"])
        needed = REQUIRED_PREDICTION_FILE_COLUMNS + [score_col]
        missing_files, bad_headers = [], {}
        for path in with_file[m]:
            if not os.path.exists(path):
                missing_files.append(path)
                continue
            if path not in header_cache:
                try:
                    header_cache[path] = _header_columns(path)
                except (OSError, EOFError, UnicodeDecodeError) as e:
                    header_cache[path] = e
            header = header_cache[path]
            if isinstance(header, Exception):
                bad_headers.setdefault(f"could not be read ({header})", []).append(path)
                continue
            missing = [c for c in needed if c not in header]
            if missing:
                hint = _did_you_mean(score_col, header) if score_col in missing else ""
                bad_headers.setdefault(f"are missing column(s) {missing}.{hint}", []).append(path)
        if missing_files:
            errors.append(f"{m}: {len(missing_files)} prediction file(s) listed in {pt} do not exist: "
                          f"{_list_some(missing_files)}")
        for problem, paths in bad_headers.items():
            errors.append(f"{m}: {len(paths)} prediction file(s) {problem} Expected columns "
                          f"{needed} (score_col '{score_col}' comes from {config['methodsTable']}). "
                          f"Files: {_list_some(paths)}")


def check_config_inputs(methods_config, config):
    errors = []
    if _check_methods_table(methods_config, config, errors):
        _check_predictions_table(methods_config, config, errors)
    if errors:
        listed = "\n".join(f"  {i}. {e}" for i, e in enumerate(errors, 1))
        raise WorkflowError(f"Config check found {len(errors)} problem(s); no jobs were run.\n{listed}")
