"""
Kayla Brand
Make eQTL predictions configuration table based on Synapse files
May 2025
"""

import os
# synapseclient is imported lazily inside the --sync path only. Building the
# predictions table from an already-downloaded tree needs no Synapse credentials,
# and requiring them made the generator unrunnable offline.
# import gzip
import pandas as pd
import numpy as np
# import sys
import re
import csv
from pathlib import Path
from collections import OrderedDict
import click

# Verify format of synced files matches expectations before running (dataset/cluster/predictions nesting)
# Can be run in eQTLEnv

# Dataset folders that are present in the Synapse sync but must never be benchmarked.
# "(Archived)Y2_versions" holds the superseded CharacterizationMcGinnis_* tree and
# "igvf7 (outdated prediction set)" is flagged as outdated by the submitters.
EXCLUDED_DATASETS = (
    'Legacy datasets',
    '(Archived)Y2_versions',
    'igvf7 (outdated prediction set)',
)

class DuplicateModelConflict(Exception):
    """Raised when it is unclear which of two files represents the E2G predictions for a particular cell type"""
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class DuplicateBiosampleConflict(Exception):
    """Raised when two dataset/cluster folders would collapse to the same biosample key.

    Biosample keys are dataset-prefixed precisely so this cannot happen silently; if it
    still does, the tree layout has changed and the table must not be written.
    """
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def configure_predictions(path, column_order):
    """
    Makes predictions configuration from a local copy of E2G Predictions project.
    REQUIRES ALL PREDICTION FILES TO BE NESTED INCLUDE CLUSTER NAME FOLDERS (none under a dataset folder)

    Args:
        path (str): Path to the E2G Predictions project directory.
        column_order (list): A list of model names (column names) in the desired order.

    Returns:
        OrderedDict: An OrderedDict where keys are cluster names and values are
                     OrderedDicts of model names and file paths.
    """
    predictions = OrderedDict()
    for dataset in sorted(os.listdir(path)):
        if dataset in EXCLUDED_DATASETS:
            print(f"Skipping excluded dataset {dataset}")
            continue
        dataset_path = os.path.join(path, dataset)
        # Top level also contains loose files (manifest.csv, SYNAPSE_METADATA_MANIFEST.tsv)
        if not os.path.isdir(dataset_path):
            continue
        # datasets level
        print(f"In dataset {dataset_path}...")
        for cluster_path in [str(p) for p in Path(dataset_path).glob('*/') if p.is_dir()]:
            cluster_models = OrderedDict()  # Use OrderedDict here
            print(f"In cluster {cluster_path}:")
            cluster_name = os.path.basename(cluster_path)
            # Biosample names are NOT unique across datasets (e.g. b_activated_d1 and
            # b_quiescent each appear under two datasets), so qualify every key with its
            # dataset. The eQTL pipeline keys on this value, and utils.smk builds a dict
            # from it, so a collision would silently drop one biosample's predictions.
            biosample_id = f"{dataset}_{cluster_name}"

            # Initialize all columns to empty strings
            for model in column_order:
                cluster_models[model] = ""  # Default to empty string

            for file in os.listdir(cluster_path):
                file_path = os.path.join(cluster_path, file)
                # files level
                if not is_candidate_prediction_file(file):
                    print(f"\tIgnoring file {file_path}")
                    continue
                
                # Extract model
                model = extract_model_name(os.path.basename(file_path))
                if model not in column_order:
                    print(f"\tWarning: Unknown model '{model}' found in file '{file_path}'. Skipping.")
                    continue
                # Add model to dictionary of cluster predictions
                if len(cluster_models[model]) > 0:  # Check if already populated
                    raise DuplicateModelConflict(f"""
                                                    Attempted to add {file_path} for {model} predictions on 
                                                    {cluster_name} when {cluster_models[model]} already exists.
                                                    """)
                else:
                    cluster_models[model] = file_path
                print(f"\t> Found model {model}")
            # Record models dictionary for this cluster
            if biosample_id in predictions:
                raise DuplicateBiosampleConflict(
                    f"Biosample key '{biosample_id}' was produced twice. "
                    f"Second occurrence: {cluster_path}. Dataset-prefixed keys are "
                    f"expected to be unique; refusing to silently overwrite."
                )
            predictions[biosample_id] = cluster_models

    return predictions


def generate_predictions_config_table(local_copy_dir, column_order, use_clusters = []):
    """
    Build a configuration table with paths to E2G prediction files to use
    Takes path to local directory synced with Synapse
    """
    # Get nested dictionary
    
    predictions = configure_predictions(local_copy_dir, column_order)

    # Ignore Legacy datasets
    if predictions.get("Legacy datasets", None) is not None:
        predictions.pop("Legacy datasets")
    
    # Print dictionary
    for cluster in predictions:
        print(f"{cluster}:")
        for model in predictions[cluster]:
            print(f"\t{model}: {predictions[cluster][model]}")

    # Use only the relevant cell clusters
    if len(use_clusters) > 0:
        relevant = dict((cluster, predictions[cluster]) 
                    for cluster in use_clusters
                    if cluster in predictions)
    else:
        relevant = predictions

    # Make the table
    predictions_df = pd.DataFrame.from_dict(relevant, orient='index')
    # Name the biosamples column
    predictions_df.index.name = "biosample"

    return predictions_df

def find_raw_baseline_predictions(base_directory):
    """Returns a dictionary of paths to raw scE2G Multiome predictions with cluster names as keys"""
    result = {}

    for root, dirs, files in os.walk(base_directory):
        for file in files:
            # print(f"root of predictions: {os.path.basename(root), root}" if file == 'scE2G_predictions.tsv.gz' else f"Ignoring file: {file}")
            if file == 'scE2G_predictions.tsv.gz' and os.path.basename(root) == 'multiome_powerlaw_v3':
                parent_directory = os.path.basename(os.path.dirname(root))
                file_path = os.path.join(root, file)
                result[parent_directory] = os.path.abspath(file_path)
    
    return result

def print_project_head(project):
    """
    Prints 2 layers of a Synapse Project
    """
    for folder_at_root in project.folders:
        print(f"Folder at root: {folder_at_root.name}")
        for file_in_root_folder in folder_at_root.files:
            print(f"File in {folder_at_root.name}: {file_in_root_folder.name}")
        for folder_in_folder in folder_at_root.folders:
            print(f"Folder in {folder_at_root.name}: {folder_in_folder.name}")
            for file_in_folder in folder_in_folder.files:
                print(f"File in {folder_in_folder.name}: {file_in_folder.name}")

# Method tokens as they appear at the END of a prediction filename, mapped to the
# method name used throughout the pipeline. Matching is anchored on this trailing
# token (see extract_model_name) rather than searched anywhere in the filename.
MODEL_FILENAME_SUFFIXES = {
    'scE2G_scATAC_powerlaw_v3': 'scE2G_ATAC',
    'scE2G_multiome_powerlaw_v3': 'scE2G_multiome',
    'scE2G_multiome_epcot': 'EPCOT',
    'Open4Gene': 'Open4Gene',
    'GEEES': 'GEEES',
    'SCARlink': 'SCARlink',
    'pgBoost': 'pgBoost',
    'Cicero': 'Cicero',
    'Signac': 'Signac',
    'SCENT': 'SCENT',
    'ArchR': 'ArchR',
    'EPCOT': 'EPCOT',
    'FigR': 'FigR',
}

# Longest first, so e.g. scE2G_multiome_powerlaw_v3 is tested before any shorter token
# that might also be a valid suffix of it.
_SUFFIXES_LONGEST_FIRST = sorted(MODEL_FILENAME_SUFFIXES, key=len, reverse=True)

def strip_prediction_extension(filename):
    """Return the filename stem with the prediction-file extension removed.

    Handles both `*.e2g.tsv.gz` (most models) and bare `*.tsv.gz` (Open4Gene ships
    without the `.e2g` infix). Returns None if this is not a gzipped TSV at all.
    """
    lowered = filename.lower()
    if lowered.endswith('.e2g.tsv.gz'):
        return filename[: -len('.e2g.tsv.gz')]
    if lowered.endswith('.tsv.gz'):
        return filename[: -len('.tsv.gz')]
    return None

def is_candidate_prediction_file(filename):
    """True if this file could be an unthresholded, standard-format prediction file.

    Only checks the extension. Derived sidecars (`*_thresholded`, `*_threshold0.177`,
    `*_percentile`, `*_absolute`, `*_boolean_score`, `*_element_list`, `*_gene_list`)
    are rejected downstream by extract_model_name, because anchoring the method token
    at the end of the stem means those names simply do not match any known model. That
    is safer than substring-matching markers like "list" against the whole filename,
    which would also match anything in the biosample or dataset name.
    """
    return strip_prediction_extension(filename) is not None

def extract_model_name(filename):
    """Returns name of the E2G model used given the file name.

    The method token is matched as a trailing, underscore-delimited component of the
    filename stem. Anchoring at the end matters: a plain substring search matches
    "scent" inside the biosample name "b_quiescent", which misclassified every file in
    those biosample folders as SCENT and raised DuplicateModelConflict.
    """
    stem = strip_prediction_extension(filename)
    if stem is None:
        return 'Unrecognized_Model'

    lowered_stem = stem.lower()
    for suffix in _SUFFIXES_LONGEST_FIRST:
        if lowered_stem.endswith('_' + suffix.lower()):
            return MODEL_FILENAME_SUFFIXES[suffix]

    return 'Unrecognized_Model'

def read_gtex_tissue_map(tsv_file_path):
    """Build {biosample_key: "tissue1,tissue2"} from a pairing table.

    Two layouts are accepted:

    1. Legacy wide  - columns `cluster`, `tissue`, where `tissue` is already a
       comma-separated list and the key is the bare cluster name.
    2. 2026 long    - columns `Dataset`, `Cluster_2026`, `eQTL_tissue`, one row per
       (biosample, tissue) pair. Keys are built as `{Dataset}_{Cluster_2026}` to match
       the dataset-prefixed biosample IDs emitted by configure_predictions, and the
       tissues for a biosample are collected into one comma-separated value.

    The long layout is what the pipeline ultimately needs anyway: utils.smk's
    process_biosample_tissue_maps splits GTExTissue on commas and emits one
    (biosample, tissue) pair per entry.
    """
    cluster_to_tissue = OrderedDict()

    with open(tsv_file_path, mode='r') as tsv_file:
        reader = csv.DictReader(tsv_file, delimiter='\t')
        fieldnames = reader.fieldnames or []

        is_long_format = 'Cluster_2026' in fieldnames and 'eQTL_tissue' in fieldnames

        for row in reader:
            if is_long_format:
                cluster = (row.get('Cluster_2026') or '').strip()
                dataset = (row.get('Dataset') or '').strip()
                tissue = (row.get('eQTL_tissue') or '').strip()
                if not cluster or not tissue:
                    continue  # trailing/blank rows
                key = f"{dataset}_{cluster}" if dataset else cluster
                existing = cluster_to_tissue.get(key, '')
                tissues = [t for t in existing.split(',') if t]
                if tissue not in tissues:
                    tissues.append(tissue)
                cluster_to_tissue[key] = ','.join(tissues)
            else:
                cluster = row['cluster']
                cluster_to_tissue[cluster] = row['tissue'] if 'tissue' in row else ''

    return cluster_to_tissue

DEFAULT_PREDICTIONS_DIR = "/scratch/users/kaybrand/E2G_Pillar_Project_predictions_Synapse_09_10_26"

# Baselines that are columns inside the scE2G multiome prediction file rather than
# separate files. config_methods.tsv distinguishes them purely by score_col
# (ABC.Score / distance / Kendall), so all three point at the same physical path.
# These names must match the `method` values in the methods config exactly, because
# utils.smk selects the predictions-table column by method name and a mismatch
# silently yields an all-NaN column rather than an error.
SCE2G_DERIVED_BASELINES = ('scABC', 'distanceToTSS', 'Kendall')

@click.command()
@click.option( # or generate this table automatically
    '-s', '--sync', 'sync_please',
    is_flag=True,  # This makes it a boolean flag
    default=False,
    help='If set, Sync from Synapse instead of generating a prediction file'
)
@click.option(
    '-d', '--predictions-dir', 'predictions_dir',
    type=click.Path(file_okay=False, path_type=Path),
    default=DEFAULT_PREDICTIONS_DIR,
    show_default=True,
    help='Local directory holding the downloaded E2G predictions (dataset/biosample/files).'
)
@click.option(
    '-o', '--output', 'output_table',
    type=click.Path(dir_okay=False, path_type=Path),
    default=os.path.join('config', 'predictions.tsv'),
    show_default=True,
    help='Where to write the generated predictions config table.'
)
@click.option(
    '-t', '--tissue-map', 'tissue_map_file',
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help='TSV with columns cluster,tissue used to fill GTExTissue. '
         'Defaults to config/eQTL_tissue_matches.tsv.'
)
def main(sync_please, predictions_dir, output_table, tissue_map_file):
    column_order = ['scE2G_multiome', 'scE2G_ATAC', 'pgBoost', 'ArchR', 'Signac', 'Cicero',
                    'SCENT', 'FigR', 'SCARlink', 'EPCOT', 'Open4Gene', 'GEEES']
    # Define paths at the module level for reliability
    THIS_SCRIPT_DIR = Path(__file__).resolve().parent
    CONFIG_SCRIPTS_DIR = THIS_SCRIPT_DIR.parent.parent / 'config'
    print(CONFIG_SCRIPTS_DIR)
    if tissue_map_file is None:
        tissue_map_file = os.path.join(CONFIG_SCRIPTS_DIR, "eQTL_tissue_matches.tsv")

    if sync_please:
        # Only the sync path needs Synapse credentials.
        import synapseclient
        from synapseutils import sync
        print(f" Synapse client version {synapseclient.__version__}")
        syn = synapseclient.Synapse()
        syn.login()
        sync.syncFromSynapse(syn, entity='syn53469845', path=str(predictions_dir), ifcollision="overwrite.local")
        return

    # Create predictions config
    predictions_df = generate_predictions_config_table(str(predictions_dir), column_order)

    # Add the baseline predictors. These are columns of the scE2G multiome prediction
    # file, so they resolve to that same path wherever it exists. Previously these were
    # looked up in a separate scE2G results tree keyed on old biosample names, which
    # matches nothing under the current naming and left all three columns empty.
    for baseline in SCE2G_DERIVED_BASELINES:
        predictions_df[baseline] = predictions_df['scE2G_multiome']

    # Add the gene expression tissue mapping column based on cell names
    gtex_tissue_map = read_gtex_tissue_map(tissue_map_file)
    predictions_df['GTExTissue'] = predictions_df.index.map(gtex_tissue_map)

    # Print for debug
    print(predictions_df.head())

    # Save the table
    os.makedirs(os.path.dirname(os.path.abspath(output_table)), exist_ok=True)
    predictions_df.to_csv(output_table, sep='\t', index=True)
    print(f"Wrote {len(predictions_df)} biosamples to {output_table}")

    

if __name__ ==  "__main__":
    main()

    # print("Testing main validation script")
    # file_path = Path("/scratch/users/kaybrand/Data/CharacterizationMcGinnis_Dataset4/imac_precursor/CharacterizationMcGinnis_Dataset4_imac_precursor_ArchR.e2g.tsv.gz")
    # validator_to_use = AbsCorrelationValidator(score_col_name="abs_Score", original_score_col="Score")
    # process_prediction_file(file_path, validator_to_use, False, False)

    # TESTED SCARlink: file_path = "/scratch/users/kaybrand/Data/CharacterizationMcGinnis_Dataset10/K562-CRISPRi/CharacterizationMcGinnis_Dataset10_K562_SCARlink.e2g.tsv.gz"
    # boolean_score.e2g.tsv.gz file passes all checks
    # TESTED corrFamily: file_path = "/scratch/users/kaybrand/Data/GM12878_10XMultiome/GM12878/GM12878_10XMultiome_GM12878_ArchR.e2g.tsv.gz" # Need to take abs value
    # TESTED corrFamily with finished file, passed all checks -> file_path = "/scratch/users/kaybrand/Data/CharacterizationMcGinnis_Dataset2/teloHAEC-24hr/CharacterizationMcGinnis_Dataset2_teloHAEC-24hr_ArchR_absolute.e2g.tsv.gz"
    # TESTED pgBoost Percentile file_path = "/scratch/users/kaybrand/Data/CharacterizationMcGinnis_Dataset9/D6/CharacterizationMcGinnis_Dataset9_D6_pgBoost_percentile.e2g.tsv.gz"
    # TESTED pgBoost rescue: file_path = "/scratch/users/kaybrand/Data/CharacterizationMcGinnis_Dataset6/S6_2D/CharacterizationMcGinnis_Dataset6_S6_2D_pgBoost.e2g.tsv.gz"
    # TESTED scE2G reformat to get scE2G Multiome core named columns - file_path = "/oak/stanford/groups/engreitz/Users/kaybrand/scE2G/results/CharacterizationMcGinnis_Dataset1/HUDEP2/multiome_powerlaw_v3/encode_e2g_predictions.tsv.gz"
    # TESTED - passed original - scE2G reformat on something that is already correctly reformatted (for ABC) - file_path = "/oak/stanford/groups/engreitz/Users/kaybrand/scE2G/results/CharacterizationMcGinnis_Dataset1/HUDEP2/multiome_powerlaw_v3/scE2G_predictions_for_eQTL_reformated.e2g.tsv.gz"
    # TESTED - passed reformatted - scE2G family reformat to get something else, like Kendall score - file_path = "/oak/stanford/groups/engreitz/Users/kaybrand/scE2G/results/CharacterizationMcGinnis_Dataset1/HUDEP2/multiome_powerlaw_v3/encode_e2g_predictions.tsv.gz"

