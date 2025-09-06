"""
Kayla Brand
Make eQTL predictions configuration table based on Synapse files
May 2025
"""

import os
import synapseclient
print(f" Synapse client version {synapseclient.__version__}")
from synapseclient.models import Folder, Project
from synapseutils import sync
# import gzip
import pandas as pd
import numpy as np
# import sys
import re
import csv
from pathlib import Path
from collections import OrderedDict

class DuplicateModelConflict(Exception):
    """Raised when it is unclear which of two files represents the E2G predictions for a particular cell type"""
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
    for dataset in os.listdir(path):
        if dataset == 'Legacy datasets':
            continue
        dataset_path = os.path.join(path, dataset)
        # datasets level
        print(f"In dataset {dataset_path}...")
        for cluster_path in [str(p) for p in Path(dataset_path).glob('*/') if p.is_dir()]:
            cluster_models = OrderedDict()  # Use OrderedDict here
            print(f"In cluster {cluster_path}:")
            cluster_name = os.path.basename(cluster_path)

            # Initialize all columns to empty strings
            for model in column_order:
                cluster_models[model] = ""  # Default to empty string

            for file in os.listdir(cluster_path):
                file_path = os.path.join(cluster_path, file)
                # files level
                if ("threshold" in file) or ("list" in file) or ('e2g.tsv.gz' not in file) or ("percentile" in file) or ("absolute" in file) or ('boolean' in file):
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
            predictions[cluster_name] = cluster_models

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
            # print(f"root of predictions: {os.path.basename(root), root}" if file == 'encode_e2g_predictions.tsv.gz' else f"Ignoring file: {file}")
            if file == 'encode_e2g_predictions.tsv.gz' and os.path.basename(root) == 'multiome_powerlaw_v3':
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

def extract_model_name(filename):
    """Returns name of the E2G model used given the file name"""
    filename = filename.lower()
    if "scent" in filename:
        return "SCENT"
    elif "signac" in filename:
        return "Signac"
    elif "cicero" in filename:
        return "Cicero"
    elif "archr" in filename:
        return "ArchR"
    elif "epcot" in filename:
        return "EPCOT"
    elif "figr" in filename:
        return "FigR"
    elif "scarlink" in filename:
        return "SCARlink"
    elif "pgboost" in filename:
        return "pgBoost"
    elif "sce2g" in filename:
        scE2G_pattern = r'_(scE2G_[^_]+_[^_]+_v\d)\.e2g\.tsv'
        scE2G_match = re.search(scE2G_pattern, filename)
        if "scatac" in filename:
            return "scE2G_ATAC"
        elif "multiome" in filename:
            return "scE2G_multiome"
        elif scE2G_match:
            return scE2G_match.group(1)
        else:
            return 'Unrecognized_Model'
    else:
        return 'Unrecognized_Model'

def read_gtex_tissue_map(tsv_file_path):
    cluster_to_tissue = {}

    with open(tsv_file_path, mode='r') as tsv_file:
        reader = csv.DictReader(tsv_file, delimiter='\t')
        
        for row in reader:
            cluster = row['cluster']
            tissue = row['tissue'] if 'tissue' in row else ''
            cluster_to_tissue[cluster] = tissue
    
    return cluster_to_tissue


def main():
    # Log into Synapse
    syn = synapseclient.Synapse()
    syn.login()

    # Sync Synapse E2G predictions Project locally
    sync_please = False
    directory_to_sync_project_to = os.path.join("/", "scratch", "users", "kaybrand", "Data")
    scE2G_results_dir = os.path.join("/", "oak", "stanford", "groups", "engreitz", "Users", "kaybrand", "scE2G", "results")
    column_order = ['scE2G_multiome', 'scE2G_ATAC', 'pgBoost', 'ArchR', 'Signac', 'Cicero', 'SCENT', 'FigR', 'SCARlink', 'EPCOT']
    # Define paths at the module level for reliability
    THIS_SCRIPT_DIR = Path(__file__).resolve().parent
    CONFIG_SCRIPTS_DIR = THIS_SCRIPT_DIR.parent.parent.parent / 'config'
    print(CONFIG_SCRIPTS_DIR)
    tissue_mapping_config = os.path.join(CONFIG_SCRIPTS_DIR, "eQTL_tissue_matches.tsv")

    if sync_please:
        # project = synapseclient.models.Project(name="E2G predictions", id='syn53469845')
        # Set the `if_collision` to `keep.local` so that we don't overwrite any files
        # project.sync_from_synapse(path=DIRECTORY_TO_SYNC_PROJECT_TO, if_collision="keep.local")
        sync.syncFromSynapse(syn, entity='syn53469845', path=directory_to_sync_project_to, ifcollision="overwrite.local")
        # print_project_head(project)

    # Create predictions config
    predictions_df = generate_predictions_config_table(directory_to_sync_project_to, column_order)

    # Add the baseline predictors
    raw_scE2G_prediction_paths = find_raw_baseline_predictions(scE2G_results_dir)
    print(raw_scE2G_prediction_paths)
    for baseline in ['scABC', 'Kendall', 'ABC_distanceToTSS']:
        predictions_df[baseline] = predictions_df.index.map(raw_scE2G_prediction_paths)

    # Add the gene expression tissue mapping column based on cell names
    gtex_tissue_map = read_gtex_tissue_map(tissue_mapping_config)
    predictions_df['GTExTissue'] = predictions_df.index.map(gtex_tissue_map)

    # Print for debug
    print(predictions_df.head())

    # Save the table
    predictions_df.to_csv(os.path.join('config', 'predictions.tsv'), sep='\t', index=True)

    

if __name__ ==  "__main__":
    main()

