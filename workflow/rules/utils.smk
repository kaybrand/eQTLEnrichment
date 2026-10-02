import os
import pandas as pd

MAX_MEM_MB = 250 * 1000  # 250GB

def determine_mem_mb(wildcards, input, attempt, min_gb=8):
    # Memory resource calculator for snakemake rules
    input_size_mb = input.size_mb
    if ".gz" in str(input):
        input_size_mb *= 8  # assume gz compressesed the file <= 8x
    attempt_multiplier = 2 ** (attempt - 1)  # Double memory for each retry
    mem_to_use_mb = attempt_multiplier *  max(4 * input_size_mb, min_gb * 1000)
    return min(mem_to_use_mb, MAX_MEM_MB)

def add_biosamples_and_files_to_config_old():
    biosampleKeys = methods_config['sampleKey']
    methods = methods_config['method']

    samples = []
    files = []
    for i in range(len(methods)):
        if biosampleKeys.iloc[i]=="None":
            samples.append("")
            files.append("")
        else:
            key = pd.read_csv(biosampleKeys.iloc[i], sep='\t')
            sampleList = key["biosample"].tolist()
            fileList =  key["predictionFile"].tolist()
            fileDict = dict(zip(sampleList, fileList))
            samples.append(sampleList)
            files.append(fileDict)

    methods_config['biosamples'] = samples # add samples to config
    methods_config['predFiles'] = files 

def process_biosample_tissue_maps_old():
    biosampleKeys = methods_config['sampleKey']
    methods = methods_config['method']

    biosample_maps = []
    GTExTissue_maps = []
    for i in range(len(methods)):
        tissues_this = []
        biosamples_this = []
        if not biosampleKeys.iloc[i]=="None":
            key = pd.read_csv(biosampleKeys.iloc[i], sep='\t')
            key = key.dropna(subset='GTExTissue')
            # make lists, accounting for replicate matches
            for index, row in key.iterrows():
                tissues_this_row = [x.strip() for x in row["GTExTissue"].split(',')]
                biosamples_this_row = [row['biosample'].strip() for x in tissues_this_row]
                tissues_this.extend(tissues_this_row)
                biosamples_this.extend(biosamples_this_row)

        biosample_maps.append(biosamples_this)
        GTExTissue_maps.append(tissues_this)

    methods_config['GTExTissue_map'] = GTExTissue_maps
    methods_config['biosample_map'] = biosample_maps

# def validate_prediction_table(config):
#     try:
#         result = subprocess.run(
#             ["conda", "run", "-n", "eQTLEnv", "python",
#             'workflow/scripts/preprocessing/validate_prediction_formats.py',
#             config['predictionsTable'],
#             config['methodsTable']],
#             check=True,
#             capture_output=True,
#             text=True,  # Crucial for text output
#         )

#         print("Prediction Validation Alerts:")
#         print(result.stdout)
#         if len(result.stdout) < 2:
#             print("All prediction files look good!")

#         return os.path.join(os.path.dirname(config['predictionsTable']), "validated_" + os.path.basename(config['predictionsTable']))

#     except subprocess.CalledProcessError as e:
#         print(f"Script failed with error code {e.returncode}:")
#         print(e.stderr)
#         sys.exit(1)  # Exit with a non-zero code on failure


def add_biosamples_and_files_to_config(methods_config, config):
    key = pd.read_csv(config["predictionsTable"], sep="\t").dropna(subset=["biosample"])
    methods_config = methods_config[methods_config["method"].isin(config["methods"])] # filter to relevant methods
    
    samples = []
    files = []
    
    # iterate over the rows of the df
    for i in range(len(methods_config)):
        this_method = methods_config.iloc[i]["method"]
        
        # Select the relevant columns for this method
        this_key = key[["biosample", this_method]].dropna()
        
        sampleList = this_key["biosample"].tolist()
        fileList = this_key[this_method].tolist()

        fileDict = dict(zip(sampleList, fileList))
        samples.append(sampleList)
        files.append(fileDict)
    
    # Add new columns to `methods_config` using Pandas
    methods_config.loc[:, 'biosamples'] = samples # add samples to config
    methods_config.loc[:, 'predFiles'] = files 
    print((methods_config.head()))
    
    return methods_config

def process_biosample_tissue_maps(methods_config, config):
    print("Called process_biosample_tissue_maps")
    print(f'config["predictionsTable"] is: {config["predictionsTable"]}')
    key = pd.read_csv(config["predictionsTable"], sep="\t").dropna(subset=["biosample", "GTExTissue"])
    methods_config = methods_config[methods_config["method"].isin(config["methods"])]

    biosample_maps = []
    GTExTissue_maps = []
    # iterate over the rows of the df
    for i in range(len(methods_config)):
        tissues_this = []
        biosamples_this = []
        this_method = methods_config.iloc[i]["method"]
        
        # select the relevant columns for this method
        this_key = key[["biosample", this_method, "GTExTissue"]].dropna()
        for index, row in this_key.iterrows():
                tissues_this_row = [x.strip() for x in row["GTExTissue"].split(',')]
                biosamples_this_row = [row['biosample'].strip() for x in tissues_this_row]
                tissues_this.extend(tissues_this_row)
                biosamples_this.extend(biosamples_this_row)

        biosample_maps.append(biosamples_this)
        GTExTissue_maps.append(tissues_this)

    methods_config.loc[:, 'GTExTissue_map'] = GTExTissue_maps
    methods_config.loc[:, 'biosample_map'] = biosample_maps
    print((methods_config.head()))

    return methods_config

def flatten(list_of_lists):
    return [x for xs in list_of_lists for x in xs]

# These two always share a candidate element-gene universe even though they
# read different prediction files.
ALWAYS_SHARED_UNIVERSE = [("scE2G_multiome", "scE2G_ATAC")]

def determine_candidate_universe_groups(methods_config, config):
    """Group config["methods"] by candidate element-gene universe.

    Two methods share a universe when the predictions table points them at the
    same file for every biosample they both have (methods that read one file,
    e.g. scABC and distanceToTSS reading the scE2G multiome file, are certain to
    share it). Different files mean different universes, except for the pairs in
    ALWAYS_SHARED_UNIVERSE. Sharing is transitive.

    Returns: dict method -> group_label, where group_label is "-".join(sorted
    member methods).
    """
    methods = [m for m in config["methods"] if m in methods_config.index]
    parent = {m: m for m in methods}

    def find(m):
        while parent[m] != m:
            parent[m] = parent[parent[m]]
            m = parent[m]
        return m

    def join(a, b):
        parent[find(a)] = find(b)

    pred_files = {m: methods_config.loc[m, "predFiles"] for m in methods}
    for i, a in enumerate(methods):
        for b in methods[i + 1:]:
            shared = set(pred_files[a]) & set(pred_files[b])
            if shared and all(pred_files[a][bio] == pred_files[b][bio] for bio in shared):
                join(a, b)
    for a, b in ALWAYS_SHARED_UNIVERSE:
        if a in parent and b in parent:
            join(a, b)

    groups = {}
    for m in methods:
        groups.setdefault(find(m), []).append(m)
    return {m: "-".join(sorted(members)) for members in groups.values() for m in members}
