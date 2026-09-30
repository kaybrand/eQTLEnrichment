import os
import pandas as pd
import subprocess

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

# Methods that always score the same underlying scE2G candidate element-gene
# list, differing only in which score column they extract from it -- confirmed
# empirically, all 10 biosamples, md5-identical (see analysis/
# enrichment_denominator_fix/compute_universe_fingerprint.R, 2026-09-30). ABC is
# included on the same domain-knowledge basis even when a run doesn't include it.
# Skipping computation for these is Tier 1 of determine_candidate_universe_groups.
KNOWN_SHARED_UNIVERSE_METHODS = {"ABC", "distanceToTSS", "scE2G_multiome", "scE2G_ATAC", "scABC"}

def _raw_row_count(predFile):
    cmd = f"zcat '{predFile}' | awk '!/^#/' | tail -n +2 | wc -l"
    return int(subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout.strip())

def _raw_universe_hash(predFile):
    # Column positions, not names: avoids a csvtk dependency at Snakefile parse
    # time (before any rule's `conda:` environment is active). ElementChr,
    # ElementStart, ElementName, GeneEnsemblID -- ElementName already encodes
    # chr:start-end, so this is a sufficient identity key without ElementEnd.
    header_cmd = f"zcat '{predFile}' | awk '!/^#/' | head -1"
    header = subprocess.run(header_cmd, shell=True, check=True, capture_output=True, text=True).stdout.strip().split("\t")
    wanted = ["ElementChr", "ElementStart", "ElementName", "GeneEnsemblID"]
    idx = [str(header.index(col) + 1) for col in wanted]  # cut is 1-indexed
    cmd = (f"zcat '{predFile}' | awk '!/^#/' | tail -n +2 | "
           f"cut -f{','.join(idx)} | sort | md5sum | cut -d' ' -f1")
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout.strip()

def determine_candidate_universe_groups(methods_config, config):
    """Group config["methods"] by candidate element-gene universe, cheaply.

    Four tiers, cheapest first, each only run when the previous one couldn't
    resolve the answer:
      Tier 0 (free)  -- identical raw predFile path for a shared biosample.
      Tier 1 (free)  -- method in KNOWN_SHARED_UNIVERSE_METHODS.
      Tier 2 (cheap) -- raw row count (no sort); a necessary condition.
      Tier 3 (only on a Tier-2 collision) -- md5 of the sorted, unique
        (ElementChr, ElementStart, ElementName, GeneEnsemblID) tuples from the
        RAW standard-format file -- before this pipeline's own gene-universe
        filter, sign handling, or renaming, so the answer reflects what the
        submitter tested, not an artifact of our processing.

    Verifies the result is identical across every biosample a method is used
    on in this run; a method whose biosamples disagree is given its own solo
    group (its own plot) rather than merged into a group that doesn't hold
    everywhere.

    Caches every (method, biosample) fingerprint under
    {outDir}/intermediate/universeFingerprints/ -- a fingerprint already on
    disk is read, never recomputed, across pipeline reruns of the same outDir.

    Returns: dict method -> group_label, where group_label is "-".join(sorted
    member methods) -- readable, and stable regardless of which member happens
    to be the representative used to read the shared enrichment table.
    """
    cache_dir = os.path.join(config["outDir"], "intermediate", "universeFingerprints")
    os.makedirs(cache_dir, exist_ok=True)

    methods = [m for m in config["methods"] if m in methods_config.index]
    biosamples_by_method = {m: methods_config.loc[m, "biosamples"] for m in methods}
    predfiles_by_method = {m: methods_config.loc[m, "predFiles"] for m in methods}
    all_biosamples = sorted({b for bios in biosamples_by_method.values() for b in bios})

    def cache_file(method, biosample):
        return os.path.join(cache_dir, f"{method}__{biosample}.tsv")

    def read_cache(method, biosample):
        f = cache_file(method, biosample)
        if not os.path.exists(f):
            return None
        with open(f) as fh:
            tier, row_count, h = fh.read().strip().split("\t")
        return tier, (int(row_count) if row_count != "NA" else None), (h if h != "NA" else None)

    def write_cache(method, biosample, tier, row_count, h):
        with open(cache_file(method, biosample), "w") as fh:
            fh.write(f"{tier}\t{row_count if row_count is not None else 'NA'}\t{h if h is not None else 'NA'}\n")

    method_tag_per_biosample = {}  # method -> {biosample: tag}

    for biosample in all_biosamples:
        methods_here = [m for m in methods if biosample in biosamples_by_method[m]]
        known_here = [m for m in methods_here if m in KNOWN_SHARED_UNIVERSE_METHODS]
        for m in known_here:
            method_tag_per_biosample.setdefault(m, {})[biosample] = "KNOWN:scE2G_candidate_universe"

        remaining = [m for m in methods_here if m not in KNOWN_SHARED_UNIVERSE_METHODS]
        if not remaining:
            continue

        pred_files = {m: predfiles_by_method[m][biosample] for m in remaining}
        row_counts, hashes = {}, {}
        for m in remaining:
            cached = read_cache(m, biosample)
            if cached is not None:
                _, rc, h = cached
            else:
                rc = _raw_row_count(pred_files[m])
                h = None
                write_cache(m, biosample, "2_rowcount", rc, None)
            row_counts[m], hashes[m] = rc, h

        by_count = {}
        for m in remaining:
            by_count.setdefault(row_counts[m], []).append(m)

        for rc, ms in by_count.items():
            if len(ms) == 1:
                m = ms[0]
                method_tag_per_biosample.setdefault(m, {})[biosample] = f"SOLO:{m}:{rc}"
                continue
            # Row-count collision -- escalate to Tier 3, once per distinct
            # predFile path (methods reading the same file need one hash, not
            # one each), and persist the result so it is never recomputed.
            path_hash = {}
            for m in ms:
                if hashes[m] is not None:
                    path_hash[pred_files[m]] = hashes[m]
            for m in ms:
                p = pred_files[m]
                if p not in path_hash:
                    path_hash[p] = _raw_universe_hash(p)
                write_cache(m, biosample, "3_hash", rc, path_hash[p])
                method_tag_per_biosample.setdefault(m, {})[biosample] = f"HASH:{path_hash[p]}"

    # Reconcile per-method: a tag (especially a HASH or SOLO one) is itself
    # biosample-specific -- two different biosamples' candidate sets never hash
    # equal, even for the same method, since the coordinates differ. So the tag
    # STRING is never the same across a method's own biosamples, by
    # construction, and comparing it directly (as an earlier version of this
    # function did) would call every non-KNOWN method "inconsistent". What must
    # be stable is WHICH OTHER METHODS a method groups with, per biosample --
    # its partner set -- not the tag value itself.
    partners_per_biosample = {}  # biosample -> {tag: [methods]}
    for m, per_bio in method_tag_per_biosample.items():
        for biosample, tag in per_bio.items():
            partners_per_biosample.setdefault(biosample, {}).setdefault(tag, []).append(m)

    method_partner_sets = {}  # method -> set of frozenset(partners), across its biosamples
    for m, per_bio in method_tag_per_biosample.items():
        for biosample, tag in per_bio.items():
            partners = frozenset(partners_per_biosample[biosample][tag])
            method_partner_sets.setdefault(m, set()).add(partners)

    method_tag = {}
    for m in methods:
        partner_sets = method_partner_sets.get(m, set())
        if len(partner_sets) == 1:
            method_tag[m] = "GROUP:" + "-".join(sorted(partner_sets.pop()))
        else:
            method_tag[m] = f"SOLO:{m}:inconsistent"
            print(f"WARNING: candidate universe for method '{m}' does not group with the "
                  f"same other methods on every biosample it appears on "
                  f"({sorted(sorted(p) for p in partner_sets)}) -- giving it its own "
                  f"candidate-universe plot instead of merging it into a group that "
                  f"does not hold on every biosample.")

    groups = {}
    for m, tag in method_tag.items():
        groups.setdefault(tag, []).append(m)

    method_to_label = {}
    for members in groups.values():
        label = "-".join(sorted(members))
        for m in members:
            method_to_label[m] = label
    return method_to_label