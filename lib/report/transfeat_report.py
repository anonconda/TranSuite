import os
import sys
import time
import numpy as np
from collections import defaultdict, namedtuple
from itertools import chain, combinations
from lib.parsing.gtf_object_tools import create_gtf_object


def get_per(n, tot):

    try:
        return n, round((n/tot)*100, 1)
    except ZeroDivisionError:
        return n, 0.0


def get_percentiles(values, perc_range=None):

    if not values:
        return []

    if perc_range is None:
        perc_range = [0, 1, 5, 10, 15, 20, 25, 33, 50, 66, 75, 80, 85, 90, 95, 99, 100]

    assert all(0 <= p <= 100 for p in perc_range)

    res = []
    for quantile in perc_range:
        quant_val = np.percentile(values, quantile)
        data = (quantile, quant_val)
        res.append(data)

    return res


def get_transfeat_data(transfeat_table):

    print(time.asctime(), f"Uploading information from TransFeat table: {transfeat_table}")

    trans_by_feature_dt = defaultdict(set)
    gene_coding_pot_dt, trans_coding_pot_dt, trans_by_gene_coding_pot_dt, gene_trans_dt = [defaultdict(set) for _ in range(4)]
    with open(transfeat_table) as fh:
        # Header:
        # [0] Gene_ID, [1] Transcript_ID, [2] Protein_coding, [3] Coding_features, [4] Alternative_ORF,
        # [5] NMD_features, [6] Transcript_coordinates, [7] CDS_coordinates, [8] Translation
        next(fh)
        for row in fh:
            # Avoid empty rows
            if not row.replace(",", ""):
                continue

            row_list = row.strip("\n").split(",")

            gene_id, trans_id, coding_pot = row_list[0], row_list[1], row_list[2].upper()  # Use upper to sanitize input
            code_bool, code_feat, nmd_feat, alt_orf = row_list[2], row_list[3], row_list[5], row_list[4]
            pep_len = len(row_list[8].replace("-", ""))  # Remove the placeholder symbol for empty sequences "-"

            # The code use only upper case letters to check for features, thus use upper() to sanitize table row
            t_desc = f"{code_bool},{code_feat},{nmd_feat},{alt_orf}".replace(";", ",").upper()

            for feat in t_desc.split(","):
                trans_by_feature_dt[feat].add(trans_id)

            # Track the Gene IDs by their Coding potentiality to include them into the per-group analysis
            gene_coding_pot_dt[coding_pot].add(gene_id)
            trans_coding_pot_dt[coding_pot].add(trans_id)

            gene_trans_dt[gene_id].add(trans_id)

    # Previous versions of TransFeat table "UNPRODUCTIVE" transcripts were tag as "NON_CODING". This accounted for that
    gene_coding_pot_dt["NON_CODING"] = gene_coding_pot_dt["NON_CODING"].difference(gene_coding_pot_dt["CODING"])

    # Track transcripts belonging to Protein-coding genes and Non-coding genes
    # Note: Protein-coding genes will contain Non-coding isoforms
    for coding_pot_key, gene_set in gene_coding_pot_dt.items():
        for g_id in gene_set:
            trans_by_gene_coding_pot_dt[coding_pot_key].update(gene_trans_dt[g_id])

    coding_categories_dt = {}
    coding_categories_dt["GENES"] = gene_coding_pot_dt
    coding_categories_dt["TRANS"] = trans_coding_pot_dt

    return trans_by_feature_dt, coding_categories_dt


def group_models_into_categories(gtf_obj, coding_categories_dt):

    print(time.asctime(), f"Classifying Genes and Transcripts into diverse categories")

    multiso_genes, singleiso_monoexon_genes, singleiso_multiexon_genes = [set() for _ in range(3)]
    multiso_trans, singleiso_monoexon_trans, singleiso_multiexon_trans = [set() for _ in range(3)]

    for g_id, g_transcripts in gtf_obj.gene_trans_dt.items():
        # Multiple isoforms
        if len(g_transcripts) > 1:
            multiso_genes.add(g_id)
            multiso_trans.update(g_transcripts)
        # Single isoform
        else:
            t_id = sorted(g_transcripts)[0]
            t_exons = gtf_obj.trans_exons_dt[t_id]
            # Single-Isoform Intron-containing
            if len(t_exons) > 1:
                singleiso_multiexon_genes.add(g_id)
                singleiso_multiexon_trans.add(t_id)
            # Single-Isoform Monoexonic
            else:
                singleiso_monoexon_genes.add(g_id)
                singleiso_monoexon_trans.add(t_id)

    # Intron-containing Genes. No need to track it as it is redundant with other subcategories
    multiexonic_genes = multiso_genes | singleiso_multiexon_genes
    multiexonic_trans = multiso_trans | singleiso_multiexon_trans

    # Total
    all_genes = multiso_genes | singleiso_multiexon_genes | singleiso_monoexon_genes
    all_trans = multiso_trans | singleiso_multiexon_trans | singleiso_monoexon_trans

    # Genes / Transcripts by their coding potentiality
    coding_genes = coding_categories_dt["GENES"]["CODING"]
    coding_trans = coding_categories_dt["TRANS"]["CODING"]
    unprod_trans = coding_categories_dt["TRANS"]["UNPRODUCTIVE"]
    coding_trans = coding_trans - unprod_trans  # Remove unproductive transcripts from the coding transcripts category

    noncoding_genes = coding_categories_dt["GENES"]["NON_CODING"]
    noncoding_trans = coding_categories_dt["TRANS"]["NON_CODING"]

    # Track the genes and transcripts in the GTF not in the TransFeat table
    notfound_genes = all_genes.difference(coding_genes | noncoding_genes)
    notfound_trans = all_trans.difference(coding_trans | unprod_trans | noncoding_trans)

    # Note: KEY names are important. Do not change!
    gene_categories_dt = {}
    gene_categories_dt["TOTAL"] = (all_genes, all_trans)

    gene_categories_dt["NOT_FOUND"] = (notfound_genes, notfound_trans)

    gene_categories_dt["PROTEIN_CODING"] = (coding_genes, coding_trans)
    gene_categories_dt["UNPRODUCTIVE"] = (coding_genes, unprod_trans)
    gene_categories_dt["NON_CODING"] = (noncoding_genes, noncoding_trans)

    gene_categories_dt["SINGLEISO_MONOEXONIC"] = (singleiso_monoexon_genes, singleiso_monoexon_trans)
    gene_categories_dt["SINGLEISO_MULTIEXONIC"] = (singleiso_multiexon_genes, singleiso_multiexon_trans)
    gene_categories_dt["MULTIISO_MULTIEXONIC"] = (multiso_genes, multiso_trans)

    # Track the number of possibly unstranded transcripts, these are not processed by TranSuite
    unstranded_genes, unstranded_trans = set(), set()
    for t_id, t_strand in gtf_obj.trans_sense_dt.items():
        if t_strand not in {"+", "-"}:
            unstranded_trans.add(t_id)
            g_id = gtf_obj.trans_gene_dt[t_id]
            unstranded_genes.add(g_id)

    if unstranded_genes:
        gene_categories_dt["UNSTRANDED"] = (unstranded_genes, unstranded_trans)

    return gene_categories_dt


def get_features_numbers(trans_group, trans_by_feature_dt):

    # 0) Define the Keys / Tags to analyze and report
    # 1) Classify the transcripts into "major" single-feature categories
    # 2) Classify the transcripts into the actual output categories (which contain combination of features)

    # 1) Macro-classification dictionary
    # 1.1) Coding features: "CODING", "NAGNAG", "AS_5UTR", "AS_3UTR", "AS_5UTR_3UTR"
    # 1.2) Non-coding features: "NOT_FOUND", "NON_CODING", "UNPRODUCTIVE", "NO_ORF", "SHORT_ORF", "PTC"
    # 1.3) NMD-related features: "NMD", "DS_SJ", "LONG_3UTR", "UORF", "OUORF", "NON_INFRAME", "INFRAME"

    # NOTE: order is important! The analysis check for substring (feature) inside a longer string (trans description)
    # Feature name are subnames of other features, ex: string CODING is contain between NON_CODING, UORF within OVUORF
    # For this reason, the analysis must be done on separate lists using break statements to don't test everything
    coding_pot_tags = ["NON_CODING", "CODING", "UNPRODUCTIVE", "NOT_FOUND"]
    noncoding_feats = ["NO_ORF", "SHORT_ORF"]
    coding_feats = ["NAGNAG", "AS_5UTR", "AS_3UTR", "AS_5UTR_3UTR"]
    unprod_feats = ["PTC", "NMD", "DS_SJ", "LONG_3UTR"]
    altorf_feats = ["OUORF", "UORF", "NOT_INFRAME", "INFRAME", "LDORF"]

    # Note: Add placeholder tag "-" to avoid warning message
    feat_func_keys = coding_pot_tags + noncoding_feats + coding_feats + unprod_feats + altorf_feats
    feat_table_keys = [k for k in trans_by_feature_dt.keys()]

    # Remove feature tags that are not used or should be ignored during this check
    ignore_feats = {"AS_5UTR_3UTR", "NOT_FOUND", "-"}
    feat_func_keys = set([f for f in feat_func_keys if f not in ignore_feats])
    feat_table_keys = set([f for f in feat_table_keys if f not in ignore_feats])

    # Check if there are feature in the tables that are not analyzed by the report pipeline
    inter_keys_func = feat_table_keys.difference(feat_func_keys)
    inter_keys_table = feat_func_keys.difference(feat_table_keys)

    if inter_keys_func:
        print(time.asctime(), f"WARNING: There are {len(inter_keys_func)} features in the table that are not analyzed "
        f"by the report function. This features will be ignored: {' '.join([e for e in sorted(inter_keys_func)])}\n")

    if inter_keys_table:
        print(time.asctime(), f"WARNING: There are {len(inter_keys_table)} of the analyzed features that seems to be "
        f"missing from the TransFeat table: {' '.join([e for e in sorted(inter_keys_table)])}")

    tag_dt = defaultdict(set)
    for t_id in sorted(trans_group):

        # A) Check to which coding functionality the transcript belongs to
        for feat_tag in coding_pot_tags:
            if feat_tag not in trans_by_feature_dt.keys():
                # print(time.asctime(), f"WARNING: Feature {feat_tag} not found in the TransFeat table.")
                continue
            if t_id in trans_by_feature_dt[feat_tag]:
                tag_dt[feat_tag].add(t_id)
                break  # break to avoid multiple classification (CODING in NON_CODING)

        # B) Check to which Non-coding feature
        for feat_tag in noncoding_feats:
            if feat_tag not in trans_by_feature_dt.keys():
                continue
            if t_id in trans_by_feature_dt[feat_tag]:
                tag_dt[feat_tag].add(t_id)
                break  # break to avoid multiple classification

        for feat_tag in altorf_feats:
            if feat_tag not in trans_by_feature_dt.keys():
                continue
            if t_id in trans_by_feature_dt[feat_tag]:
                tag_dt[feat_tag].add(t_id)
                break  # break to avoid multiple classification

        # From te other list DO NOT break, as the transcripts can contain multiple of these features
        for feat_list in [coding_feats, unprod_feats]:
            for feat_tag in feat_list:
                if feat_tag not in trans_by_feature_dt.keys():
                    continue
                if t_id in trans_by_feature_dt[feat_tag]:
                    tag_dt[feat_tag].add(t_id)

    # 2.1) Major classifications according to their coding potentiality:
    # "CODING", "UNPRODUCTIVE", "NON_CODING", "NOT_FOUND"

    # 2.2) Subgroups of "CODING" classifications:
    # "PROTEIN_VARIANT", "AS_5UTR", "AS_3UTR", "NAGNAG", "AS_5UTR_3UTR", "AS_5UTR_NAGNAG", "AS_3UTR_NAGNAG", "AS_5UTR_3UTR_NAGNAG"

    # 2.3) Subgroups of "UNPRODUCTIVE" classifications
    # "PTC", "PTC_NMD", "PTC_DSSJ", "PTC_LONG3UTR", "PTC_OVUORF",
    # "PTC_DSSJ_OVUORF", "PTC_DSSJ_LONG3UTR", "PTC_LONG3UTR_OVUORF", "PTC_DSSJ_LONG3UTR_OVUORF"

    # 2.4) Subgroups of "NON_CODING" classifications
    # "NON_CODING", "NO_ORF", "SHORT_ORF"

    # 3) Classify transcripts into their category
    # To be able to write each subcategory in a separate part of the dictionary, the function must return multi-dict
    coding_feats_dt, noncoding_feats_dt, unprod_feats_dt = [defaultdict(set) for _ in range(3)]
    for t_id in sorted(trans_group):

        # 3.1) Classify NON_CODING transcripts into their subcategory
        if t_id in tag_dt["NON_CODING"]:
            if t_id in tag_dt["NO_ORF"]:
                noncoding_feats_dt["NO_ORF"].add(t_id)
            elif t_id in tag_dt["SHORT_ORF"]:
                noncoding_feats_dt["SHORT_ORF"].add(t_id)
            elif t_id in tag_dt["PTC"]:
                noncoding_feats_dt["PTC_NC"].add(t_id)
            else:
                noncoding_feats_dt["NONCODING_UNCLASSIFIED"].add(t_id)

        # 3.2) Classify CODING transcripts into their subcategory:
        # "PROTEIN_VARIANT", "AS_5UTR", "AS_3UTR", "NAGNAG",
        # "AS_5UTR_3UTR_NAGNAG", "AS_5UTR_3UTR", "AS_5UTR_NAGNAG", "AS_3UTR_NAGNAG"
        # NOTE: Go for more complex (multiple features) to single-feature to allow for exclusive classification
        if t_id in tag_dt["CODING"]:
            if t_id in tag_dt["AS_5UTR"] and t_id in tag_dt["AS_3UTR"] and t_id in tag_dt["NAGNAG"]:
                coding_feats_dt["AS_5UTR_3UTR_NAGNAG"].add(t_id)
            elif t_id in tag_dt["AS_5UTR"] and t_id in tag_dt["AS_3UTR"]:
                coding_feats_dt["AS_5UTR_3UTR"].add(t_id)
            elif t_id in tag_dt["AS_5UTR"] and t_id in tag_dt["NAGNAG"]:
                coding_feats_dt["AS_5UTR_NAGNAG"].add(t_id)
            elif t_id in tag_dt["AS_3UTR"] and t_id in tag_dt["NAGNAG"]:
                coding_feats_dt["AS_3UTR_NAGNAG"].add(t_id)
            elif t_id in tag_dt["AS_5UTR"]:
                coding_feats_dt["AS_5UTR"].add(t_id)
            elif t_id in tag_dt["AS_3UTR"]:
                coding_feats_dt["AS_3UTR"].add(t_id)
            elif t_id in tag_dt["NAGNAG"]:
                coding_feats_dt["NAGNAG"].add(t_id)
            else:
                coding_feats_dt["PROTEIN_VARIANT"].add(t_id)

        # 3.3) Classify UNPRODUCTIVE transcripts into their subcategory
        # "PTC", "PTC_NMD", "PTC_DSSJ", "PTC_LONG3UTR", "PTC_OVUORF",
        # "PTC_DSSJ_LONG3UTR_OVUORF", "PTC_DSSJ_OVUORF", "PTC_DSSJ_LONG3UTR", "PTC_LONG3UTR_OVUORF"
        if t_id in tag_dt["UNPRODUCTIVE"]:
            # PTC is a fundamentally required signal
            if t_id in tag_dt["PTC"]:
                if t_id in tag_dt["DS_SJ"] and t_id in tag_dt["LONG_3UTR"] and t_id in tag_dt["OUORF"]:
                    unprod_feats_dt["PTC_DSSJ_LONG3UTR_OUORF"].add(t_id)
                elif t_id in tag_dt["DS_SJ"] and t_id in tag_dt["LONG_3UTR"]:
                    unprod_feats_dt["PTC_DSSJ_LONG3UTR"].add(t_id)
                elif t_id in tag_dt["DS_SJ"] and t_id in tag_dt["OUORF"]:
                    unprod_feats_dt["PTC_DSSJ_OUORF"].add(t_id)
                elif t_id in tag_dt["LONG_3UTR"] and t_id in tag_dt["OUORF"]:
                    unprod_feats_dt["PTC_LONG3UTR_OUORF"].add(t_id)
                elif t_id in tag_dt["DS_SJ"]:
                    unprod_feats_dt["PTC_DSSJ"].add(t_id)
                elif t_id in tag_dt["LONG_3UTR"]:
                    unprod_feats_dt["PTC_LONG3UTR"].add(t_id)
                elif t_id in tag_dt["OUORF"]:
                    unprod_feats_dt["PTC_OUORF"].add(t_id)
                else:
                    unprod_feats_dt["PTC"].add(t_id)
            else:
                # This shouldn't happen
                unprod_feats_dt["UNPRODUCTIVE_UNCLASSIFIED"].add(t_id)

    return coding_feats_dt, unprod_feats_dt, noncoding_feats_dt


def get_categories_numbers(models_categories_dt, trans_by_feature_dt):

    categories_dt = {}

    all_genes, all_trans = models_categories_dt["TOTAL"]
    tot_genes, tot_trans = len(all_genes), len(all_trans)

    notfound_genes, notfound_trans = models_categories_dt["NOT_FOUND"]
    noncoding_genes, noncoding_trans = models_categories_dt["NON_CODING"]
    coding_genes, coding_trans = models_categories_dt["PROTEIN_CODING"]
    _, unprod_trans = models_categories_dt["UNPRODUCTIVE"]

    for cat_name, (cat_genes, cat_trans) in models_categories_dt.items():
        # Dict to track the number in the category
        cat_dt = {}

        # Using multiple dictionaries for the sub-categories to facilitate their write order in the output table
        general_subcat, gene_subcat, trans_subcat, coding_trans_subcat, unprod_trans_subcat, noncoding_trans_subcat = [{} for _ in range(6)]

        # TOTAL genes in group
        n_genes, p_genes = get_per(len(cat_genes), tot_genes)
        n_trans, p_trans = get_per(len(cat_trans), tot_trans)

        general_subcat["GROUP_GENES"] = (n_genes, p_genes)
        general_subcat["GROUP_TRANS"] = (n_trans, p_trans)

        # NOT_FOUND subcategory
        group_notfound_genes = cat_genes.intersection(notfound_genes)
        group_notfound_trans = cat_trans.intersection(notfound_trans)

        n_genes_notfound, p_genes_notfound = get_per(len(group_notfound_genes), tot_genes)
        n_trans_notfound, p_trans_notfound = get_per(len(group_notfound_trans), tot_trans)

        gene_subcat["NOT_FOUND_GENES"] = (n_genes_notfound, p_genes_notfound)
        trans_subcat["NOT_FOUND_TRANS"] = (n_trans_notfound, p_trans_notfound)

        # CODING genes
        group_coding_genes = cat_genes.intersection(coding_genes)
        group_coding_trans = cat_trans.intersection(coding_trans)
        group_unprod_trans = cat_trans.intersection(unprod_trans)
        group_coding_trans = group_coding_trans - group_unprod_trans

        n_genes_coding, p_genes_coding = get_per(len(group_coding_genes), tot_genes)
        n_trans_coding, p_trans_coding = get_per(len(group_coding_trans), tot_trans)
        n_trans_unprod, p_trans_unprod = get_per(len(group_unprod_trans), tot_trans)

        gene_subcat["CODING_GENES"] = (n_genes_coding, p_genes_coding)
        trans_subcat["CODING_TRANS"] = (n_trans_coding, p_trans_coding)
        trans_subcat["UNPROD_TRANS"] = (n_trans_unprod, p_trans_unprod)

        # NON_CODING genes
        group_noncoding_genes = cat_genes.intersection(noncoding_genes)
        group_noncoding_trans = cat_trans.intersection(noncoding_trans)
        n_genes_noncoding, p_genes_noncoding = get_per(len(group_noncoding_genes), tot_genes)
        n_trans_noncoding, p_trans_noncoding = get_per(len(group_noncoding_trans), tot_trans)

        gene_subcat["NONCODING_GENES"] = (n_genes_noncoding, p_genes_noncoding)
        trans_subcat["NONCODING_TRANS"] = (n_trans_noncoding, p_trans_noncoding)

        coding_feats_dt, unprod_feats_dt, noncoding_feats_dt = get_features_numbers(cat_trans, trans_by_feature_dt)

        subcat_feats_dt = [coding_trans_subcat, unprod_trans_subcat, noncoding_trans_subcat]
        feats_dt = [coding_feats_dt, unprod_feats_dt, noncoding_feats_dt]

        for subcat_dt, feat_dt in zip(subcat_feats_dt, feats_dt):
            for feat_tag, feat_trans in feat_dt.items():
                n_feat_trans, p_feat_trans = get_per(len(feat_trans), tot_trans)
                subcat_dt[feat_tag] = (n_feat_trans, p_feat_trans)

        cat_dt["GROUP_TOTAL"] = general_subcat
        cat_dt["GENE_SUBCATEGORIES"] = gene_subcat
        cat_dt["TRANS_SUBCATEGORIES"] = trans_subcat
        cat_dt["CODING_TRANS_SUBCATEGORIES"] = coding_trans_subcat
        cat_dt["UNPROD_TRANS_SUBCATEGORIES"] = unprod_trans_subcat
        cat_dt["NONCODING_TRANS_SUBCATEGORIES"] = noncoding_trans_subcat

        categories_dt[cat_name] = cat_dt

    return categories_dt


def get_table_lines(cat_dt, brief=False):

    # Allow brief summary description (report only major categories) to create more readable tables
    major_cat = {"GROUP_GENES", "GROUP_TRANS", "CODING_GENES", "NONCODING_GENES", "CODING_TRANS", "NONCODING_TRANS", "UNPROD_TRANS"}

    lines = []
    for subcat_name, subcat_dt in cat_dt.items():
        # Sort subcat dict by length and then alphabetically
        if subcat_dt:
            lines.append(f"{subcat_name}")
            for feat_tag, (feat_n, feat_p) in sorted(subcat_dt.items(), key=lambda kv: (len(kv[0]), kv[0])):
                # Avoid printing/writing empty values
                if feat_n:
                    # Allow to generate more concise information for some tables
                    if brief:
                        if feat_tag in major_cat:
                            lines.append(f"{feat_tag},{feat_n},{feat_p}")
                    else:
                        lines.append(f"{feat_tag},{feat_n},{feat_p}")

    return lines


def write_table(lines, outpath, outname):

    outname = outname.replace(".csv", "")
    outname = f"{outname}"

    outfile = os.path.join(outpath, outname)

    with open(f"{outfile}.csv", "w+") as fh:
        for ln in lines:
            ln = ln.strip("\n")
            fh.write(f"{ln}\n")

    print(time.asctime(), f"{outfile}.csv")


def write_transfeat_summary_tables(categories_dt, outpath, outname, sep=True):

    # These are redundant with the category TOTAL, no need to write them into a table
    cat_ignore = {"PROTEIN_CODING", "UNPRODUCTIVE", "NON_CODING"}

    # These are groups that could be present but usually aren't
    cat_others = ["NOT_FOUND", "UNSTRANDED"]

    # These are groups that should be written into the same table
    table_groups = {"BY_CODING_POTENTIAL": ["TOTAL"],
                    "BY_STRUCTURE": ["SINGLEISO_MONOEXONIC", "SINGLEISO_MULTIEXONIC", "MULTIISO_MULTIEXONIC"]}

    # 1) First write the categories that should be written in the same table together and track them
    cat_processed = set()
    for group_name, cat_group in table_groups.items():

        # Make more concise table for "BY_STRUCTURE" category
        brief_flag = False
        if group_name == "BY_STRUCTURE":
            brief_flag = True

        group_lines = ["Category,Number,Percentage\n"]
        for cat_name in cat_group:
            try:
                cat_dt = categories_dt[cat_name]
            except KeyError as e:
                print(f"WARNING: Category {cat_name} not found")
                continue

            subgroup_lines = get_table_lines(cat_dt, brief=brief_flag)

            if subgroup_lines:
                # Optionally add separator (empty row) between groups
                if sep and len(group_lines) > 1:
                    group_lines.append("\n")

                group_lines.append(cat_name)
                group_lines.extend(subgroup_lines)

            # Track the processed groups to not write them again
            cat_processed.add(cat_name)

        # Write the group into a table
        write_table(group_lines, outpath, f"{outname}_numbers_{group_name.lower()}")

    # Update the categories to ignore (either because redundant or already processed) and write the rest
    cat_ignore = cat_ignore | cat_processed

    # 2) Write any other category that has not been processed yet
    for cat_name, cat_dt in categories_dt.items():

        # Ignore the processed categories
        if cat_name in cat_ignore:
            continue

        subcat_lines = []
        if cat_dt:
            subcat_lines.append(f"{cat_name}")
            for subcat_name, subcat_dt in cat_dt.items():
                # Sort subcat dict by length and then alphabetically
                if subcat_dt:
                    subcat_lines.append(f"{subcat_name}")
                    for feat_tag, (feat_n, feat_p) in sorted(subcat_dt.items(), key=lambda kv: (len(kv[0]), kv[0])):
                        # Avoid printing empty values
                        if feat_n:
                            subcat_lines.append(f"{feat_tag}\t{feat_n}\t{feat_p}")

            # Write the group into a table
            write_table(subcat_lines, outpath, f"{outname}_numbers_{cat_name.lower()}")


def generate_transfeat_summary(gtf_file, transfeat_table):

    # Check if files exist
    for fl in [gtf_file, transfeat_table]:
        if not os.path.isfile(fl):
            sys.exit(f"File {fl} does not exist.")

    # Create output file report and subfolders
    outpath = os.path.dirname(transfeat_table)
    outname = os.path.basename(transfeat_table).replace(".csv", "")

    # report_outfile = os.path.join(outpath, f"{outname}_report.txt")
    # table_subfolder = os.path.join(outpath, f"{outpath}_tables")

    print("\n")
    print(time.asctime(), f'Generating summary of TransFeat results ({transfeat_table})', flush=True)

    # 1) Get transcriptome information
    gtf_obj = create_gtf_object(gtf_file)

    # 2) Get information from TransFeat table
    trans_by_feature_dt, coding_categories_dt = get_transfeat_data(transfeat_table)

    # 3) Create a dictionary with the gene categories to analyze (Mono-exonic / Intron-containing genes, etc)
    models_by_categories_dt = group_models_into_categories(gtf_obj, coding_categories_dt)

    # 4) Get gene categories numbers
    categories_dt = get_categories_numbers(models_by_categories_dt, trans_by_feature_dt)

    # 5) Write tables
    write_transfeat_summary_tables(categories_dt, outpath, outname)
