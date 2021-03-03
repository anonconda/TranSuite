"""
Created on Wed May 15 14:25:00 2019
@author: Juan C Entizne
@email: e.entizne[at]dundee.ac.uk
"""
import os
import time
import warnings

from lib.findlorf.findlorf_tools import *
from lib.transfeat.identify_coding_features import *
from lib.transfeat.identify_non_coding_features import *

from lib.parsing.gtf_object_tools import create_gtf_object
from lib.parsing.fasta_parsing_tools import get_fasta_sequences, write_fasta_file
from lib.report.transfeat_report import generate_transfeat_summary

warnings.filterwarnings("ignore")


def transfeat_main(gtf, fasta, outpath, outname, pep_len=50, ptc_len=70, uorf_len=10, sj_dist=50, utr3_len=350,
                   orf_index=None):

    print("\n")
    print(time.asctime(), "Starting TransFeat analysis")

    # +1 AA to account for stop codons during the AA length check
    pep_len += 1

    # In case the user pass the name with a file extension, remove it
    if outname.endswith(".gtf"):
        outname = outname.replace(".gtf", "")
    outname += "_transfeat"

    # Create output folder
    outfolder = os.path.join(outpath, outname)
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)

    # Get transcriptome annotation
    gtf_obj = create_gtf_object(gtf)

    # Upload transcripts sequences from fasta file
    trans_seq_dt = get_fasta_sequences(fasta)

    # Generate transcripts sequence information
    fasta_header_dt, cds_seq_dt, pep_seq_dt = translate_transcript_cds(trans_seq_dt, gtf_obj)

    print(time.asctime(), "Retrieving ORF information")
    if orf_index:
        print(time.asctime(), "Uploading ORF information from ORF index file")
        with open(orf_index) as orf_index_fh:
            orf_dt = json.load(orf_index_fh)

    else:
        # Generate ORF index file
        _, orf_index = find_transcripts_orf_information(gtf, trans_seq_dt, gtf_obj, outfolder)

        print(time.asctime(), "Uploading ORF information from ORF index file")
        with open(orf_index) as orf_index_fh:
            orf_dt = json.load(orf_index_fh)

    if not orf_dt:
        sys.exit("No ORF information found.")

    # Get transcript start-codon relative position
    relative_start_dt = get_transcript_start_codon_relative_position(gtf_obj)

    # Select authentic stop-codon (at gene level)
    auth_stop_dt = get_genes_authentic_stop_codon_position(gtf_obj)

    print(time.asctime(), "Retrieving alternative ORFs information")

    # Identify transcripts with long downstream ORF
    is_longer_dorf_dt, ldorf_data_dt = identify_longer_dorf(gtf_obj, relative_start_dt, orf_dt, trans_seq_dt)

    # Identify transcripts with upstream ORF
    is_uorf_dt, uorf_data_dt, urof_categories = identify_uorf(gtf_obj, relative_start_dt, orf_dt, trans_seq_dt, uorf_len)

    print(time.asctime(), "Identifying Non-Coding features")

    # Identify transcripts without an annotated CDS
    is_orf_absent_dt = is_orf_absent(gtf_obj)

    # Identify transcripts with "Premature Termination Codons" (PTC)
    is_ptc_dt = is_ptc(gtf_obj, ptc_len)

    # Identify transcripts coding for short peptides
    # The identification of "short peptides" is done after the PTC check to avoid redundancy of classification
    is_orf_short_dt = is_orf_short(gtf_obj, pep_len)

    is_long_3utr_dt = is_long_3utr(gtf_obj, utr3_len)

    # Get transcripts groups (PTC transcripts, long 3' UTR transcripts, etc) to use for NMD classification
    ptc_trans = set([t_id for t_id, t_bool in is_ptc_dt.items() if t_bool is True])
    long_3utr_trans = set([t_id for t_id, t_bool in is_long_3utr_dt.items() if t_bool is True])
    ov_uorf_trans = urof_categories["overlapping"]
    uorf_trans = urof_categories["not_overlapping"]

    is_nmd_dt, is_dssj_dt = is_nmd(gtf_obj, auth_stop_dt,
                                   sj_dist_th=sj_dist, ptc_trans=ptc_trans,
                                   long_3utr_trans=long_3utr_trans, ov_uorf_trans=ov_uorf_trans, uorf_trans=uorf_trans)

    nmd_features_dt = generate_nmd_features_lines(gtf_obj, is_nmd_dt, is_ptc_dt, is_dssj_dt, is_long_3utr_dt, urof_categories)

    # Identify AS in UTR and NAGNAG features
    as_in_utr_dt, as_utr_location_dt, nagnag_dt = identify_similar_coding_features(gtf_obj)

    # Some transcripts are missing from the features_dict either because:
    # (1) Not present in the FASTA file, (2) it doesn't have a CDS, or (3) it's the only transcript in the overlap group
    # To avoid KeyError further downstream in 'generate_feature_tag'. These dictionaries return None by default

    # Dictionary of features to annotate
    feature_dicts = {}
    feature_dicts["Auto"] = gtf_obj.trans_gene_dt
    feature_dicts["No_ORF"] = is_orf_absent_dt
    feature_dicts["Short_ORF"] = is_orf_short_dt
    feature_dicts["Long_3UTR"] = is_long_3utr_dt
    feature_dicts["PTC"] = is_ptc_dt
    feature_dicts["NMD"] = is_nmd_dt
    feature_dicts["ds_SJ"] = is_dssj_dt
    feature_dicts["NMD_features"] = nmd_features_dt
    feature_dicts["uORF"] = urof_categories
    feature_dicts["ldORF"] = is_longer_dorf_dt
    feature_dicts["AS_in_UTR"] = as_in_utr_dt
    feature_dicts["AS_Location"] = as_utr_location_dt
    feature_dicts["NAGNAG"] = nagnag_dt

    # Generate the features to annotate into the output table
    coding_potentiality_dt, coding_features_dt, alternative_ORF_dt = generate_feature_tag(gtf_obj, feature_dicts)

    # These dictionaries are required to write the TransFeat table
    table_info_dicts = {}
    table_info_dicts["Coding_potentiality"] = coding_potentiality_dt
    table_info_dicts["Coding_features"] = coding_features_dt
    table_info_dicts["NMD_features"] = nmd_features_dt
    table_info_dicts["Alternative_ORF"] = alternative_ORF_dt

    # Get alternative ORFs IDs to validate classification on table
    ldorf_trans = set([t_id for t_id in is_longer_dorf_dt if is_longer_dorf_dt[t_id] is True])
    uorf_trans = set([t_id for t_id in is_uorf_dt if is_uorf_dt[t_id] is True])

    # Write TransFeat table output
    transfeat_table = write_transfeat_table(gtf_obj, table_info_dicts, pep_seq_dt, outfolder, outname,
                                            ldorf_ids=ldorf_trans, uorf_ids=uorf_trans, pep_len=pep_len)

    # Write GENERAL/TOTAL output fasta files
    fasta_outfile = os.path.join(outfolder, outname)
    write_fasta_file(cds_seq_dt, f'{fasta_outfile}_nuc.fasta', fasta_header_dt)
    write_fasta_file(pep_seq_dt, f'{fasta_outfile}_pep.fasta', fasta_header_dt)

    # These dictionaries are required to write the fasta files
    sequences_dicts = {}
    sequences_dicts["Headers"] = fasta_header_dt
    sequences_dicts["Exonic_seq"] = trans_seq_dt
    sequences_dicts["CDS_seq"] = cds_seq_dt
    sequences_dicts["Peptide_seq"] = pep_seq_dt
    sequences_dicts["ldORF_data"] = ldorf_data_dt
    sequences_dicts["uORF_data"] = uorf_data_dt

    # Write ADDITIONAL output fasta files
    write_subcategories_fasta(gtf_obj, transfeat_table, sequences_dicts, outfolder, outname)

    # Generate tables summarizing the main transfeat data
    generate_transfeat_summary(gtf, transfeat_table)

    # Generate table with additional NMD information
    write_NMD_table(gtf_obj, feature_dicts, sequences_dicts, outfolder, outname)

    # Return output file for TransAll function
    return transfeat_table
