"""
Created on Wed May 15 14:25:00 2019
@author: Juan C Entizne & Mark Spensley
@email: e.entizne[at]dundee.ac.uk
"""

import os
import sys
import time
import filecmp
from collections import defaultdict

from lib.parsing.gtf_object_tools import create_gtf_object

import lib.transfix.load_annotation as load
from lib.transfix.reannotate_to_max_orf import characterise_max_orfs
from lib.transfix.transfix_tools import remove_transcripts_without_cds, write_transfix_tables
from lib.transfix.transfix_tools import get_transcript_data_from_gff_obj, fix_chimeric_start_codon, fix_atg_position

from lib.parsing.fasta_parsing_tools import get_fasta_sequences, write_fasta_file
from lib.parsing.gtf_parsing_tools import annotate_cds_into_gtf, add_features_to_gtf


def transfix_main(gtf_file, fasta, outpath, outname, chimeric=None, iter_th=5):

    print("\n")
    print(time.asctime(), "Starting TransFix analysis")

    if not 0 < iter_th <= 5:
        sys.exit(f"ERROR: The number of iterations must be within 0 and 5, not {iter_th}.")

    # If iter_th value is 0, do iterations indefinitely. Disabled for now due to unexpected bug for indefinite iteration
    check_iter = False
    if iter_th:
        check_iter = True

    # In case the user pass the name with a file extension, remove it
    if outname.endswith(".gtf"):
        outname = outname.replace(".gtf", "")
    outname += "_transfix"

    # Create output folder
    outfolder = os.path.join(outpath, outname)
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)

    # Classification categories
    cat_dt = defaultdict(set)

    # Track the fixed start-codon in the 1st fixing cycle to correct chimeric models
    gene_atg_pos_dt = defaultdict(list)

    # Assorted dictionaries to track features
    trans_cds_dt = {}  # To save the CDS coordinates of the translations
    trans_cds_seq_dt, trans_header_dt = {}, {}  # To save the CDS sequences to write them into nucl/peptide fasta files
    cycle_trans_dt = defaultdict(set)  # To track the processed transcripts at each cycle

    # Get transcriptome information
    gtf_obj = create_gtf_object(gtf_file)

    # Remove transcripts that do NOT have an annotated CDS
    # This method introduce the tag ".transfix.temp." that marks the temporary files to be removed further downstream
    gtf_known_cds, trans_with_cds, trans_without_cds = remove_transcripts_without_cds(gtf_file, outfolder)

    cat_dt["cds_not_found"].update(trans_without_cds)  # Track those models without CDS

    # Assign the file with transcripts containing CDS as the starting GTF file
    gtf_1 = gtf_known_cds

    # Create empty file to make first files comparison
    gtf_2 = os.path.join(outfolder, "empty.transfix.temp.gtf")
    with open(gtf_2, "w+"):
        pass

    # Upload transcripts sequences from fasta file
    trans_sequences_dt = get_fasta_sequences(fasta)

    # Variables to define a maximum number of iterations
    i = 0
    # The iteration will stop either when the previous output is the same as the current one, or the iter limit is reach
    while not filecmp.cmp(gtf_1, gtf_2):

        # This check is done before increasing the counter because
        # if the user specify the value as 1, the clearer meaning of this is to stop after completing 1 iteration
        if check_iter:
            if i >= iter_th:
                break
        i += 1

        print("\n")
        print(time.asctime(), "Iteration number {}:".format(i))
        print(time.asctime(), "Processing annotation file: ", gtf_1)

        print(time.asctime(), "Loading transcriptome information")
        locus_dict = load.get_models_gtf(gtf_1)

        print(time.asctime(), "Converting transcriptome information into GFF3 format")
        gff3_file = load.gtf_to_gff3(gtf_1)

        print(time.asctime(), "Loading GFF3 information")
        gff3_models = load.get_models_gff3(gff3_file)

        print(time.asctime(), "Selecting Genes with AUG start codon")
        _, _, gff3_models = characterise_max_orfs(gff3_models)

        print(time.asctime(), "Fixing transcripts start codon")

        n_genes = len(locus_dict.keys())
        for z, locus_id in enumerate(sorted(locus_dict.keys())):

            progress_bar = True
            if progress_bar:
                print(f'Processing Gene {locus_id}, {(z/n_genes)*100:.1f}% complete ({z+1}/{n_genes})')

            if locus_id not in gff3_models:
                cat_dt["absent_gff3"].add(locus_id)
                # absent_gff3.add(locus_id)
                continue

            if gff3_models[locus_id].transposon is True:
                cat_dt["retro_transposons"].add(locus_id)
                # retro_transposons.add(locus_id)
                continue

            # Start CDS
            atg_pos = gff3_models[locus_id].rep_atg

            # Keep only the ATG positions found in the first cycle of fixing
            if not gene_atg_pos_dt[locus_id]:
                gene_atg_pos_dt[locus_id].append(atg_pos)

            trans_data_dt = get_transcript_data_from_gff_obj(locus_id, locus_dict, trans_sequences_dt)

            # Fix the start-codon position for the Gene group
            grp_output_dt, grp_cat_dt = fix_atg_position(trans_data_dt, atg_pos)

            # Update the transcripts categories
            cat_dt["cds_not_found"].update(grp_cat_dt["cds_not_found"])
            cat_dt["seq_not_present"].update(grp_cat_dt["seq_not_present"])
            cat_dt["atg_not_in_cds"].update(grp_cat_dt["atg_not_in_cds"])
            cat_dt["start_codon_not_atg"].update(grp_cat_dt["start_codon_not_atg"])

            cat_dt["cds_not_found_lines"].update(grp_cat_dt["cds_not_found_lines"])
            cat_dt["seq_not_present_lines"].update(grp_cat_dt["seq_not_present_lines"])
            cat_dt["atg_not_in_cds_lines"].update(grp_cat_dt["atg_not_in_cds_lines"])
            cat_dt["start_codon_not_atg_lines"].update(grp_cat_dt["start_codon_not_atg_lines"])

            cat_dt["rejected_start_codons"].update(grp_cat_dt["rejected_start_codons"])
            cat_dt["processed_transcripts"].update(grp_cat_dt["processed_transcripts"])

            trans_cds_dt.update(grp_output_dt["trans_cds_dt"])
            trans_cds_seq_dt.update(grp_output_dt["trans_cds_seq_dt"])
            trans_header_dt.update(grp_output_dt["trans_header_dt"])

        # Remove transcripts processed in current iteration
        cat_dt["cds_not_found"] -= cat_dt["atg_not_in_cds"] | cat_dt["processed_transcripts"]
        cat_dt["atg_not_in_cds"] -= cat_dt["processed_transcripts"]
        cat_dt["start_codon_not_atg"] -= cat_dt["atg_not_in_cds"] | cat_dt["cds_not_found"] | cat_dt["processed_transcripts"]

        # Filter processed transcripts
        new_gtf = load.filter_gtf(gtf_1, cat_dt["processed_transcripts"], cat_dt["rejected_start_codons"], i, outfolder)
        gtf_1, gtf_2 = new_gtf, gtf_1

        # Track the processed transcripts at each translation cycle
        cycle = f"Cycle_{i}"
        cycle_trans_dt[cycle].update(cat_dt["processed_transcripts"])
    print("\n")

    # Track transcripts that were not processed in the analysis
    removed_st = cat_dt["processed_transcripts"] | cat_dt["cds_not_found"] | cat_dt["atg_not_in_cds"] | \
                 cat_dt["start_codon_not_atg"] | cat_dt["seq_not_present"]

    cat_dt["unprocessed_transcripts"] = trans_with_cds - removed_st
    cat_dt["unprocessed_transcripts_lines"] = [e for e in sorted(cat_dt["unprocessed_transcripts"])]

    # TODO implement command line arguments to pass chimeric-IDs table
    # If a table specifying chimeric models is reported by the user, then TransFix can correct the ATG of these models
    if chimeric:
        chimeric_output_dt = fix_chimeric_start_codon(gtf_obj, chimeric, trans_cds_dt, trans_sequences_dt)

        trans_cds_dt.update(chimeric_output_dt["trans_cds_dt"])
        trans_cds_seq_dt.update(chimeric_output_dt["trans_cds_seq_dt"])
        trans_header_dt.update(chimeric_output_dt["trans_header_dt"])

    # Write annotation file with re-annotated CDS
    outfile = os.path.join(outfolder, outname + ".gtf")
    outfile = annotate_cds_into_gtf(gtf_obj, trans_cds_dt, outfile)

    # Write output fasta files
    outfile_fasta = os.path.join(outfolder, outname + "_nuc.fasta")
    write_fasta_file(trans_cds_seq_dt, outfile_fasta, trans_header_dt)

    trans_pep_dt = {}
    for trans, trans_seq in trans_cds_seq_dt.items():
        trans_pep_dt[trans] = trans_seq.translate(to_stop=True)

    outfile_fasta = os.path.join(outfolder, outname + "_pep.fasta")
    write_fasta_file(trans_pep_dt, outfile_fasta, trans_header_dt)

    # Write TransFix related tables
    write_transfix_tables(gtf_obj, cat_dt, cycle_trans_dt, trans_cds_dt, outfolder, outname)

    remove = True
    if remove:
        print(time.asctime(), "Removing temporary files")
        for (path, dirs, files) in os.walk(outfolder):
            for file in files:
                if "transfix.temp." in file or file == "empty.transfix.temp.gtf" or file == gff3_file:
                    os.remove(os.path.join(outfolder, file))

    # Return output file for TransAll function
    return outfile
