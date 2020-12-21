import os
import sys
import time
import json
from Bio.Seq import Seq
from collections import defaultdict


def convert_to_relative_position(exon_list):

    result, start = [], 0
    for exon in exon_list:
        exon_len = exon[-1]-exon[0]

        # Using a list-index starting with 1
        start = start+1
        end = start+exon_len

        result.append([start, end])
        start = end

    # Since list-index use a +1 coordinate system it is necessary to add +1 position at the end
    result[-1][-1] += 1

    return result


def convert_to_genomic_coord(cds, lookup_table, trans_strand, trans):

    # CDS must have a +1 index position
    cds += 1

    exons_relative = lookup_table[0]
    exons_genomic = lookup_table[1]

    # Initialize variables (to avoid warning by the EDI)
    cds_genomic = None
    relative_exon = None
    genomic_exon = None

    for i, exon in enumerate(exons_relative):

        exon_relative_st, exon_relative_end = exon

        if exon_relative_st <= cds <= exon_relative_end:
            exon_genomic_st, exon_genomic_end = exons_genomic[i]

            relative_exon = exon
            genomic_exon = exons_genomic[i]

            if trans_strand == "+":
                # exon_genomic_end > exon_genomic_st
                cds_genomic = exon_genomic_st + cds - exon_relative_st
                break

            elif trans_strand == "-":
                # exon_genomic_end < exon_genomic_st
                cds_genomic = exon_genomic_st - cds + exon_relative_st
                break

            else:
                print(f"Transcript {trans} strand must be \"+\" or \"-\", not \"{trans_strand}\"")
                cds_genomic = None
                relative_exon = None
                genomic_exon = None

        else:
            cds_genomic = None
            relative_exon = None
            genomic_exon = None

    return cds_genomic, relative_exon, genomic_exon


def convert_from_genomic_to_relative(exon_list):

    genomic_to_relative_coordinates = {}
    last_pos = 0

    # Initialize variable (to avoid warning)
    relative_pos = None
    for exon in exon_list:
        for genomic_pos, relative_pos in zip(range(exon[0], exon[-1]+1), range(last_pos, exon[-1]+1)):
            genomic_to_relative_coordinates[genomic_pos] = relative_pos

        last_pos = relative_pos+1

    return genomic_to_relative_coordinates


def create_coordinate_lookup_table(gtf_obj):

    exon_genomic_coords = gtf_obj.trans_exons_dt
    trans_sense_dt = gtf_obj.trans_sense_dt

    exon_lookup_table = {}
    for trans, exons in exon_genomic_coords.items():
        trans_strand = trans_sense_dt[trans]
        if trans_strand == "-":
            relative_exons = convert_to_relative_position(exons[::-1])
            genomic_exons = [(e[1], e[0]) for e in exons[::-1]]
        elif trans_strand == "+":
            relative_exons = convert_to_relative_position(exons)
            genomic_exons = exons
        else:
            sys.exit(f"Transcript {trans} strand must be \"+\" or \"-\", not \"{trans_strand}\"")

        exon_lookup_table[trans] = (relative_exons, genomic_exons)

    return exon_lookup_table


def get_orf_info(trans_nuc_seq, strand_known, min_protein_length):

    trans_nuc_seq = Seq(trans_nuc_seq)

    orf_list = []
    seq_len = len(trans_nuc_seq)
    for strand, nuc_seq in [(+1, trans_nuc_seq), (-1, trans_nuc_seq.reverse_complement())]:
        for frame in range(3):
            transcript_aa_seq = str(nuc_seq[frame:].translate())
            trans_len = len(transcript_aa_seq)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = transcript_aa_seq.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end-aa_start >= min_protein_length:
                    if strand == 1:
                        orf_start = frame+aa_start*3
                        orf_end = min(seq_len, frame+aa_end*3+3)
                    else:
                        orf_start = seq_len-frame-aa_end*3-3
                        orf_end = seq_len-frame-aa_start*3

                    peptide = transcript_aa_seq[aa_start:aa_end]

                    # Find CDS start
                    if "M" in set(peptide):
                        if strand == 1:
                            cds_start = orf_start + (peptide.find("M") * 3)
                            cds_end = orf_end
                            cds_seq = trans_nuc_seq[cds_start:cds_end]

                        else:
                            cds_start = orf_end - (peptide.find("M") * 3)
                            cds_end = orf_start
                            cds_seq = trans_nuc_seq[cds_end:cds_start]

                    else:
                        cds_start = -1
                        cds_end = -1
                        cds_seq = ""

                    pep_seq = transcript_aa_seq[aa_start:aa_end]
                    orf_list.append((orf_start, orf_end, frame, strand, cds_start, cds_end, cds_seq, pep_seq))

                aa_start = aa_end+1

    # If the strand of the transcript is known, remove reverse_complement results
    if strand_known:
        orf_list = [cds_info for cds_info in orf_list if cds_info[3] != -1]

    # Important! Sort by the longest CDS sequence
    orf_list = sorted(orf_list, key=lambda x: len(x[6]), reverse=True)

    return orf_list


def check_codons(cds_pair, trans_exons,  trans_seq, strand, trans):

    # The index to slice the sequence are found from left to right; therefore I must invert the sequence for "-" strands
    seq = Seq(trans_seq)
    if strand == "-":
        seq = seq.reverse_complement()

    lookup_table = convert_from_genomic_to_relative(trans_exons)

    try:
        cds_1_ix = lookup_table[cds_pair[0]]
        cds_2_ix = lookup_table[cds_pair[-1]]
    except KeyError:
        return False

    cds_seq = seq[cds_1_ix:cds_2_ix+1]
    cds_seq_first_codon = cds_seq[:3]
    cds_seq_last_codon = cds_seq[-3:]

    start_codon = "ATG"
    stop_codons = {"TAG", "TAA", "TGA"}

    start_boolean, end_boolean = False, False
    if strand == "+":
        if cds_seq_first_codon == start_codon:
            start_boolean = True
        if cds_seq_last_codon in stop_codons or cds_pair[1] == trans_exons[-1][-1]:
            end_boolean = True

    elif strand == "-":
        if cds_seq_last_codon.reverse_complement() == start_codon:
            start_boolean = True
        if cds_seq_first_codon.reverse_complement() in stop_codons or cds_pair[0] == trans_exons[0][0]:
            end_boolean = True

    else:
        sys.exit(f'Strand error for transcript {trans}')

    # assert start_boolean is True and end_boolean is True

    if not start_boolean or not end_boolean:
        print(f"ERROR: Incorrect START or STOP codon selected for transcript '{trans}'")
        print("INFO:")
        print(f"Start/Stop booleans: {start_boolean}/{end_boolean}")
        print(f"Transcript strand: {trans}/{strand}")
        print(f"Transcript exons: {trans_exons}")
        print(f"Transcript sequence: {trans_seq}")
        print(f"CDS pair: {cds_pair}")
        print("\n")
        exit()


def fix_cds_offset(cds_pair, strand, trans):

    cds_1, cds_2 = sorted(cds_pair)

    # Fix the +/- 1 coordinate off-set resulting from the convert_to_genomic() method
    if strand == "+":
        # Second CDS coordinate always show an off-set of +1 for "+" strands
        cds_1 = cds_1
        cds_2 -= 1

    elif strand == "-":
        # First CDS coordinate always show an off-set of -1 for "-" strands
        cds_1 += 1
        cds_2 = cds_2

    else:
        sys.exit(f'Strand error for transcript {trans}')

    fixed_cds = (cds_1, cds_2)

    return fixed_cds


def fix_exon_offset(cds_pair, exon_list, trans_strand, trans):

    cds_pair = sorted(cds_pair)

    # I have observed that there is a constant error to be fix for a subset of transcripts:
    # The STOP codon show an offset of one exon (+/- 1 nucl). The START codon seems to always be correct.
    # For "+" strand transcripts, the correct STOP codon is the Last-coord of the Previous-exon
    # For "-" strand transcripts, the correct STOP codon is the First-coord of the Next-exon

    flat = lambda l: [e for sub in l for e in sub]

    # Check that both CDS fall inside at least one of the exon coordinates
    is_inside_exon_list = [True if exon[0] <= cds <= exon[-1] else False for cds in cds_pair for exon in exon_list]
    n_trues = [e for e in is_inside_exon_list if e is True]
    if len(n_trues) > 1:
        # No correction necessary
        return cds_pair

    else:
        if trans_strand == "+":
            start_cds, stop_cds = cds_pair
            exons = sorted(flat(exon_list) + [stop_cds])
            cds_ix = exons.index(stop_cds)
            fixed_stop_cds = exons[cds_ix-1]
            cds_pair = (start_cds, fixed_stop_cds)

        elif trans_strand == "-":
            stop_cds, start_cds = cds_pair
            exons = sorted(flat(exon_list) + [stop_cds])
            cds_ix = exons.index(stop_cds)
            fixed_stop_cds = exons[cds_ix+1]
            cds_pair = (fixed_stop_cds, start_cds)

        else:
            sys.exit(f'Strand error for transcript {trans}')

    return cds_pair


def split_cds(cds_pair, exon_list, trans_id):

    cds_start, cds_end = sorted(cds_pair)

    # In case of mono-exonic transcripts return directly the CDS section
    if len(exon_list) < 2:
        return [(cds_start, cds_end)]

    # In cases where the CDS Start/End is located inside the same exon return directly the CDS section
    for (exon_st, exon_end) in exon_list:
        if exon_st <= cds_start <= exon_end and exon_st <= cds_end <= exon_end:
            return [(cds_start, cds_end)]

    # Check that the CDS is located inside at least one of the exon-coordinates
    is_inside_exon_list = [True if exon[0] <= cds <= exon[-1] else False for cds in cds_pair for exon in exon_list]
    n_trues = [e for e in is_inside_exon_list if e is True]
    if len(n_trues) > 1:
        pass
    else:
        verbose = True
        if verbose:
            print(f"Transcript {trans_id} CDS coordinates ({cds_start}-{cds_end})  is not within any of exon ranges")
        return None

    # First, get the exons that are inside the CDS range
    cds_exons = []
    for exon in sorted(exon_list):
        exon_start, exon_end = exon

        # Capture the exons that are inside the CDS range
        if cds_start <= exon_start <= cds_end or cds_start <= exon_end <= cds_end:
            cds_exons.append((exon_start, exon_end))

        # The above check will fail for CDS that are entirely contained in a single exon
        if exon_start <= cds_start <= exon_end and exon_start <= cds_end <= exon_end:
            cds_exons.append((exon_start, exon_end))

    # Second, remove the first and last exons, and substitute them with their respective CDS coordinates
    try:
        exon_first, exon_last = cds_exons[0], cds_exons[-1]
    except IndexError:
        return None

    cds_first = (cds_start, exon_first[-1])
    cds_last = (exon_last[0], cds_end)

    cds_list = [cds_first] + cds_exons[1:-1] + [cds_last]

    return cds_list


def find_transcripts_orf_information(gtf_file, sequences_dt, gtf_obj, outfolder):

    orf_index_filename = os.path.splitext(os.path.basename(gtf_file))[0] + "_ORF_index.json"
    orf_index_file = os.path.join(outfolder, orf_index_filename)

    print(time.asctime(), f"Generating ORF index file: {orf_index_file}")

    trans_relative_cds_dt, orf_data_dt = (defaultdict(list) for _ in range(2))
    for i, (trans, trans_seq) in enumerate(sequences_dt.items()):

        try:
            trans_strand = gtf_obj.trans_sense_dt[trans]
        except KeyError:
            continue

        if trans_strand not in {"+", "-"}:
            print(f"Unstranded transcript {trans} ({trans_strand}). Ignoring it.")
            continue

        # We want all the possible ORFs, thus min_prot_len=0. We ignore unknown-strand trans; thus strand_known=True
        orf_info = get_orf_info(trans_seq, strand_known=True, min_protein_length=0)

        # The output of orf_info is already sorted, with the first element being the one with the longest ORF
        # orf_data structure: (orf_start, orf_end, frame, strand, cds_start, cds_end, cds_seq, pep_seq)

        # Program fails/break here when then FASTA sequences don't match the models in the GTF file
        # For example because the user provided FASTA sequences representing the CDS regions instead of the exons
        # TODO find a way to handle this error

        longest_orf_data = orf_info[0]
        orf_data_dt[trans] = longest_orf_data

        # Generate an ORF Index file. This file is necessary for TransFeat to explore the alternative ORFs
        # Sort by starting coordinate of ORF
        orf_info = sorted(orf_info, key=lambda x: x[0])
        for orf_data in orf_info:
            try:
                *_, cds_start_relative, cds_end_relative, _, _ = orf_data
            except ValueError:
                print(f"Transcript {trans} not found in the annotation file.")
                continue

            cds_pair = (cds_start_relative, cds_end_relative)

            if -1 not in cds_pair:
                trans_relative_cds_dt[trans].append(cds_pair)

    with open(orf_index_file, "w+") as fh:
        json.dump(trans_relative_cds_dt, fh, sort_keys=True)

    return orf_data_dt, orf_index_file


def assign_longest_orf_as_cds(gtf_obj, sequences_dt, orf_data_dt, cds_th):

    print(time.asctime(), "Assigning longest ORF as putative coding region (CDS)")

    # Dict to save the new CDS coordinates
    trans_cds_dt = {}

    # Look up table to convert from relative transcript coordinates to genomic coordinates
    exon_lookup_table = create_coordinate_lookup_table(gtf_obj)

    cds_not_found_trans, short_cds_trans = (set() for _ in range(2))
    for i, (trans, trans_seq) in enumerate(sequences_dt.items()):

        if trans not in gtf_obj.trans_exons_dt.keys():
            continue

        # Skip transcripts which already have CDS coordinates
        if gtf_obj.trans_cds_dt[trans]:
            continue

        # Skip transcripts from unknown strands
        trans_strand = gtf_obj.trans_sense_dt[trans]
        if trans_strand not in {"+", "-"}:
            continue

        trans_lookup_table = exon_lookup_table[trans]

        # Get longest CDS
        try:
            # orf_data structure: (orf_start, orf_end, frame, strand, cds_start, cds_end, cds_seq, pep_seq)
            cds_start_relative, cds_end_relative = orf_data_dt[trans][4], orf_data_dt[trans][5]
        except KeyError:
            cds_not_found_trans.add(trans)
            continue

        if -1 in {cds_start_relative, cds_end_relative}:
            cds_not_found_trans.add(trans)
            continue

        # Convert the CDS coordinates, which are relative only to the transcript start-end, to their genomic coordinates
        cds_genomic_st, relative_exon_st, genomic_exon_st = \
            convert_to_genomic_coord(cds_start_relative, trans_lookup_table, trans_strand, trans)

        cds_genomic_end, relative_exon_end, genomic_exon_end = \
            convert_to_genomic_coord(cds_end_relative, trans_lookup_table, trans_strand, trans)

        try:
            cds_pair = (min({cds_genomic_st, cds_genomic_end}), max({cds_genomic_st, cds_genomic_end}))
        except Exception as err:
            print(err)

            print(trans, trans_strand, gtf_obj.trans_exons_dt[trans], cds_genomic_st, cds_genomic_end)
            print(orf_data_dt[trans])

            sys.exit(f"CDS coordinates not found for transcript {trans}")

        cds_pair = fix_cds_offset(cds_pair, trans_strand, trans)
        cds_pair = fix_exon_offset(cds_pair, gtf_obj.trans_exons_dt[trans], trans_strand, trans)

        if not cds_pair:
            print(f"Error while fixing off-set for transcript {trans}")
            continue

        check_codons(cds_pair, gtf_obj.trans_exons_dt[trans], trans_seq, trans_strand, trans)

        cds_list = split_cds(cds_pair, gtf_obj.trans_exons_dt[trans], trans)

        if cds_list is None:
            print(f"Error while slicing the Start/End CDS coordinates for transcript {trans}; Ignoring it.")
            continue

        cds_len = sum([max(cds)-min(cds)+1 for cds in cds_list])
        if cds_len < cds_th:
            short_cds_trans.add(trans)
            continue

        trans_cds_dt[trans] = cds_list

    return trans_cds_dt, cds_not_found_trans, short_cds_trans
