import os
import sys
import time
from collections import defaultdict
from itertools import combinations

from lib.transfeat.transfeat_tools import *


def is_as_in_utr(gtf_obj, trans_group):

    # AS_in_UTR_dt, AS_UTR_location_dt = ({} for _ in range(2))
    AS_in_UTR_dt, AS_UTR_location_dt = defaultdict(lambda: None), {}
    AS_in_UTR_transcripts = set()
    for trans_A, trans_B in combinations(trans_group, 2):

        trans_A_exons, trans_A_CDS = set(gtf_obj.trans_exons_dt[trans_A]), set(gtf_obj.trans_cds_dt[trans_A])
        trans_B_exons, trans_B_CDS = set(gtf_obj.trans_exons_dt[trans_B]), set(gtf_obj.trans_cds_dt[trans_B])

        if not trans_A_CDS or not trans_B_CDS:
            verbose = False
            if verbose:
                print(f"Either Transcript {trans_A} or {trans_B} does not have an annotated CDS. Skipping.")
            continue

        AS_UTR_flag = False
        if trans_A_CDS == trans_B_CDS and trans_A_exons != trans_B_exons:
            AS_UTR_flag = True

        # Two transcripts can be AS between them, but not with a third one. Thus, feature flag may be re-annotated.
        # This avoid the re-annotation of AS flags
        for trans in [trans_A, trans_B]:
            if trans not in AS_in_UTR_transcripts:
                AS_in_UTR_dt[trans] = AS_UTR_flag

        if AS_UTR_flag is True:
            AS_in_UTR_transcripts.add(trans_A)
            AS_in_UTR_transcripts.add(trans_B)

            # Identify in which UTR
            UTR_type = detect_AS_location(trans_A, trans_B, gtf_obj)
            AS_UTR_location_dt[trans_A] = UTR_type
            AS_UTR_location_dt[trans_B] = UTR_type

    return AS_in_UTR_dt, AS_UTR_location_dt


def detect_AS_location(trans_A, trans_B, gtf_obj):

    # This method assume the transcript under analysis present and AS in the UTR (determined by is_AS_in_UTR() method)
    AS_UTR_type = ""

    flat = lambda l: [e for sub in l for e in sub]

    # convert_to_introns = lambda exons: [(ex1[-1]+1, ex2[0]-1) for (ex1, ex2) in zip(exons[:-1], exons[1:])]

    trans_A_cds = sorted(flat(gtf_obj.trans_cds_dt[trans_A]))
    trans_A_min_cds, trans_A_max_cds = trans_A_cds[0], trans_A_cds[-1]
    try:
        trans_A_introns = flat(gtf_obj.trans_introns_dt[trans_A])
    except KeyError:
        trans_A_introns = []

    trans_A_introns_left_utr, trans_A_introns_right_utr = ([] for _ in range(2))
    if trans_A_introns:
        for intron in trans_A_introns:
            if intron < trans_A_min_cds:
                trans_A_introns_left_utr.append(intron)
            elif intron > trans_A_max_cds:
                trans_A_introns_right_utr.append(intron)
            else:
                pass

    trans_B_cds = sorted(flat(gtf_obj.trans_cds_dt[trans_B]))
    trans_B_min_cds, trans_B_max_cds = trans_B_cds[0], trans_B_cds[-1]
    try:
        trans_B_introns = flat(gtf_obj.trans_introns_dt[trans_B])
    except KeyError:
        trans_B_introns = []

    trans_B_introns_left_utr, trans_B_introns_right_utr = ([] for _ in range(2))
    if trans_B_introns:
        for intron in trans_B_introns:
            if intron < trans_B_min_cds:
                trans_B_introns_left_utr.append(intron)
            elif intron > trans_B_max_cds:
                trans_B_introns_right_utr.append(intron)
            else:
                pass

    trans_A_sense = gtf_obj.trans_sense_dt[trans_A]
    trans_B_sense = gtf_obj.trans_sense_dt[trans_B]
    if trans_A_sense == trans_B_sense:
        trans_sense = trans_A_sense
    else:
        sys.exit("Strand error during AS in UTR position identification. Aborting")

    if trans_sense != "+" and trans_sense != "-":
        sys.exit("Strand error during AS in UTR position identification. Aborting")


    if trans_A_introns_left_utr != trans_B_introns_left_utr:
        if trans_sense == "+":
            AS_UTR_type += "AS_5UTR"
        elif trans_sense == "-":
            AS_UTR_type += "AS_3UTR"
        else:
            pass

    if trans_A_introns_right_utr != trans_B_introns_right_utr:
        if trans_sense == "+":
            AS_UTR_type += "AS_3UTR"
        elif trans_sense == "-":
            AS_UTR_type += "AS_5UTR"
        else:
            pass

    # Add separator if necessary
    if AS_UTR_type:
        AS_UTR_type = AS_UTR_type.replace("UTRAS", "UTR;AS")

    return AS_UTR_type


def is_NAGNAG(gtf_obj, trans_group):

    flat = lambda l: [e for sub in l for e in sub]

    get_cds_len = lambda cds_list: sum([max(cds_pair)-min(cds_pair)+1 for cds_pair in cds_list])

    # NAGNAG_dt = {}
    NAGNAG_dt = defaultdict(lambda: None)
    NAGNAG_transcripts = set()
    for trans_A, trans_B in combinations(trans_group, 2):

        pooled_exons = set(flat(gtf_obj.trans_exons_dt[trans_A])) | set(flat(gtf_obj.trans_exons_dt[trans_B]))
        pooled_cds_set = set(flat(gtf_obj.trans_cds_dt[trans_A])) | set(flat(gtf_obj.trans_cds_dt[trans_B]))

        try:
            trans_A_len = get_cds_len(gtf_obj.trans_cds_dt[trans_A])
        except KeyError:
            trans_A_len = 0

        try:
            trans_B_len = get_cds_len(gtf_obj.trans_cds_dt[trans_B])
        except KeyError:
            trans_B_len = 0

        # trans_A_len, trans_B_len = len(trans_seq_dt[trans_A]), len(trans_seq_dt[trans_B])

        # Keep only CDS coordinates that are also Exon coordinates
        cds_coords_to_test = pooled_cds_set.intersection(pooled_exons)
        pooled_coords = sorted(cds_coords_to_test)

        nagnag_flag = False
        for i in range(len(pooled_coords)-1):
            coord_1 = pooled_coords[i]
            coord_2 = pooled_coords[i+1]

            if abs(coord_1 - coord_2) == 3:

                # # Also check if the translations differ by only by 1 AA
                # if abs(trans_A_len - trans_B_len) == 1:

                # USing nucleotide sequences instead of peptide sequences
                if abs(trans_A_len - trans_B_len) == 3:
                    nagnag_flag = True
                    break

        # Avoid re-annotation of features flags
        for trans in [trans_A, trans_B]:
            if trans not in NAGNAG_transcripts:
                NAGNAG_dt[trans] = nagnag_flag

        if nagnag_flag is True:
            NAGNAG_transcripts.add(trans_A)
            NAGNAG_transcripts.add(trans_B)

    return NAGNAG_dt


def identify_similar_coding_features(gtf_obj):

    print(time.asctime(), "Grouping transcripts by their genomic overlap")
    # The next comparison need to be perform within a transcripts group. Transcripts are grouped by their overlap.
    chrom_strand_transcript_dt = defaultdict(list)
    for chrom_strand, gene_list in gtf_obj.chrom_gene_dt.items():
        for gene_id in gene_list:
            transcritps = gtf_obj.gene_trans_dt[gene_id]
            chrom_strand_transcript_dt[chrom_strand].extend(transcritps)

    # AS identification is done at the loci level, which means working at the chromosome-strand level
    as_in_utr_dt, as_utr_location_dt, nagnag_dt = (defaultdict(lambda: None) for _ in range(3))
    for chrom_strand, strand_trans in chrom_strand_transcript_dt.items():
        overlap_groups = group_transcripts_by_overlap(strand_trans, gtf_obj.trans_exons_dt)

        print(time.asctime(),
              f"Identifying transcripts with AS in UTR region, and transcripts with NAGNAG features ({chrom_strand})")

        for i, trans_group in enumerate(overlap_groups):
            # Check if the CDS regions is the same (therefore, AS is in the UTR)
            group_as_in_utr_dt, group_as_utr_location_dt = is_as_in_utr(gtf_obj, trans_group)
            as_in_utr_dt.update(group_as_in_utr_dt)
            as_utr_location_dt.update(group_as_utr_location_dt)

            # Check if a CDS coordinate differs by only 3 nucl (and therefore, it is likely NAGNAG)
            group_nagnag_dt = is_NAGNAG(gtf_obj, trans_group)
            nagnag_dt.update(group_nagnag_dt)

    return as_in_utr_dt, as_utr_location_dt, nagnag_dt

