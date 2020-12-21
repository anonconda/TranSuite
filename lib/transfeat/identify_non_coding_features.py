import sys
import time
from collections import defaultdict
from lib.parsing.gtf_object_tools import find_utr_regions
from lib.transfeat.transfeat_tools import get_orf_seq


def identify_longer_dorf(gtf_obj, relative_start_dt, orf_dt, trans_seq_dt):

    print(time.asctime(), "Identifying transcripts with long downstream ORFs (ldORFs)")

    # Useful methods to calculate lengths, "get_orf_len" length include introns, "get_cds_len" is length after splicing
    get_cds_len = lambda cds_list: sum([max(cds_pair)-min(cds_pair)+1 for cds_pair in cds_list])
    get_orf_len = lambda orf: orf[1]-orf[0] # + 1 cause the miss-identification of equally longer ORF, thus is disabled

    # Report if there is a downstream ORF which is longer that the current CDS
    ldorf_coord_dt = {}
    longer_dorfs = set()
    for trans, trans_cds in gtf_obj.trans_cds_dt.items():
        if not trans_cds:
            continue

        try:
            trans_orfs = orf_dt[trans]
            trans_orf_start = relative_start_dt[trans]
        except KeyError:
            continue

        trans_dorfs = [orf for orf in trans_orfs if orf[0] > trans_orf_start]
        if trans_dorfs:
            longest_dorf = sorted(trans_dorfs, key=lambda orf: get_orf_len(orf), reverse=True)[0]
        else:
            longest_dorf = None

        if longest_dorf:
            if get_orf_len(longest_dorf) > get_cds_len(trans_cds):
                longer_dorfs.add(trans)
                # Track to annotate them into their own file
                ldorf_coord_dt[trans] = longest_dorf

    # This checks that the ldORFs coordinates are valid by trying extract its ORF sequence from the original sequence
    accepted_ldorf_coord_dt = {}
    for trans, ldORFs_coords in ldorf_coord_dt.items():
        try:
            ldorf_pep_seq, ldorf_nuc_seq = get_orf_seq(ldORFs_coords, trans_seq_dt[trans])
            accepted_ldorf_coord_dt[trans] = ldORFs_coords
        except KeyError:
            pass

    ldorf_coord_dt = accepted_ldorf_coord_dt

    # is_longer_dorf_dt = {}
    is_longer_dorf_dt = defaultdict(lambda: None)
    for trans in gtf_obj.trans_exons_dt.keys():

        if trans in longer_dorfs:
            longer_dorf_flag = True
        else:
            longer_dorf_flag = False

        is_longer_dorf_dt[trans] = longer_dorf_flag

    return is_longer_dorf_dt, ldorf_coord_dt


def identify_uorf(gtf_obj, relative_start_dt, orf_dt, fasta_nuc_dt, uorf_len_th):

    print(time.asctime(), "Identifying transcripts with upstream ORFs (uORFs)")

    # Convert the uORF length (given in a.a.) to nucleotide length
    uorf_len_th = uorf_len_th * 3

    # Filter out uORFs that are too short.

    # Function to check if two coordinates are in frame (default is 3 due to codons being made in groups of 2)
    is_in_frame = lambda n, m: True if (n - m) % 3 == 0 else False

    non_overlapping_uorfs_id, overlapping_uorfs_id = (set() for _ in range(2))
    inframe_uorfs_id, non_inframe_urofs_id = (set() for _ in range(2))

    trans_orf_seq_dt = defaultdict(list)
    for trans in gtf_obj.trans_cds_dt.keys():
        try:
            trans_orf_start = relative_start_dt[trans]
        except KeyError:
            continue

        trans_orf = orf_dt[trans]

        # Classifify the ORF whether they are upstream or downstream of the main start
        trans_uorfs, trans_dorfs = ([] for _ in range(2))
        for orf in trans_orf:
            # Important, by setting it to "<" and NOT "<=" we AVOID to analyze the ORF that originated the authentic ORF
            if orf[0] < trans_orf_start:
                trans_uorfs.append(orf)
            elif orf[0] > trans_orf_start:
                trans_dorfs.append(orf)
            else:
                continue

        selected_uorfs = []
        overlapping_uorf_coord, non_overlapping_uorf_coord = ([] for _ in range(2))
        inframe_uorf_coord, non_inframe_uorf_coord = ([] for _ in range(2))

        if trans_uorfs:
            selected_uorfs = [orf for orf in trans_uorfs if orf[1]-orf[0]+1 >= uorf_len_th]

            if not selected_uorfs:
                continue

            # Identify if uORFs overlaps authentic ATG or not
            for uorf in selected_uorfs:
                if uorf[0] <= trans_orf_start <= uorf[1]:
                    overlapping_uorf_coord.append(uorf)
                else:
                    non_overlapping_uorf_coord.append(uorf)

            # Check if the start-codon of the uORFs is in frame with the authentic start-codon
            for uorf in selected_uorfs:
                # Note: ORFs in index file have always the position: (start ORF, stop ORF)
                uorf_start = uorf[0]
                if is_in_frame(uorf_start, trans_orf_start):
                    inframe_uorf_coord.append(uorf)
                else:
                    non_inframe_uorf_coord.append(uorf)

        if not selected_uorfs:
            continue

        for uorf in selected_uorfs:
            try:
                trans_nucl_seq = fasta_nuc_dt[trans]
                uorf_pep_seq, uorf_nuc_seq = get_orf_seq(uorf, trans_nucl_seq)

            except KeyError:
                uorf_pep_seq, uorf_nuc_seq = "", ""

            uorf_seq_data = (uorf, uorf_pep_seq, uorf_nuc_seq)
            trans_orf_seq_dt[trans].append(uorf_seq_data)

        if overlapping_uorf_coord:
            overlapping_uorfs_id.add(trans)

        if non_overlapping_uorf_coord:
            non_overlapping_uorfs_id.add(trans)

        if inframe_uorf_coord:
            inframe_uorfs_id.add(trans)

        if non_inframe_uorf_coord:
            non_inframe_urofs_id.add(trans)

    # uorf_sets = [overlapping_uorfs_id, non_overlapping_uorfs_id, inframe_uorfs_id, non_inframe_urofs_id]

    urof_categories = defaultdict(set)

    urof_categories["overlapping"].update(overlapping_uorfs_id)
    urof_categories["not_overlapping"].update(overlapping_uorfs_id)
    urof_categories["inframe"].update(overlapping_uorfs_id)
    urof_categories["not_inframe"].update(overlapping_uorfs_id)

    return trans_orf_seq_dt, urof_categories


def is_orf_absent(gtf_obj):

    print(time.asctime(), "Identifying transcripts without an annotated coding region (CDS)")

    # cds_absent_dt = {}
    cds_absent_dt = defaultdict(lambda: None)
    for trans in gtf_obj.trans_exons_dt.keys():
        try:
            trans_cds = gtf_obj.trans_cds_dt[trans]
        except KeyError:
            trans_cds = []

        cds_absent_flag = False
        if not trans_cds:
            cds_absent_flag = True
        cds_absent_dt[trans] = cds_absent_flag

    return cds_absent_dt


def is_ptc(gtf_obj, ptc_len_th, to_ignore=None):

    print(time.asctime(), "Identifying transcripts with Premature-Termination-Codons (PTC)")

    if not to_ignore:
        to_ignore = set()

    # is_ptc_dt = {}
    is_ptc_dt = defaultdict(lambda: None)

    # Useful methods to calculate lengths, "get_cds_len" is length after splicing
    get_cds_len = lambda cds_list: sum([max(cds_pair)-min(cds_pair)+1 for cds_pair in cds_list])

    if not 0 < ptc_len_th < 100:
        sys.exit(f"The PTC length threshold must be within 1% and 100%, not {ptc_len_th}. Aborting.")

    # Group transcripts by their starting CDS
    cds_trans_dt = defaultdict(list)
    for trans, trans_cds in gtf_obj.trans_cds_dt.items():

        if trans in to_ignore:
            continue

        trans_sense = gtf_obj.trans_sense_dt[trans]

        if trans_cds:
            trans_cds_flat = [cds for cds_pair in trans_cds for cds in cds_pair]
            if trans_sense == '+':
                cds_start = min(trans_cds_flat)
            elif trans_sense == '-':
                cds_start = max(trans_cds_flat)
            else:
                print(f"Transcript {trans} sense must be + or - , not {trans_sense}!")
                cds_start = None
        else:
            cds_start = None

        if cds_start:
            cds_trans_dt[cds_start].append(trans)

    trans_cds_length_dt = {}
    for trans, trans_cds in gtf_obj.trans_cds_dt.items():
        trans_cds_len = get_cds_len(trans_cds)
        trans_cds_length_dt[trans] = trans_cds_len

    trans_len_method = False
    if trans_len_method:
        for _, grouped_trans in cds_trans_dt.items():
            group_data = sorted([(trans, trans_cds_length_dt[trans]) for trans in grouped_trans], key=lambda e: e[-1])
            longest_len = group_data[-1][-1]

            for (trans, trans_len) in group_data:
                len_ratio = (trans_len/longest_len)*100

                if len_ratio < ptc_len_th:
                    is_ptc_flag = True
                else:
                    is_ptc_flag = False
                is_ptc_dt[trans] = is_ptc_flag

    cds_dist_method = True
    if cds_dist_method:
        dist = lambda x, y: y-x+1
        perc = lambda x, y: (x/y)*100
        for _, grouped_trans in cds_trans_dt.items():
            group_data = sorted([(trans, trans_cds_length_dt[trans]) for trans in grouped_trans], key=lambda e: e[-1])
            trans_with_longest_cds = group_data[-1][0]

            selected_trans_cds = sorted([coord for cds_pair in gtf_obj.trans_cds_dt[trans_with_longest_cds] for coord in cds_pair])
            selected_trans_min_cds, selected_trans_max_cds = selected_trans_cds[0], selected_trans_cds[-1]

            selected_trans_cds_distance = dist(selected_trans_min_cds, selected_trans_max_cds)

            # Check if the CDS end is below or above the percentile coordinate
            for trans in grouped_trans:
                trx_cds_list = sorted([coord for cds_pair in gtf_obj.trans_cds_dt[trans] for coord in cds_pair])
                trx_cds_min, trx_cds_max = trx_cds_list[0], trx_cds_list[-1]

                trx_dist = dist(trx_cds_min, trx_cds_max)
                trans_cds_len_perc = perc(trx_dist, selected_trans_cds_distance)

                if trans_cds_len_perc < ptc_len_th:
                    is_ptc_flag = True
                else:
                    is_ptc_flag = False

                is_ptc_dt[trans] = is_ptc_flag

    return is_ptc_dt


def is_orf_short(gtf_obj, pep_len_th, to_ignore=None):

    print(time.asctime(), "Identifying transcripts coding for short peptides")

    if not to_ignore:
        to_ignore = set()

    # Convert CDS length (given in a.a.) to nucleotide
    pep_len_th = pep_len_th * 3

    # Useful methods to calculate lengths, "get_cds_len" is length after splicing
    get_cds_len = lambda cds_list: sum([max(cds_pair)-min(cds_pair)+1 for cds_pair in cds_list])

    # ptc_transcripts = set()
    # for trans, ptc_bool in is_ptc_dt.items():
    #     if ptc_bool:
    #         ptc_transcripts.add(trans)

    # is_orf_short_dt = {}
    is_orf_short_dt = defaultdict(lambda: None)
    for trans, trans_cds in gtf_obj.trans_cds_dt.items():

        if trans in to_ignore:
            continue

        if trans_cds:
            cds_len = get_cds_len(trans_cds)

            if cds_len < pep_len_th:
                is_short_flag = True
            else:
                is_short_flag = False
            is_orf_short_dt[trans] = is_short_flag

    return is_orf_short_dt


def is_long_3utr(gtf_obj, utr3_len_th, to_ignore=None):

    if not to_ignore:
        to_ignore = set()

    get_introns = lambda exons: [(ex1[-1]+1, ex2[0]-1) for (ex1, ex2) in zip(exons[:-1], exons[1:])]

    # is_long_3utr_dt = {}
    is_long_3utr_dt = defaultdict(lambda: None)
    for trans_id, trans_cds in gtf_obj.trans_cds_dt.items():

        if trans_id in to_ignore:
            continue

        # Extract 3' UTR region
        trans_exons = gtf_obj.trans_exons_dt[trans_id]
        trans_sense = gtf_obj.trans_sense_dt[trans_id]

        trans_5utr, trans_3utr, start_codon, stop_codon = \
            find_utr_regions(trans_id, trans_sense, trans_exons, trans_cds)

        if trans_3utr:
            introns_3utr = get_introns(trans_3utr)

            if trans_sense == '+':
                introns_starts = [intron[0] for intron in introns_3utr]
            elif trans_sense == '-':
                introns_starts = [intron[1] for intron in introns_3utr]
            else:
                sys.exit(f'Strand Error for Transcript "{trans_id}"')

            if introns_3utr:
                contains_dsj = True

            # a) Long 3' UTR
            utr3_len = sum([max(exon)-min(exon)+1 for exon in trans_3utr])
            long_3utr_flag = False
            if utr3_len >= utr3_len_th:
                long_3utr_flag = True

        else:
            long_3utr_flag = False

        is_long_3utr_dt[trans_id] = long_3utr_flag

    return is_long_3utr_dt


def is_nmd(gtf_obj, auth_stop_dt, sj_dist_th=50, check_ir_event=True,
           ptc_trans=None, ov_uorf_trans=None, uorf_trans=None, long_3utr_trans=None,
           signals_to_check={"DSSJ", "PTC", "OV_UORF" "LONG_3UTR"}, check="ANY"):

    print(time.asctime(), "Identifying transcripts with common Non-Sense-Mediated (NMD) signals")

    # Initialize the sets if they are not pass
    if not ptc_trans:
        ptc_trans = set()

    if not ov_uorf_trans:
        ov_uorf_trans = set()

    if not uorf_trans:
        uorf_trans = set()

    if not long_3utr_trans:
        long_3utr_trans = set()

    # Generate missing feature information (Coordinates for: introns, 5'-3' UTR regions, start/stop codons)
    get_introns = lambda exons: [(ex1[-1]+1, ex2[0]-1) for (ex1, ex2) in zip(exons[:-1], exons[1:])]

    # is_nmd_dt, is_dssj_dt = {}, {}
    is_nmd_dt, is_dssj_dt = [defaultdict(lambda: None) for _ in range(2)]
    for trans_id, trans_cds in gtf_obj.trans_cds_dt.items():
        gene_id = gtf_obj.trans_gene_dt[trans_id]

        # Initialization of NMD features flags to detect
        # The combination of this flags is analyzed further downstream to assign an NMD signal or not

        ir_flag, dssj_flag, ptc_flag = False, False, False
        long_3utr_flag, ov_uorf_flag, uorf_flag = False, False, False


        # Assume that the DS_SJ is False at the start, it will overwrite downstream if it the signal is present
        is_dssj_dt[trans_id] = False

        # 1) Check that the transcript doesn't contain an intron retention (IR)
        try:
            ir_flag = detect_intron_retention(gtf_obj, trans_id)
        except KeyError:
            ir_flag = False

        # 2) Check for the main NMD signal, a SJ distant at least 50 nt from the authentic stop codon
        # 2.1) Extract transcript 3' UTR region
        trans_exons = gtf_obj.trans_exons_dt[trans_id]
        trans_sense = gtf_obj.trans_sense_dt[trans_id]

        trans_5utr, trans_3utr, start_codon, stop_codon = \
            find_utr_regions(trans_id, trans_sense, trans_exons, trans_cds)

        if trans_3utr:
            introns_3utr = get_introns(trans_3utr)

            if trans_sense == '+':
                introns_starts = [intron[0] for intron in introns_3utr]
            elif trans_sense == '-':
                introns_starts = [intron[1] for intron in introns_3utr]
            else:
                sys.exit(f'Strand Error for Transcript "{trans_id}"')

            # 2.2) Check if 3' UTR contains:
            # A) a SJ, B) the distance of the SJ from the authentic stop codon (gene level)
            contains_sj, dssj_dist = False, False
            # A) It contains a SJ in the 3' UTR
            if introns_3utr:
                contains_sj = True

                # B) The SJ is distant at least N nt (default 50 nt) from the authentic stop codon (gene level)
                gene_stop = auth_stop_dt[gene_id]
                for intron_st in introns_starts:
                    dist_from_stop = abs(gene_stop - intron_st)
                    if dist_from_stop >= sj_dist_th:
                        dssj_dist = True
                        is_dssj_dt[trans_id] = True
        else:
            contains_sj = False
            dssj_dist = False

        # Strong signal of NMD: DS SJ distant from stop-codon
        if dssj_dist:
            # is_nmd_dt[trans_id] = True
            dssj_flag = True

        # 3) Check for other possible NMD signals:
        # C) PTC+, D) uORF overlapping start-codon (ouORF), E) long 3' UTR

        # C) PTC signal
        if trans_id in ptc_trans:
            ptc_flag = True

        # D) overlapping uORF signal
        if trans_id in ov_uorf_trans:
            ov_uorf_flag = True

        # E) long 3' UTR signal
        if trans_id in long_3utr_trans:
            long_3utr_flag = True

        # F) Non-overlapping uORF may indicate an NMD signal, but it is not as strong signal os the others.
        if trans_id in uorf_trans:
            uorf_flag = True

        # 4) Check with signals to consider to assign NMD classification
        # {"DSSJ", "PTC", "OV_UORF" "LONG_3UTR"}

        flags_to_check = []
        if "DSSJ" in signals_to_check:
            flags_to_check.append(dssj_flag)

        if "PTC" in signals_to_check:
            flags_to_check.append(ptc_flag)

        if "OV_UORF" in signals_to_check:
            flags_to_check.append(ov_uorf_flag)

        if "UORF" in signals_to_check:
            flags_to_check.append(uorf_flag)

        if "LONG_3UTR" in signals_to_check:
            flags_to_check.append(long_3utr_flag)

        # If there is an IR event, and the specify it as a (negative) signal, don't assign NMD signal and exit loop
        if check_ir_event and ir_flag is True:
            is_nmd_dt[trans_id] = False
            continue

        # PTC is a fundamental characteristics to assign the NMD tag
        if ptc_flag is False:
            is_nmd_dt[trans_id] = False
            continue

        else:
            if check == "ANY":
                nmd_flag = any(flags_to_check)
            elif check == "ALL":
                nmd_flag = all(flags_to_check)
            else:
                sys.exit(f'NMD "check" argument must be "ANY" or "ALL", not "{check}"')

            if nmd_flag is True:
                is_nmd_dt[trans_id] = True
            else:
                is_nmd_dt[trans_id] = False

    return is_nmd_dt, is_dssj_dt


def detect_intron_retention(gtf_obj, trans_id):

    # 1) Get CDS of the transcript under analysis to detect intron retention (trans_id)
    # 2) Get the (all) introns of all the other transcripts under the same gene
    # 3) Check if any intron is completely cover by any of the trans_id exons

    ir_flag = False

    # Get CDS region of trans_id under analysis
    try:
        trans_cds = gtf_obj.trans_cds_dt[trans_id]
    except KeyError:
        trans_cds = []

    if not trans_cds:
        return ir_flag

    # Get all the other transcripts under the same gene
    gene = gtf_obj.trans_gene_dt[trans_id]
    gene_transcripts = [trx for trx in gtf_obj.gene_trans_dt[gene] if trx != trans_id]

    if not gene_transcripts:
        return ir_flag

    # Get (all) the introns of all the other transcripts in the same gene
    pooled_introns = set()
    for other_trans in gene_transcripts:
        try:
            other_introns = gtf_obj.trans_introns_dt[other_trans]
        except KeyError:
            other_introns = []

        pooled_introns.update(other_introns)

    # Check if the intron is completely contained within an exon
    ir_found = False
    for exon in trans_cds:
        if not ir_found:
            for intron in pooled_introns:
                if exon[0] <= intron[0] <= exon[1] and exon[0] <= intron[1] <= exon[1]:
                    ir_found = True
                    ir_flag = True
                    break
        else:
            break

    return ir_flag


def generate_nmd_features_lines(gtf_obj, is_nmd_dt, is_ptc_dt, is_dssj_dt, is_long_3utr_dt, urof_categories):

    # Track the identified NMD features to annotate them into the report table
    nmd_features_dt = defaultdict(str)
    for trans_id, _ in gtf_obj.trans_gene_dt.items():
        gene_id = gtf_obj.trans_gene_dt[trans_id]

        # Track the NMD features of:
        # A) NMD transcripts; and
        # B) Unproductive (PTC) / Non-coding transcripts with uORF

        # A) NMD transcripts flag
        try:
            nmd_flag = is_nmd_dt[trans_id]
        except KeyError:
            nmd_flag = False

        # B) Unproductive (PTC) / Non-coding transcripts with uORF
        noncoding_uorf_flag = False  # Disabled for now (Unproductive tag is identified much later downstream)

        # If neither cases are present, don't report NMD features
        if nmd_flag is False and noncoding_uorf_flag is False:
            nmd_features_dt[trans_id] += "-"

        if nmd_features_dt[trans_id].startswith("-"):
            continue

        # Otherwise, track NMD features:
        if detect_intron_retention(gtf_obj, trans_id):
            nmd_features_dt[trans_id] += ";IR"

        if is_ptc_dt[trans_id] is True:
            # nmd_features_dt[trans_id] += ";PTC+"
            nmd_features_dt[trans_id] += ";PTC"

        if is_dssj_dt[trans_id] is True:
            nmd_features_dt[trans_id] += ";ds_SJ"

        if is_long_3utr_dt[trans_id] is True:
            nmd_features_dt[trans_id] += ";long_3UTR"

        if trans_id in urof_categories["overlapping"]:
            nmd_features_dt[trans_id] += ";ouORF"

        if trans_id in urof_categories["not_overlapping"]:
            nmd_features_dt[trans_id] += ";uORF"

        # Currently the program do not track to which uORF the inframe/not_inframe refers to, just that it exist
        if trans_id in urof_categories["inframe"] and ";inframe" not in nmd_features_dt[trans_id]:
            nmd_features_dt[trans_id] += ";inframe"

        if trans_id in urof_categories["not_inframe"] and ";not_inframe" not in nmd_features_dt[trans_id]:
            nmd_features_dt[trans_id] += ";not_inframe"

        # Clean the NMD feature lines
        if nmd_features_dt[trans_id].startswith(";"):
            nmd_features_dt[trans_id] = nmd_features_dt[trans_id][1:]

        while ";;" in nmd_features_dt[trans_id]:
            nmd_features_dt[trans_id] = nmd_features_dt[trans_id].replace(";;", ";")

    return nmd_features_dt
