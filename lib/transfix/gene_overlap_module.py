import os
import sys
import time
from collections import defaultdict


def flat(l):

    return [e for sub in l for e in sub]


def get_start_end(model):

    return (model[0][0], model[-1][-1])


def get_exon_model(gtf_obj, t_id):

    return gtf_obj.trans_exons_dt[t_id]


def get_intron_model(gtf_obj, t_id, force=False):

    # Mono-exonic transcripts does not contain introns
    if len(gtf_obj.trans_exons_dt[t_id]) == 1:
        # For mono-exons, if no intron is available, we may be interest on having the exon coordinates instead
        if force:
            return gtf_obj.trans_exons_dt[t_id]
        else:
            return None
    else:
        return gtf_obj.trans_introns_dt[t_id]


def get_cds_model(gtf_obj, t_id, force=False):

    try:
        trans_cds = gtf_obj.trans_cds_dt[t_id]
    except KeyError:
        trans_cds = []

    if not trans_cds:
        if force:
            return gtf_obj.trans_exons_dt[t_id]
        else:
            return None
    else:
        return trans_cds


def get_overlap_percentage(range_A, range_B):

    overlap_A = range(max(range_A[0], range_B[0]), min(range_A[1], range_B[1]) + 1)
    try:
        overlap_perc = len(overlap_A)/max(range_A[-1] - range_A[0], range_B[-1] - range_B[0])
    # In case the two ranges doesn't present any overlap
    except ZeroDivisionError:
        return 0.0

    return overlap_perc


def sort_by_start(gtf_obj, t_group):

    return sorted(t_group, key=lambda t_id: (get_start_end(get_exon_model(gtf_obj, t_id))[0], t_id))


def refine_overlap_group(gtf_obj, overlap_group, overlap_th=0.25):

    if not 0.0 < overlap_th < 1.0:
        sys.exit(f"The overlap percentage tolerance threshold must be within 0.0 and 1.0, not '{overlap_th}'.")

    # The method would return an empty list if the initial overlap_group contains only one transcript. Thus:
    if len(overlap_group) == 1:
        return [overlap_group]

    refined_groups = []

    # Identify transcripts that present only a small overlap in their boundary exons and put them into separate groups
    sorted_group = sorted(overlap_group, key=lambda t_id: get_start_end(get_exon_model(gtf_obj, t_id)))

    # Initialize a tracking set
    temp = set()

    for prev_trans, next_trans in zip(sorted_group[:-1], sorted_group[1:]):
        # Get boundary exons of the transcripts to compare if:
        # "LAST exon of the PREVIOUS transcript" present a PARTIAL overlap with the "FIRST exon of the NEXT transcript"

        # Previous transcript information
        prev_exons = get_exon_model(gtf_obj, prev_trans)
        _, prev_last_exon = prev_exons[0], prev_exons[-1]

        # Next transcript information
        next_exons = get_exon_model(gtf_obj, next_trans)
        next_first_exon, _ = next_exons[0], next_exons[-1]

        # Add the previous transcript to the group
        temp.add(prev_trans)

        # Compare partial overlap, not using get_overlap_percentage method as we want to check for PARTIAL overlap
        partial_overlap = False

        # Important! The requirements of a partial overlap, which should be put into separate overal groups, are:
        # 1st) The overlap must be partial:
        if next_first_exon[0] <= prev_last_exon[1] <= next_first_exon[1] and prev_last_exon[0] < next_first_exon[0]:

            prev_exon_coverage = get_overlap_percentage(prev_last_exon, next_first_exon)
            next_exon_coverage = get_overlap_percentage(prev_last_exon, next_first_exon)

            # 2nd) The overlap must be very small:
            # Thus, if there is a partial overlap, these values must exist (< 0) but be below our required threshold
            if 0 < prev_exon_coverage < overlap_th and 0 < next_exon_coverage < overlap_th:
                partial_overlap = True

        # If there is a partial overlap, put the next transcript in a different overlap subgroup
        if not partial_overlap:
            temp.add(next_trans)
        else:
            # Saved the overlapping subgroup
            refined_groups.append(sorted(temp))

            # Re-initialize the tracking subgroup with the latest non-overlapping (according to previous rules) trans_id
            temp = set([next_trans])

    # Add the last group to the list
    refined_groups.append(sorted(temp))

    return refined_groups


def group_transcripts_by_overlap(gtf_obj, t_group, feature="exon", strict=True, overlap_th=0.25):

    t_group = sort_by_start(gtf_obj, t_group)

    result = []
    overlap_group, prev_max_end = [], 0
    for i, trans_id in enumerate(t_group):

        if feature.upper() == "EXON":
            trans_start, trans_end = get_start_end(get_exon_model(gtf_obj, trans_id))
        elif feature.upper() == "INTRON":
            # For mono-exonic transcripts force=True will return the exon coordinates instead
            trans_start, trans_end = get_start_end(get_intron_model(gtf_obj, trans_id, force=True))
        elif feature.upper() == "CDS":
            trans_start, trans_end = get_start_end(get_cds_model(gtf_obj, trans_id, force=True))
        else:
            sys.exit(f'Overlapping feature must be "exon", "intron", or "CDS", not "{feature}".')

        if trans_start <= prev_max_end:
            overlap_group.append(trans_id)

            # Choose largest end:
            if trans_end > prev_max_end:
                prev_max_end = trans_end

        else:
            # First group is empty
            if overlap_group:
                result.append(overlap_group)
            prev_max_end = trans_end

            overlap_group = []
            overlap_group.append(trans_id)

        # Necessary to append the last group
        if i+1 == len(t_group):
            result.append(overlap_group)

    if not strict:
        refined_results = []
        for overlap_group in result:
            refined_overlap_group = refine_overlap_group(gtf_obj, overlap_group, overlap_th=overlap_th)
            refined_results.extend(refined_overlap_group)

        result = refined_results

    if len(t_group) != len(flat(result)):
        sys.exit(f'Error. Number of overlapping groups does not match for the group of transcripts:{flat(t_group)}.')

    return result


def group_transcripts_across_strands(gtf_obj, to_analyze=None):

    # TODO this function is not used anywhere, check if it could be useful in the future, otherwise delete

    chrom_transcripts = defaultdict(set)
    for chrom, strand_transcripts in sorted(gtf_obj.chrom_trans_dt.items()):

        # Use un-stranded chromosome/scaffold as key
        if not chrom[-1] in {"+", "-", "."}:
            print(f'WARNING: Chromosome/Scaffold "{chrom}" does not finish with strand tag. Please check.')

        chrom_key = chrom[:-1]
        chrom_transcripts[chrom_key].update(strand_transcripts)

    # Assign a Gene ID number to each overlapping group
    across_strand_overlapping_dt = defaultdict(set)
    for chrom_key, chrom_trans in chrom_transcripts.items():

        # Group transcripts by their overlap, use refine=True to avoid grouping transcripts with just small overlaps
        overlap_transcripts = group_transcripts_by_overlap(gtf_obj, chrom_trans, strict=False, overlap_th=0.25)

        for overlap_group in overlap_transcripts:
            # No need to check groups with a single transcript
            if len(overlap_group) > 1:

                # Check if the group contain transcript of interest
                contain_relevant_trans = False
                for t_id in overlap_group:
                    if to_analyze:
                        if t_id in to_analyze:
                            contain_relevant_trans = True
                    else:
                        contain_relevant_trans = True

                if not contain_relevant_trans:
                    continue

                # Get the leftmost and rightmost coordinates
                left_boundary, right_boundary = None, None
                for t_id in sorted(overlap_group):
                    t_left_boundary = gtf_obj.trans_exons_dt[t_id][0][0]
                    t_right_boundary = gtf_obj.trans_exons_dt[t_id][-1][-1]

                    if not left_boundary or left_boundary > t_left_boundary:
                        left_boundary = t_left_boundary

                    if not right_boundary or right_boundary < t_right_boundary:
                        right_boundary = t_right_boundary

                assert left_boundary and right_boundary

                group_key = f"{chrom_key}:{left_boundary}-{right_boundary}"

                across_strand_overlapping_dt[group_key].update(overlap_group)

    for _, overlap_group in across_strand_overlapping_dt.items():
        yield overlap_group
