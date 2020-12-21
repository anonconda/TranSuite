import os
import sys
import time
import linecache
from bisect import insort
from collections import defaultdict, namedtuple


def get_overlap_percentage(range_A, range_B):

    overlap_A = range(max(range_A[0], range_B[0]), min(range_A[1], range_B[1]) + 1)
    try:
        overlap_perc = len(overlap_A)/max(range_A[-1] - range_A[0], range_B[-1] - range_B[0])
    # In case the two ranges doesn't present any overlap
    except ZeroDivisionError:
        return 0.0

    return overlap_perc


def find_utr_regions(trans_id, trans_sense, trans_exons, trans_cds):

    if not trans_cds:
        trans_5utr, trans_3utr, start_codon, stop_codon = [], [], tuple(), tuple()
        return trans_5utr, trans_3utr, start_codon, stop_codon

    if trans_sense not in {"+", "-"}:
        sys.exit(f'Strand error for transcript "{trans_id}"')

    flat = lambda l: [e for sub in l for e in sub]
    group_by_pair = lambda l: [e for e in zip(l[::2], l[1::2])]

    trans_exons = flat(trans_exons)
    trans_exons_set = set(trans_exons)
    trans_cds = flat(trans_cds)

    if trans_sense == "+":
        # 5' UTR for FORWARD strands
        start = trans_cds[0]
        exons_subset = [exon for exon in trans_exons if exon < start]
        if exons_subset:
            if start in trans_exons_set:
                trans_5utr = exons_subset
            else:
                trans_5utr = exons_subset + [start-1]
        else:
            trans_5utr = []

        if not len(trans_5utr) % 2 == 0:
            trans_5utr = trans_5utr + [start-1]

        # 3' UTR for FORWARD strands
        stop = trans_cds[-1]
        exons_subset = [exon for exon in trans_exons if exon > stop]
        if exons_subset:
            if stop in trans_exons_set:
                trans_3utr = exons_subset
            else:
                trans_3utr = [stop+1] + exons_subset
        else:
            trans_3utr = []

        if not len(trans_3utr) % 2 == 0:
            trans_3utr = [stop+1] + trans_3utr

        start_codon = (start, start + 2)
        stop_codon = (stop - 2, stop)

    elif trans_sense == "-":
        # 5' UTR for REVERSE strands
        start = trans_cds[-1]
        exons_subset = [exon for exon in trans_exons if exon > start]
        if exons_subset:
            if start in trans_exons_set:
                trans_5utr = exons_subset
            else:
                trans_5utr = [start+1] + exons_subset
        else:
            trans_5utr = []

        if not len(trans_5utr) % 2 == 0:
            trans_5utr = [start+1] + trans_5utr

        # 3' UTR for REVERSE strands
        stop = trans_cds[0]
        exons_subset = [exon for exon in trans_exons if exon < stop]
        if exons_subset:
            if stop in trans_exons_set:
                trans_3utr = exons_subset
            else:
                trans_3utr = exons_subset + [stop-1]
        else:
            trans_3utr = []

        if not len(trans_3utr) % 2 == 0:
            trans_3utr = trans_3utr + [stop-1]

        start_codon = (start - 2, start)
        stop_codon = (stop, stop + 2)

    else:
        print(f'Transcript "{trans_id}" does not present a valid strand ("{trans_sense}").')
        return [], [], (), ()

    # TODO, answer these questions:
    #  Is it possible that a Start/Stop codon are formed after splicing, and thus its coordinates belong to two
    #  different exons? If so, should I check that all 3 bp of codon are contained inside the same exon?
    #  There is a function in transfix_tools file, is_inside_an_exon(), that could do this

    trans_5utr = group_by_pair(sorted(trans_5utr))
    trans_3utr = group_by_pair(sorted(trans_3utr))

    return trans_5utr, trans_3utr, start_codon, stop_codon


def get_id(gtf_line, tag):

    # Warning messages can introduce too much noise in the output, thus we silence them for now
    print_warning = True
    if tag not in gtf_line:
        if print_warning:
            print(f'WARNING: No instance of tag "{tag}" found on line: "{gtf_line}"')
        return None

    splitted = gtf_line.strip('\n').split(tag)

    # Some lines may contain multiple instances of the "splitting tag"; ex: it contains both "gene_id" AND "ref_gene_id"
    # Generally, in the GTF files the desired tag is the first one ("gene_id" come always before "ref_gene_id")
    if len(splitted) != 2:
        if print_warning:
            print(f'WARNING: Multiple matches of tag "{tag}" found on line: "{gtf_line}"')
        pass
    res_id = splitted[1].split('";')[0]

    while res_id.startswith(" ") or res_id.startswith('"') or res_id.startswith("'"):
        res_id = res_id[1:]

    while res_id.endswith(" ") or res_id.endswith('"') or res_id.endswith("'"):
        res_id = res_id[:-1]

    return res_id


def add_tag(gtf_line, tag=None):

    # If no tag is given, return the gtf_line unchanged
    # This is necessary when using the "split_redundancy_removal" method to avoid redundant re-tagging of GTF lines
    if not tag:
        return gtf_line

    if 'gene_id "' not in gtf_line or 'transcript_id "' not in gtf_line:
        return gtf_line

    # Add "_" separator to tag if not present
    if not tag.startswith("_"):
        tag = f"_{tag}"

    # Avoid re-annotating a line that has already been tagged
    try:
        if tag in gtf_line:
            return gtf_line
    except TypeError:
        # Previous check make this redundant, but better be safe than sorry
        sys.exit("Tag cannot be 'None'.")

    # Important, some lines have a "ref_gene_id" tag, this screw up the gene_id, this is handled by get_id method
    gene_id = get_id(gtf_line, 'gene_id "')
    trans_id = get_id(gtf_line, 'transcript_id "')

    # Important! Replace only the first instance found! This is necessary to don't re-annotate "ref_gene_id" IDs
    temp_line = gtf_line.replace(f'gene_id "{gene_id}', f'gene_id "{gene_id + tag}', 1)
    tagged_line = temp_line.replace(f'transcript_id "{trans_id}', f'transcript_id "{trans_id + tag}', 1)

    return tagged_line


def parse_gtf(gtf_file, tag=None):

    with open(gtf_file) as fh:
        for i, line in enumerate(fh):

            # Some assemblies/annotations (Ex: stringtie, TAIR10) start with headers, ignore these lines
            if not line or line.startswith("#"):
                continue

            # If tag is not given (None), add_tag() return the line unchanged
            line = add_tag(line, tag)

            if 'gene_id "' not in line or 'transcript_id "' not in line:
                print(f'File row "{i}" does not contain gene/transcript ID: "{line}"')
                continue

            # GTF line format: seqname, source, feature, start, end, score, strand, frame, attr
            gtf_row = line.strip('\n').split('\t')
            if len(gtf_row) != 9:
                print(f'File "{gtf_file}" row "{i}" does not contain a valid number of fields: "{line}"')
                continue

            yield gtf_row


def create_gtf_object(gtf_file, to_keep=None, tag=None):

    print(time.asctime(), f'Uploading information from annotation file: {gtf_file}')

    flat = lambda l: [item for sublist in l for item in sublist]

    # Structure of the GTF object
    GTF = namedtuple('GTF',
                     'gtf_path '
                     'chrom_gene_dt chrom_trans_dt '
                     'gene_coords_dt gene_trans_dt gene_name_dt '
                     'trans_chrom_dt trans_gene_dt trans_sense_dt '
                     'trans_exons_dt trans_introns_dt trans_cds_dt trans_5utr_dt trans_3utr_dt '
                     'trans_start_codon trans_stop_codon '
                     'trans_gtf_lines_index')

    # Initialize data structures to store relational-information of interest
    chrom_gene_dt, chrom_trans_dt, gene_trans_dt = [defaultdict(set) for _ in range(3)]

    trans_exons_dt, trans_introns_dt, trans_cds_dt, trans_5utr_dt, trans_3utr_dt, \
    trans_gtf_lines_index = [defaultdict(list) for _ in range(6)]

    gene_name_dt = defaultdict(str)

    gene_coords_dt, trans_chrom_dt, trans_gene_dt, \
    trans_sense_dt, trans_start_codon_dt, trans_stop_codon_dt = [{} for _ in range(6)]

    # This create a Transcript_Id to row position in a file, so as to recover the relevant lines when necessary
    # Important! As the row number must represent the number in the original file, this check must be done here
    # as further ignore invalid rows (and thus the row number would be skipped)
    with open(gtf_file) as fh:
        for line_ix, line in enumerate(fh, 1):

            # If tag is not given (None), add_tag() return the line unchanged
            line = add_tag(line, tag)

            if 'transcript_id "' not in line:
                continue
            else:
                trans_id = get_id(line, 'transcript_id "')
                trans_gtf_lines_index[trans_id].append(line_ix)

    # As it is a generator, it is not convenient to put the print statement inside the method
    print(time.asctime(), f'Parsing annotation file: {gtf_file}')

    # Line indexing must start from 1 for the linecache.getline() function
    for line_obj in parse_gtf(gtf_file, tag):
        # line_obj (GTF row) is a list with the following format:
        seqname, source, feature, start, end, score, strand, frame, attr = line_obj

        gene_id = get_id(attr, 'gene_id "')
        trans_id = get_id(attr, 'transcript_id "')

        # If there is a subset to keep, ignore the others
        if to_keep:
            if trans_id not in to_keep:
                continue

        # Uniquely ID each strand by pairing it to its scaffold. Ex: Chrom01+, Chrom01-, Chrom02+,  etc
        genomic_strand = f'{seqname}{strand}'

        # Get Gene Name and gene description (Note) if available
        if 'gene_name "' in attr:
            gene_name = get_id(attr, 'gene_name "')
            gene_name_dt[gene_id] = gene_name
        if 'Note "' in attr:
            gene_note = get_id(attr, 'Note "')
            gene_name_dt[gene_id] += " " + gene_note

        gene_trans_dt[gene_id].add(trans_id)
        chrom_gene_dt[genomic_strand].add(gene_id)
        chrom_trans_dt[genomic_strand].add(trans_id)

        trans_gene_dt[trans_id] = gene_id
        trans_sense_dt[trans_id] = strand
        trans_chrom_dt[trans_id] = genomic_strand

        # Ensure that start/ end are in order (start < end)
        start, end = int(start), int(end)
        start, end = min(start, end), max(start, end)
        coord_pair = (start, end)

        # Many functions downstream assume these list of coordinates are sorted. Thus the use of insort to add elements
        if feature.upper() == 'EXON':
            # Important! Newly assembled GTF may contain repeated lines/rows. Thus, check if info is already present!
            if coord_pair not in trans_exons_dt[trans_id]:
                insort(trans_exons_dt[trans_id], coord_pair)

        if feature.upper() == 'CDS':
            if coord_pair not in trans_cds_dt[trans_id]:
                insort(trans_cds_dt[trans_id], coord_pair)

        if feature.upper() == 'FIVE_PRIME_UTR':
            if coord_pair not in trans_5utr_dt[trans_id]:
                insort(trans_5utr_dt[trans_id], coord_pair)

        if feature.upper() == 'THREE_PRIME_UTR':
            if coord_pair not in trans_3utr_dt[trans_id]:
                insort(trans_3utr_dt[trans_id], coord_pair)

        if feature.upper() == 'START_CODON':
            trans_start_codon_dt[trans_id] = coord_pair

        if feature.upper() == 'STOP_CODON':
            trans_stop_codon_dt[trans_id] = coord_pair

    # Important! Some annotations may report only CDS coordinates of a transcript, and doesn't report them also as exons
    # For the program to work, it is necessary to identify these coordinates both as CDS and exons!
    for trans_id in trans_gene_dt.keys():
        if not trans_exons_dt[trans_id] and trans_cds_dt[trans_id]:
            trans_exons_dt[trans_id] = trans_cds_dt[trans_id]

    # Generate missing feature information (Coordinates for: introns, 5'-3' UTR regions, start/stop codons)
    get_introns = lambda exons: [(ex1[-1]+1, ex2[0]-1) for (ex1, ex2) in zip(exons[:-1], exons[1:])]

    for trans_id, trans_exons in trans_exons_dt.items():
        trans_sense = trans_sense_dt[trans_id]

        # Generate Intron information
        trans_introns = get_introns(trans_exons)
        trans_introns_dt[trans_id] = trans_introns

        # Generate UTR information if possible
        trans_cds = trans_cds_dt[trans_id]
        if trans_cds:
            trans_5utr, trans_3utr, start_codon, stop_codon = \
                find_utr_regions(trans_id, trans_sense, trans_exons, trans_cds)

            if not trans_5utr_dt[trans_id] and trans_5utr:
                trans_5utr_dt[trans_id] = trans_5utr

            if not trans_3utr_dt[trans_id] and trans_3utr:
                trans_3utr_dt[trans_id] = trans_3utr

            try:
                trans_start = trans_start_codon_dt[trans_id]
            except KeyError:
                if start_codon:
                    trans_start_codon_dt[trans_id] = start_codon
                else:
                    trans_start_codon_dt[trans_id] = tuple()

            try:
                trans_stop = trans_stop_codon_dt[trans_id]
            except KeyError:
                if stop_codon:
                    trans_stop_codon_dt[trans_id] = stop_codon
                else:
                    trans_stop_codon_dt[trans_id] = tuple()

    # It is posssible for incorrectly annotated transcriptomes to have transcript without any annotated exon-coordinates
    # Also, due to mis-assembly, it is possibly to encounter transcripts models with overlapping exons, remove them
    missing_exons, overlap_exons = [set() for _ in range(2)]
    for t_id, t_exons in trans_exons_dt.items():
        if not t_exons:
            missing_exons.add(t_id)
        else:
            for cur_exon, next_exon in zip(t_exons, t_exons[1:]):
                if get_overlap_percentage(cur_exon, next_exon):
                    overlap_exons.add(t_id)

    if overlap_exons:
        print(f'WARNING! The annotation file "{os.path.basename(gtf_file)}" contains "{len(overlap_exons)}" '
              f'transcripts without annotated exon-coordinates. These models will be ignored.')

        over_exons = os.path.join(os.path.dirname(gtf_file), os.path.basename(gtf_file).replace(".gtf", "_overlapping_exons.csv"))
        with open(over_exons, "a+") as fh:
            fh.write(f"Location,Transcript_ID\n")
            for t_id in sorted(overlap_exons):
                t_exons = flat(trans_exons_dt[t_id])
                t_start, t_end = min(t_exons), max(t_exons)
                t_chrom = trans_chrom_dt[t_id][:-1]
                t_loc = f"{t_chrom}:{t_start}-{t_end}"
                fh.write(f"{t_loc},{t_id}\n")

    if missing_exons:
        print(f'WARNING! The annotation file "{os.path.basename(gtf_file)}" contains "{len(missing_exons)}" '
              f'transcripts without annotated exon-coordinates. These models will be ignored.')

        no_exons = os.path.join(os.path.dirname(gtf_file), os.path.basename(gtf_file).replace(".gtf", "_missing_exons.csv"))
        with open(no_exons, "a+") as fh:
            fh.write(f"Transcript_ID\n")
            for t_id in sorted(missing_exons):
                fh.write(f"{t_id}\n")

    to_remove = missing_exons | overlap_exons

    # Annotate Genes start/end Coordinates
    for gene, trans_list in gene_trans_dt.items():
        starts, ends = (set() for _ in range(2))
        for t_id in trans_list:
            trans_exons = flat(trans_exons_dt[t_id])

            # If there are transcripts without exons, this method will run again ignoring those models,
            # thus we skip this error on the current iteration to allow this method to continue on the first iteration
            if not trans_exons:
                continue

            starts.add(trans_exons[0])
            ends.add(trans_exons[-1])

        if starts and ends:
            gene_coords_dt[gene] = (min(starts), max(ends))

    # Create a GTF object with the parsed and generated information
    gtf_object = GTF(gtf_path=gtf_file,
                     chrom_gene_dt=chrom_gene_dt, chrom_trans_dt=chrom_trans_dt,
                     gene_coords_dt=gene_coords_dt, gene_trans_dt=gene_trans_dt, gene_name_dt=gene_name_dt,
                     trans_chrom_dt=trans_chrom_dt, trans_gene_dt=trans_gene_dt, trans_sense_dt=trans_sense_dt,
                     trans_exons_dt=trans_exons_dt, trans_introns_dt=trans_introns_dt, trans_cds_dt=trans_cds_dt,
                     trans_5utr_dt=trans_5utr_dt, trans_3utr_dt=trans_3utr_dt,
                     trans_start_codon=trans_start_codon_dt, trans_stop_codon=trans_stop_codon_dt,
                     trans_gtf_lines_index=trans_gtf_lines_index)

    # If there are incorrectly annotated transcripts models, run the method again while ignoring incorrect models
    if to_remove:
        to_keep = set(gtf_object.trans_exons_dt.keys()) - to_remove
        gtf_object = create_gtf_object(gtf_file, to_keep=to_keep)

    return gtf_object


def write_gtf(gtf_obj, transcripts, outfolder, outname, w_mode="w+", all_t=False):

    # Create an output folder if it doesnt exist
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)

    outfile = os.path.join(outfolder, outname)

    # Sort transcripts by the key (Chrom, Gene_ID, Trans_ID, Leftmost_coordinate)
    get_sort_key = lambda t_id: (gtf_obj.trans_chrom_dt[t_id], gtf_obj.trans_gene_dt[t_id], t_id,
                                 gtf_obj.trans_exons_dt[t_id][0][0])

    if all_t:
        print(time.asctime(), "Including all transcripts from input annotation file into the output file")
        transcripts = set(gtf_obj.trans_exons_dt.keys())

    sorted_transcripts = sorted(transcripts, key=lambda t_id: get_sort_key(t_id))

    print(time.asctime(), f'Writing output file: {outfile}')
    with open(outfile, w_mode) as fh:
        for trans_id in sorted_transcripts:
            if trans_id in transcripts:
                for line_ix in gtf_obj.trans_gtf_lines_index[trans_id]:
                    line = linecache.getline(gtf_obj.gtf_path, line_ix)
                    fh.write(line)

    return outfile
