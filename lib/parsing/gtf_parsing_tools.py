import os
import time
import linecache
from collections import defaultdict
from lib.parsing.gtf_object_tools import create_gtf_object, get_id


def get_valid_gtf_line(gtf_file):

    with open(gtf_file) as fh:
        for i, line in enumerate(fh):
            try:
                seqname, source, feature, start, end, score, strand, frame, attr = line.split('\t')
            except ValueError:
                print(f'File "{gtf_file}" row "{i}" does not contain a valid number of fields: "{line}"')
                continue

            if feature.upper() == "EXON" and strand in {"+", "-"}:
                yield line


def filter_gtf_file(gtf_file):

    print(time.asctime(), f'Generating exon-only GTF annotation from file: {gtf_file}')

    outfile = os.path.splitext(gtf_file)[0] + '_exons.gtf'
    with open(outfile, "w+") as fh:
        for line in get_valid_gtf_line(gtf_file):
            fh.write(line)

    return outfile


def add_features_to_gtf(gtf_file):

    gtf_obj = create_gtf_object(gtf_file)
    trans_cds_dt = gtf_obj.trans_cds_dt
    trans_5utr_dt = gtf_obj.trans_5utr_dt
    trans_3utr_dt = gtf_obj.trans_3utr_dt
    trans_start_codon_dt = gtf_obj.trans_start_codon
    trans_stop_codon_dt = gtf_obj.trans_stop_codon

    trans_gene_coords_dt = {}
    for gene, trans_list in gtf_obj.gene_trans_dt.items():
        gene_coords = gtf_obj.gene_coords_dt[gene]
        for trans in trans_list:
            trans_gene_coords_dt[trans] = gene_coords

    transcripts_lines_dt = defaultdict(list)
    for trans, lines_ix_list in gtf_obj.trans_gtf_lines_index.items():

        line_1 = linecache.getline(gtf_obj.gtf_path, lines_ix_list[0])
        seqname, source, _, _, _, score, strand, frame, attr = line_1.strip('\n').split('\t')

        gene_coords = trans_gene_coords_dt[trans]
        start, end = gene_coords[0], gene_coords[-1]

        g_row = f'{seqname}\t{source}\tgene\t{start}\t{end}\t"."\t{strand}\t{frame}\t{attr}'
        t_row = f'{seqname}\t{source}\ttranscript\t{start}\t{end}\t"."\t{strand}\t{frame}\t{attr}'

        for line in [g_row, t_row]:
            if line not in transcripts_lines_dt[trans]:
                # Transcripts visualization on IGV is better without these lines; thus I disable it for the moment
                # transcripts_lines_dt[trans].append(line.strip('\n'))
                continue

        for line_ix in lines_ix_list:
            line = linecache.getline(gtf_obj.gtf_path, line_ix)

            # Important! Ignore any line that is not an exon coordinate so as the CDS re-annotation is completely new
            # Other features will be re-added by the function write_gtf_with_features further downstream
            _, _, line_feature, *_ = line.split('\t')

            if line_feature != "exon":
                continue

            if line not in transcripts_lines_dt[trans]:
                transcripts_lines_dt[trans].append(line.strip('\n'))

        feature_dicts_list = [trans_cds_dt, trans_5utr_dt, trans_3utr_dt, trans_start_codon_dt, trans_stop_codon_dt]
        feature_tags_list = ["CDS", "five_prime_utr", "three_prime_utr", "start_codon", "stop_codon"]
        for feature_dt, feature_tag in zip(feature_dicts_list, feature_tags_list):
            score = "."
            try:
                coord_list = feature_dt[trans]
                # The value of "start_codon", "stop_codon" is a tuple; thus, it must be converted to a list
                if feature_tag in {"start_codon", "stop_codon"}:
                    coord_list = [coord_list]

                for coord in coord_list:
                    start, end = sorted(coord)
                    line = f"{seqname}\t{source}\t{feature_tag}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attr}"
                    if line not in transcripts_lines_dt[trans]:
                        transcripts_lines_dt[trans].append(line.strip('\n'))
            except Exception:
                continue

    gtf_lines = []
    for trans, trans_lines in transcripts_lines_dt.items():
        for line in trans_lines:
            gtf_lines.append(line.strip('\n'))

    sorted_lines = sorted(gtf_lines, key=lambda l: (get_id(l, 'gene_id'), get_id(l, 'transcript_id'),
                                                    int(l.split('\t')[3])))

    # Write output file, overwrite input file
    with open(gtf_file, "w+") as fh:
        for line in sorted_lines:
            fh.write(line+'\n')

    return gtf_file


def annotate_cds_into_gtf(gtf_obj, trans_cds_dt, outfile):

    print(time.asctime(), f'Re-annotating CDS coordinates into file: {outfile}')

    gtf_path = gtf_obj.gtf_path

    new_gtf_lines = defaultdict(list)

    # First, Upload the transcripts for which a CDS was found
    trans_with_cds = set()
    for trans, trans_cds in trans_cds_dt.items():

        trans_with_cds.add(trans)

        # GTF line format: seqname, source, feature, start, end, score, strand, frame, attr
        # Pick a representative GTF line from for the transcript ID to fill the other fields
        trans_line_ix = gtf_obj.trans_gtf_lines_index[trans][0]
        trans_gtf_line = linecache.getline(gtf_path, trans_line_ix)
        seqname, source, _, _, _, score, _, frame, attr = trans_gtf_line.strip('\n').split('\t')

        feature_cds = "CDS"
        score_cds = "."

        strand = gtf_obj.trans_sense_dt[trans]

        # Create a new, updated, Trans_id to GTF_line dict
        for trans_line_ix in gtf_obj.trans_gtf_lines_index[trans]:

            # Important! Ignore any line that is not an exon coordinate so as the CDS re-annotation is completely new
            # Other features will be re-added by the function write_gtf_with_features further downstream
            # GTF file row format: seqname, source, feature, start, end, score, strand, frame, attr
            trans_gtf_line = linecache.getline(gtf_path, trans_line_ix)

            # Do NOT rewrite the feature variabl "line_feature", it cause conflict with same variable on top!
            _, _, line_feature, *_ = trans_gtf_line.split('\t')

            if line_feature.upper() == "EXON":
                new_gtf_lines[trans].append(trans_gtf_line)

        if trans_cds:
            for cds in trans_cds:
                start, end = sorted(cds)
                new_line = [seqname, source, feature_cds, str(start), str(end), score_cds, strand, frame, attr]
                new_line = "\t".join(new_line)+"\n"
                new_gtf_lines[trans].append(new_line)

    cds_not_found_trans = set(gtf_obj.trans_exons_dt.keys()) - trans_with_cds

    # Second, upload the lines of transcript for which a CDS could not be found
    for trans in sorted(cds_not_found_trans):
        for trans_line_ix in gtf_obj.trans_gtf_lines_index[trans]:
            trans_gtf_line = linecache.getline(gtf_path, trans_line_ix)
            new_gtf_lines[trans].append(trans_gtf_line)

    # # Third, upload the lines for the transcript for which a CDS was already known
    # for trans in sorted(trans_with_cds):
    #     for line_obj in gtf_obj.trans_gtf_lines_dt[trans]:
    #         new_gtf_lines[trans].append(line_obj)

    missing_trans = set(gtf_obj.trans_gene_dt.keys()).symmetric_difference(set(new_gtf_lines.keys()))
    for trans in missing_trans:
        for trans_line_ix in gtf_obj.trans_gtf_lines_index[trans]:
            trans_gtf_line = linecache.getline(gtf_path, trans_line_ix)
            new_gtf_lines[trans].append(trans_gtf_line)

    # To conclude, sort all the lines and write output
    all_lines = []
    for trans, trans_lines in new_gtf_lines.items():
        for line in trans_lines:
            # GTF line format: seqname, source, feature, start, end, score, strand, frame, attr
            seqname, source, feature, start, end, score, strand, frame, attr = line.split('\t')

            gene = get_id(line, 'gene_id')
            trans = get_id(line, 'transcript_id')

            sort_tag = f"{seqname}-{gene}-{trans}-{feature}-{start}"
            all_lines.append((line, sort_tag))

    # Write Output
    with open(outfile, "w+") as fh:
        for (line, sort_tag) in sorted(all_lines, key=lambda x: x[1]):
            fh.write(line)

    return outfile
