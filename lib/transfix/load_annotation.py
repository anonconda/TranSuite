# Functions to read GTF and GFF3 files and extract details to populate the Gene and Transcript models in gene_models.py
import os
from collections import defaultdict
from lib.transfix.gene_models import Gene, Transcript
from lib.parsing.gtf_object_tools import get_id


def add_gff3_feature(locus_model, line_fields):
    '''
    Called by get_models_gff3()
    instatiantiates Transcript() and adds
    features to the model
    '''

    feature = line_fields[8]

    if 'retrotransposon' in feature:
        locus_model.transposon = True

    # Transcript_ID is indicated by ID= when using the internal method gtf_to_gff3
    get_transcript_id = lambda s: s.split(";GID=")[0].split("ID=")[-1]

    transcript_ID_list = []
    transcript_ID_list.append(get_transcript_id(feature))

    for transcript_ID in transcript_ID_list:
        if transcript_ID not in locus_model.transcript_dict:
            transcript_model = Transcript(transcript_ID)
            transcript_model.sense = line_fields[6]
            locus_model.transcript_dict.update({transcript_ID: transcript_model})
        transcript_model = locus_model.transcript_dict[transcript_ID]

        if line_fields[2].upper() == 'MRNA':
            transcript_ID = get_transcript_id(feature)
            mRNA_bounds = [line_fields[3], line_fields[4]]
            transcript_model.mRNA = mRNA_bounds

        if line_fields[2].upper() == 'CDS':
            cds_coords = [int(line_fields[3]), int(line_fields[4])]
            transcript_model.add_CDS(cds_coords)

        if line_fields[2].upper() == 'EXON':
            exon_coords = [int(line_fields[3]), int(line_fields[4])]
            transcript_model.add_exon(exon_coords)
            locus_model.transcript_dict.update({transcript_ID: transcript_model})

    return locus_model


def get_models_gff3(filename):
    '''
    Opens GFF3 file with <filename>, reads and
    instantiates the Gene model at the first occurrence of
    each locus ID. Additional locus features are extracted
    by calling add_feature_gff3()
    '''

    infile = open(filename, 'rU')
    gff_lines = infile.readlines()
    infile.close()

    locus_dict = {}
    for line in gff_lines:
        line = line.strip()
        line_fields = line.split('\t')

        if line.startswith('#'):
            continue

        # Second index allows for 'genes' and 'transposable_element_gene' etc
        # Code doesn't work if GFF doesn't contain a "gene" line at the start of each new gene block
        if line_fields[2].upper() == 'GENE':
            # The locus_model (Gene ID) is always indicated by "GID=" due to using the internal method gtf_to_gff3
            locus_ID = line_fields[8].split('GID=')[-1].replace("\n", "")

            # assert locus_ID not in locus_dict

            if locus_ID not in locus_dict:
                gene_model = Gene(locus_ID)
                gene_coords = [int(line_fields[3]), int(line_fields[4])]
                gene_model.add_bounds(gene_coords)
                gene_model.add_sense(line_fields[6])
                locus_dict.update({locus_ID: gene_model})

        else:
            locus_model = locus_dict[locus_ID]
            locus_model = add_gff3_feature(locus_model, line_fields)

    return locus_dict


def get_models_gtf(gtf_filename):
    '''Build gene models from gtf file'''
    infile = open(gtf_filename, 'rU')
    gff_lines = infile.readlines()
    infile.close()

    locus_dict = {}
    for line in gff_lines:
        line = line.strip()
        if line.startswith('#'):
            continue

        line_fields = line.split('\t')

        transcript_ID = get_id(line, 'transcript_id')
        locus_ID = get_id(line, 'gene_id')

        chrom = line_fields[0]

        if locus_ID not in locus_dict:
            gene_model = Gene(locus_ID)
            gene_model.add_sense(line_fields[6])
            locus_dict.update({locus_ID: gene_model})

        else:
            gene_model = locus_dict[locus_ID]

        if transcript_ID not in gene_model.transcript_dict:
            transcript_model = Transcript(transcript_ID)
            transcript_model.sense = line_fields[6]
            gene_model.transcript_dict.update({transcript_ID: transcript_model})

            transcript_model.gene = locus_ID
            transcript_model.chrom = chrom

        transcript_model = gene_model.transcript_dict[transcript_ID]

        if line_fields[2].upper() == 'CDS':  # and proteinCoding == True:
            cds_coords = [int(line_fields[3]), int(line_fields[4])]
            transcript_model.add_CDS(cds_coords)

        if line_fields[2].upper() == 'EXON':  # and proteinCoding == True:
            exon_coords = [int(line_fields[3]), int(line_fields[4])]
            transcript_model.add_exon(exon_coords)

        gene_model.transcript_dict.update({transcript_ID: transcript_model})
        locus_dict.update({locus_ID: gene_model})

    return locus_dict


def get_fasta_sequence(fasta):

    transcript_seq_dt = defaultdict(str)
    with open(fasta) as fh:
        row_id, seq = None, ""
        for row in fh:
            if row.startswith(">"):
                seq = ""
                row_id = row.replace(" ", "#").replace("\t", "#").split("#")[0].strip(">")
            else:
                seq = row.strip("\n")
                transcript_seq_dt[row_id] += seq

    return transcript_seq_dt


def add_gtf_seqs(fasta, locus_dict):

    transcript_seq_dt = get_fasta_sequence(fasta)
    for gene in locus_dict:
        gene_model = locus_dict[gene]
        for transcript in gene_model.transcript_dict:
            try:
                transcript_model = gene_model.transcript_dict[transcript]
                seq = transcript_seq_dt[transcript]
                transcript_model.seq = seq
            except Exception as e:
                print(e)
                pass

    return locus_dict


def define_gff3_genes(gtf):

    # GFF/GTF conversion and differences article: http://blog.nextgenetics.net/?e=27

    exons_dt, gene_features_dt, gene_coord_dt = (defaultdict(list) for _ in range(3))
    with open(gtf) as fh:
        for ix, line in enumerate(fh):
            # Skip comment lines
            if line[0] != '#':
                data = line.strip().split('\t')

                trans_id = get_id(line, 'transcript_id')
                gene_id = get_id(line, 'gene_id')

                chrom, source, feat_type, exon_st, exon_end, score, strand, frame, *_ = data

                # Replace last column with a GFF formatted attributes columns
                # Added a GID attribute to conserve all the GTF data
                data[-1] = "ID=" + trans_id + ";GID=" + gene_id

                exons_dt[gene_id].extend([int(exon_st), int(exon_end)])
                gene_features_dt[gene_id] = "\t".join(
                    [chrom, source, "gene", "EXON_ST", "EXON_END", score, strand, frame, data[-1]])

    for gene in gene_features_dt:
        gene_line = gene_features_dt[gene]\
                        .replace("EXON_ST", str(min(exons_dt[gene])))\
                        .replace("EXON_END", str(max(exons_dt[gene]))) + "\n"
        gene_coord_dt[gene] = gene_line

    return gene_coord_dt


def gtf_to_gff3(gtf):

    gff3_lines = defaultdict(list)
    with open(gtf) as fh:
        for ix, line in enumerate(fh):
            # Skip comment lines
            if line[0] != '#':
                data = line.strip().split('\t')

                trans_id = get_id(line, 'transcript_id')
                gene_id = get_id(line, 'gene_id')

                chrom, source, feat_type, exon_st, exon_end, score, strand, frame, *_ = data

                # Replace last column with a GFF formatted attributes columns
                # Added a GID attribute to conserve all the GTF data
                data[-1] = "ID=" + trans_id + ";GID=" + gene_id

                # New GFF line
                gff3_lines[gene_id].append('\t'.join(data) + "\n")

    gene_coord_dt = define_gff3_genes(gtf)

    # Necessary to add a gene feature line to the GFF file in the correct order
    new_lines = []
    for gene_id in sorted(gene_coord_dt.keys()):
        new_lines.append(gene_coord_dt[gene_id])
        for ln in sorted(gff3_lines[gene_id]):
            new_lines.append(ln)

    outfile = os.path.splitext(gtf)[0] + ".gff3"
    with open(outfile, "w+") as fh:
        for line in new_lines:
            fh.write(line)

    return outfile


def filter_gtf(gtf, processed_transcripts, rejected_start_codons, iter_n, outfolder):

    lines = []
    with open(gtf) as fh:
        for ln in fh:
            # GTF line format:
            seqname, source, feature, start, end, score, strand, frame, attr = ln.split("\t")
            trans_id = ln.split("transcript_id \"")[-1].split("\"")[0]

            # Ignore transcripts translated in the current cycle
            if trans_id in processed_transcripts:
                continue

            if feature != "exon" and feature != "CDS":
                continue

            # Remove rejected CDS
            if feature.upper() == "CDS":
                if int(start) in rejected_start_codons or int(end) in rejected_start_codons:
                    continue

            lines.append(ln)

    # outname = os.path.splitext(gtf)[0].split(".iter.")[0] + f".iter.{iter_n}.gtf"
    # outfile = os.path.join(outfolder, outname)
    # print("TEST 1: ", outfile)
    # print("TEST 2: ", outname)

    outfile = os.path.splitext(gtf)[0].split(".iter.")[0] + f".iter.{iter_n}.gtf"
    with open(outfile, "w+") as fh:
        for ln in lines:
            fh.write(ln)

    return outfile
