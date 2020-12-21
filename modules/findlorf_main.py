import time
import warnings

from lib.findlorf.findlorf_tools import *
from lib.parsing.gtf_object_tools import create_gtf_object
from lib.parsing.fasta_parsing_tools import write_fasta_file, get_fasta_sequences
from lib.parsing.gtf_parsing_tools import filter_gtf_file, add_features_to_gtf, annotate_cds_into_gtf

warnings.simplefilter("ignore")


def findlorf_main(gtf_file, fasta_file, outpath, outname, cds_th=50, filter_gtf=True):

    print("\n")
    print(time.asctime(), "Starting FindLORF analysis")

    # In case the user pass the name with a file extension, remove it
    if outname.endswith(".gtf"):
        outname = outname.replace(".gtf", "")
    outname += "_transfind"

    # Create output folder
    outfolder = os.path.join(outpath, outname)
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)

    # Condver CDS length from AA to bp
    cds_th = cds_th * 3

    # Generate a GTF file only with "valid" fields (only rows containing "exon" coord, and with known (+ or -) strand)
    if filter_gtf:
        gtf_file_filtered = filter_gtf_file(gtf_file)
    else:
        gtf_file_filtered = gtf_file

    # Get information from annotation file
    gtf_obj = create_gtf_object(gtf_file_filtered)

    # Get transcripts nucleotide sequence
    sequences_dt = get_fasta_sequences(fasta_file)

    # Get ORF information of the transcripts
    orf_data_dt, orf_index_file = find_transcripts_orf_information(gtf_file, sequences_dt, gtf_obj, outfolder)

    trans_cds_dt, cds_not_found_trans, short_cds_trans = \
        assign_longest_orf_as_cds(gtf_obj, sequences_dt, orf_data_dt, cds_th)

    print(time.asctime(), "Writing re-annotated transcriptome annotation")
    outfile = os.path.join(outfolder, f"{outname}.gtf")

    # Annotate the identified/selected CDS into the output file
    outfile = annotate_cds_into_gtf(gtf_obj, trans_cds_dt, outfile)

    print(time.asctime(), "Generating headers for output fasta files")
    trans_header_dt, cds_seq_dt, pep_seq_dt = [{} for _ in range(3)]
    # orf_data structure: (orf_start, orf_end, frame, strand, cds_start, cds_end, cds_seq, pep_seq)
    for trans, orf_data in orf_data_dt.items():
        trans_cds_seq = orf_data[6]

        # Ignore transcripts without sequences
        if not trans_cds_seq:
            continue

        cds_seq_dt[trans] = str(trans_cds_seq)
        pep_seq_dt[trans] = str(trans_cds_seq.translate())

        trans_chrom = gtf_obj.trans_chrom_dt[trans][:-1]

        try:
            cds_start, cds_end = trans_cds_dt[trans][0][0], trans_cds_dt[trans][-1][-1]
            trans_header = f">{trans} | {gtf_obj.trans_gene_dt[trans]} | {trans_chrom}:{cds_start}-{cds_end}"
        except KeyError:
            trans_exons = gtf_obj.trans_exons_dt[trans]
            exon_start, exon_end = trans_exons[0][0], trans_exons[-1][-1]
            trans_header = f">{trans} | {gtf_obj.trans_gene_dt[trans]} | {trans_chrom}:{exon_start}-{exon_end}"

        trans_header_dt[trans] = trans_header

    nucl_fasta = os.path.join(outfolder, f"{outname}_nuc.fasta")
    write_fasta_file(cds_seq_dt, nucl_fasta, trans_header_dt)

    pep_fasta = os.path.join(outfolder, f"{outname}_pep.fasta")
    write_fasta_file(pep_seq_dt, pep_fasta, trans_header_dt)

    print(time.asctime(), "Writing output-related tables")
    # Identify transcripts in the GTF file that are not present in the FASTA file
    trans_seq_absent, trans_seq_absent_lines = (set() for _ in range(2))
    fasta_trans = set(sequences_dt.keys())
    for trans in gtf_obj.trans_exons_dt.keys():
        if trans not in fasta_trans:
            trans_seq_absent.add(trans)
            line = f"{gtf_obj.trans_chrom_dt[trans][:-1]}\t{gtf_obj.trans_sense_dt[trans]}\t" \
                   f"{gtf_obj.trans_gene_dt[trans]}\t{trans}\n"
            trans_seq_absent_lines.add(line)

    if trans_seq_absent_lines:
        absent_outfile = os.path.join(outfolder, outname + "_sequence_not_found.csv")
        with open(absent_outfile, "w+")as fh:
            fh.write(f"Chromosome,Strand,Gene_ID,Transcript_ID\n")
            for line in sorted(trans_seq_absent_lines):
                fh.write(line)

    # Keep track of removed transcripts due to short CDS
    if short_cds_trans:
        short_cds_table = os.path.join(outfolder, outname + "_short_CDS_transcripts.csv")
        with open(short_cds_table, "w+") as fh:
            fh.write(f"Chromosome,Strand,Gene_ID,Transcript_ID\n")
            for trans in sorted(short_cds_trans):
                line = f"{gtf_obj.trans_chrom_dt[trans][:-1]},{gtf_obj.trans_sense_dt[trans]},{gtf_obj.trans_gene_dt[trans]},{trans}\n"
                fh.write(line)

    # Keep track of the transcript for which a CDS was not found
    if cds_not_found_trans:
        no_cds_table = os.path.join(outfolder, outname + "_no_CDS_transcripts.csv")
        with open(no_cds_table, "w+") as fh:
            fh.write(f"Chromosome,Strand,Gene_ID,Transcript_ID\n")
            for trans in sorted(cds_not_found_trans):
                line = f"{gtf_obj.trans_chrom_dt[trans][:-1]},{gtf_obj.trans_sense_dt[trans]},{gtf_obj.trans_gene_dt[trans]},{trans}\n"
                fh.write(line)

    # Return output file for TransAll function
    return outfile, orf_index_file
