import os
import sys
import time
from collections import namedtuple, defaultdict
from Bio.Seq import Seq


from lib.parsing.gtf_object_tools import create_gtf_object
from lib.transfix.gene_overlap_module import get_overlap_percentage, group_transcripts_by_overlap
from lib.findlorf.findlorf_tools import convert_from_genomic_to_relative, split_cds


def write_table(my_list, to_keep, filename, sep=","):

    outfile = open(filename, 'w')
    outfile.writelines(f"Chromosome{sep}Strand{sep}Gene_ID{sep}Transcript_ID\n")
    for item in sorted(my_list):
        # Table row format: chromosome	strand	gene_ID	transcript_ID
        *_, trans_id = item.strip("\n").split(sep)
        if trans_id in to_keep:
            outfile.writelines(item+'\n')
    outfile.close()


def is_inside_an_exon(coordinate, exon_list):

    # TODO check if I can use this method to check the generate start-stop codons in the "gtf_object_tools" file

    # Check that the CDS fall inside at least one of the exon coordinates
    is_inside_exon_list = [True if exon[0] <= coordinate <= exon[-1] else False for exon in exon_list]
    n_trues = [e for e in is_inside_exon_list if e is True]
    if len(n_trues) >= 1:
        return True
    else:
        return False


def find_cds_end(nuc_seq, cds_start, trans_id, trans_sense, trans_exons):

    # Important! The to_stop=True is super duper important to get the correct length, and thus correct CDS end coord
    aa_seq = nuc_seq.translate(to_stop=True)

    res_dt = {}
    flat = lambda l: [e for sub in l for e in sub]

    if trans_sense == "+":
        # In the FASTA file the AA sequence doesn't contain the stop codon character,
        # therefore we sum the missing +3 nucleotides to its length
        aa_len = len(aa_seq)*3 + 3
        max_exon = max(flat(trans_exons))

        aa_ix = 0
        for genomic_ix in range(cds_start, max_exon+1):
            if is_inside_an_exon(genomic_ix, trans_exons):
                if aa_ix < aa_len:
                    res_dt[aa_ix] = genomic_ix
                    aa_ix += 1
                else:
                    break

        max_relative = max(res_dt.keys())
        cds_end = res_dt[max_relative]

        # If CDS end is outside the transcript model, set transcript ends as the CDS end
        if cds_end > max_exon:
            cds_end = max_exon

        # In some cases the CDS end is inside and intron, this is due to the stop codon being created by the
        # nucleotides at the edge of two exons; This code corrects the cds end position for these cases
        if not is_inside_an_exon(cds_end, trans_exons):
            for i in range(0, len(trans_exons)-1):
                exon1_st, exon1_end = trans_exons[i]
                exon2_st, exon2_end = trans_exons[i+1]

                # Check if CDS end fall inside and intron
                if exon1_end < cds_end < exon2_st:
                    # Calculate the amount by which the cds exceed the boundary and calculate the corrected CDS end pos
                    cds_offset = cds_end - exon1_end
                    cds_end = exon1_end + cds_offset
                    break

    elif trans_sense == "-":
        aa_len = len(aa_seq)*3 + 3
        min_exon = min(flat(trans_exons))

        aa_ix = 0
        # For "-" strand transcripts the sequence is read right to left, therefore we decrease the genomic_mix value
        for genomic_ix in range(cds_start, min_exon-1, -1):
            if is_inside_an_exon(genomic_ix, trans_exons):
                if aa_ix < aa_len:
                    res_dt[aa_ix] = genomic_ix
                    aa_ix += 1
                else:
                    break

        max_relative = max(res_dt.keys())
        cds_end = res_dt[max_relative]

        if cds_end < min_exon:
            cds_end = min_exon

        if not is_inside_an_exon(cds_end, trans_exons):
            for i in range(0, len(trans_exons)-1):
                exon1_st, exon1_end = trans_exons[i]
                exon2_st, exon2_end = trans_exons[i+1]

                if exon1_end < cds_end < exon2_st:
                    # Due to the "-" strand transcripts being right to left, the offset calculation have to be modified
                    cds_offset = exon2_st - cds_end
                    cds_end = exon1_end - cds_offset
                    break

    else:
        sys.exit(f"Transcript {trans_id} sense must be either + or -, not {trans_sense}")

    return cds_end


def remove_transcripts_without_cds(gtf_file, outfolder):

    print(time.asctime(), "Removing invalid transcripts models from annotation file")

    gtf_obj = create_gtf_object(gtf_file)
    trans_with_cds, trans_without_cds = (set() for _ in range(2))
    for trans, trans_cds in gtf_obj.trans_cds_dt.items():
        if not trans_cds:
            trans_without_cds.add(trans)
        else:
            trans_with_cds.add(trans)

    gtf_rows = []
    with open(gtf_file) as fh:
        for line in fh:
            row = line.strip('\n').split('\t')
            gtf_rows.append(row)

    filtered_rows = []
    for row in gtf_rows:
        attr = row[-1]
        trans_id = attr.strip("\n").split("transcript_id \"")[-1].split("\";")[0]
        if trans_id in trans_without_cds:
            continue
        else:
            filtered_rows.append(row)

    # Adding the ".transfix.temp." to the name assure that this file will be removed with the other temporary files
    gtf_name = os.path.basename(gtf_file).replace(".gtf", ".transfix.temp.gtf")

    gtf_path = os.path.join(outfolder, gtf_name)
    with open(gtf_path, "w+") as fh:
        for row in gtf_rows:
            line = "\t".join(row)+"\n"
            fh.write(line)

    return gtf_path, trans_with_cds, trans_without_cds


def get_transcript_data_from_gff_obj(gene_id, locus_dict, trans_sequences_dt):

    gene_model = locus_dict[gene_id]
    transcript_dict = gene_model.transcript_dict

    TransData = namedtuple('TransData', 'chrom sense gene id exons start end seq')

    transdata_dt = {}
    for transcript_id in transcript_dict:
        transcript_model = transcript_dict[transcript_id]

        try:
            trans_seq = trans_sequences_dt[transcript_id]
        except KeyError:
            trans_seq = ""

        trans_chrom = transcript_model.chrom
        trans_sense = transcript_model.sense
        trans_gene = transcript_model.gene
        trans_id = transcript_id

        trans_exons = transcript_model.exon_list
        trans_start = transcript_model.exon_list[0][0]
        trans_end = transcript_model.exon_list[-1][-1]

        trans_data = TransData(chrom=trans_chrom, sense=trans_sense, gene=trans_gene, id=trans_id,
                               exons=trans_exons, start=trans_start, end=trans_end, seq=trans_seq)

        transdata_dt[transcript_id] = trans_data

    return transdata_dt


def get_transcript_data_from_gtf_obj(gene_id, gtf_obj, trans_sequences_dt):

    TransData = namedtuple('TransData', 'chrom sense gene id exons start end seq')

    transdata_dt = {}
    for transcript_id in gtf_obj.gene_trans_dt[gene_id]:

        try:
            trans_seq = trans_sequences_dt[transcript_id]
        except KeyError:
            trans_seq = ""

        trans_chrom = gtf_obj.trans_chrom_dt[transcript_id]
        trans_sense = gtf_obj.trans_sense_dt[transcript_id]
        trans_gene = gtf_obj.trans_gene_dt[transcript_id]
        trans_id = transcript_id

        trans_exons = gtf_obj.trans_exons_dt[transcript_id]
        trans_start = trans_exons[0][0]
        trans_end = trans_exons[-1][-1]

        trans_data = TransData(chrom=trans_chrom, sense=trans_sense, gene=trans_gene, id=trans_id,
                               exons=trans_exons, start=trans_start, end=trans_end, seq=trans_seq)

        transdata_dt[transcript_id] = trans_data

    return transdata_dt


def fix_atg_position(transcript_data_dt, atg_pos):

    # Classification categories
    cat_dt = defaultdict(set)
    output_dt = defaultdict(dict)

    trans_cds_dt, trans_cds_seq_dt, trans_header_dt = [{} for _ in range(3)]

    for transcript_id in transcript_data_dt:

        # trans_data is a namedtuple
        trans_data = transcript_data_dt[transcript_id]

        trans_chrom = trans_data.chrom
        trans_sense = trans_data.sense
        trans_gene = trans_data.gene
        trans_id = trans_data.id

        trans_exons = trans_data.exons
        trans_start = trans_data.start
        trans_end = trans_data.end
        trans_seq = trans_data.seq

        if atg_pos is None:
            cat_dt["cds_not_found"].add(transcript_id)
            # cds_not_found.add(transcript_id)

            line = f"{trans_chrom},{trans_sense},{trans_gene},{trans_id}"
            cat_dt["cds_not_found_lines"].add(line)
            # cds_not_found_lines.add(line)

            continue

        if not trans_seq:
            cat_dt["seq_not_present"].add(transcript_id)
            # seq_not_present.add(transcript_id)

            line = f"{trans_chrom},{trans_sense},{trans_gene},{trans_id}"
            cat_dt["seq_not_present_lines"].add(line)

            # seq_not_present_lines.add(line)
            continue

        if not trans_start <= atg_pos <= trans_end:
            cat_dt["atg_not_in_cds"].add(transcript_id)
            # atg_not_in_cds.add(transcript_id)

            line = f"{trans_chrom},{trans_sense},{trans_gene},{trans_id}"
            cat_dt["atg_not_in_cds_lines"].add(line)
            # atg_not_in_cds_lines.add(line)

            continue

        lookup_table = convert_from_genomic_to_relative(trans_exons)

        # Check if the start CDS is contain within an exon
        if atg_pos not in set(lookup_table.keys()):
            cat_dt["atg_not_in_cds"].add(transcript_id)
            # atg_not_in_cds.add(transcript_id)

            line = f"{trans_chrom},{trans_sense},{trans_gene},{trans_id}"
            cat_dt["atg_not_in_cds_lines"].add(line)
            # atg_not_in_cds_lines.add(line)

            continue

        cds_index = lookup_table[atg_pos]

        seq = Seq(trans_seq)

        # Get CDS section of the sequence (from start)
        cds_seq = seq[cds_index:]

        # The sequences in the FASTA file are already reverse_complemented for transcripts in the - strand;
        # That is, sequences in FASTA file are annotated as read by the molecular machinery, right-to-left, and
        # with the complementary nucleotides;
        # However, the cds_index represent the position as read from left-to-right;
        # Therefore, the following sequence transformation are necessary to cut the sequence in the proper posit
        if trans_sense == '-':
            seq = seq.reverse_complement()
            cds_seq = seq[:cds_index+1]
            cds_seq = cds_seq.reverse_complement()

        if cds_seq.startswith('ATG'):
            trans_cds_seq_dt[transcript_id] = cds_seq
            # Translate only transcripts starting with 'ATG'
            peptide = True
        else:
            peptide = False

        if not peptide:
            cat_dt["rejected_start_codons"].add(atg_pos)
            # rejected_start_codons.add(atg_pos)

            cat_dt["start_codon_not_atg"].add(transcript_id)
            # start_codon_not_atg.add(transcript_id)

            line = f"{trans_chrom},{trans_sense},{trans_gene},{trans_id}"
            cat_dt["start_codon_not_atg_lines"].add(line)
            # start_codon_not_atg_lines.add(line)

            continue

        else:
            # Get information to annotate the updated CDS coordinates
            cds_end = find_cds_end(cds_seq, atg_pos, transcript_id, trans_sense, trans_exons)

            # cds_pair order is important for re-annotation (cds_1 < cds_2)
            cds_1, cds_2 = min({atg_pos, cds_end}), max({atg_pos, cds_end})

            cds_pair = (cds_1, cds_2)
            trans_cds_list = split_cds(cds_pair, trans_exons, transcript_id)

            trans_cds_dt[transcript_id] = trans_cds_list

            # Information for the transcripts header in fasta file
            trans_header = f">{transcript_id} | {trans_gene} | {trans_chrom}:{atg_pos}-{cds_end}"
            trans_header_dt[transcript_id] = trans_header

            # Track processed transcripts to ignore them in the next iteration
            cat_dt["processed_transcripts"].add(transcript_id)
            # processed_transcripts.add(transcript_id)

    output_dt["trans_cds_dt"] = trans_cds_dt
    output_dt["trans_cds_seq_dt"] = trans_cds_seq_dt
    output_dt["trans_header_dt"] = trans_header_dt

    return output_dt, cat_dt


def get_gene_longest_CDS(gene_id, gtf_obj):

    get_cds_len = lambda cds_list: sum([max(cds_pair)-min(cds_pair)+1 for cds_pair in cds_list])

    longest_cds, longest_cds_len = [], 0
    for trans_id in gtf_obj.gene_trans_dt[gene_id]:
        try:
            trans_cds = gtf_obj.trans_cds_dt[trans_id]
        except KeyError as err:
            trans_cds = []

        trans_cds_len = get_cds_len(trans_cds)

        if trans_cds_len >= longest_cds_len:
            longest_cds_len = trans_cds_len
            longest_cds = trans_cds

    return longest_cds


def get_gene_groups(gtf_obj, chimeric_genes=None):

    if not chimeric_genes:
        print("WARNING: No set of chimeric genes specified!")
        chimeric_genes = set()

    get_cds_len = lambda cds_list: sum([max(cds_pair)-min(cds_pair)+1 for cds_pair in cds_list])
    get_gene_ids = lambda t_group: set([gtf_obj.trans_gene_dt[t_id] for t_id in t_group])

    gene_groups_dt = {}
    for chrom, strand_transcripts in gtf_obj.chrom_trans_dt.items():

        # Strand is important to identify what is the "first" gene that makes up a chimeric
        chrom_strand = chrom[-1]

        if chrom_strand not in {"+", "-"}:
            print(f'WARNING: Skipping Scaffold "{chrom}" due to unrecognized strand.')

        overlapping_transcripts = group_transcripts_by_overlap(gtf_obj, strand_transcripts)

        for overlap_group in overlapping_transcripts:
            gene_ids = get_gene_ids(overlap_group)

            # Remove chimeric genes from the results (values) but not the keys!
            gene_group = [g_id for g_id in gene_ids if g_id not in chimeric_genes]

            # Sort by length of the CDS
            gene_ids = sorted(gene_ids, key=lambda g_id: get_cds_len(get_gene_longest_CDS(g_id, gtf_obj)))

            for g_id in gene_ids:
                # Remove self-references (Gene ID don't need to be in its group, as it is the dict key
                gene_subgroup = [over_g for over_g in gene_group if over_g != g_id]

                # If the gene is a chimeric one, keep only the genes that overlap with it
                if g_id in chimeric_genes:
                    g_coord = gtf_obj.gene_coords_dt[g_id]

                    accepted_subgroup = []
                    for sub_gene in gene_group:
                        subg_coord = gtf_obj.gene_coords_dt[sub_gene]

                        if get_overlap_percentage(g_coord, subg_coord):
                            accepted_subgroup.append(sub_gene)

                    gene_subgroup = accepted_subgroup

                gene_groups_dt[g_id] = gene_subgroup

    return gene_groups_dt


def assign_start_codon_position_to_chimeric(chimeric_gene, gene_groups_dt, trans_cds_dt, gtf_obj):

    # 1) Get the genes that are overlapped by the chimeric gene
    try:
        overlaped_genes = gene_groups_dt[chimeric_gene]
    except KeyError as err:
        sys.exit(f'ERROR: No overlapping group found for Gene {chimeric_gene}')

    if not overlaped_genes:
        # if ref_pos is None, the transcripts are assigned into the "No CDS" group
        return None

    # Sort again just in case
    overlaped_genes = sorted(overlaped_genes, key=lambda g_id: gtf_obj.gene_coords_dt[g_id])

    try:
        _ = gtf_obj.gene_trans_dt.keys()
    except KeyError as err:
        sys.exit(f'ERROR: Gene {chimeric_gene} not present in the specified GTF file')

    # 2) Get chimeric gene information
    rep_trans = sorted(gtf_obj.gene_trans_dt[chimeric_gene])[0]
    chromosome = gtf_obj.trans_chrom_dt[rep_trans]
    chrom, strand = chromosome[:-1], chromosome[-1]
    gene_coords = gtf_obj.gene_coords_dt[chimeric_gene]

    # 3) Select the "first" gene that makes the chimeric gene
    if strand == "+":
        ref_gene = overlaped_genes[0]
    elif strand == "-":
        ref_gene = overlaped_genes[-1]
    else:
        sys.exit(f'ERROR: Strand {strand} not recognized. It must be + or - strand.')

    # 4) Select the reference start-codon position and assign it to the chimeric gene
    try:
        # Select the first (by exon position) transcript of the gene to select the representative start-codon
        rep_trans = sorted(gtf_obj.gene_trans_dt[ref_gene], key=lambda t_id: gtf_obj.trans_exons_dt[t_id][0][0])[0]
        ref_pos = trans_cds_dt[rep_trans][0][0]
    except KeyError as err:
        print(f'ERROR: CDS start for gene {ref_gene} not found.')
        return None

    return ref_pos


def fix_chimeric_start_codon(gtf_obj, chimeric_table, trans_cds_dt, trans_sequences_dt):

    # 0) Identify chimeric-genes (separate script)
    # 1) Identify the non-chimeric gene ("ref_gene") at the 5' UTR region of the chimeric model
    # 2) Use the start-codon of the "ref-ref_gene" to re-fix the start-codon coordinates of the chimeric model

    print(time.asctime(), f"Re-fixing start-codon of chimeric genes: {chimeric_table}")

    chimeric_genes = set()
    with open(chimeric_table) as fh:
        next(fh)
        for row in fh:
            gene = row.strip("\n")
            chimeric_genes.add(gene)

    gene_groups_dt = get_gene_groups(gtf_obj, chimeric_genes)

    n_genes = len(chimeric_genes)

    chimeric_output_dt = defaultdict(dict)
    for z, chm_gene in enumerate(sorted(chimeric_genes)):

        print(f'Processing Gene {chm_gene}, {(z / n_genes) * 100:.1f}% complete ({z + 1}/{n_genes})')

        ref_pos = assign_start_codon_position_to_chimeric(chm_gene, gene_groups_dt, trans_cds_dt, gtf_obj)

        chm_trans_data_dt = get_transcript_data_from_gtf_obj(chm_gene, gtf_obj, trans_sequences_dt)

        gene_output_dt, gene_cat_dt = fix_atg_position(chm_trans_data_dt, ref_pos)

        chimeric_output_dt["trans_cds_dt"].update(gene_output_dt["trans_cds_dt"])
        chimeric_output_dt["trans_cds_seq_dt"].update(gene_output_dt["trans_cds_seq_dt"])
        chimeric_output_dt["trans_header_dt"].update(gene_output_dt["trans_header_dt"])

    return chimeric_output_dt


def write_transfix_tables(gtf_obj, cat_dt, cycle_trans_dt, trans_cds_dt, outfolder, outname):

    # Write table tracking the untranslatable transcripts and the transcripts translation cycle
    print(time.asctime(), "Writing transcript fixing-cycle table")
    row_list, previous_set, non_redundant_cycle_dt = [], set(), {}
    for cycle, trans_set in sorted(cycle_trans_dt.items()):
        non_redundant_set = trans_set - previous_set
        previous_set.update(trans_set)

        # Create a non-redundant cycle dictionary to visualize the number of transcripts per translation cycle
        non_redundant_cycle_dt[cycle] = non_redundant_set

        for trans in sorted(non_redundant_set):

            trans_strand = gtf_obj.trans_sense_dt[trans]

            try:
                trans_cds = trans_cds_dt[trans]
            except KeyError as err:
                trans_cds = []

            if trans_cds:
                if trans_strand == "+":
                    trans_start = trans_cds[0][0]
                elif trans_strand == "-":
                    trans_start = trans_cds[-1][-1]
                else:
                    trans_start = "None"
            else:
                trans_start = "None"

            row = f"{cycle},{trans},{trans_start}\n"
            row_list.append(row)

    cycle_table = os.path.join(outfolder, outname + "_cycle_table.csv")
    with open(cycle_table, "w+") as fh:
        fh.write(f"Translation_cycle,Transcript_ID,Start_position\n")
        for row in row_list:
            fh.write(row)

    verbose = True
    if verbose:

        tot_trans = sum([len(cat_dt[k]) for k in ["processed_transcripts", "cds_not_found", "start_codon_not_atg", "seq_not_present", "unprocessed_transcripts"]])

        # Function to calculate percentage
        get_perc = lambda n, tot: round((n/tot)*100, 1)

        n = len(cat_dt["processed_transcripts"])
        perc_n = get_perc(n, tot_trans)
        print(f"\nTranslated transcripts: {n} ({perc_n}%)")

        # Print out information of each start-codon fixing cycle
        for cycle, trans_set in non_redundant_cycle_dt.items():
            n = len(trans_set)
            perc_n = get_perc(n, tot_trans)
            print(f"{cycle}: {n} ({perc_n}%)")
        print("\n")

        print(time.asctime(), "Writing TransFix group categories tables")

        # Print out information for each removed category
        cat_lst = [("CDS_not_found", "cds_not_found"), ("ATG_not_within_CDS", "atg_not_in_cds"), ("Start-codon_is_not_ATG", "atg_not_in_cds"),
                   ("Sequence_not_found", "seq_not_present"), ("Unprocessed_transcripts", "unprocessed_transcripts")]

        for cat_name, cat_key in cat_lst:
            cat_lines_key = f"{cat_key}_lines"

            if cat_dt[cat_key]:
                n = len(cat_dt[cat_key])
                perc_n = get_perc(n, tot_trans)
                print(f"{cat_name}: {n} ({perc_n}%)")
                table_name = os.path.join(outfolder, outname + f"_{cat_name}.csv")
                write_table(cat_dt[cat_lines_key], cat_dt[cat_key], table_name)
        print("\n")
