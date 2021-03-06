import os
import sys
import time
from collections import defaultdict

from Bio.Seq import Seq

from lib.parsing.fasta_parsing_tools import write_fasta_file
from lib.report.transfeat_report import get_transfeat_data, group_models_into_categories
from lib.findlorf.findlorf_tools import convert_from_genomic_to_relative, create_coordinate_lookup_table


def get_location(f_id, gtf_obj, feature="TRANS"):

    flat = lambda l: [e for s in l for e in s]

    if feature.upper() == "TRANS":
        try:
            chrom = gtf_obj.trans_chrom_dt[f_id][:-1]
            t_exons = flat(gtf_obj.trans_exons_dt[f_id])
            coord = (min(t_exons), max(t_exons))
        except KeyError:
            return "-"

    elif feature.upper() == "GENE":
        try:
            t_id = sorted(gtf_obj.gene_trans_dt[f_id])[0]
            chrom = gtf_obj.trans_chrom_dt[t_id][:-1]
            coord = gtf_obj.gene_coords_dt[f_id]
        except KeyError:
            return "-"

    elif feature.upper() == "CDS":
        try:
            chrom = gtf_obj.trans_chrom_dt[f_id][:-1]
            t_cds = flat(gtf_obj.trans_cds_dt[f_id])
            coord = (min(t_cds), max(t_cds))
        except KeyError:
            return "-"

    else:
        sys.exit(f'Feature type must be either "trans" or "gene", not "{feature}"')

    coord_txt = f'{chrom}:{coord[0]}-{coord[1]}'

    return coord_txt


def get_orf_seq(orf, trans_seq):

    try:
        orf_nucl_seq = trans_seq[orf[0]:orf[1]+1]
    except:
        orf_nucl_seq = ""

    orf_pep_seq = str(Seq(orf_nucl_seq).translate(to_stop=True))

    return orf_pep_seq, orf_nucl_seq


def group_transcripts_by_overlap(t_group, trans_exons_dt):

    flat = lambda l: [e for sub in l for e in sub]
    get_model = lambda t_id: trans_exons_dt[t_id]
    get_transcript_start_end = lambda exon_model: (exon_model[0][0], exon_model[-1][-1])
    get_start = lambda t_id: get_transcript_start_end(get_model(t_id))[0]

    t_group = sorted(t_group, key=lambda t_id: (get_start(t_id), t_id))

    result = []
    overlap_group, prev_max_end = [], 0
    for i, trans in enumerate(t_group):
        trans_start, trans_end = get_transcript_start_end(get_model(trans))

        if trans_start <= prev_max_end:
            overlap_group.append(trans)

            # Choose largest end
            if trans_end > prev_max_end:
                prev_max_end = trans_end

        else:
            # First group is empty, this check avoid introducing empty groups
            if overlap_group:
                result.append(overlap_group)
            prev_max_end = trans_end

            overlap_group = []
            overlap_group.append(trans)

        # Necessary to append the last group
        if i+1 == len(t_group):
            result.append(overlap_group)

    assert len(t_group) == len(flat(result))

    return result


def generate_feature_tag(gtf_obj, feature_dicts):

    # Note: This method assumes that all dicts have the same keys.

    trans_features_dt = {}
    for trans in feature_dicts["Auto"].keys():
        trans_features = []
        coding_flag = False

        # The order of these IF checks is important!
        # Annotate the NON_CODING features first!
        if feature_dicts["PTC"][trans] and "PTC" not in trans_features:
            trans_features.append("PTC")

        if feature_dicts["NMD"][trans] and "NMD" not in trans_features:
            trans_features.append("NMD")

        # Avoid the redundant annotation of PTC as also short_peptide (most PTCs are also short_peptides)
        # trans_features list will still be empty if it's not a PTC
        if not trans_features:
            if feature_dicts["Short_ORF"][trans] and "Short_ORF" not in trans_features:
                trans_features.append("Short_ORF")

            if feature_dicts["No_ORF"][trans] and "No_ORF" not in trans_features:
                trans_features.append("No_ORF")

        # If trans_features list is still empty (NO not-coding features), then we tag the transcript as "Coding"
        if not trans_features:
            trans_features.append("Coding")
            coding_flag = True

        # Avoid redundant identification of PTC and AS_in_UTR (Most PTCs will have the same "AS_in_UTR" region)
        if coding_flag:
            if feature_dicts["AS_in_UTR"][trans]:
                AS_location = feature_dicts["AS_Location"][trans]
                trans_features.append(AS_location)

            if feature_dicts["NAGNAG"][trans] and "NAGNAG" not in trans_features:
                trans_features.append("NAGNAG")

        # Report ldORF only for Non-coding transcripts
        if feature_dicts["Long_3UTR"][trans] and "ldORF" not in trans_features and "Coding" not in trans_features:
            trans_features.append("ldORF")

        # If trans_features is still empty, tag as -
        if not trans_features and "-" not in trans_features:
            trans_features.append("-")

        if feature_dicts["NMD"][trans] and "NMD" not in trans_features:
            trans_features.append("NMD")

        trans_features_dt[trans] = ";".join(trans_features)

    # 2) Top-level information (aka: is transcript Coding, Non-coding, or Unproductive)
    coding_potentiality_dt, coding_features_dt, alternative_ORF_dt = [defaultdict(str) for _ in range(3)]

    # 2.1) Initial classification into Coding vs Non_coding
    for trans, trans_features in trans_features_dt.items():

        if "CODING" in trans_features.upper():
            coding_potentiality_dt[trans] += "Coding"

            # Remove redundant tags
            trans_features = trans_features.replace("Coding", "")

            if trans_features:
                coding_features_dt[trans] += trans_features
            # In case the string is empty (only Coding feature tag reported)
            else:
                coding_features_dt[trans] += "-"

        else:
            coding_potentiality_dt[trans] += "Non_Coding"
            coding_features_dt[trans] += trans_features

        try:
            trans_nmd_features = feature_dicts["NMD_features"][trans]
        except KeyError:
            trans_nmd_features = ""

        if "UORF" in trans_nmd_features.upper():
            if not alternative_ORF_dt[trans]:
                alternative_ORF_dt[trans] += "uORF"
            else:
                alternative_ORF_dt[trans] += ";uORF"

        # Move the ldORF information to the alternative ORF column
        if "LDORF" in trans_features.upper():
            if not alternative_ORF_dt[trans]:
                alternative_ORF_dt[trans] += "ldORF"
            else:
                alternative_ORF_dt[trans] += ";ldORF"

            coding_features_dt[trans] = coding_features_dt[trans].replace(";ldORF", "")

        for tag_dt in [coding_features_dt, alternative_ORF_dt]:
            tag_line = tag_dt[trans]

            # Clean the coding_features tag
            while ";;" in tag_line:
                tag_line = tag_line.replace(";;", ";")

            if tag_line:
                if tag_line[0] == ";":
                        tag_line = tag_line[1:]
            if tag_line:
                if tag_line[-1] == ";":
                        tag_line = tag_line[:-1]

            if not tag_line:
                tag_line = "-"

            tag_dt[trans] = tag_line

        # Add character for empty features
        if not coding_potentiality_dt[trans]:
            coding_potentiality_dt[trans] += "-"

        if not coding_features_dt[trans]:
            coding_features_dt[trans] += "-"

        if not alternative_ORF_dt[trans]:
            alternative_ORF_dt[trans] += "-"

    # 3) If there is at least one coding isoform in the gene, then use the terminology "Unpropductive" instead of "Non_Coding"
    noncoding_genes_set = identify_non_coding_genes(gtf_obj, coding_potentiality_dt)

    # Maintain the tag Non_Coding only for exclusivley Non-Coding genes, otherwise
    # If the gene contain at least a coding isoform, use the tag "Unproductive" for "Non-coding" isoforms
    for gene, gene_transcripts in sorted(gtf_obj.gene_trans_dt.items()):
        # Sort transcripts by their start-codon genomic position
        for trans in sorted(gene_transcripts):
            trans_coding_potent = coding_potentiality_dt[trans]
            if trans_coding_potent.upper() == "NON_CODING" and gene not in noncoding_genes_set:
                coding_potentiality_dt[trans] = trans_coding_potent.replace("Non_Coding", "Unproductive")

    return coding_potentiality_dt, coding_features_dt, alternative_ORF_dt


def get_transcript_start_codon_relative_position(gtf_obj):

    print(time.asctime(), "Obtaining start-codons relative position within transcript")

    flat = lambda l: [e for sub in l for e in sub]

    # Get transcript start-codon genomic position
    trans_genomic_start_dt = {}
    for trans, trans_cds in gtf_obj.trans_cds_dt.items():

        trans_cds = flat(trans_cds)

        # Important!, Ignore transcripts without an annotated CDS
        if not trans_cds:
            print(f'Transcript "{trans}" does not present an annotated CDS')
            continue

        trans_sense = gtf_obj.trans_sense_dt[trans]
        if trans_sense == "+":
            genomic_start = min(trans_cds)
        elif trans_sense == "-":
            genomic_start = max(trans_cds)
        else:
            print(f"Transcript {trans} strand is not valid (\"{trans_sense}\").")
            continue

        trans_genomic_start_dt[trans] = genomic_start

    # Format: trans_ID = (relative_exons, genomic_exons)
    # Beware that the relative_exons coordinates start from 1, not 0
    exon_lookup_table = create_coordinate_lookup_table(gtf_obj)

    trans_relative_start_dt = {}
    for trans, (relative_positions, genomic_positions) in exon_lookup_table.items():

        trans_sense = gtf_obj.trans_sense_dt[trans]

        try:
            trans_genomic_st = trans_genomic_start_dt[trans]
        except KeyError:
            # print(f"Transcript {trans} does not contain a start-codon. Skipping.")
            continue

        for relative_coord_pair, genomic_coord_pair in zip(relative_positions, genomic_positions):

            # Necessary because the order of the genomic coordinates for "-" strand is inverted
            genomic_coord_pair = sorted(genomic_coord_pair)
            if genomic_coord_pair[0] <= trans_genomic_st <= genomic_coord_pair[1]:

                if trans_sense == "+":
                    for i, genomic_pos in enumerate(range(genomic_coord_pair[0], genomic_coord_pair[1]+1)):
                        relative_pos = relative_coord_pair[0] + i

                        if genomic_pos == trans_genomic_st:
                            trans_relative_st = relative_pos

                            # Important! The CDS coordinates from the ORF index file (JSON) are 0 based,
                            # but this position is +1 based. Thus, we must -1 this value

                            trans_relative_st -= 1
                            assert trans_relative_st >= 0

                            trans_relative_start_dt[trans] = trans_relative_st

                elif trans_sense == "-":
                    # Since it is the negative strand, given the order of the ORF in the ORF index file,
                    # we have to use a decremental loop with an step of -1
                    for i, genomic_pos in enumerate(range(genomic_coord_pair[1], genomic_coord_pair[0]-1, -1)):
                        relative_pos = relative_coord_pair[0] + i

                        if genomic_pos == trans_genomic_st:
                            trans_relative_st = relative_pos

                            trans_relative_st -= 1
                            assert trans_relative_st >= 0

                            trans_relative_start_dt[trans] = trans_relative_st

                else:
                    sys.exit("Strand error")

    return trans_relative_start_dt


def get_genes_authentic_stop_codon_position(gtf_obj):

    print(time.asctime(), "Selecting Genes stop-codon position")

    flat = lambda l: [e for sub in l for e in sub]

    # Get transcript start-codon genomic position
    auth_stop_dt = {}
    for gene, gene_transcripts in gtf_obj.gene_trans_dt.items():

        # Initialize variables to avoid warning by EDI
        gene_strand, trans = None, "PLACEHOLDER"

        gene_cds = []
        for trans in sorted(gene_transcripts):
            try:
                trans_cds = gtf_obj.trans_cds_dt[trans]
            except KeyError:
                trans_cds = []

            trans_cds = flat(trans_cds)
            gene_cds.extend(trans_cds)

            # Since all transcript in a gene have the same strand, we only need to pick any trans to select gene strand
            gene_strand = gtf_obj.trans_sense_dt[trans]

        # Some genes contain only transcripts without CDS
        if not gene_cds:
            continue

        if gene_strand == "+":
            auth_stop = max(gene_cds)
        elif gene_strand == "-":
            auth_stop = min(gene_cds)
        else:
            sys.exit(f"Strand error for transcript ({trans})")

        auth_stop_dt[gene] = auth_stop

    return auth_stop_dt


def translate_transcript_cds(trans_seq_dt, gtf_obj):

    print(time.asctime(), "Translating transcripts sequences")

    cds_seq_dt, pep_seq_dt, fasta_header_dt = ({} for _ in range(3))
    for trans, trans_cds in gtf_obj.trans_cds_dt.items():

        trans_sense = gtf_obj.trans_sense_dt[trans]

        trans_exons = gtf_obj.trans_exons_dt[trans]
        genomic_to_relative_coords_dt = convert_from_genomic_to_relative(trans_exons)

        try:
            trans_seq = trans_seq_dt[trans]
        except KeyError:
            continue

        if not trans_cds:
            continue

        seq = Seq(trans_seq)

        if trans_sense == "+":
            genomic_cds_st = trans_cds[0][0]
            relative_cds_st = genomic_to_relative_coords_dt[genomic_cds_st]
            cds_seq = seq[relative_cds_st:]

        elif trans_sense == "-":
            genomic_cds_st = trans_cds[-1][-1]
            relative_cds_st = genomic_to_relative_coords_dt[genomic_cds_st]

            seq = seq.reverse_complement()
            cds_seq = seq[:relative_cds_st + 1]
            cds_seq = cds_seq.reverse_complement()

        else:
            sys.exit(f"Strand error found on transcript {trans}")

        if cds_seq.startswith('ATG'):
            peptide = cds_seq.translate(to_stop=True)
        else:
            peptide = "-"

        pep_seq_dt[trans] = str(peptide)
        cds_seq_dt[trans] = str(cds_seq)

        # Get information to create header
        chrom = gtf_obj.trans_chrom_dt[trans][:-1]
        trans_coord = f'{chrom}:{trans_cds[0][0]}-{trans_cds[-1][-1]}'

        gene = gtf_obj.trans_gene_dt[trans]
        gene_name = gtf_obj.gene_name_dt[gene]

        # Clean gene_name from redundant information/spaces
        gene_name = gene_name.replace(gene, '')
        while '  ' in gene_name:
            gene_name = gene_name.replace('  ', ' ')

        if gene_name:
            trans_header = f'>{trans} | {gene} | {gene_name} | {trans_coord}'
        else:
            trans_header = f'>{trans} | {gene} | {trans_coord}'

        fasta_header_dt[trans] = trans_header

    return fasta_header_dt, cds_seq_dt, pep_seq_dt


def identify_non_coding_genes(gtf_obj, coding_potentiality_dt):

    noncoding_genes_set = set()

    for gene_id, gene_transcripts in gtf_obj.gene_trans_dt.items():

        coding_gene_flag = False
        for t_id in sorted(gene_transcripts):
            try:
                t_cod_pont = coding_potentiality_dt[t_id]
            except KeyError:
                t_cod_pont = "UNKNOWN"

            if t_cod_pont.upper() == "CODING":
                coding_gene_flag = True

        if not coding_gene_flag:
            noncoding_genes_set.add(gene_id)

    return noncoding_genes_set


def write_transfeat_table(gtf_obj, features_info_dicts, pep_seq_dt, outfolder, outname,
                          ldorf_ids=None, uorf_ids=None, pep_len=50):

    outfile = os.path.join(outfolder, f"{outname}.csv")

    print(time.asctime(), f"Generating TransFeat table: {outfile}")

    # Proper way to initialize 'default empty containers' in python functions to avoid bugs; however:
    # Initializing these sets is most likely not necessary as an "if check" is done before accessing them TODO check
    if not ldorf_ids:
        ldorf_ids = set()

    if not uorf_ids:
        uorf_ids = set()

    # Generate rows for output table
    header = f"Gene_ID,Transcript_ID,Coding_potentiality,Features,Alternative_ORF,NMD_features,Transcript_coordinates,CDS_coordinates,Translation\n"
    table_rows = [header]

    # If True, it adds an empty row between genes in the output table
    add_spacer = False

    # Sort transcripts by their Gene ID
    for gene, gene_transcripts in sorted(gtf_obj.gene_trans_dt.items()):
        # Sort transcripts by their start-codon genomic position
        for trans in sorted(gene_transcripts, key=lambda t: gtf_obj.trans_exons_dt[t][0]):

            trans_chrom = gtf_obj.trans_chrom_dt[trans][:-1]
            trans_exons_flat = [exon for exon_pair in gtf_obj.trans_exons_dt[trans] for exon in exon_pair]
            trans_coords = f"{trans_chrom}:{min(trans_exons_flat)}-{max(trans_exons_flat)}"

            if gtf_obj.trans_cds_dt[trans]:
                trans_sense = gtf_obj.trans_sense_dt[trans]
                trans_cds_flat = [cds for cds_pair in gtf_obj.trans_cds_dt[trans] for cds in cds_pair]
                if trans_sense == '+':
                    cds_coords = f"{trans_chrom}:{min(trans_cds_flat)}-{max(trans_cds_flat)}"
                elif trans_sense == '-':
                    cds_coords = f"{trans_chrom}:{max(trans_cds_flat)}-{min(trans_cds_flat)}"
                else:
                    sys.exit(f"Transcript {trans} sense must be + or -, not \"{trans_sense}\". Aborting.")
            else:
                cds_coords = "-"

            try:
                aa_seq = pep_seq_dt[trans]
            except KeyError as err:
                aa_seq = "-"

            if not aa_seq:
                aa_seq = "-"

            coding_potential = features_info_dicts["Coding_potentiality"][trans]
            features = features_info_dicts["Coding_features"][trans]

            # Report alt_ORF only for Unproductive transcripts
            if coding_potential.upper() == "UNPRODUCTIVE":
                alt_ORF = features_info_dicts["Alternative_ORF"][trans]

                try:
                    NMD_features = features_info_dicts["NMD_features"][trans]
                    NMD_features = NMD_features.replace(";ldORF", "")
                    if not NMD_features:
                        NMD_features = "-"
                except KeyError:
                    NMD_features = "-"

            else:
                alt_ORF = "-"
                NMD_features = "-"

            # In some rare cases a "PTC" flag may be reported for "Non_Coding" genes if a genomic region contains
            # two different annotated gene IDs. Fix this mis-tag by replacing "PTC" to "Short_ORF" instead
            if coding_potential.upper() == "NON_CODING" and "PTC" in features.upper():
                features = "Short_ORF"

            # Fix: some sequences of sufficient length are being incorrectly assigned as Non_Coding due to
            # the 'spurious PTC flag check' stated above. Thus, perform length check again before writing output
            if coding_potential.upper() == "NON_CODING" and features.upper() == "SHORT_ORF":
                if len(aa_seq) >= pep_len:
                    coding_potential = "Coding"
                    features = "-"

            if ldorf_ids:
                if "LDORF" in alt_ORF.upper() and trans not in ldorf_ids:
                    alt_ORF = alt_ORF.replace("ldORF", "").replace(";", "")

            if uorf_ids:
                if "UORF" in alt_ORF.upper() and trans not in ldorf_ids:
                    alt_ORF = alt_ORF.replace("UORF", "").replace(";", "")

            if not alt_ORF:
                alt_ORF = "-"

            row = f"{gene},{trans},{coding_potential},{features},{alt_ORF},{NMD_features},{trans_coords},{cds_coords},{aa_seq}\n"
            table_rows.append(row)

        if add_spacer:
            # Empty row to separate genes in the table
            n_fields = 8
            spacer_row = "," * n_fields + "\n"
            table_rows.append(spacer_row)

    # TransFeat output table
    with open(outfile, "w+") as fh:
        for row in table_rows:
            fh.write(row)

    return outfile


def write_subcategories_fasta(gtf_obj, transfeat_table, sequences_dicts, outfolder, outname):

    print(time.asctime(), "Generating additional fasta files")

    # Get information from TransFeat table
    trans_by_feature_dt, coding_categories_dt = get_transfeat_data(transfeat_table)
    models_categories_dt = group_models_into_categories(gtf_obj, coding_categories_dt)

    # Main categories of interest
    noncoding_genes, noncoding_trans = models_categories_dt["NON_CODING"]
    coding_genes, coding_trans = models_categories_dt["PROTEIN_CODING"]

    # 0) Non-coding transcripts
    noncoding_seq_dt, noncoding_header_dt = {}, {}
    for t_id in sorted(noncoding_trans):
        t_loc = get_location(t_id, gtf_obj)
        try:
            t_seq = sequences_dicts["Exonic_seq"][t_id]
        except KeyError:
            t_seq = "-"

        noncoding_seq_dt[t_id] = t_seq
        noncoding_header_dt[t_id] = f'>{t_id} | {t_loc}'

    if noncoding_seq_dt:
        cat_outfile = os.path.join(outfolder, f'{outname}_non_coding_nuc.fa')
        write_fasta_file(noncoding_seq_dt, cat_outfile, noncoding_header_dt)

    # 1) Protein-coding transcripts
    coding_nuc_seq_dt, coding_pep_seq_dt, coding_header_dt = [{} for _ in range(3)]
    for t_id in sorted(coding_trans):
        t_loc = get_location(t_id, gtf_obj)
        try:
            t_nuc_seq = sequences_dicts["CDS_seq"][t_id]
            t_pep_seq = sequences_dicts["Peptide_seq"][t_id]
        except KeyError:
            t_nuc_seq, t_pep_seq = "-", "-"

        coding_nuc_seq_dt[t_id] = t_nuc_seq
        coding_pep_seq_dt[t_id] = t_pep_seq
        coding_header_dt[t_id] = f'>{t_id} | {t_loc}'

    if coding_nuc_seq_dt and coding_pep_seq_dt:
        nuc_outfile = os.path.join(outfolder, f'{outname}_coding_nuc.fa')
        write_fasta_file(coding_nuc_seq_dt, nuc_outfile, coding_header_dt)

        pep_outfile = os.path.join(outfolder, f'{outname}_coding_pep.fa')
        write_fasta_file(coding_pep_seq_dt, pep_outfile, coding_header_dt)

    # 2) Alt ORFs sequences
    # Both of this dict has the same structure:  t_id: list([orf_coords, orf_pep_seq, orf_nuc_seq])
    ldorf_data_dt = sequences_dicts["ldORF_data"]
    uorf_data_dt = sequences_dicts["uORF_data"]

    for altORF_dt, orf_tag in [(ldorf_data_dt, "ldORF"), (uorf_data_dt, "uORF")]:

        orf_nuc_seq_dt, orf_pep_seq_dt, orf_header_dt = [{} for _ in range(3)]

        for t_id in sorted(altORF_dt):
            t_chrom = gtf_obj.trans_chrom_dt[t_id][:-1]

            for i, t_orf_data in enumerate(altORF_dt[t_id], 1):
                ORF_coord, ORF_pep_seq, ORF_nuc_seq = t_orf_data

                tagged_id = f"{t_id}_{orf_tag}_{i}"
                t_loc = f"{t_chrom}:{ORF_coord[0]}-{ORF_coord[1]}"

                orf_nuc_seq_dt[tagged_id] = ORF_nuc_seq
                orf_pep_seq_dt[tagged_id] = ORF_pep_seq
                orf_header_dt[tagged_id] = f'>{tagged_id} | {t_loc}'

        if orf_nuc_seq_dt and orf_pep_seq_dt:
            orf_nuc_outfile = os.path.join(outfolder, f'{outname}_{orf_tag}_nuc.fa')
            write_fasta_file(orf_nuc_seq_dt, orf_nuc_outfile, orf_header_dt)

            orf_pep_outfile = os.path.join(outfolder, f'{outname}_{orf_tag}_pep.fa')
            write_fasta_file(orf_pep_seq_dt, orf_pep_outfile, orf_header_dt)


def write_NMD_table(gtf_obj, feature_dicts, sequences_dicts, outfolder, outname):

    print(time.asctime(), "Generating table with additional NMD feature information")

    nmd_outfile = os.path.join(outfolder, f"{outname}_NMD_features.csv")
    with open(nmd_outfile, "w+") as fh:
        # Write header
        fh.write("Transcript_ID,NMD_features,uORF_coordinates,uORF_sequence\n")
        for trans_id, trans_nmd_features in sorted(feature_dicts["NMD_features"].items()):

            # Ignore transcripts that doesn't report any NMD feature
            if trans_nmd_features == "-":
                continue

            trans_chrom = gtf_obj.trans_chrom_dt[trans_id][:-1]

            try:
                uORF_data_list = sequences_dicts["uORF_data"][trans_id]
            except KeyError:
                uORF_data_list = []

            if uORF_data_list:
                for i, (uORF, uORF_pep_seq, uORF_nuc_seq) in enumerate(uORF_data_list, 1):
                    line = f"{trans_id},{trans_nmd_features},{trans_chrom}:{uORF[0]}-{uORF[1]},{uORF_pep_seq}\n"
                    fh.write(line)
            else:
                line = f"{trans_id},{trans_nmd_features},-,-\n"
                fh.write(line)
