import sys
import time
import Bio.Seq
from collections import defaultdict


def get_fasta_sequences(fasta_file):

    print(time.asctime(), f"Uploading transcripts sequences from fasta file: {fasta_file}")

    # TODO make function return error if file is not in FASTA format

    # Upload Transcripts sequences from FASTA file
    trans_seq_dt = defaultdict(str)
    with open(fasta_file) as fh:
        row_id, seq = None, ''
        for row in fh:
            if row.startswith(">"):
                seq = ''
                row_id = row.replace(" ", "#").replace("\t", "#").split("#")[0].strip(">").strip("\n")
            else:
                # Important! Sanitize sequence so it is all in capital letters!
                seq = row.strip("\n")
                trans_seq_dt[row_id] += seq.upper()

    return trans_seq_dt


def write_fasta_file(seq_dt, outname, header_dt=None, line_len=70):

    print(time.asctime(), f"Writing output fasta file: {outname}")

    if not header_dt:
        header_dt = {t_id: t_id if t_id.startswith(">") else f">{t_id}" for t_id in seq_dt.keys()}

    # Sanitize dictionary with transcripts fasta headers to ensure that all headers starts with ">"
    # I cant use isinstance() because isinstance() of dict is the same as that of defaultdict
    if type(header_dt) == type(dict()):
        temp_dt = {k: v if v.startswith(">") else f">{v}" for k, v in header_dt.items()}
    elif type(header_dt) == type(defaultdict()):
        temp_dt = {}
        for t_id, t_headers in header_dt.items():
            t_headers_temp = []
            for t_head in t_headers:
                if not t_head.startswith(">"):
                    t_headers_temp.append(f">{t_head}")
                else:
                    t_headers_temp.append(t_head)
            temp_dt[t_id] = t_headers_temp
    else:
        sys.exit(f"ERROR. The object containing the FASTA headers must be either dict of defaultdict, not '{type(header_dt)}'")

    header_dt = temp_dt

    slice_line = lambda line, len_th: [line[i: i + len_th] for i in range(0, len(line), len_th)]

    with open(outname, 'w+') as fh:
        for i, (trans, trans_seq) in enumerate(sorted(seq_dt.items())):

            if isinstance(trans_seq, (str, Bio.Seq.Seq)):
                try:
                    trans_header = header_dt[trans]
                except KeyError:
                    trans_header = f'>{trans}'
                fh.write(trans_header + '\n')

                if not trans_seq:
                    fh.write("-" + '\n')
                    continue

                for line_segment in slice_line(trans_seq, line_len):
                    fh.write(str(line_segment) + '\n')

            elif isinstance(trans_seq, list):
                try:
                    trans_header = header_dt[trans]
                except KeyError:
                    trans_header = f'>{trans}'

                # These lists were populated together, thus their length must be the same
                assert len(trans_seq) == len(trans_header)

                for t_seq, t_header in zip(trans_seq, trans_header):
                    fh.write(t_header + '\n')

                    if not trans_seq:
                        fh.write("-" + '\n')
                        continue

                    for line_segment in slice_line(t_seq, line_len):
                        fh.write(str(line_segment) + '\n')

            else:
                sys.exit(f'Error found on transcript {trans} ({trans_seq}) while writing output fasta file: {outname} '
                         f'(type: {type(trans_seq)})')
