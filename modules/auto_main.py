import os
import time
from modules.findlorf_main import findlorf_main
from modules.transfix_main import transfix_main
from modules.transfeat_main import transfeat_main
from lib.parsing.gtf_parsing_tools import filter_gtf_file, add_features_to_gtf


def file_exist(outfile, skip_message=True):

    if os.path.exists(outfile):
        print(time.asctime(), f'File already exist: {outfile}')
        if skip_message:
            print(time.asctime(), f'Keeping current file')
        else:
            print(time.asctime(), f'Overwriting file')
        return True
    else:
        return False


def run_transuite(gtf, fasta, outpath, outname, iter_th=5, cds_th=30, pep_th=50, ptc_th=70, chimeric=None):

    # 1) Remove CDS information from input annotation file
    filtered_gtf = os.path.splitext(gtf)[0] + '_exons.gtf'
    if not file_exist(filtered_gtf):
        filtered_gtf = filter_gtf_file(gtf)

    # 2) Run FindLORF
    tfind_name = outname.replace(".gtf", "") + "_longorf"
    tfind_folder = os.path.join(outpath, tfind_name)
    transfind_gtf = os.path.join(tfind_folder, f"{tfind_name}.gtf")

    orf_index_filename = os.path.splitext(os.path.basename(filtered_gtf))[0] + "_ORF_index.json"
    orf_index_file = os.path.join(tfind_folder, orf_index_filename)

    if not file_exist(transfind_gtf) or not file_exist(orf_index_file):
        transfind_gtf, orf_index_file = findlorf_main(filtered_gtf, fasta, outpath, outname, cds_th=cds_th,
                                                      filter_gtf=False)

    # 3) Run TransFix
    tfix_name = outname.replace(".gtf", "") + "_transfix"
    tfix_folder = os.path.join(outpath, tfix_name)
    transfix_gtf = os.path.join(tfix_folder, f"{tfix_name}.gtf")

    if not file_exist(transfix_gtf):
        transfix_gtf = transfix_main(transfind_gtf, fasta, outpath, outname, iter_th=iter_th, chimeric=chimeric)

        # 3.5) Add extra features to the annotation
        transfix_gtf = add_features_to_gtf(transfix_gtf)

    # 4) Run TransFeat
    tfeat_name = outname.replace(".gtf", "") + "_transfeat"
    tfeat_folder = os.path.join(outpath, tfeat_name)
    transfeat_table = os.path.join(tfeat_folder, f"{tfeat_name}.csv")

    if not file_exist(transfeat_table):
        transfeat_table = transfeat_main(transfix_gtf, fasta, outpath, outname, pep_len=pep_th, ptc_len=ptc_th,
                                         uorf_len=10, sj_dist=50, utr3_len=350, orf_index=orf_index_file)

    print(time.asctime(), f'Annotation file with fixed start-codon coordinates: {transfix_gtf}')
    print(time.asctime(), f'Table with coding characterization of transcripts: {transfeat_table}')
