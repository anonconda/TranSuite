"""
Created on Wed May 15 14:25:00 2019
@author: Juan C Entizne
@email: e.entizne[at]dundee.ac.uk
"""

import os
import sys
import time
import traceback
from lib.logger.logger import logger, clean_log
from argparse import ArgumentParser, RawTextHelpFormatter
from modules.transfix_main import transfix_main
from lib.parsing.gtf_parsing_tools import add_features_to_gtf


description = \
    "Description:\n" + \
    "TransFix fix the same start codon for all transcripts in a gene, translate them, and annotates the resulting CDS.\n"

parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=False)

parser.add_argument('--gtf',
                    dest="gtf",
                    help="Transcriptome annotation file in GTF format.")

parser.add_argument('--fasta',
                    dest="fasta",
                    help='Transcripts fasta file (nucleotide sequence of exonic regions).')

parser.add_argument("--iter",
                    dest="iter_th", type=int, default=5,
                    help="Maximum number of 'start-fixing & translation' cycles to identify alternative start sites. "
                         "Default: 5")

parser.add_argument('--outpath',
                    dest="outpath",
                    help="(Absolute path of the output folder.")

parser.add_argument('--outname',
                    dest="outname",
                    help="Name of the output file (without file extension).")


def main():
    args = parser.parse_args()

    # Create logfile to track the analysis, overwrite it if it exist ("w+" mode)
    time_stamp = time.strftime("%Y%m%d-%H%M%S")
    logfile = os.path.join(args.outpath, f"{time_stamp}_{args.outname}_logfile_temp.out")
    logger(logfile, w_mode="w+")

    # Check input arguments
    for arg_val, arg_name in zip([args.gtf, args.fasta, args.outpath, args.outname],
                                 ["--gtf", "--fasta", "--outpath", "--outname"]):
        if not arg_val:
            sys.exit(f'Error: No value specified for argument "{arg_name}"')

    for fl, arg_name in zip([args.gtf, args.fasta], ["--gtf", "--fasta"]):
        if not os.path.exists(fl):
            sys.exit(f'Error: File "{fl}" specified for "{arg_name}" does not exist.')

    if args.iter_th < 0:
        sys.exit(f'The number of TransFix iterations ("--iter") must be positive')

    # Create output folder if it doesn't exist
    if not os.path.isdir(args.outpath):
        os.makedirs(args.outpath)

    # Record executed command
    command = " ".join(sys.argv)
    print(f"\n{command}\n", flush=True)

    # Run analysis
    try:
        transfix_gtf = transfix_main(args.gtf, args.fasta, args.outpath, args.outname, chimeric=None, iter_th=args.iter_th)

        # Annotate additional features (UTR, start/stop codons)
        _ = add_features_to_gtf(transfix_gtf)

    except SystemExit as err:
        # Valid for python 3.5+
        print("".join(traceback.TracebackException.from_exception(err).format()))
    except Exception as err:
        print("".join(traceback.TracebackException.from_exception(err).format()))
        sys.exit(f"{err}")

    clean_log(logfile)
