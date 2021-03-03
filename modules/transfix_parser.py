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
from lib.tools.input_tools import check_input


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
                    help="Path of the output folder.")

parser.add_argument('--outname',
                    dest="outname",
                    help="Prefix for the output files.")

parser.add_argument('--chimeric',
                    dest="chimeric", default=None,
                    help="Table indicating chimeric genes in the annotation.")

def main():
    args = parser.parse_args()

    # Check and sanitize user input
    args = check_input(args)

    # Create logfile to track the analysis, overwrite it if it exist ("w+" mode)
    time_stamp = time.strftime("%Y%m%d-%H%M%S")
    logfile = os.path.join(args.outpath, f"{time_stamp}_{args.outname}_logfile_temp.out")
    logger(logfile, w_mode="w+")

    # Record executed command
    command = " ".join(sys.argv)
    print(f"\n{command}\n", flush=True)

    # Run analysis
    try:
        transfix_gtf = transfix_main(args.gtf, args.fasta, args.outpath, args.outname,
                                     iter_th=args.iter_th, chimeric=args.chimeric)

        # Annotate additional features (UTR, start/stop codons)
        _ = add_features_to_gtf(transfix_gtf)

    except SystemExit as err:
        # Valid for python 3.5+
        print("".join(traceback.TracebackException.from_exception(err).format()))
    except Exception as err:
        print("".join(traceback.TracebackException.from_exception(err).format()))
        sys.exit(f"{err}")

    clean_log(logfile)
