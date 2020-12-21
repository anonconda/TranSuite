"""
Created on Wed May 15 14:25:00 2019
@author: Juan C Entizne
@email: e.entizne[at]dundee.ac.uk
"""

import os
import sys
import time
import traceback
from argparse import ArgumentParser, RawTextHelpFormatter
from modules.auto_main import run_transuite
from lib.logger.logger import logger, clean_log


description = \
    "Description:\n" + \
    "This module perform the FindLORF, TransFix, and TransFeat analysis in tandem.\n"

parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=False)

parser.add_argument('--gtf',
                    dest="gtf", default=None,
                    help="Transcriptome annotation file in GTF format.")

parser.add_argument('--fasta',
                    dest="fasta", default=None,
                    help='Fasta file of the transcripts (exonic) nucleotide sequences.')

parser.add_argument('--outpath',
                    dest="outpath", default=None,
                    help="Path of the output folder.")

parser.add_argument('--outname',
                    dest="outname", default=None,
                    help="Prefix for the output files.")


parser.add_argument("--cds",
                    dest="cds_th", type=int, default=30,
                    help="Minimum number of amino-acids an ORF must have to be considered as a potential CDS. "
                         "Default: 30 AA.")

parser.add_argument("--iter",
                    dest="iter_th", type=int, default=5,
                    help="Maximum number of 'start-fixing & translation' cycles to identify alternative start sites. "
                         "Default: 5")

parser.add_argument("--pep",
                    dest="pep_th", type=int, default=100,
                    help="Minimum number of amino-acids a translation must have to be consider a peptide. "
                         "Default: 100 AA.")

parser.add_argument("--ptc",
                    dest="ptc_th", type=int, default=70,
                    help="Minimum CDS length percentage below which a transcript is considered "
                         "prematurely terminated (PTC). Default: 70%%.")


def main():

    args = parser.parse_args()

    # Check input arguments
    for arg_val, arg_name in zip([args.gtf, args.fasta, args.outpath, args.outname],
                                 ["--gtf", "--fasta", "--outpath", "--outname"]):
        if arg_val is None:
            sys.exit(f'ERROR: No information specified for argument "{arg_name}"')

    # Create logfile to track the analysis, overwrite it if it exist ("w+" mode)
    time_stamp = time.strftime("%Y%m%d-%H%M%S")
    logfile = os.path.join(args.outpath, f"{time_stamp}_{args.outname}_logfile_temp.out")
    logger(logfile, w_mode="w+")

    # Normalize paths
    args.gtf = os.path.abspath(os.path.normpath(args.gtf))
    args.fasta = os.path.abspath(os.path.normpath(args.fasta))
    args.outpath = os.path.abspath(os.path.normpath(args.outpath))

    for fl, arg_name in zip([args.gtf, args.fasta], ["--gtf", "--fasta"]):
        if not os.path.exists(fl):
            sys.exit(f'ERROR: File "{fl}" specified for "{arg_name}" does not exist')

    if not 0 <= args.ptc_th <= 100:
        sys.exit(f'The % value specified for "--ptc" must be a number between 0 and 100')

    if args.cds_th < 0:
        sys.exit(f'The minimum length of the CDS ("--cds") must be positive')

    if args.pep_th < 0:
        sys.exit(f'The minimum length of the peptide ("--pep") must be positive')

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
        run_transuite(args.gtf, args.fasta, args.outpath, args.outname,
                      args.iter_th, args.cds_th, args.pep_th, args.ptc_th)
    except SystemExit as err:
        # Valid for python 3.5+
        print("".join(traceback.TracebackException.from_exception(err).format()))
    except Exception as err:
        print("".join(traceback.TracebackException.from_exception(err).format()))
        sys.exit(f"{err}")

    clean_log(logfile)
