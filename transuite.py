"""
Created on Wed May 15 14:25:00 2019
@author: Juan C Entizne
@email: e.entizne[at]dundee.ac.uk
"""

import sys
import logging
import argparse
import modules.findlorf_parser as FindLORF
import modules.transfix_parser as TransFix
import modules.transfeat_parser as TransFeat
import modules.auto_parser as Auto


description = "Description:\n\n" + \
              "TranSuite is a suite of software logger for the identification, annotation, translation, and feature " \
              "characterization of annotated transcripts."

parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter,
                                 prog='TranSuite')

# Version
parser.add_argument('-v', '--version', action='version', version='0.2.2')
subparsers = parser.add_subparsers()


# FindLORF parser
FindLORFSubparser = subparsers.add_parser("FindLORF",
                                          parents=[FindLORF.parser],
                                          help="FindLORF find the longest ORF of a transcript "
                                                "and annotates it as its putative CDS.")
FindLORFSubparser.set_defaults(which="FindLORF")


# TransFix parser
TransFixSubparser = subparsers.add_parser("TransFix",
                                          parents=[TransFix.parser],
                                          help="TransFix fix the same start codon for all transcripts in a gene, "
                                               "translate them, and annotates the resulting CDS.")
TransFixSubparser.set_defaults(which="TransFix")


# TransFeat parser
TransFeatSubparser = subparsers.add_parser("TransFeat",
                                           parents=[TransFeat.parser],
                                           help="TransFeat infer coding-related characteristics from the annotate "
                                                "transcript features.")
TransFeatSubparser.set_defaults(which="TransFeat")


# Auto parser
AutoSubparser = subparsers.add_parser("Auto",
                                      parents=[Auto.parser],
                                      help="This module executes FindLORF, TransFix, and TransFeat in tandem.")
AutoSubparser.set_defaults(which="Auto")


# Setting logging preferences
logger = logging.getLogger(__name__)


def main():
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    try:
        args = parser.parse_args()
        if args.which == "FindLORF":
            FindLORF.parser = parser
            FindLORF.main()
        elif args.which == "TransFix":
            TransFix.parser = parser
            TransFix.main()
        # Setting the modules parsers
        elif args.which == "TransFeat":
            TransFeat.parser = parser
            TransFeat.main()
        elif args.which == "Auto":
            Auto.parser = parser
            Auto.main()

    except Exception:
        logger.error("Unknown error: {}".format(sys.exc_info()))
        sys.exit(1)


if __name__ == '__main__':
    main()
