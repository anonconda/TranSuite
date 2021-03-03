import os
import sys
import shutil


def check_input(args):

    # 0) Check that required arguments are present
    for arg_val, arg_name in zip([args.gtf, args.fasta, args.outpath, args.outname],
                                 ["--gtf", "--fasta", "--outpath", "--outname"]):
        if arg_val is None:
            sys.exit(f'ERROR: No information specified for argument "{arg_name}"')

    # Sanitize paths
    paths_list = [args.gtf, args.fasta, args.outpath, args.outname, args.chimeric]
    for ix, dir_path in enumerate(paths_list):
        if dir_path:
            # Remove characters from path that may cause issues, such as quotation marks (i.e. ")
            for quote_mark in ['\'', '"', '`']:
                if dir_path.startswith(quote_mark) and dir_path.endswith(quote_mark):
                    dir_path = dir_path.strip(quote_mark)

            # If the path contains whitespaces, add quotes to the path
            if ' ' in dir_path:
                paths_list[ix] = f'"{dir_path}"'

    # 1) Normalize and check paths
    args.gtf = os.path.abspath(os.path.normpath(args.gtf))
    args.fasta = os.path.abspath(os.path.normpath(args.fasta))
    args.outpath = os.path.abspath(os.path.normpath(args.outpath))

    try:
        args.chimeric = os.path.abspath(os.path.normpath(args.chimeric))
    except:
        pass

    # Create output folder if it doesn't exist
    dir_out = os.path.join(args.outpath, f"{args.outname}_TranSuite_output")
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    args.outpath = dir_out

    for fl, arg_name in zip([args.gtf, args.fasta], ["--gtf", "--fasta"]):
        if not os.path.exists(fl):
            sys.exit(f'ERROR: File {fl} ({arg_name} argument) does not exist')

    # 2) Check if threshold values are within valid ranges
    try:
        if not 0 <= args.ptc_th <= 100:
            sys.exit(f'The % value specified for "--ptc" must be a number between 0 and 100')
    except AttributeError:
        pass

    try:
        if args.cds_th < 0:
            sys.exit(f'The minimum length of the CDS ("--cds") must be positive')
    except AttributeError:
        pass

    try:
        if args.pep_th < 0:
            sys.exit(f'The minimum length of the peptide ("--pep") must be positive')
    except AttributeError:
        pass

    try:
        if not 0 < args.iter_th < 6:
            sys.exit(f'The number of TransFix iterations ("--iter") must be between 1 and 5')
    except AttributeError:
        pass

    return args


def create_project_dir(outpath, outname, overwrite):

    # Unused function TODO check before deleting it
    paths = []

    # Main output folder
    dir_out = os.path.abspath(os.path.normpath(os.path.join(outpath, f"{outname}_RTDmaker_output")))

    # Check if the output folder has to be re-initialize
    if overwrite and os.path.isdir(dir_out):
        # Force the user to confirm selection
        while True:
            user_ans = input(f"Are you sure you want to delete the folder: {dir_out}? [Y/N]")
            if user_ans.upper() == "Y" or user_ans.upper() == "YES":
                shutil.rmtree(dir_out)
                break
            elif user_ans.upper() == "N" or user_ans.upper() == "NO":
                print(f"Analysis will continue with the existing files.")
                break
            else:
                print(f"'{user_ans}' is not among the valid options. Please answer YES or NO.")
                continue

    # Generate the whole project directory structure and create a dictionary to handle the multiple subfolder paths
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)

    paths_dt = {}
    for (dir_name, dir_path) in paths:
        if not os.path.isdir(dir_path):
            os.makedirs(dir_path)
        paths_dt[dir_name] = dir_path

    return paths_dt
