import os
import sys


def short_path(gtf_path):

    return os.path.join("..", os.path.join(os.path.basename(os.path.dirname(gtf_path)), os.path.basename(gtf_path)))


def sizeof_fmt(num, suffix='B'):

    # Source: https://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size

    for unit in ['','K','M','G','T','P','E','Z']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Y', suffix)


class Tee(object):

    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)

    def flush(self):
        for f in self.files:
            f.flush()


def logger(logfile, w_mode="a+"):

    if not os.path.isfile(logfile):
        # If folder tree doesn't exist, create it
        if not os.path.isdir(os.path.dirname(logfile)):
            os.makedirs(os.path.dirname(logfile))
        # Create logfile if it doesn't exist
        with open(logfile, "w+") as fh:
            pass

    logfile = open(logfile, w_mode)
    backup = sys.stdout
    sys.stdout = Tee(sys.stdout, logfile)

    # To print only to stdout
    stdout_only = False
    if stdout_only:
        sys.stdout = backup


def clean_log(logfile):

    lines_seen = set()

    new_log = logfile.replace("_temp.out", ".out")

    with open(new_log, "w+") as fh_out, open(logfile, "r") as log_fh:
        for line in log_fh:
            # Keep the the new lines that separate block of analysis
            if line == "\n":
                fh_out.write(line)

            if line not in lines_seen:  # not a duplicate
                fh_out.write(line)
                lines_seen.add(line)
