# define a function to download and merge fastq files
import os
from pathlib import Path
import re
import sys
import argparse
import shutil
import gzip

from .condaSRA2FastERq import condaSRA2FastERq
from .condaSRA2Fastq import condaSRA2Fastq
from .getAccType import getAccType
from .getDdbjFastq import getDdbjFastq
from .getEnaBrowserTools import getEnaBrowserTools
from .getEnaFastq import getEnaFastq
from .getReadInfo import getReadInfo
from .getSraFastq import getSraFastq
from .getSraFastqToolkits import getSraFastqToolkits
from .pigZip import pigZip


def downloadAndMergeFastq(acc: str, args: argparse.Namespace, finished: Path):
    OUTDIR = args.outdir
    RUNS_RESTRICT = args.runs_restrict
    PLAT_RESTRICT = args.platform_restrict
    SIZE_RESTRICT = args.filesize_restrict
    # get read information from NCBI-SRA / EBI-ENA
    read_info = getReadInfo(acc, {}, getAccType(acc))
    if not read_info.get(acc):
        sys.exit(f"ERROR: No sequence found. Please check if {acc} is a valid SRA/ERA/DRA number or your internet connection.\n")

    # initialize variables for tracking the total size and number of runs downloaded
    total_size = 0
    total_runs = 0

    # iterate over the runs associated with the given accession number
    for run_acc in read_info[acc]:
        # check if the number of runs exceeds the limit
        total_runs += 1
        if RUNS_RESTRICT and RUNS_RESTRICT < total_runs:
            sys.exit(f"ERROR: Run(s) exceed the limit ({RUNS_RESTRICT} MB).\n")

        # check if the platform is allowed
        platform = read_info[acc][run_acc]['platform']
        if PLAT_RESTRICT and not re.search(PLAT_RESTRICT, platform, re.IGNORECASE):
            sys.stderr.write(f"WARN: {platform} platform detected. Only {PLAT_RESTRICT} is allowed.\n")
            continue

        # download the fastq files
        download_functions = [condaSRA2Fastq, condaSRA2FastERq, getSraFastqToolkits, getSraFastq, getDdbjFastq, getEnaBrowserTools, getEnaFastq]
        for download_function in download_functions:
            dl_status = download_function(read_info[acc][run_acc], run_acc, args)
            if dl_status == 'success':
                break

        if dl_status == 'failed':
            sys.exit("ERROR: Please check your internet connection.\n")
        
        if SIZE_RESTRICT and SIZE_RESTRICT < total_size/1024/1024:
            raise ValueError(f"ERROR: downloaded file size exceed limitation ({SIZE_RESTRICT} MB).\n")
        else:
            sys.stderr.write(f"Succesfully downloaded {run_acc}.\n")
        
        pigZip(Path(OUTDIR, 'sra2fastq_temp'))


    for run_acc in read_info[acc]:
        # merging fastqs in multiple runs into a single file
        sys.stderr.write("Merging fastq files.\n")

        for file_suffix in [ ".fastq.gz", "_1.fastq.gz", "_2.fastq.gz"]:
            input_file = Path(OUTDIR, "sra2fastq_temp", f"{run_acc}{file_suffix}")
            output_file = Path(OUTDIR, "sra2fastq_temp", "merged", f"{acc}{file_suffix}")

            if input_file.is_file():
                with gzip.open(input_file, "r") as input_fastq, gzip.open(output_file, "a") as merged_file:
                        shutil.copyfileobj(input_fastq, merged_file)
                total_size += input_file.stat().st_size
                os.unlink(input_file)
        with open(finished, 'w') as fp:
            pass

    sys.stderr.write(f"Finished downloading acc# {acc}.\n")
