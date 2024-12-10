import argparse
from pathlib import Path
import subprocess
import sys

from .pigZip import pigZip

def condaSRA2Fastq(info: dict, run_acc: str, args: argparse.Namespace):
    OUTDIR = args.outdir
    platform = info["platform"].lower()
    fastq_dir = Path(OUTDIR, 'sra2fastq_temp')

    # fastq dump from SRA file
    options = ["--split-3", "--skip-technical"]
    if "illu" in platform:
        options.append("--split-files")
    elif "solid" in platform:
        options.extend(["--split-files", "-B"])
    sys.stderr.write(f"Running fastq-dump with options {' '.join(options)}...\n")

    cmd = ["fastq-dump", run_acc, *options, "-O", fastq_dir]

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        sys.stderr.write(f"Failed to run fasterq-dump from {run_acc}.\n")
        return "failed"
    

    sys.stderr.write("Done with fastq-dump.\n")

    return "success"