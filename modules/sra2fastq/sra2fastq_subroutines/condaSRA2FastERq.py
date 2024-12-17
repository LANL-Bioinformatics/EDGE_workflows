import argparse
from pathlib import Path
import shutil
import subprocess
import sys

def condaSRA2FastERq(info: dict, run_acc: str, args: argparse.Namespace):
    OUTDIR = args.outdir
    SIZE_RESTRICT = args.filesize_restrict
    platform = info["platform"].lower()
    fastq_dir = Path(OUTDIR, 'sra2fastq_temp')

    # prefetch SRA file 
    # needs to be in OUTDIR 
    cmd = ["prefetch", run_acc, "-O", fastq_dir]
    if SIZE_RESTRICT != None:
        cmd.extend(["--max-size", int(SIZE_RESTRICT)])
    subprocess.run(cmd)
    
    # fasterq dump from SRA file
    options = []
    if "illu" in platform:
        options.append("--split-files")
    elif "solid" in platform:
        options.extend(["--split-files", "-B"])
    sys.stderr.write(f"Running fasterq-dump with options {' '.join(options)}...\n")

    cmd = ["fasterq-dump", run_acc, *options, "-O", fastq_dir]


    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        sys.stderr.write(f"Failed to run fasterq-dump from {run_acc}.\n")
        return "failed"
    
    # remove pre-fetch file
    shutil.rmtree(Path(fastq_dir, run_acc))


    sys.stderr.write("Done with fasterq-dump.\n")

    return "success"