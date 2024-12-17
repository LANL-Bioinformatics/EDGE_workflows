import argparse
from pathlib import Path
import subprocess
import sys

def getEnaBrowserTools(info:dict, run_acc: str, args: argparse.Namespace):
    OUTDIR = args.outdir
    fastq_dir = Path(OUTDIR, 'sra2fastq_temp')

    sys.stderr.write(f"Retrieving FASTQ for {run_acc} from ENA Browser Tools...\n")

    cmd = ["enaDataGet", "-f", "fastq", "-d", fastq_dir, run_acc]

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        sys.stderr.write(f"Failed to run fasterq-dump from {run_acc}.\n")
        return "failed"
    

    sys.stderr.write("Done with fastq-dump.\n")

    return "success"
