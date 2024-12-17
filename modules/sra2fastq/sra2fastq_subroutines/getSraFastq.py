import shutil
import argparse
from pathlib import Path
import sys
from .deinterleaveFastq import deinterleaveFastq
import requests
import bz2

def getSraFastq(info: dict, run_acc: str, args: argparse.Namespace):
    OUTDIR = args.outdir
    sys.stderr.write(f"Retrieving FASTQ for {run_acc} from NCBI SRA (online converting)...")
    platform = info['platform']
    library = info['library']
    url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc={run_acc}"
    output_bz = Path(OUTDIR, "sra2fastq_temp", f"{run_acc}.fastq.bz2")
    output_file = Path(OUTDIR, "sra2fastq_temp", f"{run_acc}.fastq")

    sys.stderr.write(f"Downloading {url}...")
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for failed requests
        with open(output_bz, 'wb') as file:
            file.write(response.content)
        
        # Read the compressed data and decode it to text
        with bz2.open(output_bz, "rt") as bz2_file:
            # Write the uncompressed text to a file
            with open(output_file, "w") as file:
                file.write(bz2_file.read())

    except requests.exceptions.RequestException as e:
        sys.stderr.write(f"Failed to download SRA file from {url}. Error: {str(e)}\n")
        return "failed"
    
    # Deinterleaving if paired-end reads
    if "illu" in platform.lower() and "pair" in library.lower():
        sys.stderr.write("Paired-end reads found. Deinterleaving...\n")
        output_file1 = Path(OUTDIR, "sra2fastq_temp", f"{run_acc}_1.fastq")
        output_file2 = Path(OUTDIR, "sra2fastq_temp", f"{run_acc}_2.fastq")
        di_flag = deinterleaveFastq(output_file, output_file1, output_file2, compress=True)
        
        if di_flag > 0:
            return "failed"
        sys.stderr.write("Done.\n")
        shutil.rm(output_file)
    
    sys.stderr.write("Done with NCBI SRA.\n")
    
    return "success"