# getDdbjFastq.py
from pathlib import Path
import requests
import sys
import argparse

def getDdbjFastq(info: dict, run_acc: str, args: argparse.Namespace):
    OUTDIR = args.outdir
    sys.stderr.write(f"Retrieving FASTQ for {run_acc} from DDBJ...")

    platform = info['platform']
    library = info['library']
    exp_acc = info['exp_acc']
    sub_acc = info['sub_acc']
    sra_acc_first6 = sub_acc[:6]

    url = f"ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/{sra_acc_first6}/{sub_acc}/{exp_acc}"

    file_names = [f"{run_acc}.fastq"]
    if "illu" in platform.lower() and "pair" in library.lower():
        file_names.extend([f"{run_acc}_1.fastq", f"{run_acc}_2.fastq"])

    total_size = 0
    for file_name in file_names:
        sys.stderr.write(f"Downloading {url}/{file_name}...\n")
        filepath = Path(OUTDIR, "sra2fastq_temp", file_name)
        try:
            with requests.get(f"{url}/{file_name}", stream=True) as response:
                response.raise_for_status()
                with open(filepath, "wb") as file:
                    for chunk in response.iter_content(chunk_size=8192):
                        file.write(chunk)
            sys.stderr.write(f"{file_name} downloaded successfully.\n")

            total_size += filepath.getsize()
        except requests.exceptions.RequestException as e:
            sys.stderr.write(f"Failed to download {file_name}. Error: {str(e)}\n")
            return "failed"

    if total_size < 50:
        for file_name in file_names:
            filepath.unlink()
        sys.stderr.write("Failed to download FASTQ files from DDBJ.\n")
        return "failed"

    return "success"
