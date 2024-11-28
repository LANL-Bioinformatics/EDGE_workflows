import argparse
from pathlib import Path
import subprocess
import os
import sys
import requests

def getSraFastqToolkits(info: dict, run_acc: str, args: argparse.Namespace):
    OUTDIR = args.outdir
    HTTP_PROXY = args.http_proxy

    sys.stderr.write(f"Retrieving FASTQ for {run_acc} with NCBI SRA Toolkit...\n")
    platform = info["platform"]
    url = info["url"]
    filename = f"{run_acc}.fastq"

    sys.stderr.write(f"Downloading {url}...\n")
    out_path = Path(OUTDIR, "sra2fastq_temp", filename)
    headers = {"User-Agent": "Mozilla/5.0"}

    if HTTP_PROXY is not None:
        proxies = {"http": HTTP_PROXY, "https": HTTP_PROXY}
    else:
        proxies = None

    try:
        response = requests.get(url, headers=headers, stream=True, proxies=proxies)
        response.raise_for_status()

        with open(out_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)

    except requests.exceptions.RequestException as e:
        sys.stderr.write(f"Failed to download SRA file from {url}. Error: {str(e)} \n")
        return "failed"

    sys.stderr.write("Done.")

    # check downloaded file
    filesize = out_path.getsize()
    if not filesize:
        sys.stderr.write(f"Failed to download SRA file from {url}.\n")
        return "failed"

    # dump fastq from SRA file
    options = []
    if "illu" in platform.lower():
        options.append("--split-files")
    elif "solid" in platform.lower():
        options.extend(["--split-files", "-B"])
    sys.stderr.write(f"Running fastq-dump with options {' '.join(options)}...\n")

    try:
        subprocess.run(
            ["fastq-dump", *options, "--outdir", Path(OUTDIR, "sra2fastq_temp"), out_path],
            check=True,
        )
    except subprocess.CalledProcessError:
        sys.stderr.write(f"Failed to run fastq-dump from {out_path}.\n")
        return "failed"

    sys.stderr.write("Done using fastq-dump\n")

    return "success"
