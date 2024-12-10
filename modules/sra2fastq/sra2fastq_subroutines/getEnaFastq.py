import hashlib
import argparse
from pathlib import Path
import sys
import requests
import os
        
def getEnaFastq(info:dict, run_acc: str, args: argparse.Namespace):
    OUTDIR = args.outdir
    
    sys.stderr.write(f"Retrieving FASTQ for {run_acc} from EBI-ENA...\n")

    for i in info["fastq_ftp"]:
        url = info["fastq_ftp"][i]['url']
        md5 = info["fastq_ftp"][i]['md5']
        size = info["fastq_ftp"][i]['size']
        filename = os.path.basename(url)
    
        # Download file using requests library
        sys.stderr.write(f"Downloading {url} ...\n")
        filepath = Path(OUTDIR, 'sra2fastq_temp', filename)
        
        response = requests.get(f"ftp://{url}", stream=True)
        response.raise_for_status()
        
        with open(filepath, 'wb') as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)

        # Check downloaded file size
        filesize = filepath.getsize()
        if not filesize or filesize != size:
            sys.stderr.write(f"{filepath} incomplete/corrupted -- file sizes mismatch.\n")
            return "failed"
        
        # Check md5 checksum
        with open(filepath, 'rb') as file:
            md5sum = hashlib.md5(file.read()).hexdigest()
        
        if md5sum != md5:
            sys.stderr.write(f"{filepath} corrupted -- md5 checksum mismatch.n")
            return "failed"
        
        sys.stderr.write("Done with EBI-ENA.\n")
    
    return "success"
