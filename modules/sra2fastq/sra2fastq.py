#!/opt/conda/bin/python

import os
from pathlib import Path
import shutil
import subprocess
import sra2fastq_subroutines as s
import sys
import argparse
import glob

def main():
    parser = argparse.ArgumentParser(description=f'''A script retrieves sequence project in FASTQ files from 
        NCBI-SRA/EBI-ENA/DDBJ database using `curl` or `wget`. Input accession number
        supports studies (SRP*/ERP*/DRP*), experiments (SRX*/ERX*/DRX*), 
        samples (SRS*/ERS*/DRS*), runs (SRR*/ERR*/DRR*), or submissions 
        (SRA*/ERA*/DRA*).

        [USAGE]
            {sys.argv[0]} [OPTIONS] <Accession#> (<Accession# 2> <Accession# 3>...)
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('accessions', metavar='Accession#', nargs='*', help='accession number')

    parser.add_argument('--outdir', '-o',  type=str, default='.', help='Output directory (default: %(default)s)')
    parser.add_argument('--clean', '-c', type=bool, default=False, help='Clean up temp directory')
    parser.add_argument('--platform_restrict', '-pr', type=str, default=None, help='Only allow a specific platform (default: %(default)s)')
    parser.add_argument('--filesize_restrict', '-fr', type=int, default=None, help='(in MB) Only allow to download less than a specific total size of files. (default: %(default)s)')
    parser.add_argument('--runs_restrict', '-r',  type=int, default=None, help='Only allow download less than a specific number of runs. (default: %(default)s)')

    args = parser.parse_args()

    # defaults
    parser.add_argument("--user_proxy", default=None)
    parser.add_argument("--no_proxy", default=False)
    parser.add_argument("--http_proxy", default=os.environ.get('HTTP_PROXY') or os.environ.get('HTTP_PROXY'))
    parser.add_argument("--ftp_proxy", default=os.environ.get('FTP_PROXY') or os.environ.get('FTP_PROXY'))
    
    args = parser.parse_args()

    args.http_proxy = "" if args.no_proxy else args.http_proxy
    args.ftp_proxy = "" if args.no_proxy else args.ftp_proxy
    args.http_proxy = f"--proxy '{args.user_proxy}' " if args.user_proxy else args.http_proxy
    args.ftp_proxy = f"--proxy '{args.user_proxy}' " if args.user_proxy else args.ftp_proxy

    args.outdir = Path(args.outdir)
    outdir = args.outdir   


    print(args.accessions)
    if len(args.accessions) == 0:
        print("No accessions run")
        quit()

    for accession in args.accessions:
        # run tool if accession number is valid
        if s.isValidAcc(accession):
            args.outdir = Path(outdir, accession)
            temp_dir = Path(args.outdir, "sra2fastq_temp")
            merged_dir = Path(args.outdir, "sra2fastq_temp", "merged")
            meta_file = Path(args.outdir, f"{accession}_metadata.txt")
            finished = Path(args.outdir, ".finished")

            if os.path.isfile(finished):
                print("Accession already downloaded")
                quit()

            # init out directories
            args.outdir.mkdir(parents=True, exist_ok=True)
            temp_dir.mkdir(exist_ok=True)
            merged_dir.mkdir(exist_ok=True)

            # get metadata
            sys.stderr.write("Downloading metadata from pysradb\n")
            subprocess.run(f"pysradb metadata --detailed {accession} > {meta_file}", shell=True)
            
            # run tool
            s.downloadAndMergeFastq(accession, args, finished)
            files = glob.glob(str(Path(args.outdir, "sra2fastq_temp", "merged", "*")))
            # move files into accession number folder and clean directories

            for file in files:
                shutil.move(file, args.outdir)
            if args.clean:
                sys.stderr.write("Cleaning up temporary directories \n")
                shutil.rmtree(Path(args.outdir, "sra2fastq_temp"))
            
            sys.stderr.write(f"Done with {accession} \n")
        
        else:
            sys.stderr.write(f"Error. {accession} not a valid accession number. \n")

    sys.stderr.write(f"Completed sra2fastq for {args.accessions} \n")

if __name__ == '__main__':
    main()