from pathlib import Path
import glob
import subprocess
import sys

def pigZip(source_dir):
    files = glob.glob(str(Path(source_dir, '*.fastq')))
    for file in files:
        try:
            subprocess.run(['pigz', file])
        except subprocess.CalledProcessError:
            sys.stderr.write("Failed to zip using pigz \n")