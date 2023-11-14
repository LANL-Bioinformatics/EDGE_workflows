import subprocess
import os
import shutil

def index_reference(ref):
    subprocess.run(['bwa', 'index', ref])

def run_bwa_mem_paired(ref, bwamem_opts, num_cpu, paired_read):
    command = ['bwa', 'mem', '-t', str(num_cpu), bwamem_opts, ref, paired_read]
    return subprocess.run(command, capture_output=True, text=True).stdout

def run_bwa_mem_unpaired(ref, bwamem_opts, num_cpu, unpaired_read):
    command = ['bwa', 'mem', '-t', str(num_cpu), bwamem_opts, ref, unpaired_read]
    return subprocess.run(command, capture_output=True, text=True).stdout

def check_files(ref, out_dir, paired_reads_files, unpaired_reads_files):
    paired_files = []
    unpaired_files = []
    for paired_read in paired_reads_files:
        # Check if the file exists and is not empty
        if os.path.exists(paired_read) and os.path.getsize(paired_read) > 0:
            paired_files.append(paired_read)
        else:
            print(f"Warning: Paired read file {paired_read} is missing or empty.")

    for unpaired_read in unpaired_reads_files:
        # Check if the file exists and is not empty
        if os.path.exists(unpaired_read) and os.path.getsize(unpaired_read) > 0:
            unpaired_files.append(unpaired_read)
        else:
            print(f"Warning: Unpaired read file {unpaired_read} is missing or empty.")

    return paired_files, unpaired_files

def process_alignment(result):
    # Implement logic to process the alignment result
    # Filter reads based on similarity cutoff and generate separate output files
    pass

def generate_statistics(host_reads, non_host_reads):
    # Implement logic to generate statistics based on alignment results
    # Count reads, calculate percentages, etc.
    pass

# Input parameters (replace with actual file paths)
paired_reads_files = ['/path/to/paired_read1.fastq', '/path/to/paired_read2.fastq']
unpaired_reads_files = ['/path/to/unpaired_read1.fastq', '/path/to/unpaired_read2.fastq']
prefix = "host_clean"
ref = "/users/218819/scratch/data/databases/human_chromosomes/all_chromosome.fasta"
out_dir = "/path/to/output_directory"
bwamem_opts = "-T 50"
num_cpu = 4
output_host = False
out_fasta = False
similarity_cutoff = 90

# Index reference genome
index_reference(ref)

# Check files and get valid paired and unpaired files
paired_reads, unpaired_reads = check_files(ref, out_dir, paired_reads_files, unpaired_reads_files)

# Process paired reads
for paired_read in paired_reads:
    alignment_output = run_bwa_mem_paired(ref, bwamem_opts, num_cpu, paired_read)
    process_alignment(alignment_output)

# Process unpaired reads
for unpaired_read in unpaired_reads:
    alignment_output = run_bwa_mem_unpaired(ref, bwamem_opts, num_cpu, unpaired_read)
    process_alignment(alignment_output)

# Generate statistics
generate_statistics(host_reads, non_host_reads)
