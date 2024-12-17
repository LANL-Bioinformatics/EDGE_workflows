import re
import sys
import requests

def getReadInfo(acc: str, read_info: dict, sra_type: str):
    sys.stderr.write("Retrieving run(s) information from NCBI-SRA...")

    # get info from NCBI-SRA
    url0 = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={acc}&usehistory=y"
    response = requests.get(url0)
    web_result = response.text
    webenv = re.findall('<WebEnv>(\S+)<\/WebEnv>', web_result)[0]
    key = re.findall('<QueryKey>(\S+)<\/QueryKey>', web_result)[0]

    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&rettype=runinfo&query_key={key}&WebEnv={webenv}&retmode=text"
    sys.stderr.write(f"Retrieving run acc# from NCBI-SRA {url0} {url}...")
    response = requests.get(url)
    lines = response.text.strip().split("\n")[1:]
    sra_num_runs = len(lines)
    sys.stderr.write(f"{sra_num_runs} run(s) found from NCBI-SRA.")

    for line in lines:
        fields = line.split(',')
        sub_acc, exp_acc, run_acc, size_MB, platform, library, url = fields[42], fields[10], fields[0], fields[7], fields[18], fields[15], fields[9]
        if sra_type == "ByRun" and not re.search(acc, run_acc, re.IGNORECASE):
            continue
        if not size_MB:
            sys.stderr.write(f"Run {run_acc} has size 0 MB")
        read_info[acc] = {}
        read_info[acc][run_acc] = {}
        read_info[acc][run_acc] = {
            "exp_acc": exp_acc,
            "sub_acc": sub_acc,
            "platform": platform,
            "library": library,
            "url": url
        }
        
    # get info from EBI-ENA when NCBI-SRA fails
    sys.stderr.write("Retrieving run(s) information from EBI-ENA...\n")

    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run&fields=run_accession,submission_accession,study_accession,experiment_accession,instrument_platform,library_layout,fastq_ftp,fastq_md5,fastq_bytes"
    sys.stderr.write(f"Retrieving run acc# from EBI-ENA {url}...\n")
    response = requests.get(url)
    web_result = response.text

    if not re.search(r"^study_accession|^run_accession", web_result) and not sra_num_runs:
        raise Exception(f"ERROR: Failed to retrieve sequence information for {acc} from both SRA and ENA database.")
    elif not re.search(r"^study_accession|^run_accession", web_result) and sra_num_runs:
        sys.stderr.write(f"WARNING: {acc} only found in SRA database. The data may be not synchronized among INSDC yet.\n")
    elif re.search(r"^study_accession|^run_accession", web_result) and not sra_num_runs:
        sys.stderr.write(f"WARNING: {acc} only found in ENA database. The data may be not synchronized among INSDC yet.\n")

    #run_accession	fastq_ftp	fastq_bytes	fastq_md5	submitted_ftp	submitted_bytes	submitted_md5	sra_ftp	sra_bytes	sra_md5

    # split the string only once
    lines = web_result.split('\n')[:-1]
    # sys.stderr.write(f"{lines=}")
    # extract the number of runs found and sys.stderr.write it to stderr
    sys.stderr.write(f"{len(lines)-1} run(s) found from EBI-ENA.\n")

    if len(lines)-1 > 0:
        # compile the regular expression outside the loop
        pattern = re.compile(r"^study_accession|^run_accession")

        # extract the fields once before entering the loop
        # fields_indices = [1, 3, 0, 4, 5, 6, 7, 8]

        # process each line of the result
        for line in lines:
            if pattern.match(line):
                continue
            fields = line.strip().split('\t')

            if len(fields) >= 8:
                # extract the relevant fields
                sub_acc, exp_acc, run_acc, platform, library = [fields[i] for i in [1, 3, 0, 4, 5]] 

                url, md5, size = [fields[i].split() for i in [6,7,8]]

                if sra_type == "ByRun" and acc.lower() != str(run_acc).lower():
                    continue
                
                read_info[acc][run_acc] = {
                    "platform": platform,
                    "exp_acc": exp_acc,
                    "sub_acc": sub_acc,
                    "library": library,
                    "url": url,
                    "md5": md5,
                    "size": size

                }
                read_info[acc][run_acc]["fastq_ftp"] = {}
                for i in range(len(url)):
                    read_info[acc][run_acc]["fastq_ftp"][i] = {
                        "url": url[i],
                        "md5": md5[i],
                        "size": size[i]
                    }
    return read_info




