#!/usr/bin/env nextflow
//to run: nextflow run sra2fastq.nf -params-file [config file]


//these defaults are overriden by any config files or command-line options
params.clean = ""
params.platform_restrict = ""
params.filesize_restrict = ""
params.runs_restrict = ""
params.outdir = ""
params.accessions = "" 

accessions_ch = Channel.of(params.accessions)

process SRA2FASTQ {
    publishDir "$params.outdir"

    input: 
    val accessions //space-separated string of accessions
    //TODO: go from JSON array input to this input 

    output:
    path "*/*.fastq.gz", emit: fastq_files //TODO: check for the appropriate number (and names?) of output files
    path "*/*_metadata.txt", emit: metadata_files
    path "*/sra2fastq_temp/*", emit: temp_files, optional: true
    //TODO: allow for discovery of hidden "finished" file
    //path "*/.finished", emit: finished_files

    script: 
    //conditionally create command-line options based on non-empty parameters, for use in the command below
    def clean = params.clean != "" ? "--clean True" : "" 
    def platform_restrict = params.platform_restrict != "" ? "--platform_restrict $params.platform_restrict" : ""
    def filesize_restrict = params.filesize_restrict != "" ? "--filesize_restrict $params.filesize_restrict" : ""
    def runs_restrict = params.runs_restrict != "" ? "--runs_restrict $params.runs_restrict" : ""

    //invoke sra2fastq.py with those options
    """
    sra2fastq.py $accessions \
    $clean \
    $platform_restrict \
    $filesize_restrict \
    $runs_restrict
    """
}


workflow {
    fastq_ch = SRA2FASTQ(accessions_ch.flatten())
    }