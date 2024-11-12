#!/usr/bin/env nextflow
//to run: nextflow run sra2fastq.nf -params-file [JSON parameter file]
//not supporting filesize or run count restrictions



process sraDownload {
    label "sra2fastq"
    tag "$accession"
    publishDir "${settings["outDir"]}/SRA_Download", mode: 'copy'

    //retries download in case of transient failure, then completes any downloads that didn't fail
    maxRetries 3
    errorStrategy { (task.attempt <= maxRetries) ? 'retry' : 'finish' }

    input: 

    val accession //single accession string
    val settings

    output:
    path "$accession/${accession}*.fastq.gz", emit: files
    path "$accession/${accession}_metadata.txt"
    path "$accession/sra2fastq_temp/*", optional: true //needed output?

    script: 
    //conditionally create command-line options based on non-empty parameters, for use in the command below
    def clean = settings["clean"] ? "--clean True" : "" 
    def platform_restrict = settings["fastqSource"] != null ? "--platform_restrict ${settings["fastqSource"]}" : ""

    //invoke sra2fastq.py with those options
    """
    sra2fastq.py $accession \
    $clean \
    $platform_restrict \
    """
}


workflow SRA2FASTQ {
    take:
    settings

    main:
    accessions_ch = channel.of(settings["accessions"])
    sraDownload(accessions_ch.flatten().unique(), settings)

    fastq = sraDownload.out.files

    emit:
    fastq
}