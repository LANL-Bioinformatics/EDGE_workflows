#!/usr/bin/env nextflow
//to run: nextflow [OPT: -log /path/to/log file] run sra2fastq.nf -params-file [JSON parameter file]
//not supporting filesize or run count restrictions



process sraDownload {
    label "sra2fastq"
    tag "$accession"
    publishDir "${settings["outDir"]}/SRA_Download", mode: 'copy'

    //retries download in case of transient failure, then completes any processes that didn't fail
    maxRetries 3
    errorStrategy { (task.attempt <= maxRetries) ? 'retry' : 'finish' }

    input: 

    val accession //single accession string
    val settings

    output:
    path "$accession/${accession}.fastq.gz", emit: unpairedSRA, optional:true
    path "$accession/${accession}_{1,2}.fastq.gz", emit: pairedSRA, optional:true
    path "$accession/${accession}_metadata.txt"
    path "$accession/sra2fastq_temp/*", optional: true //needed output?

    script: 
    //conditionally create command-line options based on non-empty parameters, for use in the command below
    def clean = settings["clean"] != null ? "--clean True" : "" 
    def platform_restrict = settings["platformRestrict"] != null ? "--platform_restrict ${settings["platformRestrict"]}" : ""

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

    paired = sraDownload.out.pairedSRA
    unpaired = sraDownload.out.unpairedSRA

    emit:
    paired
    unpaired
}