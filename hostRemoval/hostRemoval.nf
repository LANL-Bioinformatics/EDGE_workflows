#!/usr/bin/env nextflow
//to run: nextflow run hostRemoval.nf -params-file [JSON parameter file]

//overwriteable parameters
params.inputFiles = null
params.ref = "NO_FILE" //setup for NF optional input file pattern
params.outDir = null
params.prefix = "host_clean" //perl host_removal script default
params.fasta = null
params.host = null
params.similarity = null
params.bwaMemOptions = null
params.cpu = null
params.h = null



process hostRemoval {
    debug true
    publishDir ".", mode: 'copy'

    input:
    path "reads"
    path ref

    output:
    //-fasta option does work, but files still have .fastq extensions 
    path "$params.outDir/${params.prefix}.{1,2,unpaired}.fastq"
    path "$params.outDir/${params.prefix}.mapping?E.log"
    path "$params.outDir/${params.prefix}.stats.txt"
    path "$params.outDir/host.fast{a,q}", optional: true

    script:
    
    def outDir = params.outDir != null ? "-o $params.outDir " : ""
    def files = (params.inputFiles != null && params.inputFiles.size() > 1) ? "-p reads1 reads2 " : "-u reads "

    def refFile = ref.name != "NO_FILE" ? "-ref $ref " : ""
    def prefix = params.prefix != "host_clean" ? "-prefix $params.prefix " : ""
    def fasta = params.fasta != null ? "-fasta $params.fasta " : ""
    def host = params.host != null ? "-host $params.host " : ""
    def similarity = params.similarity != null ? "-similarity $params.similarity " : ""
    def bwaMemOptions = params.bwaMemOptions != null ? "-bwaMemOptions \"$params.bwaMemOptions\" " : ""
    def cpu = params.cpu != null ? "-cpu $params.cpu " : ""
    
    """
    host_reads_removal_by_mapping.pl\
    $refFile\
    $prefix\
    $fasta\
    $host\
    $bwaMemOptions\
    $files\
    $outDir
    """
}

workflow {
    if (params.h != null) {
        "perl host_reads_removal_by_mapping.pl -help".execute().text
    }
    else {
        "touch NO_FILE".execute().text
        providedRef = file(params.ref, checkIfExists:true)
        hostRemoval(channel.fromPath(params.inputFiles, relative:true).collect(), providedRef)
        "rm NO_FILE".execute().text
    }
}