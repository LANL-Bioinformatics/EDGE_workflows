#!/usr/bin/env nextflow

process r2c {
    debug true
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis/readsMappingToContig",
        mode: 'copy'
    )
    input:
    path paired
    path unpaired
    path contigs

    output:
    path "*.sort.bam", emit: sortedBam
    path "*.alnstats.txt"
    path "*_coverage.table", emit: cov_table
    path "*_plots.pdf"
    path "Final_contigs.fasta", emit: contig_file
    path "mapping.log", emit: logFile

    script:
    def outPrefix = params.prefix!=null ? "$params.prefix" : "readsToContigs"
    def paired = paired.name != "NO_FILE" ? "-p \'${paired[0]} ${paired[1]}\' " : ""
    def unpaired = unpaired.name != "NO_FILE2" ? "-u $unpaired " : ""
    def cutoff = params.assembledContigs ? "-c 0 " : "-c 0.1 "
    def cpu = params.cpus != null ? "-cpu $params.cpus " : ""
    def max_clip = params.r2g_max_clip != null ? "-max_clip $params.r2g_max_clip " : ""


    def ont_flag = (params.fastq_source != null && params.fastq_source.equalsIgnoreCase("nanopore")) ? "-x ont2d " : ""
    def pb_flag = (params.fastq_source != null && params.fastq_source.equalsIgnoreCase("pacbio")) ? "-x pacbio " : ""

    def aligner_options = ""
    if(params.r2c_aligner =~ "bowtie") {
        def bowtie_options = params.r2c_aligner_options.replaceAll("-p\\s*\\d+","")
        if(!(bowtie_options =~ /-k/)) {
            bowtie_options += " -k 10 "
        }
        aligner_options = "-aligner bowtie -bowtie_options \'$bowtie_options\'"
    }
    else if(params.r2c_aligner =~ "bwa") {
        def bwa_options = params.r2c_aligner_options.replaceAll("-t\\s*\\d+","")
        if (ont_flag != "") {
            unpaired = unpaired.replaceAll("-u ","-long ")
            bwa_options += ont_flag
        }
        if (pb_flag != "") {
            unpaired = unpaired.replaceAll("-u ","-long ")
            bwa_options += pb_flag
        }
        aligner_options = "-aligner bwa -bwa_options \'$bwa_options\'"
    }
    else if (params.r2c_aligner =~ "minimap") { 
        def minimap_options = params.r2c_aligner_options.replaceAll("-t\\s*\\d+","")
        if(ont_flag != "" || pb_flag != "") {
            unpaired = unpaired.replaceAll("-u ","-long ")
        }
        if(pb_flag != "") {
            minimap_options += " -x map-pb "
        }
        aligner_options = "-aligner minimap2 -minimap2_options \'$minimap_options\'"
    }


    """
    runReadsToContig.pl \
    $cutoff\
    $cpu\
    $paired\
    $unpaired\
    -d . -pre $outPrefix\
    -ref $contigs \
    $max_clip\
    $aligner_options &> mapping.log


    awk \'{print \$1\"\\t\"\$4}\' ${outPrefix}_coverage.table > magnitudes.txt
    """

}

process r2c_jsonTable {
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis/readsMappingToContig",
        mode: 'copy',
        pattern: "*_coverage.table.json"
    )
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis",
        mode: 'copy',
        pattern: "*stats.{pdf,txt}"
    )
    input:
    path cov_table
    path contigFile

    output:
    path "contigs_stats.txt"
    path "contigs_stats.pdf"
    path "*_coverage.table.json"

    script:
    def rowLimit = params.rowLimit != null ? "$params.rowLimit" : "3000"
    def outPrefix = params.prefix!=null ? "$params.prefix" : "readsToContigs"
    
    """
    tab2Json_for_dataTable.pl -project_dir . -mode contig -limit $rowLimit  \
    ${outPrefix}_coverage.table > ${outPrefix}_coverage.table.json

    contig_stats.pl -p $contigFile > contigs_stats.txt
    """
}

process extractUnmapped {
    publishDir(
        path:"$params.outDir/AssemblyBasedAnalysis/readsMappingToContig/",
        mode: 'copy',
        overwrite: true
    )
    input:
    path bamFile
    path logFile


    output:
    path "mapping.log"
    path "Unmapped*.fastq"

    script:
    """
    echo "Extract unmapped reads" >> $logFile
    bam_to_fastq.pl -unmapped -prefix Unmapped ${bamFile} >> $logFile
    """

}

workflow {
    "mkdir nf_assets".execute().text
    "touch nf_assets/NO_FILE".execute().text
    "touch nf_assets/NO_FILE2".execute().text
    paired_ch = channel.fromPath(params.pairFile, checkIfExists:true).collect()
    unpaired_ch = channel.fromPath(params.unpairFile, checkIfExists:true)
    contig_ch = channel.fromPath(params.contigFile, checkIfExists:true)

    r2c(paired_ch, unpaired_ch, contig_ch)
    r2c_jsonTable(r2c.out.cov_table, r2c.out.contig_file)
    if(params.extractUnmapped) {
        extractUnmapped(r2c.out.sortedBam, r2c.out.logFile)
    }

}
