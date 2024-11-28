#!/usr/bin/env nextflow
//to run: nextflow run hostRemoval.nf -params-file [JSON parameter file]


//clean input FASTQ files of reads that map to provided host references, one process per given host
process hostRemoval {
    label "hostRemoval"
    publishDir(
        path: "${settings["outDir"]}/HostRemoval",
        mode: 'copy'
    )
    
    tag "${ref.name.take(ref.name.lastIndexOf('.'))}"

    input:
    val settings
    path paired
    path unpaired
    each path(ref)

    output:
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/${ref.name.take(ref.name.lastIndexOf('.'))}.clean.1.fastq", emit: cleaned1, optional:true
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/${ref.name.take(ref.name.lastIndexOf('.'))}.clean.2.fastq", emit: cleaned2, optional:true
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/${ref.name.take(ref.name.lastIndexOf('.'))}.clean.unpaired.fastq", emit: cleanedSingleton
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/${ref.name.take(ref.name.lastIndexOf('.'))}.clean.mapping?E.log"
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/${ref.name.take(ref.name.lastIndexOf('.'))}.clean.stats.txt", emit: cleanstats
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/host.fastq", emit:hostReads

    script:
    
    def pairedFiles = paired.name != "NO_FILE" ? "-p $paired" : ""
    def unpairedFiles = unpaired.name != "NO_FILE2" ? "-u $unpaired" : ""

    def refFile = ref.name != "NO_FILE" ? "-ref $ref " : ""
    def prefix = "-prefix ${ref.name.take(ref.name.lastIndexOf('.'))}.clean "
    def similarity = settings["similarity"] != null ? "-similarity ${settings["similarity"]} " : ""
    def minScore = settings["bwaMemOptions"] != null ? "${settings["bwaMemOptions"]} " : "-T 50 "
    def ontFlag = settings["fastqSource"].equalsIgnoreCase("nanopore") ? "-x ont2d " : ""
    ontFlag = settings["fastqSource"].equalsIgnoreCase("pacbio") ? "-x pacbio " : ontFlag
    minScore = ontFlag != "" ? "-T ${settings["minLen"]} " : minScore
    def bwaMemOptions = "-bwaMemOptions \"${ontFlag} ${minScore}\" "
    def cpu = settings["cpus"] != null ? "-cpu ${settings["cpus"]} " : ""
    
    """
    host_reads_removal_by_mapping.pl\
    $refFile\
    $prefix\
    $cpu\
    -host \
    $bwaMemOptions\
    $pairedFiles\
    $unpairedFiles\
    -o .

    mkdir ${ref.name.take(ref.name.lastIndexOf('.'))}
    mv host.fastq ./${ref.name.take(ref.name.lastIndexOf('.'))}
    mv ${ref.name.take(ref.name.lastIndexOf('.'))}.* ./${ref.name.take(ref.name.lastIndexOf('.'))}
    """
}

//merge cleaned FASTQ files into files cleaned of ALL reads mapping to ANY provided host reference
process collectCleanPairedReads {
    label "hostRemoval"
    publishDir(
        path: "${settings["outDir"]}/HostRemoval",
        mode: 'copy'
    )

    input:
    val settings
    path cleanedFiles1
    path cleanedFiles2
    path(hostFiles, stageAs: 'host?.fastq')

    output:
    path "hostclean.{1,2}.fastq", emit: paired
    path "merged_host_unique.fastq", emit: hostMerged
    
    script:
    """
    seqkit common $cleanedFiles1 -n > hostclean.1.fastq
    seqkit common $cleanedFiles2 -n > hostclean.2.fastq
    cat $hostFiles > merged_host.fastq
    seqkit rmdup -n merged_host.fastq > merged_host_unique.fastq
    """
}

process collectCleanPairedReadsOneHost {
    label "hostRemoval"
    publishDir(
        path: "${settings["outDir"]}/HostRemoval",
        mode: 'copy'
    )

    input:
    val settings
    path cleanedFiles

    output:
    path "hostclean.{1,2}.fastq", emit:paired

    
    script:
    """
    mv $cleanedFiles hostclean.${cleanedFiles[0].name.tokenize('.')[-2]}.fastq
    """
}

//Concatenate leftover unpaired reads that didn't map to a host reference, and remove any name duplicates (i.e., leftovers appearing in multiple cleanings)
process collectCleanSingleReads {
    label "hostRemoval"
    publishDir(
        path: "${settings["outDir"]}/HostRemoval",
        mode: 'copy'
    )
    
    input:
    val settings
    path remainingUnpairedReads

    output:
    path "hostclean.unpaired.fastq", emit:unpaired

    script:
    """
    cat $remainingUnpairedReads > hostclean.unpaired_dups.fastq
    seqkit rmdup -n hostclean.unpaired_dups.fastq > hostclean.unpaired.fastq
    """
}

process hostRemovalStats {
    label "hostRemoval"
    publishDir "${settings["outDir"]}/HostRemoval", mode: 'copy'

    input:
    val settings
    path stats
    path hostReads

    output:
    path "hostclean.stats.txt" //issue with counting here
    path "HostRemovalStats.pdf"

    script:
    """
    HOSTREADS=\$((\$(cat $hostReads | wc -l)/4))
    removal_stats.pl\
    -stats $stats\
    -hostReads \$HOSTREADS
    """
}


workflow HOSTREMOVAL{
    take:
    settings
    paired
    unpaired

    main:
    providedRef = channel.fromPath(settings["host"], checkIfExists:true)

    //remove host reads in parallel
    hostRemoval(settings, paired, unpaired, providedRef.collect())

    cleaned1_ch = hostRemoval.out.cleaned1.collect()
    cleaned2_ch = hostRemoval.out.cleaned2.collect()
    print("at host branch\n")
    //more than one host
    if (([] + settings["host"]).size() > 1) {
        print("in multihost branch\n")
        print(settings["host"] + "\n")
        print(settings["host"].size())
        print("\n")
        //merge clean paired-end reads (intersection)
        collectCleanPairedReads(settings, cleaned1_ch, cleaned2_ch, hostRemoval.out.hostReads.collect())
        paired = collectCleanPairedReads.out.paired
        //calculate overall stats and create PDF
        hostRemovalStats(settings, hostRemoval.out.cleanstats.collect(), collectCleanPairedReads.out.hostMerged)
        paired.view()
    } 
    else {
        print("in single host branch\n")
        //no need to merge if only reads from one host were removed
        paired = collectCleanPairedReadsOneHost(settings, cleaned1_ch.concat(cleaned2_ch)).collect()
        //calculate overall stats and create PDF
        hostRemovalStats(settings, hostRemoval.out.cleanstats.collect(), hostRemoval.out.hostReads)
        paired.view()
    }
    print("after host branch \n")
    //merge clean unpaired reads (removing any duplicates by read name)
    unpaired = collectCleanSingleReads(settings, hostRemoval.out.cleanedSingleton.collect())

    emit:
    paired
    unpaired

    
}