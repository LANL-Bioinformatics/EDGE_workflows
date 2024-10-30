#!/usr/bin/env nextflow

include {SRA2FASTQ} from './modules/sra2fastq/sra2fastq.nf'
include {COUNTFASTQ} from './modules/countFastq/countFastq.nf'
include {FAQCS} from './modules/runFaQCs/runFaQCs.nf'

workflow {

    //input specification

    pairedFiles = channel.fromPath(params.pairedFiles, checkIfExists:true)
    unpairedFiles = channel.fromPath(params.unpairedFiles, checkIfExists:true)

    if(params.modules.sra2fastq) {
        SRA2FASTQ(params.sra2fastq.plus(params.shared))
        pairedFiles = pairedFiles.concat(SRA2FASTQ.out.paired).flatten()
        unpairedFiles = unpairedFiles.concat(SRA2FASTQ.out.unpaired).flatten()
    }
    
    COUNTFASTQ(pairedFiles.collect(), unpairedFiles.collect())

    avgLen = COUNTFASTQ.out.avgReadLen
    paired = COUNTFASTQ.out.paired.ifEmpty(params.pairedFiles)
    unpaired = COUNTFASTQ.out.unpaired.ifEmpty(params.unpairedFiles)


    if(params.modules.faqcs) {
        FAQCS(params.faqcs.plus(params.shared),paired,unpaired,avgLen)
    }

}