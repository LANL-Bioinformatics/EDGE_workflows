#!/usr/bin/env nextflow

//to run: nextflow run runAssembly.nf -params-file [JSON parameter file]

params.unpaired = "NO_FILE"
params.paired = "NO_FILE"
params.longReads = "NO_FILE"
params.pacBio = null
params.ref = "NO_FILE"
params.pre = "ReadsMapping"
params.outDir = null
params.consensus = null
params.aligner = null
params.bwaOptions = null
params.bowtieOptions = null
params.snapOptions = null
params.minimap2Options = null
params.maq = null
params.maxClip = null
params.alignTrimBedFile = null
params.alignTrimStrand = null
params.skipAln = null
params.noPlot = null
params.nosnp = null
params.noIndels = null
params.debugFlag = null
params.cpu = null
params.ploidy = null
params.disableBAQ = null
params.maxDepth = null
params.minIndelCandidateDepth = null
params.minAltBases = null
params.minAltRatio = null
params.minDepth = null
params.snpGapFilter = null


def paramConstructor(flag, value) {
    return "$flag $value"
}

process runAssembly {
    
    debug true
    publishDir ".", mode: 'copy'

    input:
    
    path ('singleton', arity:1)
    path ('paired', arity:2)
    path ('longReads', arity:1)
    path (ref, arity:1)

    output:
    path "$params.outDir/*" //catchall for produced files.

    script:

    def unpaired = params.unpaired != "NO_FILE" ? "-u singleton" : ""
    def paired = params.paired != "NO_FILE" ? "-p \"paired1 paired2\"" : ""
    def longReads = params.longReads != "NO_FILE" ? "-long longReads" : ""
    def pacBio = params.pacBio != null ? "--pacbio $params.pacBio" : ""
    def ref = params.ref != "NO_FILE" ? "-ref $ref" : ""
    def pre = params.pre != "ReadsMapping" ? "-pre $params.pre" : ""
    def outDir = params.outDir != null ? "-d $params.outDir" : ""
    def consensus = params.consensus != null ? "-consensus $params.consensus" : ""
    def aligner = params.aligner != null ? "-aligner $params.aligner" : ""
    def bwaOptions = params.bwaOptions != null ? "-bwa_options $params.bwaOptions" : ""
    def bowtieOptions = params.bowtieOptions != null ? "-bowtie_options $params.bowtieOptions" : ""
    def snapOptions = params.snapOptions != null ? "-snap_options $params.snapOptions" : ""
    def minimap2Options = params.minimap2Options != null ? "-minimap2_options $params.minimap2Options" : ""
    def maq = params.maq != null ? "-maq $params.maq" : ""
    def maxClip = params.maxClip != null ? "-max_clip $params.maxClip" : ""
    def alignTrimBedFile = params.alignTrimBedFile != null ? "-align_trim_bed_file $params.alignTrimBedFile" : ""
    def alignTrimStrand = params.alignTrimStrand != null ? "-align_trim_strand $params.alignTrimStrand" : ""
    def skipAln = params.skipAln != null ? "-skip_aln $params.skipAln" : ""
    def noPlot = params.noPlot != null ? "-no_plot $params.noPlot" : ""
    def nosnp = params.nosnp != null ? "-no_snp $params.nosnp" : ""
    def noIndels = params.noIndels != null ? "-no_indels $params.noIndels" : ""
    def debugFlag = params.debugFlag != null ? "-debug $params.debugFlag" : ""
    def cpu = params.cpu != null ? "-cpu $params.cpu" : ""
    def ploidy = params.ploidy != null ? "-ploidy $params.ploidy" : ""
    def disableBAQ = params.disableBAQ != null ? "-disableBAQ $params.disableBAQ" : ""
    def maxDepth = params.maxDepth != null ? "-max_depth $params.maxDepth" : ""
    def minIndelCandidateDepth = params.minIndelCandidateDepth != null ? "-min_indel_candidate_depth $params.minIndelCandidateDepth" : ""
    def minAltBases = params.minAltBases != null ? "-min_alt_bases $params.minAltBases" : ""
    def minAltRatio = params.minAltRatio != null ? "-min_alt_ratio $params.minAltRatio" : ""
    def minDepth = params.minDepth != null ? "-min_depth $params.minDepth" : ""
    def snpGapFilter = params.snpGapFilter != null ? "-snp_gap_filter $params.snpGapFilter" : ""

    """
    runReadsToGenome.pl \
    $paired\
    $unpaired\
    $longReads\
    $pacBio\
    $ref\
    $pre\
    $outDir\
    $consensus\
    $aligner\
    $bwaOptions\
    $bowtieOptions\
    $snapOptions\
    $minimap2Options\
    $maq\
    $maxClip\
    $alignTrimBedFile\
    $alignTrimStrand\
    $skipAln\
    $noPlot\
    $nosnp\
    $noIndels\
    $debugFlag\
    $cpu\
    $ploidy\
    $disableBAQ\
    $maxDepth\
    $minIndelCandidateDepth\
    $minAltBases\
    $minAltRatio\
    $minDepth\
    $snpGapFilter
    """

}

workflow {
    "touch NO_FILE".execute().text //later removed by workflow completion handler in nextflow.config
    singleton_ch = Channel.fromPath(params.unpaired, checkIfExists: true, relative:true)
    paired_ch = Channel.fromPath(params.paired, checkIfExists: true, relative:true).collect()
    longReads_ch = Channel.fromPath(params.longReads, checkIfExists: true, relative:true)
    providedRef = file(params.ref, checkIfExists: true)
    runAssembly(singleton_ch, paired_ch, longReads_ch, providedRef)

}