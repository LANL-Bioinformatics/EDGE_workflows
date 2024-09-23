#!/usr/bin/env nextflow


//prepares input for phage finder (creates appropriate files)
process phageFinderPrep {

    input:
    path gff
    path fna

    output:
    path "id_map.txt", emit:idMap //separate output declaration for post-PF processing
    path "*", emit:allPFoutput //all output files will go into the next process


    script:
    """
    phageFinder_prepare.pl -o . $gff $fna
    """
} 

//calls phage_finder
process phageFinder {
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis/Prophage",
        mode: 'copy',
        pattern: "log.txt"
    )

    input:
    path prepOut
    path faa, stageAs: "Assembly.pep"

    output:
    path "PFPR_tab.txt", emit: phageTable
    path "log.txt"

    //must be on PATH
    script:
    """
    phage_finder_v2.1.sh Assembly $params.numCPU 1>log.txt 2>&1
    """

}


//creates text summary of results
process summary {
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis/Prophage",
        mode: 'copy'
    )

    input:
    path idMap
    path pfprTab

    output:
    path "phageFinder_summary.txt"

    script:
    """
    phageFinder_summary.pl -t $pfprTab -i $idMap
    """
}


workflow {
    gff_ch = channel.fromPath(params.gffFile, checkIfExists:true)
    faa_ch = channel.fromPath(params.faaFile, checkIfExists:true).filter{ it.size()>0 }
    fna_ch = channel.fromPath(params.fnaFile, checkIfExists:true)

    phageFinderPrep(gff_ch, fna_ch)
    phageFinder(phageFinderPrep.out.allPFoutput, faa_ch)
    summary(phageFinderPrep.out.idMap,phageFinder.out.phageTable)

}