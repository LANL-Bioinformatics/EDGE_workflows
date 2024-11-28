#!/usr/bin/env nextflow


process phageFinderPrep {
    container "apwat/phage_finder:noWrite"

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

process phageFinder {
    container "apwat/phage_finder:noWrite"
    publishDir(
        path: "${settings["outDir"]}/AssemblyBasedAnalysis/Prophage",
        mode: 'copy',
        pattern: "log.txt"
    )

    input:
    path prepOut
    path faa, stageAs: "Assembly.pep"
    val settings

    output:
    path "PFPR_tab.txt", emit: phageTable
    path "log.txt"

    //must be on PATH
    script:
    """
    phage_finder_v2.1.sh Assembly ${settings["cpus"]} 1>log.txt 2>&1
    """

}

process summary {
    container "apwat/phage_finder:noWrite"
    publishDir(
        path: "${settings["outDir"]}/AssemblyBasedAnalysis/Prophage",
        mode: 'copy'
    )

    input:
    path idMap
    path pfprTab
    val settings

    output:
    path "phageFinder_summary.txt"

    script:
    """
    phageFinder_summary.pl -t $pfprTab -i $idMap
    """
}


workflow PHAGE_FINDER {
    take:
    settings
    gff
    faa
    fna
    
    main:

    phageFinderPrep(gff, fna)
    phageFinder(phageFinderPrep.out.allPFoutput, faa, settings)
    summary(phageFinderPrep.out.idMap,phageFinder.out.phageTable, settings)

}