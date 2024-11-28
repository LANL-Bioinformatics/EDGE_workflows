
//runs ANTISMASH on provided .fa or .gbk input.
process antismash {
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis/AntiSmash",
        mode: 'copy'
    )

    input:
    path input

    output:
    path "output/*"
    path "antismashLog.txt"

    script:
    def taxon = params.taxon.equalsIgnoreCase("fungi") ? "--taxon fungi" : "--taxon bacteria"
    def clusterblast = params.clusterblast == true ? "--cb-general" : ""
    def subclusterblast = params.subclusterblast == true ? "--cb-subclusters" : ""
    def knownclusterblast = params.knownclusterblast == true ? "--cb-knownclusters" : ""
    def mibig = params.mibig == true ? "--cc-mibig" : ""
    def smcogs = params.smcogs == true ? "--smcog-trees" : ""
    def asf = params.asf == true ? "--asf" : ""
    def rre = params.rre == true ? "--rre" : ""
    def fullhmmer = params.fullhmm == true ? "--fullhmmer" : ""
    def tigrfam = params.tigrfam == true ? "--tigrfam" : ""
    def pfam2go = params.pfam2go == true ? "--pfam2go" : ""
    def genefinding = params.taxon.equalsIgnoreCase("fungi") ? "--genefinding-tool glimmerhmm" : "--genefinding-tool prodigal-m"
    def cassis = (params.taxon.equalsIgnoreCase("fungi") && params.cassis == true) ? "--cassis" : ""

    """
    antismash -c $params.numCPU $taxon \
    --logfile antismashLog.txt --output-dir ./output \
    --html-title $params.projName --database $params.database \
    $clusterblast \
    $subclusterblast \
    $knownclusterblast \
    $mibig \
    $smcogs \
    $asf \
    $rre \
    $fullhmmer \
    $tigrfam \
    $pfam2go \
    $cassis \
    $genefinding \
    $input
    """

}


workflow {
    antismash(channel.fromPath(params.input,checkIfExists:true))
}