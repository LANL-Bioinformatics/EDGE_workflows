#!/usr/bin/env nextflow

//sets taxonomic kingdom for analysis if none provided
process unsetKingdom {
    container 'apwat/annotation:0.4'
    containerOptions '--compat --bind .:/venv/bin/ec_info'
    input:
    path contigs

    output:
    stdout

    shell:
    //get total contig length; smaller fastas with no set kingdom are set to viral annotation if going through prokka 
    '''
    a=$(grep -v ">" !{contigs} | wc -m)
    b=$(grep -v ">" !{contigs} | wc -l)
    c=$((a-b))
    if [ c -gt 580000 ]; then echo "Bacteria"; else echo "Viruses"; fi
    '''
}

//process to invocate prokka
process prokkaAnnotate {
    container 'apwat/annotation:0.4'
    containerOptions '--compat --bind .:/venv/bin/ec_info'
    publishDir(
        path: "${settings["outDir"]}/AssemblyBasedAnalysis/Annotation",
        mode: 'copy'
    )

    input:
    path contigs
    path protein
    path hmm
    val kingdom
    val settings

    output:
    path "Annotation.log"
    path "${settings["projName"]}.gff", emit: gff
    path "${settings["projName"]}.err"
    path "${settings["projName"]}.faa", emit:faa
    path "${settings["projName"]}.ffn"
    path "${settings["projName"]}.fna", emit:fna
    path "${settings["projName"]}.fsa"
    path "${settings["projName"]}.gbk"
    path "${settings["projName"]}.log"
    path "${settings["projName"]}.sqn"
    path "${settings["projName"]}.tbl"
    path "${settings["projName"]}.tsv"
    path "${settings["projName"]}.txt"
    path "*.hmm.*", optional:true

    script:
    def kingdom = kingdom.trim()
    def protein = protein.name == "NO_FILE" ? "" : "--protein ${settings["customProtein"]}"
    def hmmPrep = hmm.name == "NO_FILE2" ? "" : "hmmpress $hmm"
    def hmm = hmm.name == "NO_FILE2" ? "" : "--hmms $hmm"
    def evalue = settings["evalue"] == null ? "" : "--evalue ${settings["evalue"]}"
    def gcode = settings["gcode"] == null ? "" : "--gcode ${settings["gcode"]}"
    def locustag = settings["projName"] == null ? "" : "--locustag ${settings["projName"]}"
    def prefix = settings["projName"] == null ? "" : "--prefix ${settings["projName"]}"
    def cpu = settings["cpus"] == null ? "" : "--cpus ${settings["cpus"]}"
    def taxKingdom = kingdom.equalsIgnoreCase("metagenome") ? "--kingdom Bacteria --metagenome" : "--kingdom $kingdom"

    """
    $hmmPrep

    prokka --quiet --force \
    $protein \
    $hmm \
    $evalue \
    $gcode \
    $locustag \
    $prefix \
    $cpu \
    --outdir . \
    $taxKingdom \
    $contigs 2>>Annotation.log 

    cat ${settings["projName"]}.log >> Annotation.log
    """
}

//process to invocate RATT
process rattAnnotate {
    container 'apwat/annotation:0.4'
    containerOptions '--compat --bind .:/venv/bin/ec_info'
    publishDir(
        path: "${settings["outDir"]}/AssemblyBasedAnalysis/Annotation",
        mode: 'copy'
    )

    input:
    path contigs
    path gbk
    val settings

    output:
    path "Annotation.log"
    path "${settings["projName"]}.gff", emit: gff
    path "${settings["projName"]}.faa", emit: faa
    path "${settings["projName"]}.fna", emit: fna
    path "${settings["projName"]}.gbk"

    shell:
    //happens in work directory
    //RATT version needs to be custom EDGE version
    '''
    mkdir -p ./RATT/source
    cp !{gbk} ./RATT/source/source.gbk
    genbank2embl.pl ./RATT/source/source.gbk
    cd RATT
    runRATT.sh $PWD/source ../!{contigs} !{settings["projName"]} Species 1>>../Annotation.log 2>&1
    cd ..
    cat ./RATT/*final.embl | fix_RATT_embl_feature.pl - > RATT/all.embl && embl2genbank.pl RATT/all.embl !{settings["projName"]}.gbk 
    genbank2fasta.pl -translation !{settings["projName"]}.gbk > !{settings["projName"]}.faa
    genbank2fasta.pl -genome !{settings["projName"]}.gbk > !{settings["projName"]}.fna
    genbank2gff3.pl -e 3 --outdir stdout --DEBUG --typesource contig $PWD/!{settings["projName"]}.gbk >!{settings["projName"]}.gff
    '''
}

//plots feature count, protein size distribution, etc.
process annPlot {
    container 'apwat/annotation:0.4'
    containerOptions '--compat --bind .:/venv/bin/ec_info'
    publishDir(
        path: "${settings["outDir"]}/AssemblyBasedAnalysis/Annotation",
        mode: 'copy'
    )
    
    input:
    path gff
    val settings

    output:
    path "plot_gff3.log"
    path "annotation_stats_plots.pdf"
    path "${settings["projName"]}.txt", optional:true

    script:
    def rattReport = settings["annotateProgram"].equalsIgnoreCase("ratt") ? "awk \'\$1 ~ /CDS|RNA/ {print \$1\": \"\$2}' plot_gff3.log > ${settings["projName"]}.txt" : ""
    """
    plot_gff3_stats.pl --input $gff --title ${settings["projName"]} --prefix ./annotation_stats --outfmt PDF 1>plot_gff3.log 2>&1
    $rattReport
    """
}

//generates KEGG pathway plots
process keggPlots {
    container 'apwat/annotation:0.4'
    containerOptions '--compat --bind .:/venv/bin/ec_info'
    publishDir(
        path: "${settings["outDir"]}/AssemblyBasedAnalysis/Annotation",
        mode: 'copy'
    )

    input:
    path gff
    val settings
    
    output:
    path "kegg_map/*"
    path "kegg_map.log"
    
    script:
    //TODO: will eventually need UI integration.
    """
    check_server_up.pl --url "http://rest.kegg.jp" && \
    opaver_anno.pl -g $gff -o ./kegg_map -p ${settings["projName"]} > kegg_map.log 2>&1
    """
    }


workflow ANNOTATION {
    take:
    settings
    data

    main:
    //optional file pattern setup
    "mkdir nf_assets".execute().text
    "touch nf_assets/NO_FILE".execute().text
    "touch nf_assets/NO_FILE2".execute().text

    //inputs
    kingdom_ch = channel.of(settings["taxKingdom"])
    hmm_ch = channel.fromPath(settings["customHMM"], checkIfExists:true)
    prot_ch = channel.fromPath(settings["customProtein"], checkIfExists:true)

    //output setup
    gff = channel.empty()
    faa = channel.empty()
    fna = channel.empty()

    //workflow logic
    if (settings["annotateProgram"].equalsIgnoreCase("prokka")) {
        if (settings["taxKingdom"] == null) {
        kingdom_ch = unsetKingdom(data)
        }
        prokkaAnnotate(data, prot_ch, hmm_ch, kingdom_ch, settings)
        annPlot(prokkaAnnotate.out.gff, settings)
        if(settings["keggView"] == true) {
            keggPlots(prokkaAnnotate.out.gff, settings)
        }
        gff = prokkaAnnotate.out.gff
        faa = prokkaAnnotate.out.faa
        fna = prokkaAnnotate.out.fna
    }
    else if (settings["annotateProgram"].equalsIgnoreCase("ratt")) {
        gb_ch = channel.fromPath(settings["sourceGBK"], checkIfExists:true)
        rattAnnotate(data, gb_ch, settings)
        annPlot(rattAnnotate.out.gff, settings)
        if(settings["keggView"] == true) {
            keggPlots(rattAnnotate.out.gff, settings)
        }
        gff = rattAnnotate.out.gff
        faa = rattAnnotate.out.faa
        fna = rattAnnotate.out.fna
    }

    emit:
    gff
    faa
    fna
}