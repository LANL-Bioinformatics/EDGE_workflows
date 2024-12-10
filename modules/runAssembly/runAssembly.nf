#!/usr/bin/env nextflow
//to run: nextflow run runAssembly.nf -params-file [JSON parameter file]

//this workflow is unable to set memory limits (used in idba and spades assemblies) by itself, 
//but a limit (in KB) can be provided as a parameter.

//main process for assembly with IDBA
process idbaUD {
    label "assembly"
    publishDir (
    path:"${settings["outDir"]}/AssemblyBasedAnalysis",
    mode: 'copy',
    saveAs: {
        filename ->
        if(filename ==~ /log/) {
            "assembly.log"
        }
        else if(filename ==~ /scaffold(s|-level-2|)\.fa(sta)?/) {
            "scaffold.fa"
        }
        else{
            null //do not publish contig-max.fa or intermediate contigs yet, but use them for downstream processes
        }
    }
    )

    input:
    val settings
    path short_paired
    path short_single
    path long_reads
    val avg_len

    output:
    path "log"
    path "scaffold*.{fa,fasta}", optional:true 
    path "contig-max.fa", emit: contigs, optional:true
    path "{contig-*,*contigs.fa,K*/final_contigs.fasta}", emit: intContigs

    script:
    def avg_len = avg_len as Integer
    def runFlag = ""
    if(short_paired.name != "NO_FILE" && short_single.name != "NO_FILE2") {
        runFlag = "-r $short_single --read_level_2 $short_paired "
    }
    else if(short_paired.name != "NO_FILE") {
        runFlag = "-r $short_paired "
    }
    else if(short_single.name != "NO_FILE2") {
        runFlag = "-r $short_single "
    }
    def longReadsFile = long_reads.name != "NO_FILE3" ? "-l $long_reads" : ""
    def maxK = 121
    def maxK_option = "--maxk $maxK "
    if (settings["idba"]["maxK"] == null || settings["idba"]["maxK"] > avg_len) {
        if(avg_len > 0 && avg_len <= 151) {
            maxK = avg_len - 1
            maxK_option = "--maxk ${avg_len - 1}"
        }
    }
    minK = settings["idba"]["minK"] != null ? "--mink ${settings["idba"]["minK"]} " : ""
    step = settings["idba"]["step"] != null ? "--step ${settings["idba"]["step"]} " : ""
    minLen = settings["minContigSize"] != null ? "--min_contig ${settings["minContigSize"]} " : ""

    memLimit = settings["memLimit"] != null ? "ulimit -v ${settings["memLimit"]} 2>/dev/null;" : ""
    """
    ${memLimit}idba_ud --pre_correction -o . --num_threads ${settings["cpus"]}\
    $runFlag\
    $longReadsFile\
    $maxK_option\
    $minK\
    $step\
    $minLen

    mv contig-${maxK}.fa contig-max.fa
    """

}

//prep for idba
process idbaExtractLong {
    label "assembly"

    input:
    path paired
    path unpaired

    output:
    path "short_paired.fa"
    path "short_single.fa"
    path "long.fa"

    script:
    def pair_file = paired.name != "NO_FILE" ? "-p $paired " : ""
    def unpaired_file = unpaired.name != "NO_FILE2" ? "-u $unpaired " : ""
    """
    extractLongReads.pl\
    $pair_file\
    $unpaired_file\
    -d .
    """
}

//prep for idba
process idbaPrepReads {
    label "assembly"
    input:
    path paired
    path unpaired

    output:
    path "pairedForAssembly.fasta", emit:idba_prep_paired
    path "unpairedForAssembly.fasta", optional:true

    script:
    def pair_process = paired.name != "NO_FILE" ? "fq2fa --filter --merge ${paired[0]} ${paired[1]} pairedForAssembly.fasta;" : "" 
    def unpair_process = unpaired.name != "NO_FILE2" ? "fq2fa --filter $unpaired unpairedForAssembly.fasta;" : "" 
    
    """
    $pair_process
    $unpair_process
    """

}


//assemble using spades
process spades {
    label "assembly"

    publishDir (
    path: "${settings["outDir"]}/AssemblyBasedAnalysis", 
    mode: 'copy',
    saveAs: {
        filename ->
        if(filename ==~ /spades\.log/) {
            "assembly.log"
        }
        else if(filename ==~ /scaffold(s|-level-2|)\.fa(sta)?/) {
            "scaffold.fa"
        }
        else if(filename ==~ /assembly_graph\.fastg/) {
            "${settings["projName"]}_contigs.fastg"
        }
        else if(filename ==~ /assembly_graph_with_scaffolds\.gfa/) {
            "assembly_graph_with_scaffolds.gfa"
        }
        else if(filename ==~ /(contigs|transcripts)\.fasta/) {
            null //do not publish contigs, but use downstream
        }
        else if(filename ==~ /((contig-.*)|(.*contigs\.fa)|(K.*\/final_contigs\.fasta))/) {
            null //do not publish intermediate contigs, but use downstream
        }
        else{
            filename
        }
    }
    )

    input:
    val settings
    path paired
    path unpaired
    path pacbio
    path nanopore

    output:
    path "scaffold*.{fa,fasta}", optional:true 
    path "spades.log" 
    path "contigs.paths"
    path "{contigs,transcripts}.fasta", emit:contigs, optional:true
    path "assembly_graph.fastg", optional:true
    path "assembly_graph_with_scaffolds.gfa", optional:true
    path "{contig-*,*contigs.fa,K*/final_contigs.fasta}", emit: intContigs


    script:
    def paired = paired.name != "NO_FILE" ? "--pe1-1 ${paired[0]} --pe1-2 ${paired[1]} " : ""
    def unpaired = unpaired.name != "NO_FILE2" ? "--s1 $unpaired " : ""
    def pacbio_file = pacbio.name != "NO_FILE3" ? "--pacbio $pacbio " : ""
    def nanopore_file = nanopore.name != "NO_FILE4" ? "--nanopore $nanopore " : ""
    def meta_flag = (paired != "" && settings["spades"]["algorithm"] == "metagenome") ? "--meta " : ""
    def sc_flag = settings["spades"]["algorithm"] == "singlecell" ? "--sc " : ""
    def rna_flag = settings["spades"]["algorithm"] == "rna" ? "--rna " : ""
    def plasmid_flag = settings["spades"]["algorithm"] == "plasmid" ? "--plasmid " : ""
    def bio_flag = settings["spades"]["algorithm"] == "bio" ? "--bio " : ""
    def corona_flag = settings["spades"]["algorithm"] == "corona" ? "--corona " : ""
    def metaviral_flag = settings["spades"]["algorithm"] == "metaviral" ? "--metaviral " : ""
    def metaplasmid_flag = settings["spades"]["algorithm"] == "metaplasmid" ? "--metaplasmid " : ""
    def rnaviral_flag = settings["spades"]["algorithm"] == "rnaviral" ? "--rnaviral " : ""
    def memLimit = settings["memLimit"] != null ? "-m ${settings["memLimit"]}" : ""

    """
    spades.py -o . -t ${settings["cpus"]}\
    $paired\
    $meta_flag\
    $sc_flag\
    $rna_flag\
    $plasmid_flag\
    $bio_flag\
    $corona_flag\
    $metaviral_flag\
    $metaplasmid_flag\
    $rnaviral_flag\
    $unpaired\
    $pacbio_file\
    $nanopore_file\
    $memLimit
    """
}

//assemble using megahit
process megahit {

    label "assembly"
    publishDir(
    path: "${settings["outDir"]}/AssemblyBasedAnalysis", 
    mode: 'copy',
    saveAs: {
        filename ->
        if(filename.equals("megahit/log")) {
            "assembly.log"
        }
        else if(filename ==~ /megahit\/scaffold(.*)\.fa(sta)?/) {
            "scaffold.fa"
        }
        else if(filename ==~ /megahit\/final\.contigs\.fa/) {
            null //don't publish, but pass to downstream "renameFilterFasta" process
        }
        else if(filename ==~ /megahit\/((contig-.*)|(.*contigs\.fa)|(K.*\/final_contigs\.fasta))/) {
            null //do not publish intermediate contigs, but use downstream
        }
        else{
            filename
        }
    }
    )

    input:
    val settings
    path paired
    path unpaired

    output:
    path "megahit/log"
    path "megahit/scaffold*.{fa,fasta}", optional:true //I don't believe this is a normal output of megahit, but just in case
    path "megahit/final.contigs.fa", emit: contigs, optional:true
    path "${settings["projName"]}_contigs.fastg"
    path "megahit/{contig-*,*contigs.fa,K*/final_contigs.fasta}", emit: intContigs

    script:
    def paired = paired.name != "NO_FILE" ? "-1 ${paired[0]} -2 ${paired[1]} " : ""
    def unpaired = unpaired.name != "NO_FILE2" ? "-r $unpaired " : ""
    def megahit_preset = settings["megahit"]["preset"] != null ? "--presets ${settings["megahit"]["preset"]} " : ""

    """
    megahit -o ./megahit -t ${settings["cpus"]}\
    $megahit_preset\
    $paired\
    $unpaired\
    2>&1

    LARGESTKMER=\$(head -n 1 megahit/final.contigs.fa | perl -ne '/^>k(\\d+)\\_/; print \$1;')

    megahit_toolkit contig2fastg \$LARGESTKMER megahit/intermediate_contigs/k\${LARGESTKMER}.contigs.fa  > ${settings["projName"]}_contigs.fastg
    """

}


//assembly using unicycler
process unicycler {
    label "assembly"
    publishDir (
        path: "${settings["outDir"]}/AssemblyBasedAnalysis", 
        mode: 'copy',
        saveAs: {
        filename ->
        if(filename ==~ /unicycler\.log/) {
            "assembly.log"
        }
        else if(filename ==~ /scaffold(s|-level-2|)\.fa(sta)?/) {
            "scaffold.fa"
        }
        else if(filename.equals("assembly.gfa")) {
            "${settings["projName"]}_contigs.fastg"
        }
        else if(filename ==~ /assembly\.fasta/) {
            null //don't publish, but emit for use in downstream process "renameFilterFasta"
        }
        else{
            filename
        }
    }
    )

    input:
    val settings
    path paired
    path unpaired
    path longreads //If present, expects filtered long reads.

    output:
    path "unicycler.log"
    path "scaffold*.{fa,fasta}", optional:true //I don't believe this is a normal output of unicycler, but just in case
    path "assembly.fasta", emit: contigs
    path "assembly.gfa", optional:true
    //path "{contig-*,*contigs.fa,K*/final_contigs.fasta}", emit: intContigs | not produced by unicycler

    script:
    def paired = paired.name != "NO_FILE" ? "-1 ${paired[0]} -2 ${paired[1]} " : ""
    def unpaired = unpaired.name != "NO_FILE2" ? "-r $unpaired " : ""
    def filt_lr = longreads.name != "NO_FILE3" ? "-l $longreads " : ""
    def bridge = settings["unicycler"]["bridgingMode"] != "normal" ? "--mode ${settings["unicycler"]["bridgingMode"]} " : "--mode normal"

    """
    export _JAVA_OPTIONS='-Xmx20G'; export TERM='xterm';

    unicycler -t ${settings["cpus"]} -o .\
    $paired\
    $filt_lr\
    $bridge 2>&1 1>/dev/null
    """

}

//filter long reads for unicycler
process unicyclerPrep {
    label "assembly"

    input:
    val settings
    path longreads


    output:
    path "long_reads.fasta"

    script:

    """
    seqtk seq -A -L \
    ${settings["unicycler"]["minLongReads"]} \
    $longreads > long_reads.fasta
    """
}


//assembly using lrasm
process lrasm {
    label "assembly"

    publishDir (
        path: "${settings["outDir"]}/AssemblyBasedAnalysis", 
        mode: 'copy',
        saveAs: {
        filename ->
        if(filename ==~ /contigs\.log/) {
            "assembly.log"
        }
        else if(filename ==~ /scaffold(s|-level-2|)\.fa(sta)?/) {
            "scaffold.fa"
        }
        else if(filename.equals("Assembly/unitig.gfa")) {
            "${settings["projName"]}_contigs.fastg"
        }
        else if(filename.equals("Assembly/assembly_graph.gfa")) {
            "${settings["projName"]}_contigs.fastg"
        }
        else if(filename.equals("Assembly/assembly_graph.gv")) {
            "${settings["projName"]}_contigs.gv"
        }
        else if(filename.equals("Assembly/assembly_info.txt")) {
            "assembly_info.txt"
        }
        else if(filename ==~ /contigs\.fa/) {
            null //do not publish, but emit for use in downstream process "renameFilterFasta"
        }
        else if(filename ==~ /((contig-.*)|(.*contigs\.fa)|(K.*\/final_contigs\.fasta))/) {
            null //do not publish intermediate contigs, but use downstream
        }
        else{
            filename
        }
    }
    )

    input:
    val settings
    path unpaired

    output:
    path "contigs.log" 
    path "scaffold*.{fa,fasta}", optional:true //I don't believe this is a normal output of lrasm, but just in case
    path "contigs.fa", emit:contigs, optional:true
    path "Assembly/unitig.gfa", optional:true
    path "Assembly/assembly_graph.gfa", optional:true
    path "Assembly/assembly_graph.gv", optional:true
    path "Assembly/assembly_info.txt", optional:true
    path "{contig-*,*contigs.fa,K*/final_contigs.fasta}", emit: intContigs

    script:
    def consensus = settings["lrasm"]["numConsensus"] != null ? "-n ${settings["lrasm"]["numConsensus"]} ": ""
    def preset = settings["lrasm"]["preset"] != null ? "-x ${settings["lrasm"]["preset"]} " : ""
    def errorCorrection = settings["lrasm"]["ec"] != null ? "-e " : ""
    def algorithm = settings["lrasm"]["algorithm"] != null ? "-a ${settings["lrasm"]["algorithm"]} " : ""
    def minLenOpt = ""
    if (settings["lrasm"]["algorithm"] == "miniasm") {
        minLenOpt = "--ao \'-s ${settings["lrasm"]["minLength"]}\' "
    }
    else if (settings["lrasm"]["algorithm"] == "wtdbg2") {
        minLenOpt = "--wo \'-L ${settings["lrasm"]["minLength"]}\' "
    }
    def flyeOpt = settings["lrasm"]["algorithm"] == "metaflye" ? "--fo '--meta' ": ""

    """
    lrasm -o . -t ${settings["cpus"]} \
    $preset\
    $consensus\
    $errorCorrection\
    $algorithm\
    $minLenOpt\
    $flyeOpt\
    $unpaired\
    """
}

process renameFilterFasta {
    label "assembly"
    publishDir(
        path: "${settings["outDir"]}/AssemblyBasedAnalysis",
        mode: 'copy'
    )
    input:
    val settings
    path contigs

    output:
    path "*_contigs.fa", emit: contigs
    path "*_contigs_*up.fa", emit: annotationContigs

    script:
    def annotation = settings["annotation"] ? "-ann 1" : ""
    """
    CONTIG_NUMBER=\$(grep -c '>' ${contigs})
    
    renameFilterFasta.pl \
    -u $contigs\
    -d .\
    -filt ${settings["minContigSize"]} \
    -maxseq \$CONTIG_NUMBER\
    -ann_size ${settings["contigSizeForAnnotation"]} \
    -n ${settings["projName"]} \
    $annotation
    """

}

process bestIncompleteAssembly {
    label "assembly"
    input:

    val x
    path intContigs

    when:
    x == 'EMPTY'

    output:
    path "bestIntContig/*"

    shell:
    '''
    #!/usr/bin/env perl
    use Cwd;
    my $dir = getcwd;
    use Cwd 'abs_path';
    my @intermediate_contigs = sort { -M $a <=> -M $b} glob("!{intContigs}");
    my $best_int = $dir . "/" . $intermediate_contigs[0];
    mkdir bestIntContig;
    system("cp $best_int ./bestIntContig");
    '''

}

workflow ASSEMBLY {
    take:
    settings
    paired
    unpaired
    avgLen

    main:

    //supplementary long read files setup
    spades_pb = file(settings["spades"]["pacbio"], checkIfExists:true)
    spades_np = file(settings["spades"]["nanopore"], checkIfExists:true)
    unicycler_lr = file(settings["unicycler"]["longreads"], checkIfExists:true)


    //output channel setup
    outContigs = channel.empty()
    annotationContigs = channel.empty()

    if (settings["assembler"].equalsIgnoreCase("IDBA_UD")) {

        (c1,c2) = idbaPrepReads(paired, unpaired)
        (sp,su,l) = idbaExtractLong(c1,c2.ifEmpty({file("nf_assets/NO_FILE")}))

        idbaUD(settings, sp.filter{ it.size()>0 }.ifEmpty({file("nf_assets/NO_FILE")}),
            su.filter{ it.size()>0 }.ifEmpty({file("nf_assets/NO_FILE2")}),
            l.filter{ it.size()>0 }.ifEmpty({file("nf_assets/NO_FILE3")}),
            avgLen)
        
        bestIncompleteAssembly(idbaUD.out.contigs.ifEmpty('EMPTY'), idbaUD.out.intContigs)
        renameFilterFasta(settings, idbaUD.out.contigs.concat(bestIncompleteAssembly.out).first())
        outContigs = renameFilterFasta.out.contigs
        annotationContigs = renameFilterFasta.out.annotationContigs

    }
    else if (settings["assembler"].equalsIgnoreCase("SPAdes")) {
        spades(settings, paired, unpaired, spades_pb, spades_np)
        
        bestIncompleteAssembly(spades.out.contigs.ifEmpty('EMPTY'), spades.out.intContigs)
        renameFilterFasta(settings, spades.out.contigs.concat(bestIncompleteAssembly.out).first())
        outContigs = renameFilterFasta.out.contigs
        annotationContigs = renameFilterFasta.out.annotationContigs
    }
    else if (settings["assembler"].equalsIgnoreCase("MEGAHIT")) {
        megahit(settings, paired, unpaired)
        
        bestIncompleteAssembly(megahit.out.contigs.ifEmpty('EMPTY'), megahit.out.intContigs)
        renameFilterFasta(settings, megahit.out.contigs.concat(bestIncompleteAssembly.out).first())
        outContigs = renameFilterFasta.out.contigs
        annotationContigs = renameFilterFasta.out.annotationContigs
    }
    else if (settings["assembler"].equalsIgnoreCase("UniCycler")) {
        if (settings["unicycler"]["longreads"] != "nf_assets/NO_FILE3") {
            println("Filter long reads with ${settings["unicycler"]["minLongReads"]} (bp) cutoff")
            unicycler(
                settings,
                paired,
                unpaired,
                unicyclerPrep(settings,unicycler_lr).filter{it.size()>0}.ifEmpty({file("nf_assets/NO_FILE3")})
                )
            //unicycler produces no intermediate contigs, we let it error out above rather than try to rescue a failed assembly
            renameFilterFasta(settings, unicycler.out.contigs)
            outContigs = renameFilterFasta.out.contigs
            annotationContigs = renameFilterFasta.out.annotationContigs
        }
        else {
            unicycler(settings, paired, unpaired, unicycler_lr)
            renameFilterFasta(settings, unicycler.out.contigs)
            outContigs = renameFilterFasta.out.contigs
            annotationContigs = renameFilterFasta.out.annotationContigs
        }
    }
    else if (settings["assembler"].equalsIgnoreCase("LRASM")) {
        lrasm(settings, unpaired)

        bestIncompleteAssembly(lrasm.out.contigs.ifEmpty('EMPTY'), lrasm.out.intContigs)
        renameFilterFasta(settings, lrasm.out.contigs.concat(bestIncompleteAssembly.out).first())
        outContigs = renameFilterFasta.out.contigs
        annotationContigs = renameFilterFasta.out.annotationContigs
    }
    else {
        error "Invalid assembler: ${settings["assembler"]}"
    }

    emit:
    outContigs
    annotationContigs

}