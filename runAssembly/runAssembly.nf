#!/usr/bin/env nextflow
//to run: nextflow run runAssembly.nf -params-file [JSON parameter file]

//this workflow is unable to set memory limits (used in idba and spades assemblies) by itself, 
//but a limit (in KB) can be provided as a parameter.

//main process for assembly with IDBA
process idbaUD {
    publishDir (
    path:"$params.outDir/AssemblyBasedAnalysis",
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
    if (params.idba.maxK == null || params.idba.maxK > avg_len) {
        if(avg_len > 0 && avg_len <= 151) {
            maxK = avg_len - 1
            maxK_option = "--maxk ${avg_len - 1}"
        }
    }
    def minK = params.idba.minK != null ? "--mink $params.idba.minK " : ""
    def step = params.idba.step != null ? "--step $params.idba.step " : ""
    def minLen = params.minContigSize != null ? "--min_contig $params.minContigSize " : ""

    def memLimit = params.memLimit != null ? "ulimit -v $params.memLimit 2>/dev/null;" : ""
    """
    ${memLimit}idba_ud --pre_correction -o . --num_threads $params.threads\
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

//prep for idba
process idbaReadFastq {
    input:
    path paired
    path unpaired

    output:
    path "fastqCount.txt"

    script:
    def paired_list = paired.name != "NO_FILE" ? "-p ${paired}" : ""
    def unpaired_list = unpaired.name != "NO_FILE2" ? "-u ${unpaired}" : ""

    """
    getAvgLen.pl\
    $paired_list\
    $unpaired_list\
    -d .
    """
}

//prep for idba
process idbaAvgLen {
    input:

    path countFastq

    output:
    stdout

    shell:
    '''
    #!/usr/bin/env perl
    my $fastq_count_file = "./!{countFastq}";
    my $total_count = 0;
    my $total_len = 0;
    open (my $fh, "<", $fastq_count_file) or die "Cannot open $fastq_count_file\n";
    while(<$fh>){
        chomp;
        my ($name,$count,$len,$avg) = split /\t/,$_;
        $total_count += $count;
        $total_len += $len;
    }
    close $fh;
    my $avg_len = ($total_count > 0)? $total_len/$total_count : 0;
    print "$avg_len";
    '''
}

//assemble using spades
process spades {

    publishDir (
    path: "$params.outDir/AssemblyBasedAnalysis", 
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
            "${params.projName}_contigs.fastg"
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
    def meta_flag = (paired != "" && params.spades.algorithm == "metagenome") ? "--meta " : ""
    def sc_flag = params.spades.algorithm == "singlecell" ? "--sc " : ""
    def rna_flag = params.spades.algorithm == "rna" ? "--rna " : ""
    def plasmid_flag = params.spades.algorithm == "plasmid" ? "--plasmid " : ""
    def bio_flag = params.spades.algorithm == "bio" ? "--bio " : ""
    def corona_flag = params.spades.algorithm == "corona" ? "--corona " : ""
    def metaviral_flag = params.spades.algorithm == "metaviral" ? "--metaviral " : ""
    def metaplasmid_flag = params.spades.algorithm == "metaplasmid" ? "--metaplasmid " : ""
    def rnaviral_flag = params.spades.algorithm == "rnaviral" ? "--rnaviral " : ""
    def memLimit = params.memLimit != null ? "-m ${params.memLimit/1024*1024}" : ""

    """
    spades.py -o . -t $params.threads\
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

    publishDir(
    path: "$params.outDir/AssemblyBasedAnalysis", 
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
    path paired
    path unpaired

    output:
    path "megahit/log"
    path "megahit/scaffold*.{fa,fasta}", optional:true //I don't believe this is a normal output of megahit, but just in case
    path "megahit/final.contigs.fa", emit: contigs, optional:true
    path "${params.projName}_contigs.fastg"
    path "megahit/{contig-*,*contigs.fa,K*/final_contigs.fasta}", emit: intContigs

    script:
    def paired = paired.name != "NO_FILE" ? "-1 ${paired[0]} -2 ${paired[1]} " : ""
    def unpaired = unpaired.name != "NO_FILE2" ? "-r $unpaired " : ""
    def megahit_preset = params.megahit.preset != null ? "--presets $params.megahit.preset " : ""

    """
    megahit -o ./megahit -t $params.threads\
    $megahit_preset\
    $paired\
    $unpaired\
    2>&1

    LARGESTKMER=\$(head -n 1 megahit/final.contigs.fa | perl -ne '/^>k(\\d+)\\_/; print \$1;')

    megahit_toolkit contig2fastg \$LARGESTKMER megahit/intermediate_contigs/k\${LARGESTKMER}.contigs.fa  > ${params.projName}_contigs.fastg
    """

}


//assembly using unicycler
process unicycler {
    publishDir (
        path: "$params.outDir/AssemblyBasedAnalysis", 
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
            "${params.projName}_contigs.fastg"
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
    def bridge = params.unicycler.bridgingMode != "normal" ? "--mode $params.unicycler.bridgingMode" : "--mode normal"

    """
    export _JAVA_OPTIONS='-Xmx20G'; export TERM='xterm';

    unicycler -t $params.threads -o .\
    $paired\
    $filt_lr\
    $bridge 2>&1 1>/dev/null
    """

}

//filter long reads for unicycler
process unicyclerPrep {
    input:
    path longreads


    output:
    path "long_reads.fasta"

    script:

    """
    seqtk seq -A -L\
    $params.unicycler.minLongReads\
    $longreads > long_reads.fasta
    """
}


//assembly using lrasm
process lrasm {

    publishDir (
        path: "$params.outDir/AssemblyBasedAnalysis", 
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
            "${params.projName}_contigs.fastg"
        }
        else if(filename.equals("Assembly/assembly_graph.gfa")) {
            "${params.projName}_contigs.fastg"
        }
        else if(filename.equals("Assembly/assembly_graph.gv")) {
            "${params.projName}_contigs.gv"
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
    def consensus = params.lrasm.numConsensus != null ? "-n $params.lrasm.numConsensus ": ""
    def preset = params.lrasm.preset != null ? "-x $params.lrasm.preset " : ""
    def errorCorrection = params.lrasm.ec != null ? "-e " : ""
    def algorithm = params.lrasm.algorithm != null ? "-a $params.lrasm.algorithm " : ""
    def minLenOpt = ""
    if (params.lrasm.algorithm == "miniasm") {
        minLenOpt = "--ao \'-s $params.lrasm.minLength\' "
    }
    else if (params.lrasm.algorithm == "wtdbg2") {
        minLenOpt = "--wo \'-L $params.lrasm.minLength\' "
    }
    def flyeOpt = params.lrasm.algorithm == "metaflye" ? "--fo '--meta' ": ""

    """
    lrasm -o . -t $params.threads\
    $preset\
    $consensus\
    $errorCorrection\
    $algorithm\
    $minLenOpt\
    $flyeOpt\
    $unpaired\
    """
    //2>/dev/null
}

process renameFilterFasta {
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis",
        mode: 'copy'
    )
    input:
    path contigs

    output:
    path "*"

    script:
    def annotation = params.annotation ? "-ann 1" : ""
    """
    CONTIG_NUMBER=\$(grep -c '>' ${contigs})
    
    renameFilterFasta.pl \
    -u $contigs\
    -d .\
    -filt $params.minContigSize\
    -maxseq \$CONTIG_NUMBER\
    -ann_size $params.contigSizeForAnnotation\
    -n $params.projName\
    $annotation
    """

}

process bestIncompleteAssembly {
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

workflow {
    "mkdir nf_assets".execute().text
    "touch nf_assets/NO_FILE".execute().text
    "touch nf_assets/NO_FILE2".execute().text
    "touch nf_assets/NO_FILE3".execute().text
    "touch nf_assets/NO_FILE4".execute().text

    paired_ch = channel.fromPath(params.pairedFiles, relative:true, checkIfExists:true).collect()
    unpaired_ch = channel.fromPath(params.unpairedFile, relative:true, checkIfExists:true)
    spades_pb = file(params.spades.pacbio, checkIfExists:true)
    spades_np = file(params.spades.nanopore, checkIfExists:true)
    unicycler_lr = file(params.unicycler.longreads, checkIfExists:true)

    if (params.assembler.equalsIgnoreCase("IDBA_UD")) {
        avg_len_ch = idbaAvgLen(idbaReadFastq(paired_ch, unpaired_ch))
        (c1,c2) = idbaPrepReads(paired_ch, unpaired_ch)
        (sp,su,l) = idbaExtractLong(c1,c2.ifEmpty({file("nf_assets/NO_FILE")}))
        idbaUD(sp.filter{ it.size()>0 }.ifEmpty({file("nf_assets/NO_FILE")}),
            su.filter{ it.size()>0 }.ifEmpty({file("nf_assets/NO_FILE2")}),
            l.filter{ it.size()>0 }.ifEmpty({file("nf_assets/NO_FILE3")}),
            avg_len_ch)
        
        bestIncompleteAssembly(idbaUD.out.contigs.ifEmpty('EMPTY'), idbaUD.out.intContigs)
        renameFilterFasta(idbaUD.out.contigs.concat(bestIncompleteAssembly.out).first())

    }
    else if (params.assembler.equalsIgnoreCase("SPAdes")) {
        spades(paired_ch, unpaired_ch, spades_pb, spades_np)
        
        bestIncompleteAssembly(spades.out.contigs.ifEmpty('EMPTY'), spades.out.intContigs)
        renameFilterFasta(spades.out.contigs.concat(bestIncompleteAssembly.out).first())
    }
    else if (params.assembler.equalsIgnoreCase("MEGAHIT")) {
        megahit(paired_ch, unpaired_ch)
        
        bestIncompleteAssembly(megahit.out.contigs.ifEmpty('EMPTY'), megahit.out.intContigs)
        renameFilterFasta(megahit.out.contigs.concat(bestIncompleteAssembly.out).first())
    }
    else if (params.assembler.equalsIgnoreCase("UniCycler")) {
        if (params.unicycler.longreads != "nf_assets/NO_FILE3") {
            println("Filter long reads with $params.unicycler.minLongReads (bp) cutoff")
            unicycler(paired_ch,
                unpaired_ch,
                unicyclerPrep(unicycler_lr).filter{it.size()>0}.ifEmpty({file("nf_assets/NO_FILE3")}))
            //unicycler produces no intermediate contigs, we let it error out above rather than try to rescue a failed assembly
            renameFilterFasta(unicycler.out.contigs)
        }
        else {
            unicycler(paired_ch, unpaired_ch, unicycler_lr)
            renameFilterFasta(unicycler.out.contigs)
        }
    }
    else if (params.assembler.equalsIgnoreCase("LRASM")) {
        lrasm(unpaired_ch)

        bestIncompleteAssembly(lrasm.out.contigs.ifEmpty('EMPTY'), lrasm.out.intContigs)
        renameFilterFasta(lrasm.out.contigs.concat(bestIncompleteAssembly.out).first())
    }
    else {
        error "Invalid assembler: $params.assembler"
    }

}