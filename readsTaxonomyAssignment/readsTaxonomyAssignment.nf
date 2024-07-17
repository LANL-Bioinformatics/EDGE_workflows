#!/usr/bin/env nextflow

process lenFile {
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

process avgLen {
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

process readsTaxonomy {
    input:
    path paired
    path unpaired
    path settings

    output:

    script:
    def debugging = (params.debugFlag != null && params.debugFlag) ? "--debug " : ""
    def numCPU = params.cpus != null ? params.cpus : 8
    """
    cat $paired $unpaired > allReads.fastq
    microbial_profiling.pl $debugging -o ./RTA \
    -s microbial_profiling.settings.ini \
    -c $numCPU \
    allReads.fastq 
    """
    //2>error.log
}

process readsTaxonomyConfig {
    publishDir(
        path: "$params.outDir",
        mode: 'copy'
    )

    input:
    val avgLen

    output:
    path "error.log"
    path "microbial_profiling.settings.ini", emit: settings 

    script:
    def bwaScoreCut = 30
    if (params.fastq_source != null && (params.fastq_source.equalsIgnoreCase("nanopore") || params.fastq_source.equalsIgnoreCase("pacbio"))) {
        if (params.minLen > 1000) {
            bwaScoreCut = params.minLen
        } 
        else {
            bwaScoreCut=1000
        }
    }
    else{
        bwaScoreCut = (avgLen as Integer)*0.8
    }
    bwaScoreCut = bwaScoreCut as Integer
    def tools = params.enabledTools != null ? "-tools \'$params.enabledTools\'" : ""
    def template = params.template != null ? "-template $params.template " : ""

    def splittrim_minq = params.split_trim_minq != null ? "-splitrim-minq $params.split_trim_minq " : ""
    def bwa = params.custom_bwa_db != null ? "-bwa-db $params.custom_bwa_db " : ""
    def metaphlan = params.custom_metaphlan_db != null ? "-metaphlan-db $params.custom_metaphlan_db " : ""
    def kraken = params.custom_kraken_db != null ? "-kraken-db $params.custom_kraken_db " : ""
    def centrifuge = params.custom_centrifuge_db != null ? "-centrifuge-db $params.custom_centrifuge_db " : ""
    def pangia = params.custom_pangia_db != null ? "-pangia-db $params.custom_pangia_db " : ""
    def diamond = params.custom_diamond_db != null ? "-diamond-db $params.custom_diamond_db " : ""

    def gottcha_speDB_v = params.custom_gottcha_speDB_v != null ? "-gottcha-speDB-v $params.custom_gottcha_speDB_v " : ""
    def gottcha_speDB_b = params.custom_gottcha_speDB_b != null ? "-gottcha-speDB-b $params.custom_gottcha_speDB_b " : ""
    def gottcha_strDB_v = params.custom_gottcha_strDB_v != null ? "-gottcha-strDB-v $params.custom_gottcha_strDB_v " : ""
    def gottcha_strDB_b = params.custom_gottcha_strDB_b != null ? "-gottcha-strDB-b $params.custom_gottcha_strDB_b " : ""
    def gottcha_genDB_v = params.custom_gottcha_genDB_v != null ? "-gottcha-genDB-v $params.custom_gottcha_genDB_v " : ""
    def gottcha_genDB_b = params.custom_gottcha_genDB_b != null ? "-gottcha-genDB-b $params.custom_gottcha_genDB_b " : ""

    def gottcha2_genDB_v = params.custom_gottcha2_genDB_v != null ? "-gottcha2-genDB-v $params.custom_gottcha2_genDB_v " : ""
    def gottcha2_speDB_v = params.custom_gottcha2_speDB_v != null ? "-gottcha2-speDB-v $params.custom_gottcha2_speDB_v " : ""
    def gottcha2_speDB_b = params.custom_gottcha2_speDB_b != null ? "-gottcha2-speDB-b $params.custom_gottcha2_speDB_b " : ""

    def np = (params.fastq_source != null && params.fastq_source.equalsIgnoreCase("nanopore")) ? "--nanopore " : ""

    """
    microbial_profiling_configure.pl $template \
    $tools -bwaScoreCut $bwaScoreCut\
    $bwa\
    $metaphlan\
    $kraken\
    $centrifuge\
    $pangia\
    $diamond\
    $gottcha_speDB_v\
    $gottcha_speDB_b\
    $gottcha_strDB_v\
    $gottcha_strDB_b\
    $gottcha_genDB_v\
    $gottcha_genDB_b\
    $gottcha2_genDB_v\
    $gottcha2_speDB_v\
    $gottcha2_speDB_b\
    $np >microbial_profiling.settings.ini 2>error.log
    """

}

workflow {
    "mkdir nf_assets".execute().text
    "touch nf_assets/NO_FILE".execute().text
    "touch nf_assets/NO_FILE2".execute().text 
    paired_ch = channel.fromPath(params.pairFile, checkIfExists:true).collect()
    unpaired_ch = channel.fromPath(params.unpairFile, checkIfExists:true)
    avg_len_ch = avgLen(lenFile(paired_ch, unpaired_ch))
    readsTaxonomyConfig(avg_len_ch)
    readsTaxonomy(paired_ch, unpaired_ch, readsTaxonomyConfig.out.settings)
    //TODO: allow unmapped read input
    //TODO: access EDGE databases and troubleshoot main profiling script.
    //TODO: clean up output
}
