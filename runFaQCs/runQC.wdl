# RunQC workflow
# Declare WDL version 1.0 if working in Terra
version 1.0


workflow runQC{
  input{
    String outDir
    Array[File] inputFastq
    Boolean pairedFile

    String? trimMode
    Int? trimQual
    Int? trim5end
    Int? trim3end
    Boolean? trimAdapter
    Float? trimRate
    Boolean? trimPolyA
    File? artifactFile

    Int? minLen
    Int? avgQual
    Int? numN
    Float? filtLC
    Boolean? filtPhiX

    Int? ascii
    Int? outAscii

    String? outPrefix
    File? outStats

    Int? numCPU
    Int? splitSize
    Boolean? qcOnly
    Boolean? kmerCalc
    Int? kmerNum
    Int? splitSubset
    Boolean? discard
    Boolean? substitute
    Boolean? trimOnly
    Int? replaceGN
    Boolean? trim5off
    Boolean? debug
  }

  call faqcs{
    input:
    outDir = outDir,
    inputFastq = inputFastq,
    pairedFile = pairedFile,

    trimMode = trimMode,
    trimQual = trimQual,
    trim5end = trim5end,
    trim3end = trim3end,
    trimAdapter = trimAdapter,
    trimRate = trimRate,
    trimPolyA = trimPolyA,
    artifactFile = artifactFile,

    minLen = minLen,
    avgQual = avgQual,
    numN = numN,
    filtLC = filtLC,
    filtPhiX = filtPhiX,

    ascii = ascii,
    outAscii = outAscii,

    outPrefix = outPrefix,
    outStats = outStats,

    numCPU = numCPU,
    splitSize = splitSize,
    qcOnly = qcOnly,
    kmerCalc = kmerCalc,
    kmerNum = kmerNum,
    splitSubset = splitSubset,
    discard = discard,
    substitute = substitute,
    trimOnly = trimOnly,
    replaceGN = replaceGN,
    trim5off = trim5off,
    debug = debug
  }
  output {
    Array[File?] outputFiles = faqcs.outputFiles
    Array[File?] trimmedFastq = faqcs.trimmedFastq
    File? stats = faqcs.stats
    File? outPDF = faqcs.outPDF  
  }
}




task faqcs {
  input{
    String outDir
    Array[File] inputFastq
    Boolean pairedFile

    String? trimMode
    Int? trimQual
    Int? trim5end
    Int? trim3end
    Boolean? trimAdapter
    Float? trimRate
    Boolean? trimPolyA
    File? artifactFile

    Int? minLen
    Int? avgQual
    Int? numN
    Float? filtLC
    Boolean? filtPhiX

    Int? ascii
    Int? outAscii

    String? outPrefix
    File? outStats

    Int? numCPU
    Int? splitSize
    Boolean? qcOnly
    Boolean? kmerCalc
    Int? kmerNum
    Int? splitSubset
    Boolean? discard
    Boolean? substitute
    Boolean? trimOnly
    Int? replaceGN
    Boolean? trim5off
    Boolean? debug
  }
    # FaQCs [options] [-u unpaired.fastq] -p reads1.fastq reads2.fastq -d out_directory
  command <<<

    # ln /localization/path/file.txt /analysis/path/file.txt

    FaQCs ~{"--mode" + trimMode} \
    ~{"-q" + trimQual} \
    ~{"--5end" + trim5end} \
    ~{"--3end" + trim3end} \
    ~{true="--adapter True" false="" trimAdapter} \
    ~{"--rate" + trimRate} \
    ~{true="--polyA True" false="" trimPolyA} \
    ~{"--artifactFile" + artifactFile} \
    ~{"--min_L" + minLen} \
    ~{"--avg_q" + avgQual} \
    ~{"-n" + numN} \
    ~{"--lc" + filtLC} \
    ~{true="--phiX True" false="" filtPhiX} \
    ~{"--ascii" + ascii} \
    ~{"--out_ascii" + outAscii} \
    ~{"--prefix" + outPrefix} \
    ~{"--stats" + outStats} \
    ~{"-t" + numCPU} \
    ~{"--split_size" + splitSize} \
    ~{true="--qc_only True" false="" qcOnly} \
    ~{true="--kmer_rarefaction True" false="" kmerCalc} \
    ~{"-m" + kmerNum} \
    ~{"--subset" + splitSubset} \
    ~{true="--discard True" false="" discard} \
    ~{true="--substitute True" false="" substitute} \
    ~{true="--trim_only True" false="" trimOnly} \
    ~{"--replace_to_N_q" + replaceGN} \
    ~{true="--5trim_off True" false="" trim5off} \
    ~{true="--debug True" false="" debug} \
    ~{if pairedFile then "-1 "+ inputFastq[0]+ " -2 "+inputFastq[1] else "-u "+inputFastq[0] } \
    ~{"-d " + outDir}
  >>>

  output {
    Array[File] outputFiles = glob("${outDir}/*")
    Array[File] trimmedFastq = glob("${outDir}/*.trimmed.fastq")
    File? stats = if (defined(outPrefix)) then "${outDir}/${outPrefix}.stats.txt" else  "${outDir}/QC.stats.txt"
    File? outPDF =  if (defined(outPrefix)) then "${outDir}/${outPrefix}_qc_report.pdf" else "${outDir}/QC_qc_report.pdf"
  }

  runtime {
        docker: "kaijli/runqc:1.0"
        continueOnReturnCode: true
    }
}
