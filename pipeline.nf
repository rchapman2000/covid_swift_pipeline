#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    covid pipeline with primerclip :)

    Usage: 

    An example command for running the pipeline is as follows:
    nextflow run vpeddu/lava \\
        --INPUT         Input folder where all fastqs are located.
                        ./ can be used for current directory.
                        Fastqs should all be gzipped. This can be done with the command gzip *.fastq. [REQUIRED]
        
        --OUTDIR        Output directory. [REQUIRED]
        
        --SINGLE_END    Optional flag for single end reads. By default, this pipeline does 
                        paired-end reads.
        
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

params.INPUT = false
params.OUTDIR= false
params.SINGLE_END = false

// if INPUT not set
if (params.INPUT == false) {
    println( "Must provide an input directory with --INPUT") 
    exit(1)
}
// Make sure INPUT ends with trailing slash
if (!params.INPUT.endsWith("/")){
   params.INPUT = "${params.INPUT}/"
}
// if OUTDIR not set
if (params.OUTDIR == false) {
    println( "Must provide an output directory with --OUTDIR") 
    exit(1)
}
// Make sure OUTDIR ends with trailing slash
if (!params.OUTDIR.endsWith("/")){
   params.OUTDIR = "${params.OUTDIR}/"
}

if(params.SINGLE_END == false){ 
    input_read_ch = Channel
        .fromFilePairs("${params.INPUT}*_R{1,2}*.gz")
        .ifEmpty { error "Cannot find any FASTQ pairs in ${params.INPUT} ending with .gz" }
        .map { it -> [it[0], it[1][0], it[1][1]]}
    } else { 
    input_read_ch = Channel
        .fromPath("${params.INPUT}*.gz")
        .map { it -> [ file(it)]}
}

REFERENCE_FASTA = file("${baseDir}/NC_045512.2.fasta")
MASTERFILE = file("${baseDir}/sarscov2_masterfile.txt")



process Trimming { 
    container "quay.io/biocontainers/trimmomatic:0.35--6"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
      tuple val(base), file(R1), file(R2) from input_read_ch
    output: 
      tuple val(base), file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz") into Trim_out_ch

    script:
    """
    #!/bin/bash

    trimmomatic PE ${R1} ${R2} ${base}.R1.paired.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}.R2.paired.fastq.gz ${base}.R2.unpaired.fastq.gz \
    ILLUMINACLIP:\$HOME/downloads/trimmomatic-0.38/adapters/All_adapters.fa:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:75

    """
}

process Aligning {
     container "quay.io/biocontainers/bbmap:38.86--h1296035_0"
    //container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7	"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input: 
      tuple val(base), file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz") from Trim_out_ch
      file REFERENCE_FASTA
    output:
      tuple val (base), file("${base}.bam") into Aligned_bam_ch

    script:
    """
    #!/bin/bash

    /usr/local/bin/bbmap.sh in1=${base}.R1.paired.fastq.gz in2=${base}.R2.paired.fastq.gz out=${base}.bam ref=${REFERENCE_FASTA}

    """

    // bwa?
    // script:
    // """
    // #!/bin/bash
    // bwa mem $workflow.projectDir/NC_045512.2.fasta ${base}.R1.paired.fastq.gz ${base}.R2.paired.fastq.gz > aln.sam
    // """
}

process NameSorting { 
    container "quay.io/biocontainers/samtools:1.3--h0592bc0_3"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
      tuple val (base), file("${base}.bam") from Aligned_bam_ch
    output:
      tuple val (base), file("${base}.sorted.sam") into Sorted_sam_ch

    script:
    """
    #!/bin/bash
    samtools sort -n -O sam ${base}.bam > ${base}.sorted.sam

    """
}

process Clipping { 
    container "quay.io/greninger-lab/swift-pipeline"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
      tuple val (base), file("${base}.sorted.sam") from Sorted_sam_ch
      file MASTERFILE
    output:
      tuple val (base), file("${base}.clipped.bam") into Clipped_bam_ch

    script:
    """
    #!/bin/bash
    /./root/.local/bin/primerclip ${MASTERFILE} ${base}.sorted.sam ${base}.clipped.sam
    /usr/local/miniconda/bin/samtools view -S -b ${base}.clipped.sam > ${base}.clipped.bam

    """
}
