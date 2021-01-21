#!/usr/bin/env nextflow
nextflow.preview.dsl=2

def helpMessage() {
    log.info"""
    Usage: 
    An example command for running the pipeline is as follows:
    nextflow run greninger-lab/covid_swift_pipeline -resume -with-docker ubuntu:18.04 --INPUT example/ --OUTDIR example/output/ --PRIMERS v2
    
    Parameters:
        --INPUT         Input folder where all fastqs are located.
                        ./ can be used for current directory.
                        Fastqs should all be gzipped. This can be done with the command gzip *.fastq. [REQUIRED]
        --OUTDIR        Output directory. [REQUIRED]
        --PRIMERS       Primer masterfile to run. By default, this pipeline uses the original Swift V1 primers. --PRIMERS V2 can be specified for the Swift V2 primer set.
        --SINGLE_END    Optional flag for single end reads. By default, this pipeline does 
                        paired-end reads.

        -with-docker ubuntu:18.04   [REQUIRED]
        -resume [RECOMMENDED]
        
    """.stripIndent()
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*          SET UP CONFIGURATION VARIABLES            */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

params.INPUT = false
params.OUTDIR= false
params.SINGLE_END = false
params.PRIMERS = false

TRIM_ENDS=file("${baseDir}/trim_ends.py")
VCFUTILS=file("${baseDir}/vcfutils.pl")
SPLITCHR=file("${baseDir}/splitchr.txt")

// if INPUT not set
if (params.INPUT == false) {
    println( "Must provide an input directory with --INPUT") 
    exit(1)
}
// Make sure INPUT ends with trailing slash
if (!params.INPUT.endsWith("/")){
    println("Make sure your input directory ends with trailing slash.")
   exit(1)
}
// if OUTDIR not set
if (params.OUTDIR == false) {
    println( "Must provide an output directory with --OUTDIR") 
    exit(1)
}
// Make sure OUTDIR ends with trailing slash
if (!params.OUTDIR.endsWith("/")){
   println("Make sure your output directory ends with trailing slash.")
   exit(1)
}

//// 
//Setting up files 
////

REFERENCE_FASTA = file("${baseDir}/NC_045512.2.fasta")
REFERENCE_FASTA_FAI = file("${baseDir}/NC_045512.2.fasta.fai")
if (params.PRIMERS == "V2" | params.PRIMERS == "v2") {
    MASTERFILE = file("${baseDir}/sarscov2_v2_masterfile.txt")
    println("Using Swift V2 primerset...")
}
else {
    MASTERFILE = file("${baseDir}/sarscov2_masterfile.txt")
    println("Using Swift V1 primerset [default]...")
}
ADAPTERS = file("${baseDir}/All_adapters.fa")
FIX_COVERAGE = file("${baseDir}/fix_coverage.py")
PROTEINS = file("${baseDir}/NC_045512_proteins.txt")
AT_REFGENE = file("${baseDir}/annotation/AT_refGene.txt")
AT_REFGENE_MRNA = file("${baseDir}/annotation/AT_refGeneMrna.fa")
LAVA_GFF = file("${baseDir}/annotation/lava_ref.gff")
MAT_PEPTIDES = file("${baseDir}/annotation/mat_peptides.txt")
MAT_PEPTIDE_ADDITION = file("${baseDir}/annotation/mat_peptide_addition.py")
RIBOSOMAL_START = file("${baseDir}/annotation/ribosomal_start.txt")
RIBOSOMAL_SLIPPAGE = file("${baseDir}/annotation/ribosomal_slippage.py")
PROTEINS = file("${baseDir}/annotation/proteins.csv")
CORRECT_AF = file("${baseDir}/annotation/correct_AF.py")
CORRECT_AF_BCFTOOLS = file("${baseDir}/annotation/correct_AF_bcftools.py")

//
// Import processes
// 

include Trimming from './modules'
include Cutadapt from './modules'
include Fastqc from './modules'
include Aligning from './modules'
include Trimming_SE from './modules'
include Fastqc_SE from './modules'
include Clipping from './modules'
include NameSorting from './modules'
include GenerateConsensus from './modules'
include Lofreq from './modules'
include AnnotateVariants_Bcftools from './modules'
include AnnotateVariants_Lofreq from './modules'
include Varscan2 from './modules'

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

workflow {
    // Paired end trimming
    if(params.SINGLE_END == false){ 
        input_read_ch = Channel
            .fromFilePairs("${params.INPUT}*_R{1,2}*.gz")
            .ifEmpty { error "Cannot find any FASTQ pairs in ${params.INPUT} ending with .gz" }
            .map { it -> [it[0], it[1][0], it[1][1]]}
        
        Trimming (
            input_read_ch,
            ADAPTERS
        )
        Cutadapt (
            Trimming.out[0],
            REFERENCE_FASTA
        )
        Fastqc (
            Trimming.out[1]
        )
    }
    // Single end first three steps
    else {
        input_read_ch = Channel
            .fromPath("${params.INPUT}*.gz")
            //.map { it -> [ file(it)]}
            .map { it -> file(it)}
        
        Trimming_SE (
            input_read_ch,
            ADAPTERS
        )
        Cutadapt (
            Trimming_SE.out[0],
            REFERENCE_FASTA
        )
        Fastqc_SE (
            Trimming_SE.out[1]
        )
    }

    Aligning (
        Cutadapt.out[0],
        REFERENCE_FASTA
    )
    
    NameSorting (
        Aligning.out[0]
    )

    Clipping (
        NameSorting.out[0],
        MASTERFILE
    )

    // GenerateConsensus (
    //     Clipping.out[0],
    //     REFERENCE_FASTA,
    //     TRIM_ENDS,
    //     FIX_COVERAGE,
    //     VCFUTILS,
    //     REFERENCE_FASTA_FAI,
    //     SPLITCHR
    // )

    Lofreq (
        Clipping.out[1],
        REFERENCE_FASTA
    )
    // AnnotateVariants_Bcftools (
    //     GenerateConsensus.out[5],
    //     MAT_PEPTIDES,
    //     MAT_PEPTIDE_ADDITION,
    //     RIBOSOMAL_SLIPPAGE,
    //     RIBOSOMAL_START,
    //     PROTEINS,
    //     AT_REFGENE,
    //     AT_REFGENE_MRNA,
    //     CORRECT_AF_BCFTOOLS
    // )
    
    // AnnotateVariants_Lofreq (
    //     Lofreq.out[0],
    //     MAT_PEPTIDES,
    //     MAT_PEPTIDE_ADDITION,
    //     RIBOSOMAL_SLIPPAGE,
    //     RIBOSOMAL_START,
    //     PROTEINS,
    //     AT_REFGENE,
    //     AT_REFGENE_MRNA,
    //     CORRECT_AF
    // )

    Varscan2 ( 
        Clipping.out[2],
        REFERENCE_FASTA,
        REFERENCE_FASTA_FAI,
        SPLITCHR
    )

    AnnotateVariants (
        Varscan.out[0].groupTuple(
            ).join(
                Lofreq.out[1].groupTuple()
            )
    )

}
