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

//files 
REFERENCE_FASTA = file("${baseDir}/NC_045512.2.fasta")
MASTERFILE = file("${baseDir}/sarscov2_masterfile.txt")
ADAPTERS = file("${baseDir}/All_adapters.fa")

if(params.SINGLE_END == false){ 
    input_read_ch = Channel
        .fromFilePairs("${params.INPUT}*_R{1,2}*.gz")
        .ifEmpty { error "Cannot find any FASTQ pairs in ${params.INPUT} ending with .gz" }
        .map { it -> [it[0], it[1][0], it[1][1]]}
    




process Trimming { 
    container "quay.io/biocontainers/trimmomatic:0.35--6"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
      tuple val(base), file(R1), file(R2) from input_read_ch
      file ADAPTERS
    output: 
      tuple val(base), file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz") into Trim_out_ch

    script:
    """
    #!/bin/bash

    trimmomatic PE ${R1} ${R2} ${base}.R1.paired.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}.R2.paired.fastq.gz ${base}.R2.unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:75

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

    cpus 4 
    memory '6 GB'

    script:
    """
    #!/bin/bash

    /usr/local/bin/bbmap.sh in1=${base}.R1.paired.fastq.gz in2=${base}.R2.paired.fastq.gz outm=${base}.bam ref=${REFERENCE_FASTA} -Xmx6g

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
    #/usr/local/miniconda/bin/samtools sort -n -O sam ${base}.clipped.sam > ${base}.clipped.sorted.sam
    #/usr/local/miniconda/bin/samtools view -Sb ${base}.clipped.sorted.sam > ${base}.clipped.unsorted.bam
    #/usr/local/miniconda/bin/samtools sort -o ${base}.clipped.unsorted.bam ${base}.clipped.bam
     /usr/local/miniconda/bin/samtools sort ${base}.clipped.sam -o ${base}.clipped.bam

    """
}

process generateConsensus {
    container "quay.io/greninger-lab/swift-pipeline"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val (base), file(BAMFILE) from Clipped_bam_ch
        file REFERENCE_FASTA
    output:
        file("${base}_swift.fasta")
        file("${base}.clipped.bam")

    publishDir params.OUTDIR, mode: 'copy'

    shell:
    '''
    #!/bin/bash
    /usr/local/miniconda/bin/bcftools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 500000 \\
        --fasta-ref !{REFERENCE_FASTA} \\
        --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
        !{BAMFILE} \\
        | /usr/local/miniconda/bin/bcftools call --output-type v --ploidy 1 --keep-alts --keep-masked-ref --multiallelic-caller --variants-only -P 0 \\
        | /usr/local/miniconda/bin/bcftools reheader --samples sample_name.list \\
        | /usr/local/miniconda/bin/bcftools view --output-file !{base}.vcf.gz --output-type z


    /usr/local/miniconda/bin/tabix -p vcf -f !{base}.vcf.gz


    cat !{REFERENCE_FASTA} | /usr/local/miniconda/bin/bcftools consensus !{base}.vcf.gz > !{base}.consensus.fa

    /usr/local/miniconda/bin/bedtools genomecov \\
        -bga \\
        -ibam !{BAMFILE} \\
        -g !{REFERENCE_FASTA} \\
        | awk '\$4 < 2' | /usr/local/miniconda/bin/bedtools merge > !{base}.mask.bed
    
    /usr/local/miniconda/bin/bedtools maskfasta \\
        -fi !{base}.consensus.fa \\
        -bed !{base}.mask.bed \\
        -fo !{base}.consensus.masked.fa

    cat !{REFERENCE_FASTA} !{base}.consensus.masked.fa > align_input.fasta
    /usr/local/miniconda/bin/mafft --auto align_input.fasta > repositioned.fasta
    awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' repositioned.fasta > repositioned_unwrap.fasta
    
    python3 !{baseDir}/trim_ends.py !{base}



    '''
}


} else { 
    input_read_ch = Channel
        .fromPath("${params.INPUT}*.gz")
        //.map { it -> [ file(it)]}
        .map { it -> file(it)}


        

    // # This approach gives a fasta with all N's (work/84/7e5c92/), vcf with lines starting from 20000s
    // #/usr/local/miniconda/bin/samtools mpileup --max-depth 500000 -uf !{REFERENCE_FASTA} !{base}.clipped.bam | \
    // #/usr/local/miniconda/bin/bcftools call -c -o !{base}.consensus.vcf
    // #/usr/local/miniconda/bin/vcfutils.pl vcf2fq !{base}.consensus.vcf > !{base}.consensus.fastq 
    // #/usr/local/miniconda/bin/seqtk seq -aQ64 -q20 -n N !{base}.consensus.fastq > !{base}.consensus.fasta

    // # This approach gives a fasta identical to ref, blank vcf
    // #/usr/local/miniconda/bin/bcftools mpileup -Ou --max-depth 500000 -f !{REFERENCE_FASTA} !{BAMFILE} | \
    // #/usr/local/miniconda/bin/bcftools call -c -o !{base}.vcf
    // ##/usr/local/miniconda/bin/bcftools call -mv -Oz -o !{base}.vcf.gz
    // #/usr/local/miniconda/bin/bgzip !{base}.vcf
    // #/usr/local/miniconda/bin/bcftools index !{base}.vcf.gz
    // #cat !{REFERENCE_FASTA} | /usr/local/miniconda/bin/bcftools consensus !{base}.vcf.gz > !{base}.consensus.fasta

    // This approach gave a bunch of ns in between
    // /usr/local/miniconda/bin/samtools mpileup -uf !{REFERENCE_FASTA} !{BAMFILE} | /usr/local/miniconda/bin/bcftools call -c | /usr/local/miniconda/bin/vcfutils.pl vcf2fq > out.fastq
    // /usr/local/miniconda/bin/samtools mpileup -uf !{REFERENCE_FASTA} !{BAMFILE} | /usr/local/miniconda/bin/bcftools call -c | /usr/local/miniconda/bin/bcftools view 

    // /usr/local/miniconda/bin/seqtk seq -aQ64 -q20 -n N out.fastq > !{base}.consensus.fasta
    
// Single end



process Trimming_SE { 
    container "quay.io/biocontainers/trimmomatic:0.35--6"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
      file R1 from input_read_ch
      file ADAPTERS
    output: 
      file "*.trimmed.fastq.gz" into Trim_out_ch_SE

    script:
    """
    #!/bin/bash

    base=`basename ${R1} ".fastq.gz"`

    echo \$base

    trimmomatic SE ${R1} \$base.trimmed.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:75

    """

}


process Aligning_SE {
     container "quay.io/biocontainers/bbmap:38.86--h1296035_0"
    //container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7	"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input: 
      file R1 from Trim_out_ch_SE
      file REFERENCE_FASTA
    output:
       file "*.bam"  into Aligned_bam_ch_SE

    cpus 4 
    memory '6 GB'

    script:
    """
    #!/bin/bash

    base=`basename ${R1} ".trimmed.fastq.gz"`
    /usr/local/bin/bbmap.sh in1="\$base".trimmed.fastq.gz  outm="\$base".bam ref=${REFERENCE_FASTA} -Xmx6g

    """
}

process NameSorting_SE { 
    container "quay.io/biocontainers/samtools:1.3--h0592bc0_3"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        file(ALIGNED) from Aligned_bam_ch_SE
    output:
        file "*.sorted.sam" into Sorted_sam_ch_SE

    script:
    """
    #!/bin/bash
    base=`basename ${ALIGNED} ".bam"`
    samtools sort -n -O sam ${ALIGNED} > \${base}.sorted.sam

    """
}

process Clipping_SE { 
    container "quay.io/greninger-lab/swift-pipeline"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
      file SORTED_SAM  from Sorted_sam_ch_SE
      file MASTERFILE
    output:
      file "*.clipped.bam"  into Clipped_bam_ch_SE

    script:
    """
    #!/bin/bash
    R1=`basename ${SORTED_SAM} ".sorted.sam"`
    /./root/.local/bin/primerclip -s ${MASTERFILE} ${SORTED_SAM} \${R1}.clipped.sam
    #/usr/local/miniconda/bin/samtools sort -n -O sam \${R1}.clipped.sam > \${R1}.clipped.sorted.sam
    #/usr/local/miniconda/bin/samtools view -Sb \${R1}.clipped.sorted.sam > \${R1}.clipped.unsorted.bam
    #/usr/local/miniconda/bin/samtools sort -o \${R1}.clipped.unsorted.bam \${R1}.clipped.bam
     /usr/local/miniconda/bin/samtools sort \${R1}.clipped.sam -o \${R1}.clipped.bam

    """
}

process generateConsensus_SE {
    container "quay.io/greninger-lab/swift-pipeline"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        file(BAMFILE) from Clipped_bam_ch_SE
        file REFERENCE_FASTA
    output:
        file("*_swift.fasta")
        file(BAMFILE)

    publishDir params.OUTDIR, mode: 'copy'

    shell:
    '''
    #!/bin/bash

    R1=`basename !{BAMFILE} .clipped.bam`
    echo \${R1}

    /usr/local/miniconda/bin/bcftools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 500000 \\
        --fasta-ref !{REFERENCE_FASTA} \\
        --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
        !{BAMFILE} \\
        | /usr/local/miniconda/bin/bcftools call --output-type v --ploidy 1 --keep-alts --keep-masked-ref --multiallelic-caller --variants-only -P 0 \\
        | /usr/local/miniconda/bin/bcftools reheader --samples sample_name.list \\
        | /usr/local/miniconda/bin/bcftools view --output-file \${R1}.vcf.gz --output-type z


    /usr/local/miniconda/bin/tabix -p vcf -f \${R1}.vcf.gz


    cat !{REFERENCE_FASTA} | /usr/local/miniconda/bin/bcftools consensus \${R1}.vcf.gz > \${R1}.consensus.fa

    /usr/local/miniconda/bin/bedtools genomecov \\
        -bga \\
        -ibam !{BAMFILE} \\
        -g !{REFERENCE_FASTA} \\
        | awk '\$4 < 2' | /usr/local/miniconda/bin/bedtools merge > \${R1}.mask.bed
    
    /usr/local/miniconda/bin/bedtools maskfasta \\
        -fi \${R1}.consensus.fa \\
        -bed \${R1}.mask.bed \\
        -fo \${R1}.consensus.masked.fa

    cat !{REFERENCE_FASTA} \${R1}.consensus.masked.fa > align_input.fasta
    /usr/local/miniconda/bin/mafft --auto align_input.fasta > repositioned.fasta
    awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' repositioned.fasta > repositioned_unwrap.fasta
    
    python3 !{baseDir}/trim_ends.py \${R1}



    '''
    }
}
