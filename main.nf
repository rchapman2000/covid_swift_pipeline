#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Usage: 
    An example command for running the pipeline is as follows:
    nextflow run greninger-lab/covid_swift_pipeline -resume -with-docker ubuntu:18.04 --INPUT example/ --OUTDIR example/output/
    
    Parameters:
        --INPUT         Input folder where all fastqs are located.
                        ./ can be used for current directory.
                        Fastqs should all be gzipped. This can be done with the command gzip *.fastq. [REQUIRED]
        --OUTDIR        Output directory. [REQUIRED]
        --SINGLE_END    Optional flag for single end reads. By default, this pipeline does 
                        paired-end reads.

        -with-docker ubuntu:18.04   [REQUIRED]
        -resume [RECOMMENDED]
        
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
TRIM_ENDS=file("${baseDir}/trim_ends.py")

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
      tuple val(base), file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz"),file("${base}.R1.unpaired.fastq.gz"), file("${base}.R2.unpaired.fastq.gz"),file("${base}_summary.csv") into Trim_out_ch
      tuple val(base), file(R1),file(R2),file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz"),file("${base}.R1.unpaired.fastq.gz"), file("${base}.R2.unpaired.fastq.gz") into Trim_out_ch2

    publishDir "${params.OUTDIR}trimmed_fastqs", mode: 'copy'

    script:
    """
    #!/bin/bash

    trimmomatic PE -threads ${task.cpus} ${R1} ${R2} ${base}.R1.paired.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}.R2.paired.fastq.gz ${base}.R2.unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

    num_r1_untrimmed=\$(gunzip -c ${R1} | wc -l)
    num_r2_untrimmed=\$(gunzip -c ${R2} | wc -l)
    num_untrimmed=\$((\$((num_r1_untrimmed + num_r2_untrimmed))/4))

    num_r1_paired=\$(gunzip -c ${base}.R1.paired.fastq.gz | wc -l)
    num_r2_paired=\$(gunzip -c ${base}.R2.paired.fastq.gz | wc -l)
    num_paired=\$((\$((num_r1_paired + num_r2_paired))/4))

    num_r1_unpaired=\$(gunzip -c ${base}.R1.unpaired.fastq.gz | wc -l)
    num_r2_unpaired=\$(gunzip -c ${base}.R2.unpaired.fastq.gz | wc -l)
    num_unpaired=\$((\$((num_r1_unpaired + num_r2_unpaired))/4))

    num_trimmed=\$((num_paired + num_unpaired))
    
    percent_trimmed=\$((100-\$((100*num_trimmed/num_untrimmed))))
    
    echo Sample_Name,Raw_Reads,Trimmed_Paired_Reads,Trimmed_Unpaired_Reads,Total_Trimmed_Reads,Percent_Trimmed,Mapped_Reads,Clipped_Mapped_Reads,Mean_Coverage,Percent_N > ${base}_summary.csv
    printf "${base},\$num_untrimmed,\$num_paired,\$num_unpaired,\$num_trimmed,\$percent_trimmed" >> ${base}_summary.csv
    
    cp .command.log ${base}_trimmosummary.txt

    """
}

process fastQc {
    container "quay.io/biocontainers/fastqc:0.11.9--0"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
      tuple val(base), file(R1),file(R2),file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz"),file("${base}.R1.unpaired.fastq.gz"), file("${base}.R2.unpaired.fastq.gz") from Trim_out_ch2
    output: 
      file("*fastqc*") into Fastqc_ch 

    publishDir "${params.OUTDIR}fastqc", mode: 'copy'

    script:
    """
    #!/bin/bash

    /usr/local/bin/fastqc ${R1} ${R2} ${base}.R1.paired.fastq.gz ${base}.R2.paired.fastq.gz

    """
}

process Aligning {
     container "quay.io/biocontainers/bbmap:38.86--h1296035_0"
    //container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7	"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input: 
      tuple val(base), file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz"),file("${base}.R1.unpaired.fastq.gz"), file("${base}.R2.unpaired.fastq.gz"),file("${base}_summary.csv") from Trim_out_ch
      file REFERENCE_FASTA
    output:
      tuple val (base), file("${base}.bam"),file("${base}_summary.csv") into Aligned_bam_ch

    cpus 4 
    memory '6 GB'

    script:
    """
    #!/bin/bash

    cat ${base}*.fastq.gz > ${base}_cat.fastq.gz
    /usr/local/bin/bbmap.sh in=${base}_cat.fastq.gz outm=${base}.bam ref=${REFERENCE_FASTA} -Xmx6g
    reads_mapped=\$(cat .command.log | grep "mapped:" | cut -d\$'\\t' -f3)
    printf ",\$reads_mapped" >> ${base}_summary.csv

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
      tuple val (base), file("${base}.bam"),file("${base}_summary.csv") from Aligned_bam_ch
    output:
      tuple val (base), file("${base}.sorted.sam"),file("${base}_summary.csv") into Sorted_sam_ch

    script:
    """
    #!/bin/bash
    samtools sort -@ ${task.cpus} -n -O sam ${base}.bam > ${base}.sorted.sam

    """
}

process Clipping { 
    container "quay.io/greninger-lab/swift-pipeline:latest"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
      tuple val (base), file("${base}.sorted.sam"),file("${base}_summary.csv") from Sorted_sam_ch
      file MASTERFILE
    output:
      tuple val (base), file("${base}.clipped.bam"), file("*.bai"),file("${base}_summary.csv") into Clipped_bam_ch
      tuple val (base), file("${base}.clipped.bam"), file("*.bai") into Clipped_bam_ch2

    script:
    """
    #!/bin/bash
    /./root/.local/bin/primerclip -s ${MASTERFILE} ${base}.sorted.sam ${base}.clipped.sam
    #/usr/local/miniconda/bin/samtools sort -@ ${task.cpus} -n -O sam ${base}.clipped.sam > ${base}.clipped.sorted.sam
    #/usr/local/miniconda/bin/samtools view -@ ${task.cpus} -Sb ${base}.clipped.sorted.sam > ${base}.clipped.unsorted.bam
    #/usr/local/miniconda/bin/samtools sort -@ ${task.cpus} -o ${base}.clipped.unsorted.bam ${base}.clipped.bam
     /usr/local/miniconda/bin/samtools sort -@ ${task.cpus} ${base}.clipped.sam -o ${base}.clipped.bam
     /usr/local/miniconda/bin/samtools index ${base}.clipped.bam

    clipped_reads=\$(cat .command.log | grep "Total mapped alignments" | cut -d\$'\\t' -f2)

    meancoverage=\$(/usr/local/miniconda/bin/samtools depth -a ${base}.clipped.bam | awk '{sum+=\$3} END { print sum/NR}')
    printf ",\$clipped_reads,\$meancoverage" >> ${base}_summary.csv

    """
}

process lofreq {
    container "quay.io/biocontainers/lofreq:2.1.5--py38h1bd3507_3"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
      tuple val (base), file("${base}.clipped.bam"), file("*.bai") from Clipped_bam_ch2
      file REFERENCE_FASTA
    output:
      file("${base}_lofreq.vcf")
    
    publishDir params.OUTDIR, mode: 'copy'

    script:
    """
    #!/bin/bash
    /usr/local/bin/lofreq call -f ${REFERENCE_FASTA} -o ${base}_lofreq.vcf ${base}.clipped.bam

    """
}

process generateConsensus {
    container "quay.io/greninger-lab/swift-pipeline:latest"
    echo true

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val (base), file(BAMFILE),file(INDEX_FILE),file("${base}_summary.csv") from Clipped_bam_ch
        file REFERENCE_FASTA
        file TRIM_ENDS
    output:
        file("${base}_swift.fasta")
        file("${base}.clipped.bam")
        file("${base}_bcftools.vcf")
        file("${base}.pileup")
        file(INDEX_FILE)
        file("${base}_summary.csv")

    publishDir params.OUTDIR, mode: 'copy'

    shell:
    '''
    #!/bin/bash

    R1=`basename !{BAMFILE} .clipped.bam`

    # Using LAVA consensus calling!
    /usr/local/miniconda/bin/bcftools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 500000 \\
        --max-idepth 500000 \\
        --fasta-ref !{REFERENCE_FASTA} \\
        --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
        --threads 10 \\
        !{BAMFILE} > \${R1}.pileup
    
    cat \${R1}.pileup | /usr/local/miniconda/bin/bcftools call -m -Oz -o \${R1}_pre.vcf.gz
    
    /usr/local/miniconda/bin/tabix \${R1}_pre.vcf.gz
    gunzip \${R1}_pre.vcf.gz
     /usr/local/miniconda/bin/bcftools filter -i '(DP4[0]+DP4[1]) < (DP4[2]+DP4[3]) && ((DP4[2]+DP4[3]) > 0)' \${R1}_pre.vcf -o \${R1}.vcf
    /usr/local/miniconda/bin/bgzip \${R1}.vcf
    /usr/local/miniconda/bin/tabix \${R1}.vcf.gz 
    cat !{REFERENCE_FASTA} | /usr/local/miniconda/bin/bcftools consensus \${R1}.vcf.gz > \${R1}.consensus.fa

    if [ -s !{BAMFILE} ]
    then
        /usr/local/miniconda/bin/bedtools genomecov \\
            -bga \\
            -ibam !{BAMFILE} \\
            -g !{REFERENCE_FASTA} \\
            | awk '\$4 < 10' | /usr/local/miniconda/bin/bedtools merge > \${R1}.mask.bed
        /usr/local/miniconda/bin/bedtools maskfasta \\
        -fi \${R1}.consensus.fa \\
        -bed \${R1}.mask.bed \\
        -fo \${R1}.consensus.masked.fa

        cat !{REFERENCE_FASTA} \${R1}.consensus.masked.fa > align_input.fasta
        /usr/local/miniconda/bin/mafft --auto align_input.fasta > repositioned.fasta
        awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' repositioned.fasta > repositioned_unwrap.fasta
        
        python3 !{TRIM_ENDS} \${R1}

        gunzip \${R1}.vcf.gz
        mv \${R1}.vcf \${R1}_bcftools.vcf
    else
       printf '>\$R1\n' > \${R1}_swift.fasta
       printf 'n%.0s' {1..29539}
    fi
    
    # Find percent ns, doesn't work, fix later in python script
    num_bases=$(grep -v ">" \${R1}_swift.fasta | wc | awk '{print $3-$1}')
    num_ns=$(grep -v ">" \${R1}_swift.fasta | awk -F"n" '{print NF-1}')
    percent_n=$(($num_ns/$num_bases*100))
    printf ",\$percent_n" >> \${R1}_summary.csv

    cat \${R1}_summary.csv | tr -d "[:blank:]" > a.tmp
    mv a.tmp \${R1}_summary.csv

    [ -s \${R1}_swift.fasta ] || echo "WARNING: \${R1} produced blank output. Manual review may be needed."

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
    /usr/local/bin/bbmap.sh in1="\$base".trimmed.fastq.gz  outm="\$base".bam ref=${REFERENCE_FASTA} -Xmx6g sam=1.3

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
    container "quay.io/greninger-lab/swift-pipeline:latest"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
      file SORTED_SAM  from Sorted_sam_ch_SE
      file MASTERFILE
    output:
      tuple file("*.clipped.bam"), file("*.bai")  into Clipped_bam_ch_SE
      file "*."

    script:
    """
    #!/bin/bash
    R1=`basename ${SORTED_SAM} ".sorted.sam"`
    /./root/.local/bin/primerclip -s ${MASTERFILE} ${SORTED_SAM} \${R1}.clipped.sam
    #/usr/local/miniconda/bin/samtools sort -n -O sam \${R1}.clipped.sam > \${R1}.clipped.sorted.sam
    #/usr/local/miniconda/bin/samtools view -Sb \${R1}.clipped.sorted.sam > \${R1}.clipped.unsorted.bam
    #/usr/local/miniconda/bin/samtools sort -o \${R1}.clipped.unsorted.bam \${R1}.clipped.bam
     /usr/local/miniconda/bin/samtools sort \${R1}.clipped.sam -o \${R1}.clipped.bam
     /usr/local/miniconda/bin/samtools index \${R1}.clipped.bam

    """
}

process generateConsensus_SE {
    container "quay.io/greninger-lab/swift-pipeline:latest"
    echo true

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple file(BAMFILE), file("*.bai") from Clipped_bam_ch_SE
        file REFERENCE_FASTA
        file TRIM_ENDS
    output:
        file("*_swift.fasta")
        file(BAMFILE)
        file("*.bai")

    publishDir params.OUTDIR, mode: 'copy'

    shell:
    '''
    #!/bin/bash

    R1=`basename !{BAMFILE} .clipped.bam`

    # Using LAVA consensus calling!
    /usr/local/miniconda/bin/bcftools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 500000 \\
        --max-idepth 500000 \\
        --fasta-ref !{REFERENCE_FASTA} \\
        --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
        --threads 10 \\
        !{BAMFILE} \\
        | /usr/local/miniconda/bin/bcftools call -m -Oz -o \${R1}_pre.vcf.gz
    
    /usr/local/miniconda/bin/tabix \${R1}_pre.vcf.gz
    gunzip \${R1}_pre.vcf.gz
     /usr/local/miniconda/bin/bcftools filter -i '(DP4[0]+DP4[1]) < (DP4[2]+DP4[3]) && ((DP4[2]+DP4[3]) > 0)' \${R1}_pre.vcf -o \${R1}.vcf
    /usr/local/miniconda/bin/bgzip \${R1}.vcf
    /usr/local/miniconda/bin/tabix \${R1}.vcf.gz 
    cat !{REFERENCE_FASTA} | /usr/local/miniconda/bin/bcftools consensus \${R1}.vcf.gz > \${R1}.consensus.fa

    /usr/local/miniconda/bin/bedtools genomecov \\
        -bga \\
        -ibam !{BAMFILE} \\
        -g !{REFERENCE_FASTA} \\
        | awk '\$4 < 10' | /usr/local/miniconda/bin/bedtools merge > \${R1}.mask.bed
    
    /usr/local/miniconda/bin/bedtools maskfasta \\
        -fi \${R1}.consensus.fa \\
        -bed \${R1}.mask.bed \\
        -fo \${R1}.consensus.masked.fa

    cat !{REFERENCE_FASTA} \${R1}.consensus.masked.fa > align_input.fasta
    /usr/local/miniconda/bin/mafft --auto align_input.fasta > repositioned.fasta
    awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' repositioned.fasta > repositioned_unwrap.fasta
    
    python3 !{TRIM_ENDS} \${R1}

    [ -s \${R1}_swift.fasta ] || echo "WARNING: \${R1} produced blank output. Manual review may be needed."
        
    '''
    }
}

// previous consensus calling, misses non-clonal variants
// /usr/local/miniconda/bin/bcftools mpileup \\
//         --count-orphans \\
//         --no-BAQ \\
//         --max-depth 500000 \\
//         --max-idepth 500000 \\
//         --fasta-ref !{REFERENCE_FASTA} \\
//         --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
//         --threads 10 \\
//         !{BAMFILE} \\
//         | /usr/local/miniconda/bin/bcftools call --output-type v --ploidy 1 --keep-alts --keep-masked-ref --multiallelic-caller -P 0 \\
//         | /usr/local/miniconda/bin/bcftools reheader --samples sample_name.list \\
//         | /usr/local/miniconda/bin/bcftools view --output-file \${R1}_pre.vcf.gz --output-type z

//     /usr/local/miniconda/bin/bcftools norm -f !{REFERENCE_FASTA} -m +any -Oz -o \${R1}.vcf.gz \${R1}_pre.vcf.gz

//     /usr/local/miniconda/bin/tabix -p vcf -f \${R1}.vcf.gz

//     cat !{REFERENCE_FASTA} | /usr/local/miniconda/bin/bcftools consensus -H1 \${R1}.vcf.gz > \${R1}.consensus.fa





// spades that doesn't work
// process bamToFastq {
//     container "quay.io/biocontainers/samtools:1.3--h0592bc0_3"

//     input:
//         file(BAMFILE) from Clipped_bam_ch_SE
//     output:
//         file("preprocessed_reads.fastq") into Preprocessed_fastq_ch

//     script:
//     """
//     #!/bin/bash
//     samtools bam2fq ${BAMFILE} > preprocessed_reads.fastq
//     """

// }

// process deNovoAssembly {
//     container "quay.io/biocontainers/spades:3.14.0--h2d02072_0"

// // 	// Retry on fail at most three times 
// //     errorStrategy 'retry'
// //     maxRetries 3

//     input: 
//         file(PREPROCESSED_FASTQ) from Preprocessed_fastq_ch

//     output:
//         file("scaffolds.fasta") into de_novo_ch

//     script:
//     """
//     spades.py -s ${PREPROCESSED_FASTQ} -o ./ --isolate -t 6
//     """
// }
 
