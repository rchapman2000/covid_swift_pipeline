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
        --PRIMERS       Primer masterfile to run. By default, this pipeline uses the original Swift V1 primers. --PRIMERS V2 can be specified for the Swift V2 primer set.

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

//files 
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

    publishDir "${params.OUTDIR}trimmed_fastqs", mode: 'copy',pattern:'*fastq*'

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
    
    echo Sample_Name,Raw_Reads,Trimmed_Paired_Reads,Trimmed_Unpaired_Reads,Total_Trimmed_Reads,Percent_Trimmed,Mapped_Reads,Clipped_Mapped_Reads,Mean_Coverage,Spike_Mean_Coverage,Spike_100X_Cov_Percentage,Spike_200X_Cov_Percentage,Percent_N > ${base}_summary.csv
    printf "${base},\$num_untrimmed,\$num_paired,\$num_unpaired,\$num_trimmed,\$percent_trimmed" >> ${base}_summary.csv
    
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
      tuple val(base), file("${base}.bam"),file("${base}_summary2.csv") into Aligned_bam_ch
      tuple val (base), file("*") into Dump_ch

    script:
    """
    #!/bin/bash

    cat ${base}*.fastq.gz > ${base}_cat.fastq.gz
    /usr/local/bin/bbmap.sh in=${base}_cat.fastq.gz outm=${base}.bam ref=${REFERENCE_FASTA} -Xmx6g > bbmap_out.txt 2>&1
    reads_mapped=\$(cat bbmap_out.txt | grep "mapped:" | cut -d\$'\\t' -f3)

    cp ${base}_summary.csv ${base}_summary2.csv
    printf ",\$reads_mapped" >> ${base}_summary2.csv

    """

    // bwa?
    // script:
    // """
    // #!/bin/bash
    // bwa mem $workflow.projectDir/NC_045512.2.fasta ${base}.R1.paired.fastq.gz ${base}.R2.paired.fastq.gz > aln.sam
    // """
}

// Single end first three steps
} else { 
    input_read_ch = Channel
        .fromPath("${params.INPUT}*.gz")
        //.map { it -> [ file(it)]}
        .map { it -> file(it)}
    
process Trimming_SE { 
    container "quay.io/biocontainers/trimmomatic:0.35--6"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
      file R1 from input_read_ch
      file ADAPTERS
    output: 
      tuple env(base),file("*.trimmed.fastq.gz"),file("*summary.csv") into Trim_out_ch_SE
      tuple env(base),file("*.trimmed.fastq.gz") into Trim_out_ch2_SE

    publishDir "${params.OUTDIR}trimmed_fastqs", mode: 'copy',pattern:'*fastq*'

    script:
    """
    #!/bin/bash

    base=`basename ${R1} ".fastq.gz"`

    echo \$base

    trimmomatic SE -threads ${task.cpus} ${R1} \$base.trimmed.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

    num_untrimmed=\$((\$(gunzip -c ${R1} | wc -l)/4))
    num_trimmed=\$((\$(gunzip -c \$base'.trimmed.fastq.gz' | wc -l)/4))
    
    percent_trimmed=\$((100-\$((100*num_trimmed/num_untrimmed))))
    
    echo Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Mapped_Reads,Clipped_Mapped_Reads,Mean_Coverage,Spike_Mean_Coverage,Spike_100X_Cov_Percentage,Spike_200X_Cov_Percentage,Percent_N > \$base'_summary.csv'
    printf "\$base,\$num_untrimmed,\$num_trimmed,\$percent_trimmed" >> \$base'_summary.csv'
    
    ls -latr
    """
}

process fastQc_SE {
    container "quay.io/biocontainers/fastqc:0.11.9--0"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
      tuple val(base),file("${base}.trimmed.fastq.gz") from Trim_out_ch2_SE
    output: 
      file("*fastqc*") into Fastqc_ch 

    publishDir "${params.OUTDIR}fastqc", mode: 'copy'

    script:
    """
    #!/bin/bash

    /usr/local/bin/fastqc ${base}.trimmed.fastq.gz
    """
}

process Aligning_SE {
     container "quay.io/biocontainers/bbmap:38.86--h1296035_0"
    //container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7	"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input: 
      tuple val(base),file("${base}.trimmed.fastq.gz"),file("${base}_summary.csv") from Trim_out_ch_SE
      file REFERENCE_FASTA
    output:
      tuple val (base), file("${base}.bam"),file("${base}_summary2.csv") into Aligned_bam_ch

    cpus 4 
    memory '6 GB'

    script:
    """
    #!/bin/bash

    base=`basename ${base}.trimmed.fastq.gz ".trimmed.fastq.gz"`
    /usr/local/bin/bbmap.sh in1="\$base".trimmed.fastq.gz  outm="\$base".bam ref=${REFERENCE_FASTA} -Xmx6g sam=1.3 > bbmap_out.txt 2>&1
    reads_mapped=\$(cat bbmap_out.txt | grep "mapped:" | cut -d\$'\\t' -f3)
    
    cp ${base}_summary.csv ${base}_summary2.csv
    printf ",\$reads_mapped" >> ${base}_summary2.csv

    """
}
}

process NameSorting { 
    container "quay.io/biocontainers/samtools:1.3--h0592bc0_3"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
      tuple val (base), file("${base}.bam"),file("${base}_summary2.csv") from Aligned_bam_ch
    output:
      tuple val (base), file("${base}.sorted.sam"),file("${base}_summary2.csv") into Sorted_sam_ch
    
    publishDir "${params.OUTDIR}inprogress_summary", mode: 'copy', pattern: '*summary.csv'

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
      tuple val (base), file("${base}.sorted.sam"),file("${base}_summary2.csv") from Sorted_sam_ch
      file MASTERFILE
    output:
      //tuple val (base), file("${base}.clipped.bam"), file("*.bai"),file("${base}_summary3.csv"),env(bamsize) into Clipped_bam_ch
      tuple val (base), file("${base}.clipped.bam"), file("${base}.clipped.bam.bai"),file("${base}_summary3.csv"),env(bamsize) into Clipped_bam_ch2

    publishDir params.OUTDIR, mode: 'copy', pattern: '*.clipped.bam'
    publishDir "${params.OUTDIR}inprogress_summary", mode: 'copy', pattern: '*summary3.csv'

    script:
        """
        #!/bin/bash

        ls -latr
        /./root/.local/bin/primerclip -s ${MASTERFILE} ${base}.sorted.sam ${base}.clipped.sam
        #/usr/local/miniconda/bin/samtools sort -@ ${task.cpus} -n -O sam ${base}.clipped.sam > ${base}.clipped.sorted.sam
        #/usr/local/miniconda/bin/samtools view -@ ${task.cpus} -Sb ${base}.clipped.sorted.sam > ${base}.clipped.unsorted.bam
        #/usr/local/miniconda/bin/samtools sort -@ ${task.cpus} -o ${base}.clipped.unsorted.bam ${base}.clipped.bam
        /usr/local/miniconda/bin/samtools sort -@ ${task.cpus} ${base}.clipped.sam -o ${base}.clipped.bam
        /usr/local/miniconda/bin/samtools index ${base}.clipped.bam

        clipped_reads=\$(/usr/local/miniconda/bin/samtools flagstat ${base}.clipped.bam | grep "mapped (" | awk '{print \$1}')
        echo "clipped reads: \$clipped_reads"

        meancoverage=\$(/usr/local/miniconda/bin/samtools depth -m 0 -a ${base}.clipped.bam | awk '{sum+=\$3} END { print sum/NR}')

        bamsize=\$((\$(wc -c ${base}.clipped.bam | awk '{print \$1'})+0))
        echo "bamsize: \$bamsize"

        if (( \$bamsize > 92 ))
        then
            # Spike protein coverage
            /usr/local/miniconda/bin/samtools depth -a -r NC_045512.2:21563-25384 -m 0 ${base}.clipped.bam > ${base}_spike_coverage.txt
            avgcoverage=\$(cat ${base}_spike_coverage.txt | awk '{sum+=\$3} END { print sum/NR}')
            proteinlength=\$((25384-21563+1))
            cov100=\$((100*\$(cat ${base}_spike_coverage.txt | awk '\$3>=100' | wc -l)/3822))
            cov200=\$((100*\$(cat ${base}_spike_coverage.txt | awk '\$3>=200' | wc -l)/3822))

        else
            avgcoverage=0
            cov100=0
            cov200=0
        fi
        
        cp ${base}_summary2.csv ${base}_summary3.csv
        printf ",\$clipped_reads,\$meancoverage,\$avgcoverage,\$cov100,\$cov200" >> ${base}_summary3.csv

        """
    
}

process lofreq {
    container "quay.io/biocontainers/lofreq:2.1.5--py38h1bd3507_3"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
      tuple val (base), file("${base}.clipped.bam"), file("${base}.clipped.bam.bai"),file("${base}_summary3.csv"),val(bamsize) from Clipped_bam_ch2
      file REFERENCE_FASTA
    output:
      tuple val (base), file("${base}.clipped.bam"), file("${base}.clipped.bam.bai"),file("${base}_summary3.csv"),val(bamsize),file("${base}_lofreq.vcf") into Lofreq_ch
    
    script:
    """
    #!/bin/bash

    echo ${bamsize}
    if (( ${bamsize} > 92))
    then
        lofreq faidx ${REFERENCE_FASTA}
        /usr/local/bin/lofreq call-parallel --pp-threads ${task.cpus} -f ${REFERENCE_FASTA} -o ${base}_lofreq.vcf ${base}.clipped.bam
    else
        touch ${base}_lofreq.vcf
    fi

    """
}

process generateConsensus {
    container "quay.io/greninger-lab/swift-pipeline:latest"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val (base), file("${base}.clipped.bam"), file("${base}.clipped.bam.bai"),file("${base}_summary3.csv"),val(bamsize),file("${base}_lofreq.vcf") from Lofreq_ch
        file REFERENCE_FASTA
        file TRIM_ENDS
        file FIX_COVERAGE
        file VCFUTILS
        file REFERENCE_FASTA_FAI
        file SPLITCHR
    output:
        file("${base}_swift.fasta")
        file("${base}_lofreq.vcf")
        file("${base}.vcf")
        file(INDEX_FILE)
        file("${base}_summary.csv")
        file("${base}.pileup")

    publishDir params.OUTDIR, mode: 'copy'

    shell:
    '''
    #!/bin/bash
    ls -latr

    R1=`basename !{BAMFILE} .clipped.bam`

    echo "bamsize: !{bamsize}"

    #if [ -s !{BAMFILE} ]
    # More reliable way of checking bam size, because of aliases
    if (( !{bamsize} > 92 ))
    then
        /usr/local/miniconda/bin/bcftools filter -i 'AF>=0.5' --threads !{task.cpus} \${R1}_lofreq.vcf -o \${R1}.vcf
        /usr/local/miniconda/bin/bgzip \${R1}.vcf
        /usr/local/miniconda/bin/tabix \${R1}.vcf.gz
        gunzip \${R1}.vcf.gz
        
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
        /usr/local/miniconda/bin/mafft --auto --thread !{task.cpus} align_input.fasta > repositioned.fasta
        awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' repositioned.fasta > repositioned_unwrap.fasta
        
        python3 !{TRIM_ENDS} \${R1}

        # Find percent ns, doesn't work, fix later in python script
        num_bases=$(grep -v ">" \${R1}_swift.fasta | wc | awk '{print $3-$1}')
        num_ns=$(grep -v ">" \${R1}_swift.fasta | awk -F"n" '{print NF-1}')
        percent_n=$(($(($num_ns/$num_bases))*100))

        echo "num_bases=$num_bases"
        echo "num_ns=$num_ns"
        echo "percent_n=$percent_n"

        gunzip \${R1}.vcf.gz

    else
       echo "Empty bam detected. Generating empty consensus fasta file..."
       printf '>!{base}\n' > \${R1}_swift.fasta
       printf 'n%.0s' {1..29539} >> \${R1}_swift.fasta
       percent_n=100

       touch \${R1}.vcf
    fi
    
    cp \${R1}_summary3.csv \${R1}_summary.csv
    printf ",\$percent_n" >> \${R1}_summary.csv

    cat \${R1}_summary.csv | tr -d "[:blank:]" > a.tmp
    mv a.tmp \${R1}_summary.csv

    if [[ !{bamsize} > 92 ]]
    then
        python3 !{FIX_COVERAGE} \${R1}
        mv \${R1}_summary_fixed.csv \${R1}_summary.csv
    fi

    [ -s \${R1}_swift.fasta ] || echo "WARNING: \${R1} produced blank output. Manual review may be needed."

    '''
}




// // Deprecated bcftools consensus generation
// process generateConsensus {
//     container "quay.io/greninger-lab/swift-pipeline:latest"

// 	// Retry on fail at most three times 
//     errorStrategy 'retry'
//     maxRetries 3

//     input:
//         tuple val (base), file(BAMFILE),file(INDEX_FILE),file("${base}_summary3.csv"),val(bamsize) from Clipped_bam_ch
//         file REFERENCE_FASTA
//         file TRIM_ENDS
//         file FIX_COVERAGE
//         file VCFUTILS
//         file REFERENCE_FASTA_FAI
//         file SPLITCHR
//     output:
//         file("${base}_swift.fasta")
//         file("${base}_bcftools.vcf")
//         file(INDEX_FILE)
//         file("${base}_summary.csv")
//         file("${base}.pileup")

//     publishDir params.OUTDIR, mode: 'copy'

//     shell:
//     '''
//     #!/bin/bash
//     ls -latr

//     R1=`basename !{BAMFILE} .clipped.bam`

//     echo "bamsize: !{bamsize}"

//     #if [ -s !{BAMFILE} ]
//     # More reliable way of checking bam size, because of aliases
//     if (( !{bamsize} > 92 ))
//     then
//         # Parallelize pileup based on number of cores
//         splitnum=$(($((29903/!{task.cpus}))+1))
//         #perl !{VCFUTILS} splitchr -l $splitnum !{REFERENCE_FASTA_FAI} | \\
//         cat !{SPLITCHR} | \\
//             xargs -I {} -n 1 -P !{task.cpus} sh -c \\
//                 "/usr/local/miniconda/bin/bcftools mpileup \\
//                     -f !{REFERENCE_FASTA} -r {} \\
//                     --count-orphans \\
//                     --no-BAQ \\
//                     --max-depth 50000 \\
//                     --max-idepth 500000 \\
//                     --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
//                 !{BAMFILE} > \${R1}.pileup
        
//         cat \${R1}.pileup | /usr/local/miniconda/bin/bcftools call -m -Oz --ploidy 1 - > tmp.{}.vcf.gz"
        
//         cat *.vcf.gz > \${R1}_catted.vcf.gz
//         /usr/local/miniconda/bin/tabix \${R1}_catted.vcf.gz
//         gunzip \${R1}_catted.vcf.gz
//         cat \${R1}_catted.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > \${R1}_pre.vcf
        
//         /usr/local/miniconda/bin/bcftools filter -i '(DP4[0]+DP4[1]) < (DP4[2]+DP4[3]) && ((DP4[2]+DP4[3]) > 0)' --threads !{task.cpus} \${R1}_pre.vcf -o \${R1}.vcf
//         /usr/local/miniconda/bin/bgzip \${R1}.vcf
//         /usr/local/miniconda/bin/tabix \${R1}.vcf.gz 
//         cat !{REFERENCE_FASTA} | /usr/local/miniconda/bin/bcftools consensus \${R1}.vcf.gz > \${R1}.consensus.fa

//         /usr/local/miniconda/bin/bedtools genomecov \\
//             -bga \\
//             -ibam !{BAMFILE} \\
//             -g !{REFERENCE_FASTA} \\
//             | awk '\$4 < 10' | /usr/local/miniconda/bin/bedtools merge > \${R1}.mask.bed
//         /usr/local/miniconda/bin/bedtools maskfasta \\
//         -fi \${R1}.consensus.fa \\
//         -bed \${R1}.mask.bed \\
//         -fo \${R1}.consensus.masked.fa

//         cat !{REFERENCE_FASTA} \${R1}.consensus.masked.fa > align_input.fasta
//         /usr/local/miniconda/bin/mafft --auto --thread !{task.cpus} align_input.fasta > repositioned.fasta
//         awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' repositioned.fasta > repositioned_unwrap.fasta
        
//         python3 !{TRIM_ENDS} \${R1}

//         # Find percent ns, doesn't work, fix later in python script
//         num_bases=$(grep -v ">" \${R1}_swift.fasta | wc | awk '{print $3-$1}')
//         num_ns=$(grep -v ">" \${R1}_swift.fasta | awk -F"n" '{print NF-1}')
//         percent_n=$(($(($num_ns/$num_bases))*100))

//         echo "num_bases=$num_bases"
//         echo "num_ns=$num_ns"
//         echo "percent_n=$percent_n"

//         gunzip \${R1}.vcf.gz
//         mv \${R1}.vcf \${R1}_bcftools.vcf

//     else
//        echo "Empty bam detected. Generating empty consensus fasta file..."
//        printf '>!{base}\n' > \${R1}_swift.fasta
//        printf 'n%.0s' {1..29539} >> \${R1}_swift.fasta
//        percent_n=100

//        touch \${R1}_bcftools.vcf
//     fi
    
//     cp \${R1}_summary3.csv \${R1}_summary.csv
//     printf ",\$percent_n" >> \${R1}_summary.csv

//     cat \${R1}_summary.csv | tr -d "[:blank:]" > a.tmp
//     mv a.tmp \${R1}_summary.csv

//     if [[ !{bamsize} > 92 ]]
//     then
//         python3 !{FIX_COVERAGE} \${R1}
//         mv \${R1}_summary_fixed.csv \${R1}_summary.csv
//     fi

//     [ -s \${R1}_swift.fasta ] || echo "WARNING: \${R1} produced blank output. Manual review may be needed."

//     '''
// }