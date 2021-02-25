# TAYLOR (Trimming Amplicons You LOve Rapidly)
This pipeline takes gzipped fastq files and outputs .bam files aligned to NC_045512.2 as well as consensus fastas. Notably this pipeline incorporates [primerclip](https://github.com/swiftbiosciences/primerclip/tree/deltest). By default input files should be paired-end fastq files to cover the 116-255 bp amplicons produced from the Swift Amplicon SARS-CoV-2 Panel. Expected format is `*.R1.paired.fastq.gz` and `*.R2.paired.fastq.gz`. `--SINGLE_END` can also be specified for single end reads.

This pipeline can also be run without the primerclip option by specifying `--NO_CLIPPING` for consensus generation of non-Swift SARS-CoV-2 samples. 

## Installation

1. Install [nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation).
   - Make sure you move nextflow to a directory in your PATH variable.
2. Install [docker](https://docs.docker.com/get-docker/).

## Usage
| Command  | Description |
| ---      | ---         | 
| --INPUT  | __(Required)__ Input folder where gzipped fastqs are located. For current  directory, `./` can be used.
| --OUTDIR | __(Required)__ Output folder where .bams and consensus fastas will be piped into.
| --SINGLE_END | (Optional) Flag to indicate input reads are single end. By default this pipeline expects paired end reads.
| --NO_CLIPPING | (Optional) Skip primerclip option for shotgun/CovidSeq runs.
| --SGRNA_COUNT | (Optional) Add extra step to count sgRNAs.
| --VARIANTS | (Optional) Specify which Swift primerset to use. Default: v1. 
| -with-docker ubuntu:18.04 | __(Required)__ Runs command with Ubuntu docker.
| -resume  | __(Recommended)__ nextflow will pick up where it left off if the previous command was interrupted for some reason.
| -with-trace | __(Recommended)__ Outputs a trace.txt that shows which processes end up in which work/ directories. 

Example paired fastqs are provided in the example/ folder. These can be run with the command:
- Example command for example fastqs: ```nextflow run greninger-lab/covid_swift_pipeline --INPUT example/ --OUTDIR output/ --PAIRED_END -resume -with-trace with-docker ubuntu:16.04```
