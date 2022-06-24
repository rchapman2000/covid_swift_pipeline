# TAYLOR (Trimming Amplicons You LOve Rapidly)
This pipeline takes gzipped fastq files and outputs .bam files aligned to NC_045512.2 as well as consensus fastas. Notably this pipeline incorporates [primerclip](https://github.com/swiftbiosciences/primerclip/tree/deltest), and can handle both Swift and QiaSeq primersets. Without specifying any additional options, default input files are paired-end fastq files that cover the 116-255 bp amplicons produced from the Swift Amplicon SARS-CoV-2 Panel. `--SINGLE_END` can also be specified for single end reads. This pipeline can also be run without the primerclip option by specifying `--NO_CLIPPING` for consensus generation of non-Swift or non-QiaSeq SARS-CoV-2 samples. 

## Installation

1. Install [nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation).
   - Make sure you move nextflow to a directory in your PATH variable.
2. Install [docker](https://docs.docker.com/get-docker/).

## Usage
| Command  | Description |
| ---      | ---         | 
| --INPUT  | __(Required)__ Input folder where bams from main pipeline are output. For current  directory, `./` can be used.
| --OUTDIR | __(Required)__ Output folder where pileups and new consensuses will be output into.
| -with-docker ubuntu:18.04 | __(Required)__ Runs command with Ubuntu docker.
| -resume  | __(Recommended)__ nextflow will pick up where it left off if the previous command was interrupted for some reason.
| -with-trace | __(Recommended)__ Outputs a trace.txt that shows which processes end up in which work/ directories. 
| -with-report | __(Recommended)__ Outputs a report.html that gives basic stats and work directories for each process.
| -latest | __(Recommended)__ Pulls the most recent github version.

- Example command: ```nextflow run greninger-lab/covid_swift_pipeline -r iupac_edit --INPUT ./ --OUTDIR ./output/ -with-docker ubuntu:18.04 -resume -c ~/nextflow.speedy.config```