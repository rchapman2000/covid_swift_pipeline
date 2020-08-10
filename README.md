# SARS-CoV-2 Swift Pipeline
This pipeline takes gzipped fastq files and outputs .bam files aligned to NC_045512 as well as consensus fastas. Notably this pipeline incorporates [primerclip](https://github.com/swiftbiosciences/primerclip/tree/deltest).

## Installation

1. Install [nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation).
   - Make sure you move nextflow to a directory in your PATH variable.
2. Install [docker](https://docs.docker.com/get-docker/).

## Usage
- Example command for fastqs in current directory: ```nextflow run michellejlin/covid_swift_pipeline --INPUT ./ --OUTDIR output/ -resume -with-trace```


| Command  | Description |
| ---      | ---         | 
| --INPUT  | Input folder where gzipped fastqs are located. For current  directory, `./` can be used.
| --OUTDIR | Output folder where .bams and consensus fastas will be piped into.
| -resume  | nextflow will pick up where it left off if the previous command was interrupted for some reason.
| -with-docker ubuntu:18.04 | Runs command with Ubuntu docker.
| -with-trace | Outputs a trace.txt that shows which processes end up in which work/ folders. 
