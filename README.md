# TAYLOR (Trimming Amplicons You LOve Rapidly)
This short post-processing pipeline takes a folder with .bam files previously run with the main TAYLOR pipeline and generates consensuses following CDC guidelines (IUPAC calling at 25% minor allele frequency, 50X coverage). 

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