#!/bin/bash

trimmomatic PE L001_R1_001.gz L001_R2_001.gz L001.R1.paired.fastq.gz L001.R1.unpaired.fastq.gz L001.R2.paired.fastq.gz L001.R2.unpaired.fastq.gz     ILLUMINACLIP:$HOME/downloads/trimmomatic-0.38/adapters/All_adapters.fa:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:75
