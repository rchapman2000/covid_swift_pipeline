#!/bin/bash
/./root/.local/bin/primerclip sarscov2_masterfile.txt L001.sorted.sam L001.clipped.sam
samtools view -S -b L001.clipped.sam > L001.clipped.bam
