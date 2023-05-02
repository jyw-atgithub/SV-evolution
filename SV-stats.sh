#!/bin/bash

## vk tajima & vcftools --TajimaD only works on SNPs
## vcftools --TajimaD Expected at least 2 parts in FORMAT entry: ID=DR,Number=2,Type=Integer

vcftools --vcf /home/jenyuw/SV-project/result/merged_SVs/all.consensus.vcf --stdout --window-pi 1000000 > all.consensus.pi


BASECALLS="/home/jenyuw/SV-project/raw/nv107_combined.fastq"
DRAFT="/home/jenyuw/SV-project/result/assemble/nv107_Flye/assembly.fasta"
OUTDIR="/home/jenyuw/SV-project/temp3"
medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t 20 -m r941_min_high_g303

virtualenv medaka --python=python3 --prompt "(medaka)"
##activate the virtial env
. medaka/bin/activate
pip install --upgrade pip
pip install medaka

## deactivate the virtialenv
deactivate