#!/bin/bash

# in the interactive mode. 
source ~/.bashrc

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
merged_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/merged_SVs"


#parallel -j 12 --progress bcftools index -f -t {} ::: ${SVs}/*.svimASM.filtered.vcf.gz
bcftools merge -m none ${SVs}/*.svimASM.filtered.vcf.gz -O z -o - |\
bcftools sort --max-mem 2G -O z -o - > ${merged_SVs}/merge.svimASM.sort.vcf.gz
bcftools index -f -t ${merged_SVs}/merge.svimASM.sort.vcf.gz

module load python/3.10.2

truvari collapse --sizemax 200000000 -k common \
-i ${merged_SVs}/merge.svimASM.sort.vcf.gz \
-c ${merged_SVs}/truvari.svimASM.collapse.vcf.gz -f ${ref_genome} |\
bcftools sort --max-mem 2G |\
bgzip -@ 8 > ${merged_SVs}/truvari.svimASM.vcf.gz

bcftools index -f -t ${merged_SVs}/truvari.svimASM.vcf.gz

module unload python/3.10.2
echo "This is the end!"