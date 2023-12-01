#!/bin/bash

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
SNPs="/dfs7/jje/jenyuw/SV-project-temp/result/SNPs"
Merged_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/Merged_SNP"

## in an interactive mode


# Avoid using --missing-to-ref, because this does look like a good assumption
bcftools merge --threads 8 ${SNPs}/*.filtered.snps.vcf.gz -O z > ${Merged_SNP}/all.snps.vcf.gz
tabix -p vcf ${Merged_SNP}/all.snps.vcf.gz

conda activate everything
snpEff -v BDGP6.32.105 ${Merged_SNP}/all.snps.vcf.gz | bgzip -@ 8 -c  >${Merged_SNP}/all.snps.annotated.vcf.gz
conda deactivate

