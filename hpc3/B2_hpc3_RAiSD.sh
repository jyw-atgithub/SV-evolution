#!/bin/bash

#in interactive mode
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
Merged_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/Merged_SNP"
raisd="/dfs7/jje/jenyuw/SV-project-temp/result/RAiSD"

##The RAiSD program is installed locally

bcftools view -O v ${polarizing}/3corrected.polarized.asm.vcf.gz >${raisd}/SV.asm.vcf

RAiSD -f -n test -I SV.asm.vcf ##not working now