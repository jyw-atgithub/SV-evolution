#!/bin/bash

#in interactive mode
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
Merged_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/Merged_SNP"
processed_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/processed_SNP"
fst="/dfs7/jje/jenyuw/SV-project-temp/result/fst"
nT=$SLURM_CPUS_PER_TASK

#cp ${polarizing}/3corrected.polarized.asm.vcf.gz ${fst}/3corrected.polarized.asm.vcf.gz
#bgzip -d ${fst}/3corrected.polarized.asm.vcf.gz
#manually change the VCF version from 4.3 to 4.2, because v4.3 made vcftools crash. Stupid
vcftools --vcf ${fst}/3corrected.polarized.asm.vcf \
--weir-fst-pop ${fst}/AF_pop.txt \
--weir-fst-pop ${fst}/EU_pop.txt \
--out ${fst}/fst_SV_AF-EU.tsv

vcftools --vcf ${fst}/3corrected.polarized.asm.vcf \
--weir-fst-pop ${fst}/EU_pop.txt \
--weir-fst-pop ${fst}/AM_pop.txt \
--out ${fst}/fst_SV_EU-AM
########vcftools seems do NOT work on SVs########
bcftools view --threads ${nT} -m 2 -M 2 ${Merged_SNP}/all.snps.vcf.gz -O z -o ${fst}/all.snps.vcf.gz
bcftools view --threads ${nT} -m 2 -M 2 ${Merged_SNP}/all.snps.vcf.gz|sed -e 's/\.\/\./0\/0/g' >${fst}/all.snps_missingas0_.vcf
cat ${fst}/AF_pop.txt|sed 's/_mumco//' > ${fst}/AF_pop-snp.txt
cat ${fst}/AM_pop.txt|sed 's/_mumco//' > ${fst}/AM_pop-snp.txt
cat ${fst}/EU_pop.txt|sed 's/_mumco//' > ${fst}/EU_pop-snp.txt

vcftools --gzvcf ${fst}/all.snps.vcf.gz \
--weir-fst-pop ${fst}/EU_pop-snp.txt \
--weir-fst-pop ${fst}/AM_pop-snp.txt \
--fst-window-size 5000 \
--fst-window-step 5000 \
--out ${fst}/fst_snp_EU-AM

vcftools --gzvcf ${fst}/all.snps.vcf.gz \
--weir-fst-pop ${fst}/AF_pop-snp.txt \
--weir-fst-pop ${fst}/EU_pop-snp.txt \
--fst-window-size 5000 \
--fst-window-step 5000 \
--out ${fst}/fst_snp_AF-EU

vcftools --vcf ${fst}/all.snps_missingas0_.vcf \
--weir-fst-pop ${fst}/AF_pop-snp.txt \
--weir-fst-pop ${fst}/EU_pop-snp.txt \
--fst-window-size 5000 \
--fst-window-step 5000 \
--out ${fst}/fst_snp_AF-EU_2

rm *.log