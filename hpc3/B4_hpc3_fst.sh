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
########vcftools seems do NOT work on SVs########
#So, we need to convert the format to SNP#
bcftools view -h ${polarizing}/3corrected.polarized.asm.vcf.gz|sed -e 's/VCFv4\.3/VCFv4\.2/' |\
sed -r 's/##contig=<ID=2[0-9].*//g; s/##ALT=<ID=.*//g' | grep "#"> ${fst}/recoded.SV.vcf
bcftools view --threads ${nT} -r "2L,2R,3L,3R,4,X,Y"  ${polarizing}/3corrected.polarized.asm.vcf.gz |\
bcftools query -f '%CHROM\t%POS\t%ID\tA\tT\t%QUAL\t%FILTER\t%INFO\t%FORMAT'|\
##change missing data to 0/0
sed -e 's/\.\/\./0\/0/g'|tr -d " "  >> ${fst}/recoded.SV.vcf
bcftools +fill-tags ${fst}/recoded.SV.vcf -O z -o - -- -t AC,NS |\
bcftools view -i 'AC>0' -O z -o ${fst}/recoded.SV.vcf.gz
bcftools index --threads ${nT} -f -t ${fst}/recoded.SV.vcf.gz
rm ${fst}/recoded.SV.vcf

vcftools --gzvcf ${fst}/recoded.SV.vcf.gz \
--weir-fst-pop ${fst}/AF_pop.txt \
--weir-fst-pop ${fst}/EU_pop.txt \
--fst-window-size 500000 \
--fst-window-step 250000 \
--out ${fst}/fst_SV_AF-EU_2

vcftools --gzvcf ${fst}/recoded.SV.vcf.gz \
--weir-fst-pop ${fst}/EU_pop.txt \
--weir-fst-pop ${fst}/AM_pop.txt \
--fst-window-size 500000 \
--fst-window-step 250000 \
--out ${fst}/fst_SV_EU-AM_2

bcftools view --threads ${nT} -m 2 -M 2 ${Merged_SNP}/all.snps.vcf.gz -O z -o ${fst}/bi.snps.vcf.gz
##change missing data to 0/0
bcftools view --threads ${nT} -m 2 -M 2 ${Merged_SNP}/all.snps.vcf.gz|sed -e 's/\.\/\./0\/0/g' >${fst}/bi.snps_missingas0.vcf
cat ${fst}/AF_pop.txt|sed 's/_mumco//' > ${fst}/AF_pop-snp.txt
cat ${fst}/AM_pop.txt|sed 's/_mumco//' > ${fst}/AM_pop-snp.txt
cat ${fst}/EU_pop.txt|sed 's/_mumco//' > ${fst}/EU_pop-snp.txt

vcftools --gzvcf ${fst}/bi.snps.vcf.gz \
--weir-fst-pop ${fst}/EU_pop-snp.txt \
--weir-fst-pop ${fst}/AM_pop-snp.txt \
--fst-window-size 10000 \
--fst-window-step 5000 \
--out ${fst}/fst_snp_EU-AM

vcftools --vcf ${fst}/bi.snps_missingas0.vcf \
--weir-fst-pop ${fst}/EU_pop-snp.txt \
--weir-fst-pop ${fst}/AM_pop-snp.txt \
--fst-window-size 10000 \
--fst-window-step 5000 \
--out ${fst}/fst_snp_EU-AM_2

vcftools --gzvcf ${fst}/bi.snps.vcf.gz \
--weir-fst-pop ${fst}/AF_pop-snp.txt \
--weir-fst-pop ${fst}/EU_pop-snp.txt \
--fst-window-size 10000 \
--fst-window-step 5000 \
--out ${fst}/fst_snp_AF-EU

vcftools --vcf ${fst}/bi.snps_missingas0.vcf \
--weir-fst-pop ${fst}/AF_pop-snp.txt \
--weir-fst-pop ${fst}/EU_pop-snp.txt \
--fst-window-size 10000 \
--fst-window-step 5000 \
--out ${fst}/fst_snp_AF-EU_2

rm *.log