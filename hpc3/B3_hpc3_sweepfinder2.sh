#!/bin/bash

#in interactive mode
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
Merged_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/Merged_SNP"
processed_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/processed_SNP"
SF2="/dfs7/jje/jenyuw/SV-project-temp/result/sweepfinder2"

nT=$SLURM_CPUS_PER_TASK
############Working on SNPs first############
#add the number of alternate alleles and number of samples (NS; is AC divided by 2) and then create a new copy
bcftools +fill-tags ${Merged_SNP}/all.snps.vcf.gz -O z -o ${SF2}/filled.snp.vcf.gz -- -t AC,NS
bcftools index --threads ${nT} -f -t ${SF2}/filled.snp.vcf.gz

#transform the vcf files into the required format
#the SNPs are folded now, so the values at folded colum is 1
sample_n=42 ## 21 samples in AM population, 21*2=42 haplotypes. 
echo -e "position\tx\tn\tfolded" > ${SF2}/AM.snp.background.txt
bcftools view --threads ${nT} -m2 -M2 -v snps -r "2L,2R,3L,3R,4,X,Y" -S ${SF2}/AM_pop-snp.txt ${SF2}/filled.snp.vcf.gz|\
bcftools query -f "%POS\t%AC\t${sample_n}\t1\n" >> ${SF2}/AM.snp.background.txt
#Compute empirical frequency spectrum (-f)
SweepFinder2 -f ${SF2}/AM.snp.background.txt ${SF2}/AM.snp.spect.txt

##Create  empirical frequency spectrum of synonymous sites
#bcftools +fill-tags ${processed_SNP}/synSNPs.vcf.gz -O z -o ${SF2}/filled.synSNPs.vcf.gz -- -t AC,NS
#bcftools index --threads ${nT} -f -t ${SF2}/filled.synSNPs.vcf.gz
#echo -e "position\tx\tn\tfolded" > ${SF2}/snp.syn.txt
#bcftools view --threads ${nT} -m2 -M2 -v snps -r "2L,2R,3L,3R,4,X,Y" ${SF2}/filled.synSNPs.vcf.gz|\
#bcftools query -f '%POS\t%AC\t124\t1\n' >> ${SF2}/snp.syn.txt
#SweepFinder2 -f ${SF2}/snp.syn.txt ${SF2}/snp.syn.spect.txt

##Create another empirical frequency spectrum with 31 samples and only synonymous sites. For SV analysis. 
echo -e "position\tx\tn\tfolded" > ${SF2}/31smaples.syn.txt
bcftools view --threads ${nT} -m2 -M2 -v snps -r "2L,2R,3L,3R,4,X,Y"  -S ${SF2}/sample_subset.txt ${SF2}/filled.snp.vcf.gz|\
bcftools query -f '%POS\t%AC\t62\t1\n'|gawk '$2>0' >> ${SF2}/31smaples.syn.txt
SweepFinder2 -f ${SF2}/31smaples.syn.txt ${SF2}/31smaples.syn.spect.txt

#extract the SNPs in chromosome 2R, AM population
echo -e "position\tx\tn\tfolded" > ${SF2}/2R.AM.snp.input.txt
bcftools view --threads ${nT} -m2 -M2 -v snps -r "2R" -S ${SF2}/AM_pop-snp.txt ${SF2}/filled.snp.vcf.gz|\
bcftools query -f "%POS\t%AC\t${sample_n}\t1\n" >> ${SF2}/2R.AM.snp.input.txt

#Actual sweep scan
#1000 means devide the chromosome into 1000 windows
SweepFinder2 -l 1000 ${SF2}/2R.AM.snp.input.txt ${SF2}/AM.snp.spect.txt ${SF2}/2R.AM.snp.output.txt
#SweepFinder2 -l 1000 ${SF2}/2R.AM.snp.input.txt ${SF2}/snp.syn.spect.txt ${SF2}/2R.AM.snp.output-2.txt

############        Working on SV now      ############
#add the number of alternate alleles and number of samples (NS; is AC divided by 2) and then create a new copy
bcftools +fill-tags ${polarizing}/3corrected.polarized.asm.vcf.gz -O z -o ${SF2}/filled.SV.asm.vcf.gz -- -t AC,NS
bcftools index -f -t ${SF2}/filled.SV.asm.vcf.gz
#extract the SV in chromosome 2R
#transform the SV vcf files into the required format
#Actually NS is what we want for "the allele count (x)"
#now, the sample number (n) is 62
#everything is unfolded, so all the values at folded colum is 0
#remove the rows with x=0
echo -e "position\tx\tn\tfolded" > ${SF2}/2R.SV.asm.input.txt
bcftools query -r "2R" -f '%POS\t%NS\t62\t0\n' ${SF2}/filled.SV.asm.vcf.gz|\
gawk '$2>0' >> ${SF2}/2R.SV.asm.input.txt

echo -e "position\tx\tn\tfolded" > ${SF2}/SV.background.txt
bcftools view -m2 -M2 -r "2L,2R,3L,3R,4,X,Y" ${SF2}/filled.SV.asm.vcf.gz|\
#Actually NS is what we want for "the allele count (x)"
bcftools query -f '%POS\t%NS\t62\t0\n'|\
gawk '$2>0' >> ${SF2}/SV.background.txt
#Compute empirical frequency spectrum (-f)
SweepFinder2 -f ${SF2}/SV.background.txt ${SF2}/SV.spect.txt
#Actual sweep scan
SweepFinder2 -l 1000 ${SF2}/2R.SV.asm.input.txt ${SF2}/SV.spect.txt ${SF2}/2R.SV.asm.output.txt
SweepFinder2 -l 1000 ${SF2}/2R.SV.asm.input.txt ${SF2}/31smaples.syn.spect.txt ${SF2}/2R.SV.asm.output-2.txt