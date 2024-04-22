#!/bin/bash

#in interactive mode
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
Merged_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/Merged_SNP"
H12_scan="/dfs7/jje/jenyuw/SV-project-temp/result/H12_scan"

#transform the SV vcf files into the required format
bcftools query -r "2R" -f '%POS[ %GT]\n' ${polarizing}/3corrected.polarized.asm.vcf.gz|\
sed 's@1/1@A@g' |sed -e "s/\.\/\./T/g"|tr " " "," > ${H12_scan}/2R.SV.csv

##It is written in python2, not python3!!
module load python/2.7.17
module load R/4.3.3
python ${H12_scan}/scripts/H12_H2H1.py ${H12_scan}/2R.SV.csv 62 -o ${H12_scan}/2R.SV.scan.txt -w 5
python ${H12_scan}/scripts/H12peakFinder.py ${H12_scan}/2R.SV.scan.txt -t "-1" -o ${H12_scan}/2R.SV.peak.txt
Rscript ${H12_scan}/scripts/H12_viz.R ${H12_scan}/2R.SV.scan.txt ${H12_scan}/2R.SV.peak.txt 2R.SV.pdf 5


#transform the SNP vcf files into the required format
#Sed and grep use basic regular expressions, in which only a few characters ( . * ^ $ [ ] ) have special meaning. If you want to use extended regex operators ( | () {} ), you need to precede them with a backslash.
bcftools view -m2 -M2 -v snps -r "2R" ${Merged_SNP}/all.snps.vcf.gz|\
bcftools query -f '%POS[ %IUPACGT]\n' |sed -e "s/\.\/\./T/g"|\
sed -e "s/R\|Y\|S\|W\|K\|M/./g"|tr " " "," > ${H12_scan}/2R.snp.csv

python ${H12_scan}/scripts/H12_H2H1.py ${H12_scan}/2R.snp.csv 62 -o ${H12_scan}/2R.SNP.scan.txt
python ${H12_scan}/scripts/H12peakFinder.py ${H12_scan}/2R.SNP.scan.txt -t 0.02 -o ${H12_scan}/2R.SNP.peak.txt
Rscript ${H12_scan}/scripts/H12_viz.R ${H12_scan}/2R.SNP.scan.txt ${H12_scan}/2R.SNP.peak.txt 2R.SNP.pdf 10