#!/bin/bash

#in interactive mode
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
H12_scan="/dfs7/jje/jenyuw/SV-project-temp/result/H12_scan"
Merged_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/Merged_SNP"

bcftools query -r "2L" -f '%POS[ %GT]\n' ${polarizing}/3corrected.polarized.asm.vcf.gz|\
sed 's@1/1@A@g' |sed -e "s/\.\/\./T/g"|tr " " "," > ${H12_scan}/2L.SV.csv
##It is written in python2, not python3!!
module load python/2.7.17
module load R/4.3.3
python ${H12_scan}/scripts/H12_H2H1.py ${H12_scan}/2L.SV.csv 62 -o "Scan.2L.txt" -w 3
python ${H12_scan}/scripts/H12peakFinder.py ${H12_scan}/Scan.2L.txt -t "-1" -o "Peak.2L.txt"
Rscript ${H12_scan}/scripts/H12_viz.R ${H12_scan}/Scan.2L.txt ${H12_scan}/Peak.2L.txt out.2L.pdf 2


#Sed and grep use basic regular expressions, in which only a few characters ( . * ^ $ [ ] ) have special meaning. If you want to use extended regex operators ( | () {} ), you need to precede them with a backslash.
bcftools view -m2 -M2 -v snps -r "2L" ${Merged_SNP}/all.snps.vcf.gz|\
bcftools query -f '%POS[ %IUPACGT]\n' |sed -e "s/\.\/\./T/g"|\
sed -e "s/R\|Y\|S\|W\|K\|M/./g"|tr " " "," > ${H12_scan}/2L.snp.csv

python ${H12_scan}/scripts/H12_H2H1.py ${H12_scan}/2L.snp.csv 62 -o "Scan.2L.snp.txt"
python ${H12_scan}/scripts/H12peakFinder.py ${H12_scan}/Scan.2L.snp.txt -t 0.02 -o "Peak.2L.snp.txt"
Rscript ${H12_scan}/scripts/H12_viz.R ${H12_scan}/Scan.2L.snp.txt ${H12_scan}/Peak.2L.snp.txt out.2L.snp.pdf 10