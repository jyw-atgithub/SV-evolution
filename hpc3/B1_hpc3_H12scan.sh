#!/bin/bash

#in interactive mode
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
Merged_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/Merged_SNP"
H12_scan="/dfs7/jje/jenyuw/SV-project-temp/result/H12_scan"

#transform the SV vcf files into the required format
bcftools query -r "2L,2R" -f '%POS[ %GT]\n' ${polarizing}/3corrected.polarized.asm.vcf.gz|\
sed 's@1/1@A@g' |sed -e "s/\.\/\./T/g"|tr " " "," > ${H12_scan}/chr2.SV.csv

##It is written in python2, not python3!!
module load python/2.7.17
module load R/4.3.3
python ${H12_scan}/scripts/H12_H2H1.py ${H12_scan}/chr2.SV.csv 62 -o ${H12_scan}/chr2.SV.scan.txt -w 5
python ${H12_scan}/scripts/H12peakFinder.py ${H12_scan}/chr2.SV.scan.txt -t "-1" -o ${H12_scan}/chr2.SV.peak.txt
Rscript ${H12_scan}/scripts/H12_viz.R ${H12_scan}/chr2.SV.scan.txt ${H12_scan}/chr2.SV.peak.txt chr2.SV.pdf 2


#transform the SNP vcf files into the required format
#Sed and grep use basic regular expressions, in which only a few characters ( . * ^ $ [ ] ) have special meaning. If you want to use extended regex operators ( | () {} ), you need to precede them with a backslash.
bcftools view -m2 -M2 -v snps -r "2L,2R" ${Merged_SNP}/all.snps.vcf.gz|\
bcftools query -f '%POS[ %IUPACGT]\n' |sed -e "s/\.\/\./T/g"|\
sed -e "s/R\|Y\|S\|W\|K\|M/./g"|tr " " "," > ${H12_scan}/chr2.snp.csv

python ${H12_scan}/scripts/H12_H2H1.py ${H12_scan}/chr2.snp.csv 62 -o ${H12_scan}/chr2.SNP.scan.txt
python ${H12_scan}/scripts/H12peakFinder.py ${H12_scan}/chr2.SNP.scan.txt -t 0.02 -o ${H12_scan}/chr2.SNP.peak.txt
Rscript ${H12_scan}/scripts/H12_viz.R ${H12_scan}/chr2.SNP.scan.txt ${H12_scan}/chr2.SNP.peak.txt chr2.SNP.pdf 10