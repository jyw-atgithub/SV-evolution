#!/bin/bash

#in interactive mode
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
Merged_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/Merged_SNP"
H12_scan="/dfs7/jje/jenyuw/SV-project-temp/result/H12_scan"
nT=$SLURM_CPUS_PER_TASK
#transform the SV vcf files into the required format
##All samples
bcftools query -r "2R" -f '%POS[%GT]\n' ${polarizing}/3corrected.polarized.asm.vcf.gz|\
sed 's@1/1@A@g' |sed -e "s/\.\/\./T/g"|tr " " "," > ${H12_scan}/2R.SV.csv

##It is written in python2, not python3!!
module load python/2.7.17
module load R/4.3.3
python ${H12_scan}/scripts/H12_H2H1.py ${H12_scan}/2R.SV.csv 62 -o ${H12_scan}/2R.SV.scan.txt -w 5
python ${H12_scan}/scripts/H12peakFinder.py ${H12_scan}/2R.SV.scan.txt -t "-1" -o ${H12_scan}/2R.SV.peak.txt
Rscript ${H12_scan}/scripts/H12_viz.R ${H12_scan}/2R.SV.scan.txt ${H12_scan}/2R.SV.peak.txt 2R.SV.pdf 5
## It turned out H12 scan does not work on SVs. 

#transform the SNP vcf files into the required format
#Sed and grep use basic regular expressions, in which only a few characters ( . * ^ $ [ ] ) have special meaning. If you want to use extended regex operators ( | () {} ), you need to precede them with a backslash.
####All samples
bcftools view --threads ${nT} -m2 -M2 -v snps -r "2R" ${Merged_SNP}/all.snps.vcf.gz| sed -e 's/\.\/\./0\/0/g'|\
bcftools query -f '%POS[ %IUPACGT]\n' |\
sed -e "s/R\|Y\|S\|W\|K\|M/./g"|tr " " "," > ${H12_scan}/2R.snp.csv

#### Only AM samples, n=21
for i in 2L 2R 3L 3R X Y
do
bcftools view --threads ${nT} -m2 -M2 -v snps -r ${i} \
-S ${H12_scan}/AM_pop-snp.txt ${Merged_SNP}/all.snps.vcf.gz | sed -e 's/\.\/\./0\/0/g'|\
bcftools query -f '%POS[ %IUPACGT]\n' |\
sed -e "s/R\|Y\|S\|W\|K\|M/./g"|tr " " "," > ${H12_scan}/${i}.AM.snp.csv
done
#### Only EU samples, n=33
bcftools view --threads ${nT} -m2 -M2 -v snps -r "2R" \
-S ${H12_scan}/EU_pop-snp.txt ${Merged_SNP}/all.snps.vcf.gz | sed -e 's/\.\/\./0\/0/g'|\
bcftools query -f '%POS[ %IUPACGT]\n' |\
sed -e "s/R\|Y\|S\|W\|K\|M/./g"|tr " " "," > ${H12_scan}/2R.EU.snp.csv

####All samples
python ${H12_scan}/scripts/H12_H2H1.py ${H12_scan}/2R.snp.csv 62 -o ${H12_scan}/2R.SNP.scan.txt -w 400 -j 50
python ${H12_scan}/scripts/H12peakFinder.py ${H12_scan}/2R.SNP.scan.txt -t 0.01 -o ${H12_scan}/2R.SNP.peak.txt
Rscript ${H12_scan}/scripts/H12_viz.R ${H12_scan}/2R.SNP.scan.txt ${H12_scan}/2R.SNP.peak.txt 2R.SNP.pdf 10
#### Only AM samples, n=21
for i in 2L 2R 3L 3R X Y
do
python ${H12_scan}/scripts/H12_H2H1.py ${H12_scan}/${i}.AM.snp.csv 21 -o ${H12_scan}/${i}.AM.snp.scan.txt -w 250 -j 50
python ${H12_scan}/scripts/H12peakFinder.py ${H12_scan}/${i}.AM.snp.scan.txt -t 0.02 -o ${H12_scan}/${i}.AM.snp.peak.txt
Rscript ${H12_scan}/scripts/H12_viz.R ${H12_scan}/${i}.AM.snp.scan.txt ${H12_scan}/${i}.AM.snp.peak.txt ${i}.AM.snp.pdf 10
done
#### Only EU samples, n=33
python ${H12_scan}/scripts/H12_H2H1.py ${H12_scan}/2R.EU.snp.csv 33 -o ${H12_scan}/2R.EU.snp.scan.txt -w 400 -j 50
python ${H12_scan}/scripts/H12peakFinder.py ${H12_scan}/2R.EU.snp.scan.txt -t 0.02 -o ${H12_scan}/2R.EU.snp.peak.txt
Rscript ${H12_scan}/scripts/H12_viz.R ${H12_scan}/2R.EU.snp.scan.txt ${H12_scan}/2R.EU.snp.peak.txt 2R.EU.snp.pdf 10