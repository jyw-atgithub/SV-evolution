#!/bin/bash
conda activate everything
cd /home/jenyuw/DSPR_snp/snp_processing

basepath=/home/jenyuw/DSPR_snp/snp_processing
raw="/home/jenyuw/Reference_genome"
together=/home/jenyuw/DSPR_snp/raw_fq/together
processed=$basepath/processed
results=/home/jenyuw/DSPR_snp/snp_processing/results
mergedSNP=/home/jenyuw/DSPR_snp/raw_fq/together/parr_out.annotated.vcf

#getting the Tajima's D in sliding windows
## with VCFtools. Many "nan" values were generated. non-overlaping windows
#vcftools --gzvcf dmel-all-r6.46.digSNPs.vcf.gz --out tajimasd --TajimaD 10000
## with vcf-kit. non-overlaping windows
#vk tajima 10000 10000 dmel-all-r6.46.digSNPs.vcf.gz > tajimasd.vk

for i in $processed/*.vcf.gz
do
#echo $i
tempname=$(basename ${i})
echo $(basename ${i})
echo $tempname
tempname=$(sed 's/dmel-all-r6.46.//g' <(echo $tempname))
tempname=$(sed 's/SNPs.vcf.gz//g' <(echo $tempname))
name=$(sed 's/.annotated//g' <(echo $tempname))
echo $name

vk tajima 100000 100000 $i| gawk -v type=$name ' BEGIN {print "CHROM" "\t" "BIN_START" "\t" "BIN_END" "\t" "N_Sites" "\t" "N_SNPs" "\t" "TajimaD"  "\t"  "TYPE"} \
NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" type}' \
> $results/$(basename ${i}).tajimasd2.vk
done

| gawk -v type=$name ' BEGIN {print "CHROM" "\t" "BIN_START" "\t" "BIN_END" "\t" "N_Sites" "\t" "N_SNPs" "\t" "TajimaD"  "\t"  "TYPE"} \
NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" type}' \
> $results/$(basename ${i}).tajimasd2.vk

| gawk ' { NR > 1 } '



#vcftools --vcf - --freq  --out $processed/$dmel.${i}SNPs

#use --freq2 and bioawk to retain bi-allelic sites and remove allele freq

: <<'SKIP'
vcftools --gzvcf dmel-all-r6.46.synSNPs.vcf.gz --stdout --freq2 \
| bioawk -t -c header 'BEGIN {print "CHROM" "\t" "POS" "\t" "N_ALLELES" "\t" "N_CHR"} \
$3=="2" {print $1 "\t" $2 "\t" $3 "\t" $4}' > dmel-all-r6.46.synSNPs.freq2
SKIP

for i in $processed/*.vcf.gz
do
#echo $i
intm=$(basename $i)
name=$(echo $intm | sed 's/dmel-all-r6.46.//; s/.vcf.gz//')
#echo "name is $name"
vcftools --gzvcf $i --stdout --freq2 \
| bioawk -t -c header 'BEGIN {print "CHROM" "\t" "POS" "\t" "N_ALLELES" "\t" "N_CHR"} \
$3=="2" {print $1 "\t" $2 "\t" $3 "\t" $4}' > $results/dmel-all-r6.46.$name.freq2 &
done

vcftools --gzvcf DSPR.r6.SNPs.vcf.gz --stdout --freq2 \
| bioawk -t -c header 'BEGIN {print "CHROM" "\t" "POS" "\t" "N_ALLELES" "\t" "N_CHR"} \
$3=="2" {print $1 "\t" $2 "\t" $3 "\t" $4}' > DSPR.r6.SNPs.freq2

vcftools --gzvcf $together/parr_out.annotated.vcf.gz --stdout --freq2 \
| bioawk -t -c header 'BEGIN {print "CHROM" "\t" "POS" "\t" "N_ALLELES" "\t" "N_CHR"} \
$3=="2" {print $1 "\t" $2 "\t" $3 "\t" $4}' > parr_out.annotated.freq2

vcftools --gzvcf A1.annotated.vcf.gz --stdout --freq2 \
| bioawk -t -c header 'BEGIN {print "CHROM" "\t" "POS" "\t" "N_ALLELES" "\t" "N_CHR"} \
$3=="2" {print $1 "\t" $2 "\t" $3 "\t" $4}' > A1.annotated.freq2



echo "start to wait"
wait
echo "waiting ends"

for i in $processed/*.vcf.gz
do
intn=$(basename $i)
name=$(echo $intn| sed 's/dmel-all-r6.46.//; s/SNPs.vcf.gz//; s/.annotated// ; s/.vcf.gz//')
echo $name
vcftools --gzvcf $i --stdout --window-pi 10000 \
| gawk -v type=$name 'BEGIN {print "CHROM" "\t" "BIN_START" "\t" "BIN_END" "\t" "N_VARIANTS" "\t" "PI" "\t" "TYPE"} \
NR>1 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" type}' \
> $results/dmel-all-r6.46.${name}.pi
done





vcftools --gzvcf ${together}/parr_out.annotated.vcf.gz --stdout --window-pi 10000 > parr_out.annotated.pi
vcftools --gzvcf parr_out.annotated.vcf.gz --stdout --window-pi 100000 > parr_out.annotated.pi2
vcftools --gzvcf A1.annotated.vcf.gz --stdout --window-pi 10000 > A1.annotated.pi
#vcftools --gzvcf A1.annotated.vcf.gz --stdout --window-pi 100000 > A1.annotated.pi2
#vcftools --gzvcf A1.annotated.vcf.gz --stdout --window-pi 1000000 > A1.annotated.pi3
vcftools --gzvcf DSPR.r6.SNPs.vcf.gz --stdout --window-pi 10000 > r6.annotated.pi
#vcftools --gzvcf DSPR.r6.SNPs.vcf.gz --stdout --window-pi 100000 > r6.annotated.pi2