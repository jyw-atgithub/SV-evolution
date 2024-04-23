#!/bin/bash

## path
merged_SVs="/home/jenyuw/SV-project/result/merged_SVs"
merged_SNP="/home/jenyuw/SV-project/result/merged_SNP"
coordinate_bed="/home/jenyuw/SV-project/result/coordinate_bed"

input="/home/jenyuw/SV-project/result/merged_SVs/all.consensus-005.vcf"
##either %END or  %INFO/END works
##ant character in the '' will be output, even a space!!
bcftools query -f '%CHROM\t%POS\t%END\n' ${input} > ${coordinate_bed}/SV-coordinate.bed

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\n' ${input} > ${coordinate_bed}/intermediate.tsv

len=1000
printf "" > ${coordinate_bed}/outside-coordinate.bed
while read i
do
echo $i | gawk -v len=$len ' $2-len > 0 {print $1 "\t" $2-len "\t" $2 "\n" $1 "\t" $3 "\t" $3+len}' >> ${coordinate_bed}/outside-coordinate.bed
# Didn't care about the o-base or 1-base coordinates.
done <${coordinate_bed}/intermediate.tsv

bedtools subtract -a ${coordinate_bed}/outside-coordinate.bed -b ${coordinate_bed}/SV-coordinate.bed >${coordinate_bed}/pure-outside.bed

bedtools intersect -header -a ${merged_SNP}/all.snps.vcf.gz -b ${coordinate_bed}/pure-outside.bed |\
bgzip -@ 8 -c > ${merged_SNP}/all.outside-1000.snps.vcf.gz

# dmel-all-r6.46.CDSspans.bed was generated in the DSPR-snp project
declare -A excl; excl[nonsyn]='synonymous_variant'; excl[syn]='missense_variant'
declare -A incl; incl[nonsyn]='missense_variant'; incl[syn]='synonymous_variant'
for i in "${!excl[@]}"
do
  bedtools intersect \
    -header \
    -a ${merged_SNP}/all.snps.annotated.vcf.gz \
    -b ${coordinate_bed}/dmel-all-r6.46.CDSspans.bed\
  | grep -v ${excl[${i}]} \
  | awk -v m=${incl[${i}]} ' $0 ~ /^#/ || $0 ~ m { print } ' \
  | bgzip -c -@ 8 > ${merged_SNP}/all.${i}SNPs.vcf.gz
done

