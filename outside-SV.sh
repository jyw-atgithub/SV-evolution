#!/bin/bash

## path
merged_SVs="/home/jenyuw/SV-project/result/merged_SVs"
coordinate_bed="/home/jenyuw/SV-project/result/coordinate_bed"

##either %END or  %INFO/END works
##ant character in the '' will be output, even a space!!
bcftools query -f '%CHROM\t%POS\t%END\n' ${merged_SVs}/all.consensus.vcf > ${coordinate_bed}/SV-coordinate.bed

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\n' ${merged_SVs}/all.consensus.vcf > ${coordinate_bed}/intermediate.tsv

len=1000
printf "" > ${coordinate_bed}/outside-coordinate.bed
while read i
do
echo $i | gawk -v len=$len ' $2-1000 > 0 {print $1 "\t" $2-len "\t" $2 "\n" $1 "\t" $3 "\t" $3+len}' >> ${coordinate_bed}/outside-coordinate.bed
done <${coordinate_bed}/intermediate.tsv


