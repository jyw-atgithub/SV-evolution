#!/bin/bash

cd /home/jenyuw/DSPR_snp/snp_processing/processed

basepath=/home/jenyuw/DSPR_snp/snp_processing
raw="/home/jenyuw/Reference_genome"
processed=$basepath/processed
results=/home/jenyuw/DSPR_snp/snp_processing/results
mergedSNP=/home/jenyuw/DSPR_snp/raw_fq/together/parr_out.annotated.vcf

dmel=dmel-all-r6.46
dsim=dsim-r3_1
dmel_chr="Dmel6/dmel-all-chromosome-r6.46.fasta" #original
dsim_chr="Dsim3/Dsim.3.1.fasta" #original
dmel_gff="Dmel6/dmel-all-filtered-r6.46.gff" #original from flybase
#dyak=dyak-r2
#merged DSPR SNP vcf is saved here


: <<'SKIP'
view -m2 -M2 -v snps $raw/SCGB.founders-only-PASS-Q30-SNPs.vcf.gz \
| vcf-subset -c A1,A2ST,A3,A4,A5,A6,A7,B1,B2,B3,B4,B6,AB8 \
| gzip \
> $processed/DSPR.SNPs.vcf.gz
SKIP


#process before merging vcf files
for i in *annotated.vcf; 
do
bgzip -k -c -@ 59 $i > ../modified_snp/$i.gz
done
## removing a tag because "Incorrect number of FORMAT/AD values at 2L:172309, cannot merge. The tag is defined as Number=R, but found 3 values and 2 alleles"
bcftools annotate -x FORMAT/AD --threads 56 A1.annotated.vcf.gz >A1.annotated.noAD.vcf
## compressing the vcf files so bcftools can be used
bgzip -k -@ 56 A1.annotated.noAD.vcf
## index the vcf files with tabix (part of samtools)
tabix -p vcf A1.annotated.noAD.vcf.gz
## merge them with bcftools
bcftools merge --threads 56 -O v -o merged.annotated.vcf *.annotated.noAD.vcf.gz


# NON_SYNONYMOUS={missense_variant,start_lost,stop_gained,stop_lost}
# now we only compare between synonymous and missense variants
declare -A excl; excl[nonsyn]='synonymous_variant'; excl[syn]='missense_variant'
declare -A incl; incl[nonsyn]='missense_variant'; incl[syn]='synonymous_variant'
for i in "${!excl[@]}"
do
  bedtools intersect \
    -header \
    -a $mergedSNP \
    -b $processed/$dmel.CDSspans.bed\
  | grep -v ${excl[${i}]} \
  | awk -v m=${incl[${i}]} ' $0 ~ /^#/ || $0 ~ m { print } ' \
  | tee >(
      vcftools \
        --vcf - \
        --freq \
        --out $results/$dmel.${i}SNPs
    ) \
  | gzip \
  > $processed/$dmel.${i}SNPs.vcf.gz
done

for i in dig pig si li three_prime_UTR five_prime_UTR; do
  echo $i
  bedtools intersect \
    -header \
    -a $mergedSNP \
    -b $processed/$dmel.${i}spans.bed \
  | tee >(
      vcftools \
        --vcf - \
        --freq \
        --out $results/$dmel.${i}SNPs
    ) \
  | bgzip -@ 60 \
  > $processed/$dmel.${i}SNPs.vcf.gz
done


#for testing
: <<'SKIP'
for i in dig pig five_prime_UTR; do
  echo $i
  bedtools intersect \
    -header \
    -a $mergedSNP \
    -b $processed/$dmel.${i}spans.bed \
  | tee >(
      vcftools \
        --vcf - \
        --freq \
        --out $processed/$dmel.${i}SNPs
    ) \
  | bgzip -@ 56 \
  > $processed/$dmel.${i}SNPs.vcf.gz
done
SKPI