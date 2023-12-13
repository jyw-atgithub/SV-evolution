#!/bin/bash

#SBATCH --job-name=snp    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=2-62%1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=60   ## number of cores the job needs
#SBATCH --mem-per-cpu=6G     # requesting 6 GB memory per CPU, the max

#Because nanocaller produce intermediate files of the same name, we can only run one sample at a time

## path
ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
SNPs="/dfs7/jje/jenyuw/SV-project-temp/result/SNPs"
Merged_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/Merged_SNP"
processed="/dfs7/jje/jenyuw/SV-project-temp/result/processed_SNP"
coordinate_bed="/dfs7/jje/jenyuw/SV-project-temp/result/coordinate_bed"
merged_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/merged_SVs"
dmel="dmel-all-r6.46"
## prep
nT=$SLURM_CPUS_PER_TASK

source ~/.bashrc

## in an interactive mode
cd ${processed}

# Avoid using --missing-to-ref, because this does look like a good assumption
bcftools merge --threads ${nT} ${SNPs}/*.filtered.snps.vcf.gz -O z > ${Merged_SNP}/all.snps.vcf.gz
tabix -p vcf ${Merged_SNP}/all.snps.vcf.gz

conda activate everything
snpEff -v BDGP6.32.105 ${Merged_SNP}/all.snps.vcf.gz | bgzip -@ ${nT} -c  >${Merged_SNP}/all.snps.annotated.vcf.gz
conda deactivate

# NON_SYNONYMOUS={missense_variant,start_lost,stop_gained,stop_lost}
# now we only compare between synonymous and missense variants
declare -A excl; excl[nonsyn]='synonymous_variant'; excl[syn]='missense_variant'
declare -A incl; incl[nonsyn]='missense_variant'; incl[syn]='synonymous_variant'
for i in "${!excl[@]}"
do
  bedtools intersect \
    -header \
    -a ${Merged_SNP}/all.snps.annotated.vcf.gz \
    -b ${processed}/${dmel}.CDSspans.bed\
  | grep -v ${excl[${i}]} \
  | awk -v m=${incl[${i}]} ' $0 ~ /^#/ || $0 ~ m { print } ' \
  | bgzip -@ ${nT} -c \
  > ${processed}/${i}SNPs.vcf.gz
done

for i in dig pig si li three_prime_UTR five_prime_UTR; do
    echo $i
    bedtools intersect \
    -header \
    -a ${Merged_SNP}/all.snps.annotated.vcf.gz \
    -b ${processed}/${dmel}.${i}spans.bed \
    | bgzip -@ ${nT} -c \
    > $processed/${i}SNPs.vcf.gz
done

input="${merged_SVs}/truvari.svimASM.vcf.gz"
##either %END or  %INFO/END works
##ant character in the '' will be output, even a space!!
bcftools query -f '%CHROM\t%POS\t%END\n' ${input} > ${coordinate_bed}/SV-svimasm-coordinate.bed

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\n' ${input} |\
sed 's@\.@66666@g' > ${coordinate_bed}/intermediate-svimasm.tsv

len=1000
printf "" > ${coordinate_bed}/outside-svimasm-coordinate.bed
while read i
do
echo $i | gawk -v len=$len ' $2-len > 0 {print $1 "\t" $2-len "\t" $2 "\n" $1 "\t" $3 "\t" $3+len}' >> ${coordinate_bed}/outside-svimasm-coordinate.bed
# Didn't care about the o-base or 1-base coordinates.
done <${coordinate_bed}/intermediate-svimasm.tsv

bedtools subtract -a ${coordinate_bed}/outside-svimasm-coordinate.bed \
-b ${coordinate_bed}/SV-svimasm-coordinate.bed >${coordinate_bed}/pure-outside-svimasm.bed

bedtools intersect -header -a ${processed}/synSNPs.vcf.gz -b ${coordinate_bed}/pure-outside-svimasm.bed |\
bcftools sort --max-mem 2G |\
bgzip -@ 8 -c > ${processed}/syn.outside-1000-svimasm.snps.vcf.gz