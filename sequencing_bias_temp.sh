#!/bin/sh
#SBATCH --job-name=read
#SBATCH -A jje_lab
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=6G
#SBATCH --output="read_bias.out"
source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK
ref_R6="/dfs7/jje/jenyuw/Assembling_ISO1/reference/NCBI.r6.renamed.fasta.gz"
our_asm="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly/hifiasm_19/ONT_19.fasta"
alignment="/dfs7/jje/jenyuw/Assembling_ISO1/results/alignment"
filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"

cd ${alignment}
##Meryl is installed with winnowmap. Another version is installed independently. Also, another installed with canu.
##We have to specify the exact version with winnowmap2
/pub/jenyuw/Software/Winnowmap/bin/meryl k=15 count ${ref_R6} threads=${nT} output refDB
/pub/jenyuw/Software/Winnowmap/bin/meryl print greater-than distinct=0.9998 refDB > ${alignment}/r6_repetitive_15.txt

read="regular_40k.fastq.gz adaptive_8k.fastq.gz"
for i in ${read}
do
prefix=${i%.fastq.gz}
winnowmap -a --cs -W ${alignment}/r6_repetitive_15.txt -x map-ont ${ref_R6} ${filtered}/${i} |\
samtools view -@ ${nT} -bh -q 10 |samtools sort -@ ${nT} -m 2G -o ${alignment}/${prefix}.winn.bam
samtools index ${alignment}/${prefix}.winn.bam

#minimap2 -a --cs -x lr:hq ${ref_R6} ${filtered}/${i} |\
#samtools view -@ ${nT} -bh -q 10 |samtools sort -@ ${nT} -m 2G -o ${alignment}/${prefix}.min.bam
#samtools index ${alignment}/${prefix}.min.bam
done

# $4: leftmost mapping position; $5: mapping quality; $9: template length; $10: read sequence
samtools view -@ ${nT} ${alignment}/regular_40k.winn.bam | gawk -v OFS='\t' '{print $4, $4+$9-1, $5, length($10)}' |\
> ${alignment}/regular_40k.winn.read_info.tsv
#output: leftmost, rightmost, mapQ, read_length