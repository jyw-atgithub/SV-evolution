#!/bin/sh
#SBATCH --job-name=read
#SBATCH -A jje_lab
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=6G
#SBATCH --output="read_bias.out"
#sbatch --constraint=nvme
source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK
ref_R6="/dfs7/jje/jenyuw/Assembling_ISO1/reference/NCBI.r6.renamed.fasta.gz"
our_asm="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly/hifiasm_19/ONT_19.fasta"
alignment="/dfs7/jje/jenyuw/Assembling_ISO1/results/alignment"
filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"

read="regular_40k.fastq.gz adaptive_8k.fastq.gz"

read="adaptive_8k.fastq.gz"
cd ${alignment}
##Meryl is installed with winnowmap. Another version is installed independently. Also, another installed with canu.
##We have to specify the exact version with winnowmap2
/pub/jenyuw/Software/Winnowmap/bin/meryl k=15 count ${ref_R6} threads=${nT} output refDB
/pub/jenyuw/Software/Winnowmap/bin/meryl print greater-than distinct=0.9998 refDB > ${alignment}/r6_repetitive_15.txt

target=${ref_R6}
for i in ${read}
do
prefix=${i%.fastq.gz}
winnowmap -a --cs -W ${alignment}/r6_repetitive_15.txt -x map-ont ${target} ${filtered}/${i} |\
samtools view -@ ${nT} -bh -q 10 |samtools sort -@ ${nT} -m 2G -o ${alignment}/${prefix}.winn.bam
samtools index ${alignment}/${prefix}.winn.bam

minimap2 -t ${nT} -a --cs -x lr:hq ${target} ${filtered}/${i} |\
samtools view -@ ${nT} -bh -q 10 |samtools sort -@ ${nT} -m 2G -o ${alignment}/${prefix}.mini.bam
samtools index ${alignment}/${prefix}.mini.bam
done



target=${our_asm}
##Meryl is installed with winnowmap. Another version is installed independently. Also, another installed with canu.
##We have to specify the exact version with winnowmap2
/pub/jenyuw/Software/Winnowmap/bin/meryl k=15 count ${target} threads=${nT} output ONT19DB
/pub/jenyuw/Software/Winnowmap/bin/meryl print greater-than distinct=0.9998 ONT19DB > ${alignment}/ONT19_repetitive_15.txt

for i in ${read}
do
prefix=${i%.fastq.gz}
asm_name=`basename ${our_asm} .fasta`
winnowmap -t ${nT} -a --cs -W ${alignment}/ONT19_repetitive_15.txt -x map-ont ${target} ${filtered}/${i} |\
samtools view -@ ${nT} -bh -q 10 |samtools sort -@ ${nT} -m 2G -o ${alignment}/${prefix}.${asm_name}.winn.bam
samtools index ${alignment}/${prefix}.${asm_name}.winn.bam
echo "winnowmap2 finished for ${i}"

minimap2 -t ${nT} -a --cs -x lr:hq ${target} ${filtered}/${i} |\
samtools view -@ ${nT} -bh -q 10 |samtools sort -@ ${nT} -m 2G -o ${alignment}/${prefix}.${asm_name}.mini.bam
samtools index ${alignment}/${prefix}.${asm_name}.mini.bam
echo "minimap2 finished for ${i}"
done


# $2: Flag; $3: reference chr; $4: leftmost mapping position; $5: mapping quality; $9: template length; $10: read sequence
# Because this is alignment of ONT long reads and not paired-end reads, there is no "template length"
for i in `ls ${alignment}/*.bam`
do
prefix=`basename ${i} .bam`
samtools view -@ ${nT} -F 2048 ${i} | gawk -v OFS="\t" '{print $3, $2, $4, $4+length($10), $5, length($10)}' > ${alignment}/${prefix}.read_info.tsv
done
#output: reference name, flag, leftmost, rightmost, mapQ, read_length