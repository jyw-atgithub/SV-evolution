#!/bin/bash
## path
ref_genome="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"
raw="/home/jenyuw/SV-project/raw"
aligned_bam="/home/jenyuw/SV-project/result/aligned_bam"
SVs="/home/jenyuw/SV-project/result/SVs"

assemble="/home/jenyuw/SV-project/result/assemble"
polishing="/home/jenyuw/SV-project/result/polishing"
canu_proc="/home/jenyuw/SV-project/result/canu_processing"
scaffold="/home/jenyuw/SV-project/result/scaffold"
busco_out="/home/jenyuw/SV-project/result/busco_out"
## prep
source ~/.bashrc
nT=10

#on THOTH
cond activate sv-calling

## Assembly based

minimap2 -t 20 -B 5 -a -x asm5 --cs \
${ref_genome} \
/home/jenyuw/SV-project-backup/result/assembly/nv107_Flye_assembly.fasta \
|samtools view -b -h -@ 20 -o -|samtools sort -@ 20 -o nv107_mapped.sort.2.bam
#SVIM-asm
svim-asm haploid svim-asm-nv107.2 nv107_mapped.sort.2.bam /home/jenyuw/SV-project-backup/reference_genome/dmel-all-chromosome-r6.49.fasta

## Mapping-based
##9mapping-based.sh

for i in $(ls ${raw}/*_combined.fastq)
do
name=$(basename ${i}|sed s/"_combined.fastq"//g)
echo $name
echo ${aligned_bam}/${name}_ONT.sort.bam
#minimap2 -t ${nT} -B 5 -a -x map-ont \
#${ref_genome} $i |\
#samtools view -b -h -@ ${nT} -o - |\
#samtools sort -@ ${nT} -o ${aligned_bam}/${name}_ONT.sort.bam
#samtools index -@ ${nT} ${aligned_bam}/${name}_ONT.sort.bam
done

for j in $(ls ${aligned_bam}/*_ONT.sort.bam)
do
name=$(basename ${i}|sed s/"_ONT.sort.bam"//g)
echo $name

cd ${SVs}
#cuteSV
cuteSV --threads 20 \
-l 50 -L 500000 \
-r 500 -q 20 -s 3 \
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
$j ${ref_genome} nv107-cutesv.vcf .
#sniffles
sniffles --threads 20 \
--reference ${ref_genome} \
--input $j --vcf nv107-sniffles.vcf 
done


#on HYDRA
conda activate long-read

minimap2 -t 20 -B 5 -a -x map-ont \
/home/jenyuw/SV-project-backup/reference_genome/dmel-all-chromosome-r6.49.fasta \
/home/jenyuw/SV-project-backup/result/assembly/nv107_Flye_assembly.fasta |samtools view -b -h -@ 20 -o nv107_mapped.bam
samtools sort -@ 20 -o nv107_mapped.sort.bam nv107_mapped.bam

samtools index -@ 20 nv107_mapped.sort.bam nv107_mapped.bam




samtools index -@ 20 nv107_mapped.sort.bam nv107_mapped.2.bam
cuteSV --threads 20 \
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
nv107_mapped.sort.bam /home/jenyuw/SV-project-backup/reference_genome/dmel-all-chromosome-r6.49.fasta \
nv107.vcf .

cuteSV --threads 20 \
-l 50 -L 500000 \
-r 500 -q 20 -s 3 \
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
nv107_mapped.sort.bam /home/jenyuw/SV-project-backup/reference_genome/dmel-all-chromosome-r6.49.fasta \
cutesv-nv107.4.vcf .

cuteSV --threads 20 \
-l 50 -L 500000 \
-r 500 -q 20 -s 3 \
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
nv107_mapped.sort.2.bam /home/jenyuw/SV-project-backup/reference_genome/dmel-all-chromosome-r6.49.fasta \
cutesv-nv107.5.vcf .

svim-asm haploid svim-asm-nv107 nv107_mapped.sort.bam /home/jenyuw/SV-project-backup/reference_genome/dmel-all-chromosome-r6.49.fasta


sniffles --threads 20 \
--reference /home/jenyuw/SV-project-backup/reference_genome/dmel-all-chromosome-r6.49.fasta \
--input nv107_mapped.sort.bam --vcf nv107-sniffles.vcf 