#!/bin/bash
## path
ref_genome="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"

raw="/home/jenyuw/SV-project/raw"
raw_prj="/home/jenyuw/SV-project/raw/PRJNA929424"

trimmed="/home/jenyuw/SV-project/result/trimmed"
assemble="/home/jenyuw/SV-project/result/assemble"
aligned_bam="/home/jenyuw/SV-project/result/aligned_bam"
polishing="/home/jenyuw/SV-project/result/polishing"
canu_proc="/home/jenyuw/SV-project/result/canu_processing"
scaffold="/home/jenyuw/SV-project/result/scaffold"
SVs="/home/jenyuw/SV-project/result/SVs"
## prep
source ~/.bashrc
nT=8

#on THOTH
conda activate sv-calling
: <<'SKIP'
## Assembly based
for i in $(ls ${polishing}/nv*.polished.pilon.1.fasta)
do
name=$(basename $i | sed s/.polished.pilon.1.fasta//g)
echo "mapping ${name} polished assembly to reference genome"
minimap2 -t ${nT} -a -x asm5 --cs \
${ref_genome} ${i} \
|samtools view -b -h -@ ${nT} -o -|samtools sort -@ ${nT} -o ${aligned_bam}/${name}.Flye-ref.sort.bam
samtools index ${aligned_bam}/${name}.Flye-ref.sort.bam
#SVIM-asm
echo "calling SVs of ${name} wiht assembly based methods"
svim-asm haploid ${SVs}/${name}-svim-asm ${aligned_bam}/${name}.Flye-ref.sort.bam ${ref_genome}
done


## Mapping-based, so we have to map the "trimmed&filtered" ONT reads to the reference genome.
##9mapping-based.sh

for i in $(ls ${trimmed}/nv*.trimmed.fastq)
do
name=$(basename ${i}|sed s/".trimmed.fastq"//g)
echo "mapping ${name} trimmed reads to reference genome"

minimap2 -t ${nT} -a -x map-ont \
${ref_genome} $i |\
samtools view -b -h -@ ${nT} -o - |\
samtools sort -@ ${nT} -o ${aligned_bam}/${name}.trimmed-ref.sort.bam
samtools index -@ ${nT} ${aligned_bam}/${name}.trimmed-ref.sort.bam
done
SKIP

for j in $(ls ${aligned_bam}/*.trimmed-ref.sort.bam)
do
name=$(basename ${j}|sed s/".trimmed-ref.sort.bam"//g)
echo "calling SVs of ${name} wiht mapping based methods"
cd ${SVs}
#cuteSV
cuteSV --threads 20 \
-l 50 -L 500000 \
-r 500 -q 20 -s 3 \
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
$j ${ref_genome} "${name}-cutesv.vcf" .
#sniffles
sniffles --threads 20 \
--reference ${ref_genome} \
--input $j --vcf "${name}-sniffles.vcf" 
#SVIM
svim alignment --min_sv_size 50 --sample ${name} \
${SVs}/${name}-SVIM $j ${ref_genome}
done


#on HYDRA
#conda activate long-read
