#!/bin/bash
dmel_ref="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"
dsim_ref="/home/jenyuw/SV-project/reference_genome/dsim-all-chromosome-r2.02.fasta"
trimmed="/home/jenyuw/SV-project/result/trimmed"
merged_SVs="/home/jenyuw/SV-project/result/merged_SVs"
aligned_bam="/home/jenyuw/SV-project/result/aligned_bam"
polarizing="/home/jenyuw/SV-project/result/polarizing"
len=1000
nT=24

bcftools view --thread 4 -i ' SVLEN<200 && SVLEN>90 && SVTYPE="INS" ' -r 2L:9000000-10000000 \
${merged_SVs}/all.consensus-005.vcf.gz 2> /dev/null |\
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\n' |\
gawk -v len=$len ' $2-len > 0 {print $1 "\t" $2-len "\t" $3+len}' |sort|uniq >${polarizing}/part_ins.bed

for i in $(ls ${trimmed}/*.trimmed.fastq)
do
name=$(basename ${i}|sed s/".trimmed.fastq"//g)
echo "mapping ${name} trimmed reads to D.sim. reference genome"

minimap2 -t ${nT} -a -x map-ont \
${dsim_ref} $i |\
samtools view -b -h -@ ${nT} -o - |\
samtools sort -@ ${nT} -o ${aligned_bam}/${name}.trimmed-dsim.sort.bam
samtools index -@ ${nT} ${aligned_bam}/${name}.trimmed-dsim.sort.bam
done

for i in $(ls ${trimmed}/*.trimmed.rn.fastq.gz)
do
name=$(basename ${i}|sed s/".trimmed.rn.fastq.gz"//g)
echo "mapping ${name} trimmed reads to D.sim. reference genome"

minimap2 -t ${nT} -a -x map-pb \
${dsim_ref} $i |\
samtools view -b -h -@ ${nT} -o - |\
samtools sort -@ ${nT} -o ${aligned_bam}/${name}.trimmed-dsim.sort.bam
samtools index -@ ${nT} ${aligned_bam}/${name}.trimmed-dsim.sort.bam
done

#Do NOT use minimap2 to perform asm-to-asm alignemnt, because it is not good at this.
#minimap2 -t 10 -a -x asm10 ${dmel_ref} ${dsim_ref} |\
#samtools view -b -h -@ 10 -o - |\
#samtools sort -@ 10 -o ${polarizing}/dmel-dsim.sort.bam

nucmer -t 20 --sam-long=${polarizing}/dmel-dsim.sam ${dmel_ref} ${dsim_ref}


for i in $(ls ${aligned_bam}/*.trimmed-ref.sort.bam)
do
name=$(basename $i|sed s/.trimmed-ref.sort.bam//g)
bedtools intersect -a ${i} -b ${polarizing}/part_ins.bed >${polarizing}/${name}.part-ins.bam
done

for i in $(ls ${aligned_bam}/*.trimmed-dsim.sort.bam)
do
name=$(basename $i|sed s/.trimmed-dsim.sort.bam//g)
bedtools intersect -a ${i} -b ${polarizing}/part_ins.bed >${polarizing}/${name}.dsim.part-ins.bam
done


samtools merge -o ${polarizing}/all.dmel.part-ins.bam ${polarizing}/*.dmel.part-ins.bam
samtools merge -o ${polarizing}/all.dsim.part-ins.bam ${polarizing}/*.dsim.part-ins.bam
##GNU parallel keeps failing
#ls ${aligned_bam}/nv*.trimmed-ref.sort.bam| parallel -j 3 "bedtools intersect -a ${polarizing}/part_ins.bed -b {} > {/.}.bam"
#parallel -j 10 echo {#}{} ::: $(ls ${aligned_bam}/*.trimmed-ref.sort.bam)
#parallel -j 10 echo {#}{} ::: $(ls ${aligned_bam}/*.trimmed-ref.sort.bam |basename |sed s/.trimmed-ref.sort.bam//g)
#::: $(ls ${aligned_bam}/*.trimmed-ref.sort.bam) 