#!/bin/bash
dmel_ref="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"
dsim_ref="/home/jenyuw/SV-project/reference_genome/dsim-all-chromosome-r2.02.fasta"

merged_SVs="/home/jenyuw/SV-project/result/merged_SVs"
aligned_bam="/home/jenyuw/SV-project/result/aligned_bam"
polarizing="/home/jenyuw/SV-project/result/polarizing"
len=1000

bcftools view --thread 4 -i ' SVLEN<200 && SVLEN>90 && SVTYPE="INS" ' -r 2L:9000000-10000000 \
${merged_SVs}/all.consensus-005.vcf.gz 2> /dev/null |\
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\n' |\
gawk -v len=$len ' $2-len > 0 {print $1 "\t" $2-len "\t" $3+len}' |sort|uniq >${polarizing}/part_ins.bed

for i in $(ls ${aligned_bam}/*.trimmed-ref.sort.bam)
do
name=$(basename $i|sed s/.trimmed-ref.sort.bam//g)
bedtools intersect -a ${i} -b ${polarizing}/part_ins.bed >${polarizing}/${name}.part-ins.bam
done

parallel -j 3 bedtools intersect -a {} -b ${polarizing}/part_ins.bed > {#}.part_ins.bam ::: $(ls ${aligned_bam}/nv*.trimmed-ref.sort.bam)
parallel -j 10 echo {#}{} ::: $(ls ${aligned_bam}/*.trimmed-ref.sort.bam)
parallel -j 10 echo {#}{} ::: $(ls ${aligned_bam}/*.trimmed-ref.sort.bam |basename |sed s/.trimmed-ref.sort.bam//g)

::: $(ls ${aligned_bam}/*.trimmed-ref.sort.bam) 