#!/bin/bash
dmel_ref="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"
dsim_ref="/home/jenyuw/SV-project/reference_genome/dsim-all-chromosome-r2.02.fasta"
trimmed="/home/jenyuw/SV-project/result/trimmed"
merged_SVs="/home/jenyuw/SV-project/result/merged_SVs"
aligned_bam="/home/jenyuw/SV-project/result/aligned_bam"
polarizing="/home/jenyuw/SV-project/result/polarizing"
len=1000
nT=24

bcftools view --thread 4 -i ' SVLEN<2000 && SVLEN>1200 && SVTYPE="INS" ' -r 2L:9000000-10000000 \
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
## Parsing the sam file generated by Murmmer4
## origin: https://github.com/mummer4/mummer/issues/24

cat ${polarizing}/dmel-dsim.sam | sed 's/HD\ /HD/; s/1.0\ /1.0/; s/\tSO:coordinate/SO:coordinate/; s/VN1/VN:1/; s/HD/HD\t/; s/SO:unsorted/\tSO:unsorted/; s/@PG /@PG\t/; s/ PN/\tPN/; s/ VN/\tVN/; s/ CL/\tCL/' |\
samtools view -@ 10 -h --reference ${dmel_ref} - |samtools sort -@ 10 -O bam -o ${polarizing}/corrected.bam
samtools index -@ 10 ${polarizing}/corrected.bam

for i in $(ls ${aligned_bam}/*.trimmed-ref.sort.bam)
do
name=$(basename $i|sed s/.trimmed-ref.sort.bam//g)
bedtools intersect -a ${i} -b ${polarizing}/part_ins.bed >${polarizing}/${name}.part-ins.bam
done

doit () { name=$(basename $1|sed s/trimmed-ref.sort.bam/dmel/);
echo $name;
bedtools intersect -a $1 -b /home/jenyuw/SV-project/result/polarizing/part_ins.bed >/home/jenyuw/SV-project/result/polarizing/${name}.part-ins.bam; 
 }
export -f doit
#declare -f doit
ls ${aligned_bam}/*.trimmed-ref.sort.bam | parallel -j 6 --eta "doit {}"

# This is much slower. Avoid -j+0.
#ls ${aligned_bam}/*.trimmed-ref.sort.bam | parallel -j+0 --eta "bedtools intersect -a {} -b ${polarizing}/part_ins.bed >${polarizing}/{/.}.part-ins.bam"


for i in $(ls ${aligned_bam}/*.trimmed-dsim.sort.bam)
do
name=$(basename $i|sed s/.trimmed-dsim.sort.bam//g)
bedtools intersect -a ${i} -b ${polarizing}/part_ins.bed >${polarizing}/${name}.dsim.part-ins.bam
done


samtools merge -f -o ${polarizing}/all.dmel.part-ins.bam ${polarizing}/*.dmel.part-ins.bam
samtools index -@ 8 ${polarizing}/all.dmel.part-ins.bam
samtools merge -o ${polarizing}/all.dsim.part-ins.bam ${polarizing}/*.dsim.part-ins.bam
##GNU parallel keeps failing
#ls ${aligned_bam}/nv*.trimmed-ref.sort.bam| parallel -j 3 "bedtools intersect -a ${polarizing}/part_ins.bed -b {} > {/.}.bam"
#parallel -j 10 echo {#}{} ::: $(ls ${aligned_bam}/*.trimmed-ref.sort.bam)
#parallel -j 10 echo {#}{} ::: $(ls ${aligned_bam}/*.trimmed-ref.sort.bam |basename |sed s/.trimmed-ref.sort.bam//g)
#::: $(ls ${aligned_bam}/*.trimmed-ref.sort.bam) 


dmel_ref="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"
dsim_ref="/home/jenyuw/SV-project/reference_genome/dsim-all-chromosome-r2.02.fasta"

svim-asm haploid --sample sim_to_mel --min_sv_size 50 --min_mapq 0 \
/home/jenyuw/SV-project/result/polarizing /home/jenyuw/SV-project/result/polarizing/corrected.bam ${dmel_ref}

lra global -CONTIG ${dmel_ref}
lra align -CONTIG ${dmel_ref} ${dsim_ref} -t 16 -p s > lra-dsim-dmel.sam

samtools view -@ 8 -b lra-dsim-dmel.sam -o -|samtools sort -@ 8 - > lra-dsim-dmel.bam

svim-asm haploid --sample "sim_to_mel" --min_sv_size 50 \
/home/jenyuw/SV-project/result/polarizing/lra /home/jenyuw/SV-project/result/polarizing/lra-dsim-dmel.bam ${dmel_ref}

minimap2 -a -x asm5 -t 16 ${dmel_ref} ${dsim_ref} |\
samtools view -@ 8 -b lra-dsim-dmel.sam -o -|samtools sort -@ 8 - > mm2-dsim-dmel.bam

svim-asm haploid --sample "sim_to_mel" --min_sv_size 50 \
/home/jenyuw/SV-project/result/polarizing/mm2 /home/jenyuw/SV-project/result/polarizing/mm2-dsim-dmel.bam ${dmel_ref}
