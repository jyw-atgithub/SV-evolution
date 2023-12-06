#!/bin/bash

#SBATCH --job-name=syri    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=10   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU

source ~/.bashrc

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
assemble="/dfs7/jje/jenyuw/SV-project-temp/result/assemble"
scaffold="/dfs7/jje/jenyuw/SV-project-temp/result/scaffold"
nT=$SLURM_CPUS_PER_TASK

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
ls ${scaffold}/*/ragtag.scaffold.fasta > ${scaffold}/scfd_list.txt
fi

file=`head -n $SLURM_ARRAY_TASK_ID ${scaffold}/scfd_list.txt |tail -n 1`
name=$(echo ${file} | cut -d '/' -f 8)
echo "the file is ${file}"
echo "the name is ${name}"


mkdir ${SVs}/${name}_SYRI
cd ${SVs}/${name}_SYRI

samtools view -@ 2 -h ${aligned_bam}/${name}.final-ref.sort.bam|sed s/_RagTag//g |\
samtools view -@ 2 -h --bam -o ${SVs}/${name}_SYRI/${name}.final-ref.sort.bam

#cp ${aligned_bam}/${name}.final-ref.sort.bam.bai ${SVs}/${name}_SYRI/${name}.final-ref.sort.bam.bai
cp ${ref_genome} ${SVs}/${name}_SYRI/r649.ref.fasta
cat ${file}|sed s/_RagTag//g > ${SVs}/${name}_SYRI/${name}.scaffold.fasta

module load python/3.10.2

fixchr -c ${SVs}/${name}_SYRI/${name}.final-ref.sort.bam -F B \
-r ${SVs}/${name}_SYRI/r649.ref.fasta -q ${SVs}/${name}_SYRI/${name}.scaffold.fasta \
-f --contig_size 100000 --dir ${SVs}/${name}_SYRI --prefix ${name}

cat ${SVs}/${name}_SYRI/ref.filtered.fa| sed s/'>4'/'>chr4'/g > ${SVs}/${name}_SYRI/ref.filtered.rn.fa
cat ${SVs}/${name}_SYRI/qry.filtered.fa| sed s/'>4'/'>chr4'/g > ${SVs}/${name}_SYRI/qry.filtered.rn.fa

minimap2 -t ${nT} -a -x asm5 --cs --eqx \
${SVs}/${name}_SYRI/ref.filtered.rn.fa  ${SVs}/${name}_SYRI/qry.filtered.rn.fa \
|samtools view -b -h -@ ${nT} -o -|samtools sort -@ ${nT} -o ${SVs}/${name}_SYRI/sort.bam
samtools index ${SVs}/${name}_SYRI/sort.bam

syri  -F B -c ${SVs}/${name}_SYRI/sort.bam -r ${SVs}/${name}_SYRI/ref.filtered.rn.fa \
-q ${SVs}/${name}_SYRI/qry.filtered.rn.fa --dir ${SVs}/${name}_SYRI \
--nc 3 --samplename ${name}_SYRI --prefix ${name}_mm2 --unic 500

module unload python/3.10.2

rm ${SVs}/${name}_SYRI/${name}.final-ref.sort.bam
rm ${SVs}/${name}_SYRI/r649.ref.fasta
rm ${SVs}/${name}_SYRI/${name}.scaffold.fasta
rm ${SVs}/${name}_SYRI/ref.filtered.fa
rm ${SVs}/${name}_SYRI/qry.filtered.fa
##grep -v "^#" result/SVs/A1_CLR_SYRI/syri.vcf | gawk '{print $3}' | cut -c 1-3 |sort |uniq -c

echo "This is the end!!"