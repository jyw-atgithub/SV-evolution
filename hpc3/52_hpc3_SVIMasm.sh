#!/bin/bash

#SBATCH --job-name=svim-asm    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=10   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=10G     # requesting memory per CPU

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


#mv ${scaffold}/rr${SLURM_ARRAY_TASK_ID}/ragtag.scaffold.fasta ${scaffold}/${name}.scaffold.fasta

### USE asm10 instead of asm5 ###

minimap2 -t ${nT} -a -x asm10 --cs --eqx \
${ref_genome} ${scaffold}/${name}/ragtag.scaffold.fasta \
|samtools view -b -h -@ ${nT} -o -|samtools sort -@ ${nT} -o ${aligned_bam}/${name}.final-ref.sort.bam
samtools index ${aligned_bam}/${name}.final-ref.sort.bam

conda activate sv-calling
svim-asm haploid --sample ${name}_svimASM --min_sv_size 50 \
${SVs}/${name}_svimASM ${aligned_bam}/${name}.final-ref.sort.bam ${ref_genome}

mv ${SVs}/${name}_svimASM/variants.vcf ${SVs}/${name}.svimASM.vcf 
conda deactivate