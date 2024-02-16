#!/bin/bash

#SBATCH --job-name=svim-asm    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p highmem       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=8   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=10G     # requesting memory per CPU

source ~/.bashrc

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
trimmed="/dfs7/jje/jenyuw/SV-project-temp/result/trimmed"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
purge_dups="/dfs7/jje/jenyuw/SV-project-temp/result/purge_dups"

nT=$SLURM_CPUS_PER_TASK

file=`head -n $SLURM_ARRAY_TASK_ID ${trimmed}/namelist.txt |tail -n 1`
name=$(basename ${file}|sed s/".trimmed.fastq.gz"//g)
read_type=`echo ${name} | gawk -F "_" '{print $2}'`
echo "the file is ${file}"
echo "the name is ${name}"


### USE asm5, use assembly as query

minimap2 -t ${nT} -a -x asm5 --cs --eqx \
${ref_genome} ${purge_dups}/${name}.final.fasta |\
samtools view -b -h -@ ${nT} -o -|samtools sort -@ ${nT} -o ${aligned_bam}/${name}.final-ref.sort.bam
samtools index ${aligned_bam}/${name}.final-ref.sort.bam

conda activate sv-calling
svim-asm haploid --sample ${name}_svimASM --min_sv_size 50 \
${SVs}/${name}_svimASM ${aligned_bam}/${name}.final-ref.sort.bam ${ref_genome}
deactivate

mv ${SVs}/${name}_svimASM/variants.vcf ${SVs}/${name}.svimASM.vcf 

bgzip -f --keep -@ ${nT} ${SVs}/${name}.svimASM.vcf 
bcftools sort --write-index --max-mem 4G -O z -o ${SVs}/${name}.svimASM.sort.gz ${SVs}/${name}.svimASM.vcf.gz