#!/bin/bash

#SBATCH --job-name=map    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=16   ## number of cores the job needs
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
trimmed="/dfs7/jje/jenyuw/SV-project-temp/result/trimmed"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"

source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK

file=`head -n $SLURM_ARRAY_TASK_ID ${trimmed}/namelist.txt |tail -n 1`
name=$(basename ${file}|sed s/".trimmed.fastq.gz"//g)
read_type=`echo ${name} | gawk -F "_" '{print $2}'`
mapping_option=(["CLR"]="map-pb" ["hifi"]="asm20" ["ONT"]="map-ont")

module load samtools/1.15.1

minimap2 -t ${nT} -a -x ${mapping_option[$read_type]} \
${ref_genome} ${file} |\
samtools view -b -h -@ ${nT} -o - |\
samtools sort -@ ${nT} -o ${aligned_bam}/${name}.trimmed-ref.sort.bam
samtools index -@ ${nT} ${aligned_bam}/${name}.trimmed-ref.sort.bam

module unload samtools/1.15.1