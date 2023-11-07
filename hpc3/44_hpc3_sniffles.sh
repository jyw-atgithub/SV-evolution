#!/bin/bash

#SBATCH --job-name=sniffles    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=8   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
trimmed="/dfs7/jje/jenyuw/SV-project-temp/result/trimmed"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"

source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK

file=`head -n $SLURM_ARRAY_TASK_ID ${aligned_bam}/alignedlist.txt |tail -n 1`
name=$(basename ${file}|sed s/".trimmed-ref.sort.bam"//g)
read_type=`echo ${name} | gawk -F "_" '{print $2}'`


module load python/3.10.2

sniffles --threads ${nT} --allow-overwrite --sample-id ${name}-snif \
--minsupport 10 \
--minsvlen 50 --mapq 20 --min-alignment-length 500 \
--cluster-merge-pos 270 \
--max-del-seq-len 100000 \
--reference ${ref_genome} \
--input ${file} --vcf "${SVs}/${name}.sniffles.vcf"

module unload python/3.10.2

echo "This is the end"