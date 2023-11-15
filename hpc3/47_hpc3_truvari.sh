#!/bin/bash

#SBATCH --job-name=truvari    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=6   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU
source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
con_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/consensus_SVs"
merged_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/merged_SVs"


module load python/3.10.2

bcftools merge -m none --force-samples ${con_SVs}/*.tru_con.vcf | bgzip -@ ${nT} > ${merged_SVs}/merge.vcf.gz

tabix -f -p vcf ${merged_SVs}/merge.vcf.gz

truvari collapse --sizemax 200000000 \
-i ${merged_SVs}/merge.vcf.gz -o ${merged_SVs}/truvari_merge.vcf \
-c ${merged_SVs}/truvari_collapsed.vcf -f ${ref_genome}

module unload python/3.10.2
echo " This is the end!"