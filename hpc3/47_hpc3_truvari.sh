#!/bin/bash

#SBATCH --job-name=truvari    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=16   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU
source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
con_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/consensus_SVs"
merged_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/merged_SVs"


module load python/3.10.2

bcftools merge -m none --force-samples ${con_SVs}/*.tru_con.sort.vcf | bgzip -@ ${nT} > ${merged_SVs}/merge.tru_con.vcf.gz

tabix -f -p vcf ${merged_SVs}/merge.tru_con.vcf.gz

truvari collapse --sizemax 200000000 -k maxqual \
-i ${merged_SVs}/merge.tru_con.vcf.gz \
-c ${merged_SVs}/truvari_collapsed.vcf -f ${ref_genome} | bgzip -@ ${nT} > ${merged_SVs}/truvari.tru_con.vcf.gz

bcftools sort --max-mem 2G -O v -o ${merged_SVs}/truvari.tru_con.sort.vcf ${merged_SVs}/truvari.tru_con.vcf.gz

module unload python/3.10.2
echo " This is the end!"