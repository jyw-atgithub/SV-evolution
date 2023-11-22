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

#parallel -j 8 bgzip -dk -@ 2 {} ::: `ls ${con_SVs}/*.tru_con.sort.vcf.gz`
#this only need once
#SURVIVOR can NOT take vcf.gz
ls ${con_SVs}/*.tru_con.sort.vcf >${con_SVs}/sample_files.txt
SURVIVOR merge ${con_SVs}/sample_files.txt 0 1 1 1 0 50 ${merged_SVs}/merge.tru_con.vcf
## SURVIVOR requires much memorey!
#bcftools merge -m none ${con_SVs}/*.tru_con.sort.vcf.gz | bgzip -@ ${nT} > ${merged_SVs}/merge.tru_con.vcf.gz
#I don't know why bcftools merge report error, so I use SURVIVOR instead.
bgzip -k -f -@ ${nT} ${merged_SVs}/merge.tru_con.vcf
bcftools sort --max-mem 2G ${merged_SVs}/merge.tru_con.vcf.gz |bgzip -@ ${nT} > ${merged_SVs}/merge.tru_con.sort.vcf.gz
#This sort step is pretty slow and takes much temp space. HPC3 ALWAYS FAILED ON IT! Do on Thoth
bcftools index -f -t ${merged_SVs}/merge.tru_con.sort.vcf.gz

wait
module load python/3.10.2

truvari collapse --sizemax 200000000 -k maxqual \
-i ${merged_SVs}/merge.tru_con.sort.vcf.gz \
-c ${merged_SVs}/truvari_collapsed.vcf -f ${ref_genome} |\
bcftools sort --max-mem 2G |\
bgzip -@ ${nT} > ${merged_SVs}/truvari.tru_con.vcf.gz

bcftools index -f -t ${merged_SVs}/truvari.tru_con.vcf.gz

module unload python/3.10.2
echo "This is the end!"