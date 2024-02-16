#!/bin/bash

#SBATCH --job-name=asm-truvari    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p highmem       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=10   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=10G     # requesting memory per CPUe. 
source ~/.bashrc

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
merged_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/merged_SVs"
con_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/consensus_SVs"
nT=$SLURM_CPUS_PER_TASK

bcftools merge -m none --threads ${nT} ${con_SVs}/*.asm-2.tru_con.sort.vcf.gz -O z -o - |\
bcftools sort --max-mem 4G -O z -o - > ${merged_SVs}/merge.asm-2.sort.vcf.gz
bcftools index -f -t ${merged_SVs}/merge.asm-2.sort.vcf.gz

module load python/3.10.2
##Do not use "-f" in truvari collapse, because this is an old function and less accurate
truvari collapse --sizemax 5000000 -k common \
-i ${merged_SVs}/merge.asm-2.sort.vcf.gz \
-c ${merged_SVs}/truvari.asm-2.collapse.vcf.gz |\
bcftools sort --max-mem 4G |\
bgzip -@ ${nT} > ${merged_SVs}/truvari.asm-2.vcf.gz

bcftools index -f -t ${merged_SVs}/truvari.asm-2.vcf.gz

module unload python/3.10.2
echo "This is the end!"