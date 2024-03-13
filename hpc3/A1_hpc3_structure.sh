#!/bin/bash

#SBATCH --job-name=STRUCTURE    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard     ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=4   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU

merged_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/Merged_SNP"
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
STRUCTURE="/dfs7/jje/jenyuw/SV-project-temp/result/structure"
nT=$SLURM_CPUS_PER_TASK

source ~/.bashrc

bcftools query -l ${merged_SNP}/all.snps.vcf.gz > ${STRUCTURE}/snp_name.txt

bcftools view -m2 -M2 -v snps ${merged_SNP}/all.snps.vcf.gz|\
bcftools query -f '[ %GT]\n'  | sed 's@1/1@2@g ; s@0/1@1@g ; s@1/0@1@g'|sed -e "s/\.\/\./0/g" > ${STRUCTURE}/snp_temp.txt
#cat ${STRUCTURE}/snp_temp.txt | head -n 100 > ${STRUCTURE}/snp_temp_100.txt

cd ${STRUCTURE}
module load R/4.2.2
Rscript -e '
snp=read.table("snp_temp.txt", header = FALSE)
snp=t(snp)
name=read.table("snp_name.txt", header = FALSE)
t=cbind(name, snp)
write.table(t, "snp_input.txt", row.names = FALSE, col.names = FALSE)
'

Rscript -e '
library(dplyr)
snp=read.table("snp_temp.txt", header = FALSE)
snp=t(snp)
print("This is the first row of snp")
head(snp,1)
name=read.table("snp_name.txt", header = FALSE)
t=bind_cols(name, snp)
print("This is the first row of t")
head(t,1)
write.table(t, "snp_input.txt", row.names = FALSE, col.names = FALSE)
'