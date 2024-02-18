#!/bin/bash

#SBATCH --job-name=extraction    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name: standard or highmem
#SBATCH --array=1     ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=30    ## number of cores the job needs
#SBATCH --mem-per-cpu=6G

#"#SBATCH --mem=200G", "#SBATCH --nodes=1", "#SBATCH --ntasks=1"

## path
dmel_ref="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
dsim_ref="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing/GCF_016746395.2_Dsim_3.1.fasta"
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
merged_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/merged_SVs"
con_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/consensus_SVs"
processed="/dfs7/jje/jenyuw/SV-project-temp/result/processed_SNP"

nT=$SLURM_CPUS_PER_TASK
source ~/.bashrc

bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\t[ %GT] \n' ${polarizing}/polarized.asm.vcf.gz >${polarizing}/extraction.polarized.asm.tsv


bcftools query -f '%CHROM\t%POS\t SNP \t 1 \t[ %GT] \n' ${processed}/synSNPs.vcf.gz >${processed}/extraction.syn-snp.tsv

bcftools query -f '%CHROM\t%POS\t SNP \t 1 \t[ %GT] \n' ${processed}/syn.outside-1000-svimasm.snps.vcf.gz \
>${processed}/extraction.syn-snp.outside-1000-svimasm.tsv