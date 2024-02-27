#!/bin/bash

#SBATCH --job-name=TE_overlap    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=12   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=3G     # requesting memory per CPU

dmel_ref="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
ref_gff="/dfs7/jje/jenyuw/SV-project-temp/reference/r649.corrected.gff"
TE="/dfs7/jje/jenyuw/SV-project-temp/result/TE_repeat"
org_SV="/dfs7/jje/jenyuw/SV-project-temp/result/organized_SVs"
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
pro_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/processed_SNP"
source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK

#transform vcf to bed format
#bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' ${polarizing}/polarized.asm.vcf.gz >${org_SV}/polarized.asm.bed


##Remove duplicated SVs
bcftools view -O v ${polarizing}/3corrected.polarized.asm.vcf.gz |\
bedtools intersect -header -a - -b ${pro_SNP}/dmel-all-r6.49.exonspans.bed |\
tee >(grep "#">${org_SV}/DUP.exon-overlap.asm.vcf) |\
bcftools +missing2ref |bcftools +fill-tags -- -t AF |bcftools filter -i 'AF>0.8' |\
bcftools sort -O v |\
##Remove duplicated SVs
bcftools norm --rm-dup all -O z -o ${org_SV}/high_AF_exon-overlap.asm.vcf.gz

bcftools index -f -t ${org_SV}/high_AF_exon-overlap.asm.vcf.gz
bcftools query -f '%SVTYPE\n' ${org_SV}/high_AF_exon-overlap.asm.vcf.gz |sort|uniq -c

bedtools intersect -header -a ${ref_gff} -b ${org_SV}/high_AF_exon-overlap.asm.vcf.gz >${org_SV}/high_AF_exon-overlap.asm.gff