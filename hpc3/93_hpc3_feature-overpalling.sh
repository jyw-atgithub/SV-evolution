#!/bin/bash

#SBATCH --job-name=TE_overlap    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=12   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=3G     # requesting memory per CPU

dmel_ref="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
TE="/dfs7/jje/jenyuw/SV-project-temp/result/TE_repeat"
org_SV="/dfs7/jje/jenyuw/SV-project-temp/result/organized_SVs"
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
pro_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/processed_SNP"
source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK

#transform vcf to bed format
#bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' ${polarizing}/polarized.asm.vcf.gz >${org_SV}/polarized.asm.bed

bedtools intersect -header -a ${polarizing}/polarized.asm.vcf.gz -b ${TE}/dmel-r649/dmel-r6.49.fasta.out.gff |\
tee >(bcftools view -i 'SVTYPE="DEL"' -O z >${org_SV}/DEL.TE-overlap.asm.vcf.gz) |\
tee >(bcftools view -i 'SVTYPE="INS"' -O z >${org_SV}/INS.TE-overlap.asm.vcf.gz) |\
tee >(bcftools view -i 'SVTYPE="DUP"' -O z >${org_SV}/DUP.TE-overlap.asm.vcf.gz) |\
tee >(bcftools view -i 'SVTYPE="INV"' -O z >${org_SV}/INV.TE-overlap.asm.vcf.gz) |\
bcftools view -i 'SVTYPE="TRA"' -O z >${org_SV}/TRA.TE-overlap.asm.vcf.gz

bedtools intersect -header -v -a ${polarizing}/polarized.asm.vcf.gz -b ${TE}/dmel-r649/dmel-r6.49.fasta.out.gff |\
tee >(bcftools view -i 'SVTYPE="DEL"' -O z >${org_SV}/DEL.no-TE.asm.vcf.gz) |\
tee >(bcftools view -i 'SVTYPE="INS"' -O z >${org_SV}/INS.no-TE.asm.vcf.gz) |\
tee >(bcftools view -i 'SVTYPE="DUP"' -O z >${org_SV}/DUP.no-TE.asm.vcf.gz) |\
tee >(bcftools view -i 'SVTYPE="INV"' -O z >${org_SV}/INV.no-TE.asm.vcf.gz) |\
bcftools view -i 'SVTYPE="TRA"' -O z >${org_SV}/TRA.no-TE.asm.vcf.gz

##Remove duplicated SVs
bcftools view -i 'SVTYPE="DUP"' -O v ${polarizing}/polarized.asm.vcf.gz |\
bedtools intersect -header -a - -b ${pro_SNP}/dmel-all-r6.49.exonspans.bed |\
tee >(grep "#">${org_SV}/DUP.exon-overlap.asm.vcf) |\
bcftools sort -O v | grep -v "#" |sort -k1,1 -k2,2n |uniq >>${org_SV}/DUP.exon-overlap.asm.vcf
bgzip -f -k -@ 4 ${org_SV}/DUP.exon-overlap.asm.vcf
bcftools index -f -t ${org_SV}/DUP.exon-overlap.asm.vcf.gz

bcftools view -i 'SVTYPE="DEL"' -O z -o ${org_SV}/DEL.asm.vcf  ${polarizing}/polarized.asm.vcf.gz 