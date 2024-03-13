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

nT=$SLURM_CPUS_PER_TASK

source ~/.bashrc

#transform vcf to bed format
#bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' ${polarizing}/polarized.asm.vcf.gz >${org_SV}/polarized.asm.bed

###1. We want to see what common SVs are overlapping with exons

bcftools view -O v ${polarizing}/3corrected.polarized.asm.vcf.gz |\
bedtools intersect -header -a - -b ${pro_SNP}/dmel-all-r6.49.exonspans.bed |\
bcftools +missing2ref |bcftools +fill-tags -- -t AF |bcftools filter -i 'AF>0.8' |\
bcftools sort -O v |\
##Remove duplicated SVs
bcftools norm --rm-dup all -O z -o ${org_SV}/high_AF_exon-overlap.asm.vcf.gz

bcftools index -f -t ${org_SV}/high_AF_exon-overlap.asm.vcf.gz
bcftools query -f '%SVTYPE\n' ${org_SV}/high_AF_exon-overlap.asm.vcf.gz |sort|uniq -c
##The gff output is too long
#bedtools intersect -header -a ${ref_gff} -b ${org_SV}/high_AF_exon-overlap.asm.vcf.gz >${org_SV}/high_AF_exon-overlap.asm.gff


###2. We want to find the common SVs in different populations

for i in AF AM AS EU
do
cat ${org_SV}/popfile.tsv|gawk -v i=${i} ' $2==i {print $1}' >${org_SV}/temp_pop.txt
bcftools view -S ${org_SV}/temp_pop.txt -O v ${polarizing}/3corrected.polarized.asm.vcf.gz |\
bedtools intersect -header -a - -b ${pro_SNP}/dmel-all-r6.49.CDSspans.bed |\
bcftools +missing2ref |bcftools +fill-tags -- -t AF |bcftools filter -i 'AF>0.8' |\
bcftools sort -O v |\
bcftools norm --rm-dup all -O z -o ${org_SV}/high_AF.${i}.CDS.asm.vcf.gz
bcftools index -f -t ${org_SV}/high_AF.${i}.CDS.asm.vcf.gz
done
for i in AF AM AS EU
do
cat ${org_SV}/popfile.tsv|gawk -v i=${i} ' $2==i {print $1}' >${org_SV}/temp_pop.txt
bcftools view -S ${org_SV}/temp_pop.txt -O v ${polarizing}/3corrected.polarized.asm.vcf.gz |\
bedtools intersect -header -a - -b ${pro_SNP}/dmel-all-r6.49.exonspans.bed |\
bcftools +missing2ref |bcftools +fill-tags -- -t AF |bcftools filter -i 'AF>0.8' |\
bcftools sort -O v |\
bcftools norm --rm-dup all -O z -o ${org_SV}/high_AF.${i}.exon.asm.vcf.gz
bcftools index -f -t ${org_SV}/high_AF.${i}.exon.asm.vcf.gz
done

rm ${org_SV}/temp_pop.txt

bedtools intersect -header -a high_AF.AF.CDS.asm.vcf.gz -b high_AF.AM.CDS.asm.vcf.gz |\
bcftools norm --rm-dup all | bgzip -@ ${nT} >${org_SV}/AF_mi_AM.CDS.asm.vcf.gz

bedtools intersect -header -a high_AF.AF.CDS.asm.vcf.gz -b high_AF.AM.CDS.asm.vcf.gz \
-b high_AF.AS.CDS.asm.vcf.gz -b high_AF.EU.CDS.asm.vcf.gz |\
bcftools norm --rm-dup all | bgzip -@ ${nT} >${org_SV}/AF_mi_other.CDS.asm.vcf.gz

bedtools intersect -header -a high_AF.AF.exon.asm.vcf.gz -b high_AF.AM.exon.asm.vcf.gz \
-b high_AF.AS.exon.asm.vcf.gz -b high_AF.EU.exon.asm.vcf.gz |\
bcftools norm --rm-dup all | bgzip -@ ${nT} >${org_SV}/AF_mi_other.exon.asm.vcf.gz
bcftools index -f -t ${org_SV}/AF_mi_other.exon.asm.vcf.gz
#Then we check each entry in IGV

