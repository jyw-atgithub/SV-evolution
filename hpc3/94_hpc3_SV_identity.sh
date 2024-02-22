#!/bin/bash

#SBATCH --job-name=SV_identity    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p highmem      ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=20   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=10G     # requesting memory per CPU

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
trimmed="/dfs7/jje/jenyuw/SV-project-temp/result/trimmed"
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
TE="/dfs7/jje/jenyuw/SV-project-temp/result/TE_repeat"
source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK

: <<'SKIP'
########################This may cause redudant!!########################
cat ${trimmed}/namelist.txt|while read i
do
name=$(basename ${i} |sed s/".trimmed.fastq.gz"//g)
#TRA always gives error
################## bcftools consensus skips overlapping SVs, so we don't use it ##################
#bcftools consensus -i 'SVTYPE="DEL"' -s ${name}_mumco --fasta-ref ${ref_genome} --haplotype LR ${polarizing}/polarized.asm.vcf.gz >${TE}/${name}.SV.fa
#for DUP and DEL, just extract the REF allele
echo "Processing ${name}"
bcftools query -i 'SVTYPE="DUP" || SVTYPE="DEL"' -f '>%CHROM-%POS-%END-%ID\n%REF\n' \
-s ${name}_mumco ${polarizing}/polarized.asm.vcf.gz >${TE}/${name}.SV.fa
#for INS and INV, just extract the alt allele
bcftools query -i 'SVTYPE="INS" || SVTYPE="INV"' -f '>%CHROM-%POS-%END-%ID\n%ALT\n' \
-s ${name}_mumco ${polarizing}/polarized.asm.vcf.gz >>${TE}/${name}.SV.fa
done
cat ${TE}/*.SV.fa |gawk -F ";" '{print $1 }' >${TE}/SV.fa
#rm ${TE}/*.SV.fa
SKIP
echo "Processing DUP & DEL"
bcftools query -i 'SVTYPE="DUP" || SVTYPE="DEL"' -f '>%CHROM-%POS-%END-%ID\n%REF\n' \
${polarizing}/polarized.asm.vcf.gz >${TE}/SV.DUP-DEL.fa
#for INS and INV, just extract the alt allele
echo "Processing INS & INV"
bcftools query -i 'SVTYPE="INS" || SVTYPE="INV"' -f '>%CHROM-%POS-%END-%ID\n%ALT\n' \
${polarizing}/polarized.asm.vcf.gz >${TE}/SV.INS-INV.fa
cat ${TE}/SV.DUP-DEL.fa ${TE}/SV.INS-INV.fa |gawk -F ";" '{print $1 }' >${TE}/SV.fa

module load singularity/3.11.3
singularity exec -B ${TE} \
-B /dfs7/jje/jenyuw/SV-project-temp/raw/Libraries:/opt/RepeatMasker/Libraries \
/pub/jenyuw/Software/dfam-tetools-latest.sif \
RepeatMasker -gff -s -xsmall \
-species "drosophila melanogaster" -dir ${TE}/extracted_SV \
${TE}/SV.fa