#!/bin/bash

#SBATCH --job-name=pol    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=24   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU
#SBATCH --tmp=100G                ## requesting 100 GB local scratch
#SBATCH --constraint=fastscratch  ## requesting nodes with fast scratch in /tmp

dmel_ref="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
dsim_ref="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing/GCF_016746395.2_Dsim_3.1.fasta"
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
merged_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/merged_SVs"
con_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/consensus_SVs"
nT=$SLURM_CPUS_PER_TASK
source ~/.bashrc

: <<'SKIP'
########################################################################################################################
##Mapping-based 
ls ${con_SVs}/*.tru_con.sort.vcf >${polarizing}/sample_files.txt
ls ${polarizing}/mm2.vcf >>${polarizing}/sample_files.txt
SURVIVOR merge ${polarizing}/sample_files.txt 0 1 1 1 0 50 ${polarizing}/allandsim.merged.vcf
bgzip -k -f -@ ${nT} ${polarizing}/allandsim.merged.vcf
#it requires temporary space larger than 800GB
#This sort step is pretty slow
bcftools sort --max-mem 2G ${polarizing}/allandsim.merged.vcf.gz |bgzip -@ ${nT} > ${polarizing}/allandsim.merged.sort.vcf.gz
bcftools index -f -t ${polarizing}/allandsim.merged.sort.vcf.gz

module load python/3.10.2
#Failed
#truvari collapse --intra -k common --median-info --sizemax 200000000 \
#-i ${polarizing}/allandsim.merged.sort.vcf.gz -f ${dmel_ref} -o /dev/null -c ${polarizing}/all2sim.INTRA-collapsed.vcf 
truvari collapse -k common --sizemax 200000000 \
-i ${polarizing}/allandsim.merged.sort.vcf.gz -f ${dmel_ref} -o /dev/null -c ${polarizing}/all2sim.collapsed.vcf 
module unload python/3.10.2
bgzip -k ${polarizing}/all2sim.collapsed.vcf 
bcftools sort -O z -o ${polarizing}/all2sim.collapsed.sort.vcf.gz  ${polarizing}/all2sim.collapsed.vcf.gz
bcftools index --threads ${nT} -t -f ${polarizing}/all2sim.collapsed.sort.vcf.gz
########################################################################################################################
SKIP

##Assembly-based
ls ${con_SVs}/*.asm-2.tru_con.sort.vcf.gz >${polarizing}/conSV_files.txt
ls ${polarizing}/dsim.final.vcf.gz >>${polarizing}/conSV_files.txt
bcftools merge -m none --threads ${nT} `cat ${polarizing}/conSV_files.txt` | bgzip -@ ${nT} > ${polarizing}/allandsim.asm.merged.vcf.gz
bcftools sort --max-mem 4G ${polarizing}/allandsim.asm.merged.vcf.gz |bgzip -@ ${nT} > ${polarizing}/allandsim.asm.sort.vcf.gz
bcftools index --threads ${nT} -t -f ${polarizing}/allandsim.asm.sort.vcf.gz

###First way to extract overlapping SVs
module load python/3.10.2
truvari collapse -k common --sizemax 5000000 --sizemin 50 --refdist 700 --pctseq 0.85 --pctsize 0.85 \
-i ${polarizing}/allandsim.asm.sort.vcf.gz -o ${polarizing}/allandsim.asm.truvari.vcf.gz \
-c ${polarizing}/all2sim.asm.collapsed.vcf
module unload python/3.10.2
printf "dsim_mumco2">${polarizing}/filter.txt
#bcftools filter -i 'GT[@filter.txt]="hom"' ${polarizing}/all2sim.asm.collapsed.vcf.gz |\
#bcftools sort -O z >  ${polarizing}/all2sim.asm.onlysim.vcf.gz
#bcftools index -t -f ${polarizing}/all2sim.asm.onlysim.vcf.gz
# --> 601 SVs
bcftools filter --threads ${nT} -i 'GT[@filter.txt]="hom"' ${polarizing}/allandsim.asm.truvari.vcf.gz |\
#bcftools filter --threads ${nT} -i 'GT="hom"' ${polarizing}/allandsim.asm.truvari.vcf.gz
bcftools view --threads ${nT} -i 'NumCollapsed >= 1'|\
bcftools sort -O z >  ${polarizing}/all2sim.asm.onlysim.vcf.gz
bcftools index -t -f ${polarizing}/all2sim.asm.onlysim.vcf.gz
# --> zcat all2sim.asm.onlysim.vcf.gz |grep -v "#"|wc -l
# --> 2623

#bcftools query -f '%CHROM %POS % \n' {polarizing}/all2sim.asm.onlysim.vcf.gz

: << 'SKIP'
###Second way to extract overlapping SVs
bcftools isec --threads ${nT} -c none --regions-overlap pos --nfiles=2 -O z \
-p /dfs7/jje/jenyuw/SV-project-temp/result/polarizing/bcftools_isect \
${merged_SVs}/truvari.svimASM.vcf.gz ${polarizing}/mm2.vcf.gz
# --> 315 SVs 

###Third way to extract overlapping SVs
bedtools intersect -header -f 0.95 -wb -a ${merged_SVs}/truvari.svimASM.vcf.gz -b ${polarizing}/mm2.vcf.gz
# --> 4910 SVs
##-v??

###Fourth way to extract overlapping SVs
ls ${SVs}/*svimASM.filtered.vcf >${polarizing}/svimasm_vcf.txt
ls ${polarizing}/mm2.vcf >>${polarizing}/svimasm_vcf.txt
SURVIVOR merge ${polarizing}/svimasm_vcf.txt 0.05 2 1 1 0 50 ${polarizing}/allandsim.asm.SURVIVOR.vcf
bcftools filter -i 'GT[@filter.txt]="hom"' ${polarizing}/allandsim.asm.SURVIVOR.vcf > all2sim.asm.SURVIVOR.onlysim.vcf

## Now, extract the overlapping SVs
#bcftools isec -c none --regions-overlap pos -O z \
#-p /dfs7/jje/jenyuw/SV-project-temp/result/polarizing/temp \
#${merged_SVs}/truvari.svimASM.vcf.gz ${polarizing}/all2sim.asm.onlysim.vcf.gz
#this yeiled no overlapping SVs

##
##This strategy is not good because there are overlapping names, which causes many many false positives
##
$bcftools query -f '%ID\n' ${polarizing}/all2sim.asm.onlysim.vcf.gz >${polarizing}/onlysim.id
#cat ${polarizing}/onlysim.id |while read line
#do
#grep -E "${line}" ${merged_SVs}/truvari.svimASM.vcf |\
#sed s@'\.\/\.'@'hahaha'@g |sed s@'1\/1'@'\.\/\.'@g|sed s@'hahaha'@'1\/1'@g >>${polarizing}/reversed.txt
#done
#cat ${polarizing}/reversed.txt|sort|uniq >${polarizing}/reversed.uniq.txt

grep "#" ${merged_SVs}/truvari.svimASM.vcf >${polarizing}/header.txt
#cat  ${polarizing}/header.txt ${polarizing}/reversed.uniq.txt |bgzip -@ ${nT} -c |bcftools sort -O z -o - >${polarizing}/reversed.vcf.gz
#bcftools index -f -t ${polarizing}/reversed.vcf.gz
bedtools intersect -header  -a ${merged_SVs}/truvari.svimASM.vcf.gz -b ${polarizing}/reversed.vcf.gz |\
bgzip -@ ${nT} -c |bcftools sort -O z -o ->${polarizing}/nochange.vcf.gz
bcftools index -f -t ${polarizing}/nochange.vcf.gz
bcftools concat --threads ${nT} -a -O z ${polarizing}/nochange.vcf.gz ${polarizing}/reversed.vcf.gz |\
bcftools sort --max-mem 2G -O z -o polarized.asm.sort.vcf.gz
bgzip -@ ${nT} -d -k polarized.asm.sort.vcf.gz
SKIP

bedtools subtract -A -header -a ${merged_SVs}/truvari.asm-2.vcf.gz -b ${polarizing}/all2sim.asm.onlysim.vcf.gz |\
bcftools sort -O v -o ${polarizing}/nochange.vcf

bedtools subtract -A -a ${merged_SVs}/truvari.asm-2.vcf.gz -b ${polarizing}/nochange.vcf |\
sed s@'\.\/\.'@'hahaha'@g |sed s@'1\/1'@'\.\/\.'@g|sed s@'hahaha'@'1\/1'@g >${polarizing}/reversed.txt

cat ${polarizing}/nochange.vcf ${polarizing}/reversed.txt |bcftools sort --max-mem 2G -O z -o ${polarizing}/polarized.asm.vcf.gz
bcftools index -f -t ${polarizing}/polarized.asm.vcf.gz