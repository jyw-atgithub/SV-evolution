#!/bin/bash

#SBATCH --job-name=pol    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=24   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU

dmel_ref="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
dsim_ref="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing/GCF_016746395.2_Dsim_3.1.fasta"
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
merged_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/merged_SVs"
con_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/consensus_SVs"
nT=$SLURM_CPUS_PER_TASK
source ~/.bashrc
## In interactive mode


: <<'SKIP'
nucmer -t ${nT} --maxmatch -l 100 -c 500 --sam-long=${polarizing}/mumm4.sam ${dmel_ref} ${dsim_ref}
cat ${polarizing}/mumm4.sam |\
sed 's/HD\ /HD/; s/1.0\ /1.0/; s/\tSO:coordinate/SO:coordinate/; s/VN1/VN:1/; s/HD/HD\t/; s/SO:unsorted/\tSO:unsorted/; s/@PG /@PG\t/; s/ PN/\tPN/; s/ VN/\tVN/; s/ CL/\tCL/' |\
sed 's@\t10\t@\t30\t@g'|\
samtools view -@ ${nT} -h --reference ${dmel_ref} - |samtools sort -@ ${nT} -O bam -o ${polarizing}/corrected.mumm4.bam
samtools index -@ ${nT} ${polarizing}/corrected.mumm4.bam

#nucmer -t ${nT} --maxmatch --delta=${polarizing}/mumm4.delta ${dmel_ref} ${dsim_ref}
#conda activate sv-calling
#lastz ${dmel_ref}[multiple] ${dsim_ref}[multiple] --chain --format=general:name1,strand1,start1,end1,name2,strand2,start2,end2 > lastz_out.txt
#conda deactivate
#svmu ${polarizing}/mumm4.delta ${dmel_ref} ${dsim_ref} l sam_lastz.txt svmu > svmu.tsv
#Does not work well
SKIP

cd ${polarizing}

bash /pub/jenyuw/Software/MUMandCo-MUMandCov3.8/mumandco_v3.8.sh  \
-r ${dmel_ref} -q ${dsim_ref} \
-g 135000000 -o dmel-dsim_mumco -t ${nT} -ml 50

cd ${polarizing}/dmel-dsim_mumco

## minimap2-asm20
minimap2 -a -x asm20 -t ${nT} ${dmel_ref} ${dsim_ref} |\
samtools sort -@ ${nT} -O bam -o ${polarizing}/mm2-20.bam
samtools index -@ ${nT} ${polarizing}/mm2-20.bam
## minimap2-asm10
minimap2 -a -x asm10 -t ${nT} ${dmel_ref} ${dsim_ref} |\
samtools sort -@ ${nT} -O bam -o ${polarizing}/mm2-10.bam
samtools index -@ ${nT} ${polarizing}/mm2-10.bam
## minimap2-asm5
minimap2 -a -x asm10 -t ${nT} ${dmel_ref} ${dsim_ref} |\
samtools sort -@ ${nT} -O bam -o ${polarizing}/mm2-5.bam
samtools index -@ ${nT} ${polarizing}/mm2-5.bam

: <<'SKIP'
conda activate sv-calling
lra global -CONTIG ${dmel_ref}
lra align -CONTIG --refineBreakpoints ${dmel_ref} ${dsim_ref} -t 16 -p s |
samtools sort -@ ${nT} -O bam -o ${polarizing}/lra.bam
samtools index -@ ${nT} ${polarizing}/lra.bam
conda deactivate
svim-asm haploid --sample "sim2mel_lra_svimASM" --min_sv_size 50 \
${polarizing}/lra_svimASM ${polarizing}/lra.bam ${dmel_ref}

svim-asm haploid --sample "sim2mel_mumm4_svimASM" --min_sv_size 50 \
${polarizing}/mumm4_svimASM ${polarizing}/corrected.mumm4.bam ${dmel_ref}
SKIP

conda activate sv-calling
svim-asm haploid --sample "sim2mel_mm2_svimASM" --min_sv_size 50 \
${polarizing}/mm2_svimASM ${polarizing}/mm2-20.bam ${dmel_ref}

svim-asm haploid --sample "sim2mel_mm2-10_svimASM" --min_sv_size 50 \
${polarizing}/mm2-10_svimASM ${polarizing}/mm2-10.bam ${dmel_ref}

svim-asm haploid --sample "sim2mel_mm2-5_svimASM" --min_sv_size 50 \
${polarizing}/mm2-5_svimASM ${polarizing}/mm2-5.bam ${dmel_ref}

conda deactivate

## We choose the SV called by minimap2-asm20 and svim-asm
#cp ./mm2_svimASM/variants.vcf mm2.vcf
#bgzip -@ -k mm2.vcf
#bcftools index -t -f mm2.vcf.gz
## create the collapsed (overlapping) SVs between the population (all) and D. simulans
#
#Avoid collapsing the SVs from samples TWICE

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

##Assembly-based
ls ${SVs}/*svimASM.filtered.vcf.gz >${polarizing}/svimasm_files.txt
ls ${polarizing}/mm2.vcf.gz >>${polarizing}/svimasm_files.txt
bcftools merge -m none --threads ${nT} `cat ${polarizing}/svimasm_files.txt` | bgzip -@ ${nT} > ${polarizing}/allandsim.asm.merged.vcf.gz
bcftools sort --max-mem 2G ${polarizing}/allandsim.asm.merged.vcf.gz |bgzip -@ ${nT} > ${polarizing}/allandsim.asm.sort.vcf.gz
bcftools index --threads ${nT}  -t -f ${polarizing}/allandsim.asm.sort.vcf.gz

###First way to extract overlapping SVs
module load python/3.10.2
nohup time truvari collapse -k common --sizemax 200000000 \
-i ${polarizing}/allandsim.asm.sort.vcf.gz -f ${dmel_ref} -o ${polarizing}/allandsim.asm.truvari.vcf.gz \
-c ${polarizing}/all2sim.asm.collapsed.vcf &
module unload python/3.10.2
printf "sim2mel_mm2_svimASM">${polarizing}/filter.txt
#bcftools filter -i 'GT[@filter.txt]="hom"' ${polarizing}/all2sim.asm.collapsed.vcf.gz |\
#bcftools sort -O z >  ${polarizing}/all2sim.asm.onlysim.vcf.gz
#bcftools index -t -f ${polarizing}/all2sim.asm.onlysim.vcf.gz
# --> 601 SVs
bcftools filter --threads ${nT} -i 'GT[@filter.txt]="hom"' ${polarizing}/allandsim.asm.truvari.vcf.gz |\
bcftools view --threads ${nT} -i 'NumCollapsed >= 1'|\
bcftools sort -O z >  ${polarizing}/all2sim.asm.onlysim.vcf.gz
bcftools index -t -f ${polarizing}/all2sim.asm.onlysim.vcf.gz
# --> 834 SVs
#bcftools query -f '%CHROM %POS % \n' {polarizing}/all2sim.asm.onlysim.vcf.gz


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

: << 'SKIP'
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

bedtools subtract -A -header -a ${merged_SVs}/truvari.svimASM.vcf.gz -b ${polarizing}/all2sim.asm.onlysim.vcf.gz |\
bcftools sort -O v -o ${polarizing}/nochange.vcf

bedtools subtract -A -a ${merged_SVs}/truvari.svimASM.vcf.gz -b ${polarizing}/nochange.vcf.gz|\
sed s@'\.\/\.'@'hahaha'@g |sed s@'1\/1'@'\.\/\.'@g|sed s@'hahaha'@'1\/1'@g >${polarizing}/reversed.txt

cat ${polarizing}/nochange.vcf ${polarizing}/reversed.txt |bcftools sort --max-mem 2G -O z -o ${polarizing}/polarized.asm.vcf.gz
