#!/bin/bash

#SBATCH --job-name=map    ## Name of the job.
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

## In interactive mode
cd /dfs7/jje/jenyuw/SV-project-temp/result/polarizing

nucmer -t ${nT} --maxmatch -l 100 -c 500 --sam-long=${polarizing}/mumm4.sam ${dmel_ref} ${dsim_ref}
cat ${polarizing}/mumm4.sam | sed 's/HD\ /HD/; s/1.0\ /1.0/; s/\tSO:coordinate/SO:coordinate/; s/VN1/VN:1/; s/HD/HD\t/; s/SO:unsorted/\tSO:unsorted/; s/@PG /@PG\t/; s/ PN/\tPN/; s/ VN/\tVN/; s/ CL/\tCL/' |\
samtools view -@ ${nT} -h --reference ${dmel_ref} - |samtools sort -@ ${nT} -O bam -o ${polarizing}/corrected.mumm4.bam
samtools index -@ ${nT} ${polarizing}/corrected.mumm4.bam

#nucmer -t ${nT} --maxmatch --delta=${polarizing}/mumm4.delta ${dmel_ref} ${dsim_ref}
#svmu ${polarizing}/mumm4.delta ${dmel_ref} ${dsim_ref} l sam_lastz.txt prefix
#Does not work well

bash /pub/jenyuw/Software/MUMandCo-MUMandCov3.8/mumandco_v3.8.sh  \
-r ${dmel_ref} -q ${dsim_ref} \
-g 135000000 -o mumco -t ${nT} -ml 50

## minimap2-asm20
minimap2 -a -x asm20 -t ${nT} ${dmel_ref} ${dsim_ref} |\
samtools sort -@ ${nT} -O bam -o ${polarizing}/mm2.bam
samtools index -@ ${nT} ${polarizing}/mm2.bam
## minimap2-asm10
minimap2 -a -x asm10 -t ${nT} ${dmel_ref} ${dsim_ref} |\
samtools sort -@ ${nT} -O bam -o ${polarizing}/mm2-10.bam
samtools index -@ ${nT} ${polarizing}/mm2-10.bam
## minimap2-asm5
minimap2 -a -x asm10 -t ${nT} ${dmel_ref} ${dsim_ref} |\
samtools sort -@ ${nT} -O bam -o ${polarizing}/mm2-5.bam
samtools index -@ ${nT} ${polarizing}/mm2-5.bam


conda activate sv-calling
lra global -CONTIG ${dmel_ref}
lra align -CONTIG --refineBreakpoints ${dmel_ref} ${dsim_ref} -t 16 -p s |
samtools sort -@ ${nT} -O bam -o ${polarizing}/lra.bam
samtools index -@ ${nT} ${polarizing}/lra.bam
conda deactivate


conda activate sv-calling
svim-asm haploid --sample "sim2mel_mm2_svimASM" --min_sv_size 50 \
${polarizing}/mm2_svimASM ${polarizing}/mm2.bam ${dmel_ref}

svim-asm haploid --sample "sim2mel_mm2-10_svimASM" --min_sv_size 50 \
${polarizing}/mm2-10_svimASM ${polarizing}/mm2-10.bam ${dmel_ref}

svim-asm haploid --sample "sim2mel_mm2-5_svimASM" --min_sv_size 50 \
${polarizing}/mm2-5_svimASM ${polarizing}/mm2-5.bam ${dmel_ref}

svim-asm haploid --sample "sim2mel_lra_svimASM" --min_sv_size 50 \
${polarizing}/lra_svimASM ${polarizing}/lra.bam ${dmel_ref}

svim-asm haploid --sample "sim2mel_mumm4_svimASM" --min_sv_size 50 \
${polarizing}/mumm4_svimASM ${polarizing}/corrected.mumm4.bam ${dmel_ref}

mv ${SVs}/${name}_svimASM/variants.vcf ${SVs}/${name}.svimASM.vcf 
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
bcftools filter -i 'GT[@filter.txt]="hom"' ${polarizing}/all2sim.asm.collapsed.vcf.gz |\
bcftools sort -O z >  ${polarizing}/all2sim.asm.onlysim.vcf.gz
bcftools index -t -f ${polarizing}/all2sim.asm.onlysim.vcf.gz
# --> 601 SVs

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
$bcftools query -f '%ID\n' ${polarizing}/all2sim.asm.onlysim.vcf.gz >${polarizing}/onlysim.id
cat ${polarizing}/onlysim.id |while read line
do
grep -E "${line}" ${merged_SVs}/truvari.svimASM.vcf |\
sed s@'\.\/\.'@'hahaha'@g |sed s@'1\/1'@'\.\/\.'@g|sed s@'hahaha'@'1\/1'@g >>${polarizing}/reversed.txt
grep -v "${line}" ${merged_SVs}/truvari.svimASM.vcf >>${polarizing}/nochange.txt
done
grep "#" ${merged_SVs}/truvari.svimASM.vcf >${polarizing}/header.txt
cat ${polarizing}/header.txt ${polarizing}/nochange.txt ${polarizing}/reversed.txt >${polarizing}/polarized.vcf
bcftools sort --max-mem 2G ${polarizing}/polarized.vcf.gz -O z -o polarized.sort.vcf.gz
