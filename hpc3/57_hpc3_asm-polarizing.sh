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

#Avoid collapsing the SVs from samples

ls ${con_SVs}/*.tru_con.sort.vcf >${polarizing}/sample_files.txt
ls ${polarizing}/mm2.vcf >>${polarizing}/sample_files.txt

SURVIVOR merge ${polarizing}/sample_files.txt 0 1 1 1 0 50 ${polarizing}/allandsim.merged.vcf
bgzip -k -f -@ ${nT} ${polarizing}/allandsim.merged.vcf

#200G temp space
bcftools sort --max-mem 2G ${polarizing}/allandsim.merged.vcf.gz |bgzip -@ ${nT} > ${polarizing}/allandsim.merged.sort.vcf.gz

#This sort step is pretty slow and takes much temp space.
bcftools index -f -t ${polarizing}/allandsim.merged.vcf.gz

bcftools merge -m none -O z -o ${polarizing}/allandsim.merged.vcf.gz ${polarizing}/mm2.vcf.gz ${merged_SVs}/truvari.svimASM.vcf.gz
bcftools index -t -f ${polarizing}/allandsim.merged.vcf.gz
module load python/3.10.2
truvari collapse -k common --sizemax 200000000 \
-i ${polarizing}/allandsim.merged.vcf.gz -f ${dmel_ref} -o /dev/null -c ${polarizing}/all2sim.collapsed.vcf 
module unload python/3.10.2
bgzip -k ${polarizing}/all2sim.collapsed.vcf 
bcftools sort -O z -o ${polarizing}/all2sim.collapsed.sort.vcf.gz  ${polarizing}/all2sim.collapsed.vcf.gz
bcftools index -t -f ${polarizing}/all2sim.collapsed.sort.vcf.gz

## Now, extract the overlapping SVs
bcftools isec -c none --regions-overlap pos -O z \
-p /dfs7/jje/jenyuw/SV-project-temp/result/polarizing/temp \
${merged_SVs}/truvari.svimASM.vcf.gz ${polarizing}/all2sim.collapsed.sort.vcf.gz 

bedtools intersect 

-a ${polarizing}/all2sim.collapsed.sort.vcf.gz -b ${merged_SVs}/truvari.svimASM.vcf.gz -header > ${polarizing}/all2sim.collapsed.sort.overlap.vcf