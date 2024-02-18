#!/bin/bash

#SBATCH --job-name=calling    ## Name of the job.
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

##cuteSV
minimap2 -t ${nT} -a -x map-pb ${dmel_ref} ${polarizing}/dsim_pcbio.fastq.gz |\
samtools view -b -h -@ ${nT} -o - |\
samtools sort -@ ${nT} -o ${polarizing}/dsim_read-dmel.sort.bam
samtools index -@ ${nT} ${polarizing}/dsim_read-dmel.sort.bam
module load python/3.10.2
cuteSV --threads ${nT} --genotype --sample "sim2mel_cute" \
--min_support 10 \
--min_size 50 --min_mapq 20 --min_read_len 500 \
-L '-1' \
--merge_del_threshold 270 --merge_ins_threshold 270 \
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5 \
"${polarizing}/dsim_read-dmel.sort.bam" "${dmel_ref}" "dsim-dmel.cutesv.vcf" "${polarizing}"
module unload python/3.10.2

bgzip -k -@ ${nT} ${polarizing}/dsim-dmel.cutesv.vcf
bcftools index -t -f ${polarizing}/dsim-dmel.cutesv.vcf.gz
bcftools view --threads ${nT} -r 2L,2R,3L,3R,4,X,Y \
-i 'QUAL >= 20 && FILTER = "PASS"'  -O v -o - ${polarizing}/dsim-dmel.cutesv.vcf.gz |\
sed 's/DUP_TANDEM/DUP/g; s/DUP:TANDEM/DUP/g; s/DUP_INT/DUP/g; s/DUP:INT/DUP/g' |\
grep -v "NNNNNNNNNNNNNNNNNNNN" |bcftools sort --max-mem 4G -O z -o ${polarizing}/dsim-dmel.cutesv.filtered.vcf.gz

## minimap2-asm20
minimap2 -a -x asm20 -t ${nT} ${dmel_ref} ${dsim_ref} |\
samtools sort -@ ${nT} -O bam -o ${polarizing}/mm2-20.bam
samtools index -@ ${nT} ${polarizing}/mm2-20.bam

conda activate sv-calling
svim-asm haploid --sample "sim2mel_mm2_svimASM" --min_sv_size 50 \
${polarizing}/mm2-20_svimASM ${polarizing}/mm2-20.bam ${dmel_ref}
conda deactivate

cd ${polarizing}
#bash /pub/jenyuw/Software/MUMandCo-MUMandCov3.8/mumandco_v3.8.sh  \
#-r ${dmel_ref} -q ${dsim_ref} \
#-g 135000000 -o dmel-dsim_mumco -t ${nT} -ml 50

#Change the genome size to force MumandCo to use "nucmer --madmatch"
bash /pub/jenyuw/Software/MUMandCo-MUMandCov3.8/mumandco_v3.8.sh  \
-r ${dmel_ref} -q ${dsim_ref} \
-g 1000 -o dmel-dsim_mumco2 -t ${nT} -ml 50

cd ${polarizing}/dmel-dsim_mumco



: <<'SKIP'
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
svim-asm haploid --sample "sim2mel_lra_svimASM" --min_sv_size 50 \
${polarizing}/lra_svimASM ${polarizing}/lra.bam ${dmel_ref}

svim-asm haploid --sample "sim2mel_mumm4_svimASM" --min_sv_size 50 \
${polarizing}/mumm4_svimASM ${polarizing}/corrected.mumm4.bam ${dmel_ref}
SKIP

: <<'SKIP'
svim-asm haploid --sample "sim2mel_mm2-10_svimASM" --min_sv_size 50 \
${polarizing}/mm2-10_svimASM ${polarizing}/mm2-10.bam ${dmel_ref}

svim-asm haploid --sample "sim2mel_mm2-5_svimASM" --min_sv_size 50 \
${polarizing}/mm2-5_svimASM ${polarizing}/mm2-5.bam ${dmel_ref}
SKIP

## We choose the SV called by minimap2-asm20 and svim-asm
#cp ./mm2_svimASM/variants.vcf mm2.vcf
#bgzip -@ -k mm2.vcf
#bcftools index -t -f mm2.vcf.gz
## create the collapsed (overlapping) SVs between the population (all) and D. simulans
#
#Avoid collapsing the SVs from samples TWICE

