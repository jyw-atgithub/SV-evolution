#!/bin/bash

#SBATCH --job-name=svim-asm    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=16   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU

##############################################################################
###This ia a crude version. We need to feed descent assemblies afterwards. ###
##############################################################################
source ~/.bashrc

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
assemble="/dfs7/jje/jenyuw/SV-project-temp/result/assemble"
scaffold="/dfs7/jje/jenyuw/SV-project-temp/result/scaffold"
nT=$SLURM_CPUS_PER_TASK



file=`head -n $SLURM_ARRAY_TASK_ID ${scaffold}/scfd_list.txt |tail -n 1`
name=$(echo ${file} | cut -d '/' -f 8)

module load python/3.10.2

svim-asm haploid --sample ${name}_svimASM --min_sv_size 50 \
${SVs}/${name}_svimASM ${aligned_bam}/${name}.pat1-ref.sort.bam ${ref_genome}

syri  ${aligned_bam}/${name}.pat1-ref.sort.bam -r ${ref_genome} -q ${file} -F B --nc 4 --samplename ${name}_syri 
[--dir DIR] [--prefix PREFIX] [--nc NCORES] [--novcf]
            [] [--nosr] [--invgaplen INVGL] [--tdgaplen TDGL] [--tdmaxolp TDOLP] [-b BRUTERUNTIME] [--unic TRANSUNICOUNT]
            [--unip TRANSUNIPERCENT] [--inc INCREASEBY] [--no-chrmatch] [--nosv] [--nosnp] [--all] [--allow-offset OFFSET] [--cigar] [-s SSPATH] [--hdrseq]

mv ${SVs}/${name}_svimASM/variants.vcf ${SVs}/${name}.svimASM.vcf 

module unload python/3.10.2
