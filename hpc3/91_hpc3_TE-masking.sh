#!/bin/bash

#SBATCH --job-name=TE    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=30   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU
#SBATCH --tmp=100G                ## requesting 100 GB local scratch
#SBATCH --constraint=fastscratch  ## requesting nodes with fast scratch in /tmp

dmel_ref="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
dsim_ref="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing/GCF_016746395.2_Dsim_3.1.fasta"
TE="/dfs7/jje/jenyuw/SV-project-temp/result/TE_repeat"
nT=$SLURM_CPUS_PER_TASK
source ~/.bashrc

cp ${dmel_ref} ${TE}/dmel-r6.49.fasta
cd ${TE}

module load singularity/3.11.3
singularity exec -B ${TE} /pub/jenyuw/Software/dfam-tetools-latest.sif RepeatMasker -gff -s -xsmall \
-lib "${TE}/D_mel_transposon_sequence_set.fa"  -dir ${TE}/dmel-r649 \
${TE}/dmel-r6.49.fasta

#-species "drosophila"