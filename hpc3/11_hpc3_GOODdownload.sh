#!/bin/bash
#SBATCH --job-name=sra3           ## job name
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard               ## partition name
#SBATCH -N 1                      ## run on a single node, cant run across multiple nodes
#SBATCH --ntasks=12                ## CPUs to use as threads in fasterq-dump command
#SBATCH --tmp=300G                ## requesting 100 GB local scratch
#SBATCH --constraint=fastscratch  ## requesting nodes with fast scratch in /tmp


##################################################################################
###############Long read data are NOT conpatible with fasterq-dump################
##################################################################################


###Running multi-threaded (read multi-CPU) fastq-dump, fasterq-dump or parallel-fastq-dump command 
###directly in any directory on DFS filesystem results in 
###a deadlock and makes a server UNUSABLE for ALL users until the server is rebooted.###

# IMPORTANT: load the latest SRA-tools, earlier versions do not handle temporary disk
module load sra-tools/3.0.0

# TMPDIR is created automatically by SLURM
# change to your temp directory assigned by SLURM to your job
cd $TMPDIR
pwd
# here we work on just 2 sequences
for f in {37..47}
do
# generate ID to prefetch, each ID is SRR1196 plus what is contained in $f variable
ID=SRR96099${f}

# prefetch SRA file
prefetch $ID
#mv ${TMPDIR}/$ID/$ID.sra /dfs7/jje/jenyuw/SV-project-temp/result/polarizing
# convert sra format to fastq format using requested number of threads (slurm tasks)
# temp files are written to fastscratch in $TMPDIR with a 100G limit
#fasterq-dump is newer than fastq-dump
fasterq-dump /dfs7/jje/jenyuw/SV-project-temp/result/polarizing/${ID}.sra  --split-spot \
-e $SLURM_NTASKS --temp $TMPDIR --disk-limit-tmp 300G --progress

#pigz -p 12 -v  >$TMPDIR/${ID}.fastq.gz
gzip ${TMPDIR}/${ID}*fastq

# compress resulting fastq files
# move  results to desired location in DFS, directory must exists
mv ${TMPDIR}/${ID}.fastq.gz /dfs7/jje/jenyuw/SV-project-temp/result/polarizing
done

for f in {37..47}
do
ID=SRR96099${f}
fastq-dump ${ID}.sra --split-spot --stdout |\
bgzip -@ 12 -c  >${ID}.fastq.gz
done