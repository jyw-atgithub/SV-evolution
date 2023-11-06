#!/bin/bash
#SBATCH --job-name=sra3           ## job name
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard               ## partition name
#SBATCH -N 1                      ## run on a single node, cant run across multiple nodes
#SBATCH --ntasks=8                ## CPUs to use as threads in fasterq-dump command
#SBATCH --tmp=100G                ## requesting 100 GB local scratch
#SBATCH --constraint=fastscratch  ## requesting nodes with fast scratch in /tmp

###Running multi-threaded (read multi-CPU) fastq-dump, fasterq-dump or parallel-fastq-dump command 
###directly in any directory on DFS filesystem results in 
###a deadlock and makes a server UNUSABLE for ALL users until the server is rebooted.###

# IMPORTANT: load the latest SRA-tools, earlier versions do not handle temporary disk
module load sra-tools/3.0.0

# TMPDIR is created automatically by SLURM
# change to your temp directory assigned by SLURM to your job
cd $TMPDIR

# here we work on just 2 sequences
for f in {447..448}
do
# generate ID to prefetch, each ID is SRR1196 plus what is contained in $f variable
ID=SRR1196${f}

# prefetch SRA file
prefetch $ID

# convert sra format to fastq format using requested number of threads (slurm tasks)
# temp files are written to fastscratch in $TMPDIR with a 100G limit
fasterq-dump ./$ID/$ID.sra -e $SLURM_NTASKS --temp $TMPDIR --disk-limit-tmp 100G  

# compress resulting fastq files
gzip $ID*fastq
done

# move all results to desired location in DFS, directory must exists
mv *fastq.gz /dfsX/panteater_lab/SRA/results/project1/