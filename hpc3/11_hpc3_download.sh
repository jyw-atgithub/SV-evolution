#!/bin/bash

#SBATCH --job-name=d    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=2-32      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=16   ## number of cores the job needs
#SBATCH --mem-per-cpu=6G     # requesting 6 GB memory per CPU, the max



line=`head -n $SLURM_ARRAY_TASK_ID download_list.csv |tail -n 1`
echo $line
strain=`echo $line |gawk -F "," '{print $1}'`
tech=`echo $line |gawk -F "," '{print $2}'`
echo $line
echo $strain
echo $tech is tech

module load sra-tools/3.0.0
prefetch -pcv ${strain}
wait
fastq-dump --split-spot --stdout  /dfs7/jje/jenyuw/SV-project-temp/raw/${strain}/*.sra |\
pigz -p 12 -v  >/dfs7/jje/jenyuw/SV-project-temp/raw/${strain}_${tech}.fastq.gz

###Running multi-threaded (read multi-CPU) fastq-dump, fasterq-dump or parallel-fastq-dump command 
###directly in any directory on DFS filesystem results in 
###a deadlock and makes a server UNUSABLE for ALL users until the server is rebooted.###

#!/bin/bash
#SBATCH --job-name=sra3           ## job name
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard               ## partition name
#SBATCH -N 1                      ## run on a single node, cant run across multiple nodes
#SBATCH --ntasks=8                ## CPUs to use as threads in fasterq-dump command
#SBATCH --tmp=100G                ## requesting 100 GB local scratch
#SBATCH --constraint=fastscratch  ## requesting nodes with fast scratch in /tmp

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