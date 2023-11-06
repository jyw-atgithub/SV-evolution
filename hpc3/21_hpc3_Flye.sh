#!/bin/bash

#SBATCH --job-name=flye    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p highmem        ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=36   ## number of cores the job needs
#SBATCH --mem-per-cpu=10G     # requesting 6 GB memory per CPU, the max

trimmed="/dfs7/jje/jenyuw/SV-project-temp/result/trimmed"
assemble="/dfs7/jje/jenyuw/SV-project-temp/result/assemble"
source ~/.bashrc

nT=$SLURM_CPUS_PER_TASK
echo "cpu number is $SLURM_CPUS_PER_TASK"

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
echo "Yes, ARRAY_TASK_ID=1"
ls ${trimmed}/*.fastq.gz >${trimmed}/namelist.txt
else
echo "No need to list the file again"
fi

file=`head -n $SLURM_ARRAY_TASK_ID ${trimmed}/namelist.txt |tail -n 1`
name=$(basename ${file}|sed s/".trimmed.fastq.gz"//g)
read_type=`echo ${name} | gawk -F "_" '{print $2}'`

preset_option=(["CLR"]="--pacbio-raw" ["ONT"]="--nano-raw")

conda activate assemble

flye --threads $nT --genome-size 135m ${preset_option[read_type]} ${file} --out-dir ${assemble}/${name}_Flye

conda deactivate