#!/bin/bash

#SBATCH --job-name=qc-tr    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=16   ## number of cores the job needs
#SBATCH --mem-per-cpu=6G     # requesting 6 GB memory per CPU, the max

source ~/.bashrc

lqc="/pub/jenyuw/Software/LongQC-1.2.0c"
#This version of longqc finally works
# compiled minimap2-coverage with gcc module, executed by python3.8 module
# Using Conda versionof Longqc  on HPC3 kept failing
raw="/dfs7/jje/jenyuw/SV-project-temp/raw"
trimmed="/dfs7/jje/jenyuw/SV-project-temp/result/trimmed"
qc_report="/dfs7/jje/jenyuw/SV-project-temp/result/qc_report"
jellyfish=""

nT=$SLURM_CPUS_PER_TASK
echo "cpu number is $SLURM_CPUS_PER_TASK"

##make name list
##Because of the lazy way of generating file list. Always submit --array=1 at the begining and then submit the rest
if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
echo "Yes, ARRAY_TASK_ID=1"
ls ${raw}/*.fastq.gz >${raw}/namelist.txt
else
echo "No need to list the file again"
fi

file=`head -n $SLURM_ARRAY_TASK_ID ${raw}/namelist.txt |tail -n 1`
name=`basename ${file} | sed s/.fastq.gz//`
read_type=`echo ${name} | gawk -F "_" '{print $2}'`


conda activate qc 
porechop_abi -abi --threads $nT -i ${file} -o ${trimmed}/${name}.abi.fastq
cat ${trimmed}/${name}.abi.fastq |chopper -l 520 --headcrop 10 --tailcrop 10 --threads ${nT} |\
pigz -p ${nT} > ${trimmed}/${name}.trimmed.fastq.gz
rm ${trimmed}/${name}.abi.fastq
conda deactivate