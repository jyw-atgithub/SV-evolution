#!/bin/bash

#SBATCH --job-name=ndenovo    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p highmem        ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=36   ## number of cores the job needs
#SBATCH --mem-per-cpu=10G     # requesting memory per CPU

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
read_type1=`echo ${name} | gawk -F "_" '{print $2}'`
read_type=`echo "${read_type1,,}"` #convert uppercase to lowercase

echo "the file is ${file}"
echo "the read type is ${read_type}"


echo ${file} > ${assemble}/${SLURM_ARRAY_TASK_ID}.input.fofn

echo -e "
job_type = local
job_prefix = nextDenovo
task = all
rewrite = yes
deltmp = yes

parallel_jobs =9 #M gb memory, between M/64~M/32
input_type = raw
read_type = ${read_type} # clr, ont, hifi
input_fofn = ${assemble}/${SLURM_ARRAY_TASK_ID}.input.fofn
workdir = ${assemble}/${name}_nextdenovo-45

[correct_option]
read_cutoff = 1k
genome_size = 135m
seed_depth = 45 #you can try to set it 30-45 to get a better assembly result
seed_cutoff = 0
sort_options = -m 40g -t 4 #m=M/(TOTAL_INPUT_BASES * 1.2/4)
minimap2_options_raw = -t 4
pa_correction = 3 #M/(TOTAL_INPUT_BASES * 1.2/4)
correction_options = -p 4 #P cores, P/parallel_jobs

[assemble_option]
minimap2_options_cns = -t 4 -k17 -w17
minimap2_options_map = -t 4 #P cores, P/parallel_jobs
nextgraph_options = -a 1 -q 10 #usuallu best according to the authors
" >${assemble}/${SLURM_ARRAY_TASK_ID}.run.cfg

nextDenovo ${assemble}/${SLURM_ARRAY_TASK_ID}.run.cfg

echo "It is the end"