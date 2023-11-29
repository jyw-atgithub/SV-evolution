#!/bin/bash

#SBATCH --job-name=scaf    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=2-61     ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=12   ## number of cores the job needs
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU

##############################################################################
###This ia a crude version. We need to feed descent assemblies afterwards. ###
##############################################################################

purge_dups="/dfs7/jje/jenyuw/SV-project-temp/result/purge_dups"
ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
scaffold="/dfs7/jje/jenyuw/SV-project-temp/result/scaffold"
patched="/dfs7/jje/jenyuw/SV-project-temp/result/patched"

source ~/.bashrc

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
ls ${patched}/*_1/ragtag.patch.fasta > ${patched}/patched_1_list.txt
fi

file=`head -n $SLURM_ARRAY_TASK_ID ${patched}/patched_1_list.txt |tail -n 1`
name=$(echo ${file} | cut -d '/' -f 8 |cut -d '_' -f 1-2)
echo "the file is ${file}"
echo "the name is ${name}"

conda activate ragtag

ragtag.py scaffold -r -w -u --aligner 'nucmer' -o ${scaffold}/${name} ${ref_genome} ${file}
##Nnucmer canNOT be multithreaded.
conda deactivate

