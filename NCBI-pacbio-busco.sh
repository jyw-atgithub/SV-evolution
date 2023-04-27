#!/bin/bash

#SBATCH --job-name=DSPR-flye-busco    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=1-14      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=20    ## number of cores the job needs
#SBATCH --mem-per-cpu=4G     # requesting 6 GB memory per CPU, the max


##pauvre is useful to check the sequencing stats (pauvre stats)

##pacbio_pipeline for ALL NCBI sequences
### on HPC3 !!!!###
wd="/pub/jenyuw/SV-project-temp"
qc_report="/pub/jenyuw/SV-project-temp/qc_report"
busco_out=/pub/jenyuw/SV-project-temp/busco
## prep
source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK
echo "cpu number is $SLURM_CPUS_PER_TASK"

file=$(head -n $SLURM_ARRAY_TASK_ID /pub/jenyuw/SV-project-temp/pacbio_list.txt|tail -n 1)
echo $file
name=$(basename ${file}|sed s/"_combine_pacbio.fastq.gz"//g)
echo $name

#busco
conda activate BUSCO

busco -i ${wd}/${name}_Flye/assembly.fasta -l diptera_odb10 -o ${busco_out}/${name} -m genome -c ${nT}