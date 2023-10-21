#!/bin/bash

#SBATCH --job-name=purged-busco    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=32   ## number of cores the job needs
#SBATCH --mem-per-cpu=6G     # requesting 6 GB memory per CPU, the max


##pauvre is useful to check the sequencing stats (pauvre stats)

##pacbio_pipeline for ALL NCBI sequences
### on HPC3 !!!!###

purge_dups="/pub/jenyuw/SV-project-temp/result/purge_dups"
busco_out="/pub/jenyuw/SV-project-temp/result/busco_out"
## prep
source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK
echo "cpu number is $SLURM_CPUS_PER_TASK"

##make name list
##Because of the lazy way of generating file list. Always submit --array=1 at the begining and then submit the rest
if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
echo "Yes, ARRAY_TASK_ID=1"
ls ${purge_dups}/*/purged.fa >${busco_out}/namelist.txt
else
echo "No need to list the file again"
fi

n=`cat ${busco_out}/namelist.txt| wc -l`
echo "there are ${n} files"
file=`head -n $SLURM_ARRAY_TASK_ID ${busco_out}/namelist.txt |tail -n 1`
echo "file is ${file}"
name=`echo ${file}|gawk -F "/" '{print $7}'`
echo "name is ${name} "

#busco
conda activate BUSCO

busco -i ${file} -l diptera_odb10 --out_path ${busco_out} -o ${name} -m genome -c ${nT}
