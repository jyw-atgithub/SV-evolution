#!/bin/bash

#SBATCH --job-name=busco    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=24   ## number of cores the job needs
#SBATCH --mem-per-cpu=6G     # requesting 6 GB memory per CPU, the max


##pauvre is useful to check the sequencing stats (pauvre stats)

##pacbio_pipeline for ALL NCBI sequences
### on HPC3 !!!!###
trimmed="/dfs7/jje/jenyuw/SV-project-temp/result/trimmed"
assemble="/dfs7/jje/jenyuw/SV-project-temp/result/assemble"
purge_dups="/dfs7/jje/jenyuw/SV-project-temp/result/purge_dups"
patched="/dfs7/jje/jenyuw/SV-project-temp/result/patched"
busco_out="/dfs7/jje/jenyuw/SV-project-temp/result/busco_out"

## prep
source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK
echo "cpu number is $SLURM_CPUS_PER_TASK"

##make name list
##Because of the lazy way of generating file list. Always submit --array=1 at the begining and then submit the rest

n=`cat ${trimmed}/namelist.txt| wc -l`
echo "there are ${n} strains"
file=`head -n $SLURM_ARRAY_TASK_ID ${trimmed}/namelist.txt |tail -n 1`
echo "file is ${file}"
name=`echo ${file}|gawk -F "/" '{print $8}'|cut -d "." -f 1`
echo "name is ${name}"

#busco
conda activate BUSCO

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
printf "" >${busco_out}/busco_missing.txt
fi

assembler="canu flye nextdenovo-45"
for i in `echo $assembler`
do
    if [[ $i == "flye" ]]
    then 
        if [[ -f "${assemble}/${name}_${i}/assembly.fasta" ]]
        then
        echo "${assemble}/${name}_${i}/assembly.fasta"
        this="${assemble}/${name}_${i}/assembly.fasta"
        busco -i ${this} -l diptera_odb10 --out_path ${busco_out} -o ${name}_${i}_ori -m genome -c ${nT}
        else
        this="${assemble}/${name}_${i}/${name}.contigs.fasta"
        echo "${this}" " does not exist"
        echo "${this}" >> ${busco_out}/busco_missing.txt
        fi
    fi 

    if [[ $i == "canu" ]]
    then
        if [[ -f "${assemble}/${name}_${i}/${name}.contigs.fasta" ]]
        then
        echo "${assemble}/${name}_${i}/${name}.contigs.fasta"
        this="${assemble}/${name}_${i}/${name}.contigs.fasta"
        busco -i ${this} -l diptera_odb10 --out_path ${busco_out} -o ${name}_${i}_ori -m genome -c ${nT}
        else
        this="${assemble}/${name}_${i}/${name}.contigs.fasta"
        echo "${this}" " does not exist"
        echo "${this}" >> ${busco_out}/busco_missing.txt
        fi
    fi

    if [[ $i == "nextdenovo-45" ]]
    then
        if [[ -f "${assemble}/${name}_${i}/03.ctg_graph/nd.asm.fasta" ]]
        then
        echo "${assemble}/${name}_${i}/03.ctg_graph/nd.asm.fasta"
        this="${assemble}/${name}_${i}/03.ctg_graph/nd.asm.fasta"
        busco -i ${this} -l diptera_odb10 --out_path ${busco_out} -o ${name}_${i}_ori -m genome -c ${nT}
        else
        this="${assemble}/${name}_${i}/${name}.contigs.fasta"
        echo "${this}" " does not exist"
        echo "${this}" >> ${busco_out}/busco_missing.txt
        fi
    fi
done

conda deactivate

conda activate BUSCO

conda deactivate