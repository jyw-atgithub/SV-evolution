#!/bin/bash

#SBATCH --job-name=DSPR-flye-assembly    ## Name of the job.
#SBATCH -A ecoevo283       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=1-6   ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=10    ## number of cores the job needs 


##pauvre is useful to check the sequencing stats (pauvre stats)

##pacbio_pipeline for ALL NCBI sequences
### on HPC3 !!!!###
wd="/pub/jenyuw/SV-project-temp/"
qc_report="/pub/jenyuw/SV-project-temp/qc_report"
## prep
source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK
echo "cpu number is $SLURM_CPUS_PER_TASK"

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
echo "yes"
ls /pub/jenyuw/SV-project-temp/*_combine_pacbio.fastq.gz >/pub/jenyuw/SV-project-temp/pacbio_list.txt
else
echo "no need to list the file again"
fi

file=$(head -n $SLURM_ARRAY_TASK_ID /pub/jenyuw/SV-project-temp/pacbio_list.txt|tail -n 1)
echo $file
name=$(basename ${file}|sed s/"_combine_pacbio.fastq.gz"//g)
echo $name

#LongQC report
#conda activate longqc
#python /pub/jenyuw/Software/LongQC-1.2.0c/longQC.py sampleqc -p ${nT} -x pb-sequel -n 6000 \
#-o ${qc_report}/${name}_longQC -s ${name} ${file}

conda activate assemble

# only chopper is enough for PacBio
#unpigz -p ${nT} -c ${file} | chopper -l 500 --headcrop 10 --tailcrop 10 -q 7 --threads ${nT} > ${wd}/${name}.trimmed.fastq 

flye --threads $nT --genome-size 170m --pacbio-raw ${wd}/${name}.trimmed.fastq --out-dir ${wd}/${name}_Flye
#flye does NOt accept stdin (/dev/stdin)
#unpigz -p 4 -c A2_combine_pacbio.fastq.gz|head -n 10000 |grep "@" | sort | uniq -d

flye --threads 6 --genome-size 170m --pacbio-raw B1-2_combine_pacbio.fastq.gz --out-dir ./B1_Flye