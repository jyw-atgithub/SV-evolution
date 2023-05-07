#!/bin/bash

#SBATCH --job-name=DSPR-flye-assembly    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=1-14      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=40    ## number of cores the job needs
#SBATCH --mem-per-cpu=6G     # requesting 6 GB memory per CPU


##pauvre is useful to check the sequencing stats (pauvre stats)

##pacbio_pipeline for ALL NCBI sequences
### on HPC3 !!!!###
wd="/pub/jenyuw/SV-project-temp"
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

# Porechop-abi is not required for PacBio
# "The input contain reads with duplicated IDs." ==> Error From Flye
# because we merged the reads from different runs
# need to use Seqkit to rename those duplicated IDs

unpigz -p ${nT} -c ${file} | chopper -l 600 --headcrop 46 --tailcrop 46 --threads ${nT} |\
#tee ${wd}/${name}.trimmed.fastq.gz |
seqkit rename -j ${nT} |pigz -p ${nT} > ${wd}/${name}.trimmed.rn.fastq.gz 
wait
flye --threads $nT --genome-size 170m --asm-coverage 90 --pacbio-raw ${wd}/${name}.trimmed.rn.fastq.gz \
--out-dir ${wd}/${name}_Flye

#flye does NOT accept stdin (/dev/stdin)
#unpigz -p 4 -c A1_combine_pacbio.fastq.gz|head -n 10000 |grep "@" | sort | uniq -d
