#!/bin/bash

#SBATCH --job-name=canu    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=60   ## number of cores the job needs
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU
#SBATCH --time=7-00:00:00   # 7 days

###
# It looks wierd on hpc3 right now. canu detected slurm automatically but only utilized 8 cpus
###

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
echo file is $file
echo name is $name
echo read type is $read_type

##For NANOPORE flip-flop R9.4 or R10.3
#If you have over 30x coverage add the options: 'corMhapOptions=--threshold 0.8 --ordered-sketch-size 1000 --ordered-kmer-size 14 correctedErrorRate=0.105.
## useGrid=false, so canu won't detect other machine and charge my personal account
if [[ ${read_type} == "ONT" ]]
then 
canu -p ${name} -d ${assemble}/${name}_canu \
genomeSize=135m \
maxInputCoverage=90 \
minReadLength=500 \
maxThreads=${nT} \
correctedErrorRate=0.105 \
'corMhapOptions=--threshold 0.8 --ordered-sketch-size 1000 --ordered-kmer-size 14' \
useGrid=false \
-nanopore-raw ${file}

elif [[ ${read_type} == "CLR" ]]
then
##For PACBIO Sequel II
canu -p ${name} -d ${assemble}/${name}_canu \
genomeSize=135m \
maxInputCoverage=90 \
minReadLength=500 \
minOverlapLength=500 \
maxThreads=${nT} \
correctedErrorRate=0.035 utgOvlErrorRate=0.065 trimReadsCoverage=2 trimReadsOverlap=500 \
stopOnLowCoverage=2 minInputCoverage=2.5 \
useGrid=false \
-raw -pacbio ${file}

#correctedErrorRate=0.085 corMhapSensitivity=normal
fi