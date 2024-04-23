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


## On THOTH

#!/bin/bash
## path
trimmed="/home/jenyuw/SV-project/result/trimmed"
assemble="/home/jenyuw/SV-project/result/assemble"

for i in $(ls ${trimmed}/*.trimmed.rn.fastq.gz)
do
name=$(basename ${i}|sed s/".trimmed.rn.fastq.gz"//g)

echo -e "
job_type = local
job_prefix = nextDenovo
task = all
rewrite = yes
deltmp = yes

parallel_jobs =8 #M gb memory, between M/64~M/32
input_type = raw
read_type = clr # clr, ont, hifi
input_fofn = /home/jenyuw/SV-project/result/assemble/input.fofn
workdir = /home/jenyuw/SV-project/result/assemble/${name}_nextdenovo-30

[correct_option]
read_cutoff = 1k
genome_size = 135m
seed_depth = 30 #you can try to set it 30-45 to get a better assembly result
seed_cutoff = 0
sort_options = -m 100g -t 4 #m=M/(TOTAL_INPUT_BASES * 1.2/4)
minimap2_options_raw = -t 4
pa_correction = 5 #M/(TOTAL_INPUT_BASES * 1.2/4)
correction_options = -p 4 #P cores, P/parallel_jobs

[assemble_option]
minimap2_options_cns = -t 4 -k17 -w17
minimap2_options_map = -t 4 #P cores, P/parallel_jobs
nextgraph_options = -a 1
" >${assemble}/run.cfg

ls $i > ${assemble}/input.fofn
nextDenovo ${assemble}/run.cfg
done




##For PACBIO
for i in $(ls ${trimmed}/ORE.rn.trimmed.fastq.gz)
do
name=$(basename ${i}|sed s/".rn.trimmed.fastq.gz"//g)
canu -p ${name} -d ${assemble}/${name}_canu \
genomeSize=135m \
maxInputCoverage=90 \
minReadLength=300 \
minOverlapLength=300 \
maxThreads=60 \
correctedErrorRate=0.085 corMhapSensitivity=normal \
stopOnLowCoverage=2 minInputCoverage=2.5 \
-raw -pacbio ${i}
done
