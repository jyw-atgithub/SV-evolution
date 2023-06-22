#!/bin/bash
## path
ref_genome="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"
raw="/home/jenyuw/SV-project/raw"
qc_report="/home/jenyuw/SV-project/result/qc_report"
trimmed="/home/jenyuw/SV-project/result/trimmed"
assemble="/home/jenyuw/SV-project/result/assemble"
aligned_bam="/home/jenyuw/SV-project/result/aligned_bam"
polishing="/home/jenyuw/SV-project/result/polishing"
canu_proc="/home/jenyuw/SV-project/result/canu_processing"
scaffold="/home/jenyuw/SV-project/result/scaffold"
busco_out="/home/jenyuw/SV-project/result/busco_out"
## prep
source ~/.bashrc
nT=30


## The initial QC report SHORT reads
conda activate qc #this contain fastqc, porechop_abi, chopper
fastqc -t $nT -f fastq -o ${qc_report} ${raw}/nv1*_illumina_r?.fastq
#parallel "fastqc -f fastq -o ${qc_report} " ::: ${raw}/*.fastq
#fastq has built-in parallel option
wait


## The initial QC report LONG reads
conda activate longqc #this contains longQc
# Yes, it is really tricky. Install the full version of LongQC on Anaconda (for all the implicit depency) 
# and also install the longQc locally from github. 
# Then cite the path to locally installed LongQC
for i in $(ls ${raw}/*_combined.fastq)
do
name=$(basename ${i}|sed s/"_combined.fastq"//g)
python /home/jenyuw/Software/LongQC/longQC.py sampleqc -p 10 -x ont-ligation -n 9999 \
-o ${qc_report}/${name}_longQC -s ${name} ${i}
done
#The command below does not work properly bacause it cannot general reports.
#python /home/jenyuw/anaconda3/envs/longqc/bin/longQC.py sampleqc 




# pycoQC is also installed locally
# pycoQC uses the summary file generated by the basecaller, such as guppy, NOT the fastq file.
# pycoQC -f nv107_illumina_r1.fastq -o ${qc_report}/nv107_illumina_r1.html # this must fail

## Removing the adapters and trimmming. For NANOPORE

for i in $(ls ${raw}/*_combined.fastq)
do
name=$(basename ${i}|sed s/"_combined.fastq"//g)
porechop_abi -abi --threads $nT -i ${i} -o ${trimmed}/${name}.abi.fastq
cat ${trimmed}/${name}.abi.fastq | chopper -l 500 --headcrop 10 --tailcrop 10 --threads $nT  > ${trimmed}/${name}.trimmed.fastq
done
#porechop_abi -abi --threads 10 -i ${raw}/nv107_combined.fastq -o /dev/stdout |\
#chopper -l 500 --headcrop 10 --tailcrop 10 --threads 10 > ${trimmed}/nv107.trimmed.fastq
# Chopper only accept stdin, but this does NOT work.


## Long-read assembly with flye
conda activate assemble #this include flye, canu, bwa, fastp, trimmomatic, assembly-stats, mummer(4), bwa-mem2
# --nano-raw, The expected error rate is 10-15%
for i in $(ls ${trimmed}/nv*.trimmed.fastq)
do
name=$(basename ${i}|sed s/".trimmed.fastq"//g)
flye --threads $nT --genome-size 135m --nano-raw ${i} --out-dir ${assemble}/${name}_Flye
done

## Long-read assembly with Canu
# use Canu to correct, trim and assemble, all at once
for i in $(ls ${trimmed}/nv*.trimmed.fastq)
do
name=$(basename ${i}|sed s/".trimmed.fastq"//g)

canu \
-p ${name}.corrected -d ${assemble}/${name}_canu \
genomeSize=170m \
maxInputCoverage=90 \
minReadLength=500 \
-raw -nanopore ${i}
done

## Long-read assembly with nextDenovo

for i in $(ls ${trimmed}/nv*.trimmed.fastq)
do
name=$(basename ${i}|sed s/".trimmed.fastq"//g)

echo -e "
job_type = local
job_prefix = nextDenovo
task = all
rewrite = yes
deltmp = yes

parallel_jobs =8 #M gb memory, between M/64~M/32
input_type = raw
read_type = ont # clr, ont, hifi
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

: << 'SKIP'
## Long-read assembly with GoldRush
## Not working right now
conda activate goldrush
for i in $(ls ${trimmed}/nv107.trimmed.fastq)
do
name=$(basename ${i}|sed s/".trimmed.fastq"//g)
input=$(echo ${i}|sed s/.fastq//)
mkdir ${assemble}/${name}_goldpath
goldrush run reads=${input} G=1.35e8 t=22 track_time=1 p=$name
done
SKIP

##Quality control
assembly-stats
#dnadiff is a part of mummer
dnadiff –p out assembly.fasta ${ref_genome}