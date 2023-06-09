#!/bin/bash

##ONT_pipeline for ALL NCBI sequences
#!/bin/bash
## path
ref_genome="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"

raw="/home/jenyuw/SV-project/raw"
raw_prj="/home/jenyuw/SV-project/raw/PRJNA929424"
SRR_num="SRR232695"

trimmed="/home/jenyuw/SV-project/result/trimmed"
assemble="/home/jenyuw/SV-project/result/assemble"
aligned_bam="/home/jenyuw/SV-project/result/aligned_bam"
polishing="/home/jenyuw/SV-project/result/polishing"
canu_proc="/home/jenyuw/SV-project/result/canu_processing"
scaffold="/home/jenyuw/SV-project/result/scaffold"
## prep
source ~/.bashrc
nT=20

conda activate longqc

for i in $(ls ${raw_prj}/*_ONT.fastq.gz)
    do
    name=$(basename ${i}|sed s/".fastq.gz"//g)
    echo $i 
    echo $name
    python /home/jenyuw/Software/LongQC/longQC.py sampleqc -p 10 -x ont-ligation -n 6000 \
    -o ${qc_report}/${name}_longQC -s ${name} ${i}
    done

conda activate qc 

for i in $(ls ${raw_prj}/*_ONT.fastq.gz)
do
    name=$(basename ${i}|sed s/".fastq.gz"//g)
    echo $i 
    echo $name
    porechop_abi -abi --threads $nT -i ${i} -o ${trimmed}/${name}.abi.fastq
    cat ${trimmed}/${name}.abi.fastq | chopper -l 500 --headcrop 10 --tailcrop 10 --threads $nT  > ${trimmed}/${name}.trimmed.fastq
    done

conda activate assemble

for i in $(ls ${trimmed}/${SRR_num}*_ONT.trimmed.fastq)
    do
    name=$(basename ${i}|sed s/".trimmed.fastq"//g)
    echo $i 
    echo $name
    flye --threads $nT --genome-size 135m --nano-raw ${i} --out-dir ${assemble}/${name}_Flye
    done



## Long-read assembly with nextDenovo

for i in $(ls ${trimmed}/${SRR_num}*_ONT.trimmed.fastq)
do
name=$(basename ${i}|sed s/".trimmed.fastq"//g)
ls $i > ${assemble}/input.fofn

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
workdir = /home/jenyuw/SV-project/result/assemble/${name}_nextdenovo

[correct_option]
read_cutoff = 1k
genome_size = 135m
seed_depth = 45 #you can try to set it 30-45 to get a better assembly result
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

nextDenovo ${assemble}/run.cfg
done