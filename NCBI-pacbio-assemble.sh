#!/bin/bash
##pauvre is useful to check the sequencing stats (pauvre stats)


##pacbio_pipeline for ALL NCBI sequences

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
nT=8

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
# only chopper is enough for PacBio
    cat ${trimmed}/${name}.abi.fastq | chopper -l 500 --headcrop 10 --tailcrop 10 --threads $nT  > ${trimmed}/${name}.trimmed.fastq
    done

#First assembly with Flye
conda activate assemble

for i in $(ls ${trimmed}/${SRR_num}*_ONT.trimmed.fastq)
    do
    name=$(basename ${i}|sed s/".trimmed.fastq"//g)
    echo $i 
    echo $name
    flye --threads $nT --genome-size 170m --nano-raw ${i} --out-dir ${assemble}/${name}_Flye
    done


conda activate post-proc

##polishing






