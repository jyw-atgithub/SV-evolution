#!/bin/bash
## path

raw="/home/jenyuw/SV-project/raw"
qc_report="/home/jenyuw/SV-project/result/qc_report"
trimmed="/home/jenyuw/SV-project/result/trimmed"
#PRJNA574592
## prep
source ~/.bashrc
nT=4

## The initial QC report LONG reads
conda activate longqc #this contains longQc
# Yes, it is really tricky. Install the full version of LongQC on Anaconda (for all the implicit depency) 
# and also install the longQc locally from github. 
# Then cite the path to locally installed LongQC
for i in $(ls ${raw}/PRJNA574592/*_CLR.fastq.gz)
do
name=$(basename ${i}|sed s/".fastq.gz"//g)
python /home/jenyuw/Software/LongQC/longQC.py sampleqc -p ${nT} -x pb-rs2 -n 6000 \
-o ${qc_report}/${name}_longQC -s ${name} ${i}
done
conda deactivate

conda activate qc 
for i in $(ls ${raw}/PRJNA574592/*_CLR.fastq.gz)
do
name=$(basename ${i}|sed s/".fastq.gz"//g)
unpigz -p ${nT} -c ${i} | chopper -l 600 --headcrop 46 --tailcrop 46 --threads ${nT} |\
pigz -p ${nT} > ${trimmed}/${name}.trimmed.fastq.gz
done
conda deactivate