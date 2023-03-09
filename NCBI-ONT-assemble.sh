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
    porechop_abi -abi --threads $nT -i ${i} -o ${trimmed}/${name}.abi.fastq
    cat ${trimmed}/${name}.abi.fastq | chopper -l 500 --headcrop 10 --tailcrop 10 --threads $nT  > ${trimmed}/${name}.trimmed.fastq
    done

conda activate assemble

for i in $(ls ${trimmed}/${SRR_num}*_ONT.trimmed.fastq)
    do
    name=$(basename ${i}|sed s/".trimmed.fastq"//g)
    echo $i 
    echo $name
    flye --threads $nT --genome-size 170m --nano-raw ${i} --out-dir ${assemble}/${name}_Flye
    done


##polishing, without illumina reads or FAST5. --> Medaka or Racon.
##Medaka is designed to be used on Flye assembly directly!!
##Medaka was installed with "pip" because Anaconda kept failing.
conda activate post-proc


for i in $(ls ${assemble}/SRR*_ONT_Flye/assembly.fasta)
do
echo $i 
name=$(echo $i | gawk -F "\/" '{print $7}' 2>>/dev/null| sed s/_ONT_Flye//g)
echo $name

medaka_consensus -i ${raw}/PRJNA929424/${name}_ONT.fastq.gz -d ${i} -o ${polishing}/${name}-medaka-1 \
-t ${nT} -m r941_min_hac_g507
done


for i in /home/jenyuw/SV-project/result/assemble/SRR23269563_ONT_Flye/assembly.fasta
do
echo $i 
name=$(echo $i | gawk -F "\/" '{print $7}' 2>>/dev/null| sed s/_ONT_Flye//g)
echo $name

medaka_consensus -i /home/jenyuw/SV-project/raw/PRJNA929424/SRR23269563_ONT.fastq.gz \
-d /home/jenyuw/SV-project/result/assemble/SRR23269563_ONT_Flye/assembly.fasta \
-o ${polishing}/SRR23269563-medaka-1/ \
-t ${nT} -m r941_min_hac_g507
done


minimap2 -t ${nT} -B 5 -a -x map-ont \
$i ${trimmed}/${name}_ONT.trimmed.fastq |\
samtools view -b -h -@ ${nT} -o - |\
samtools sort -@ ${nT} -o ${aligned_bam}/${name}_ONT.trimmed_assembly.sort.bam
samtools index -@ ${nT} ${aligned_bam}/${name}_ONT.trimmed_assembly.sort.bam

medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${nT}\
-m r941_min_hac_g507



racon -t ${nT} \
${trimmed}/${name}_ONT.trimmed.fastq \
${aligned_bam}/${name}_ONT.trimmed_assembly.sort.bam \
$i \ 
> ${polishing}/${name}.polished.racon.1.fasta

rm ${aligned_bam}/${name}_ONT.trimmed_assembly.sort.bam
rm ${aligned_bam}/${name}_ONT.trimmed_assembly.sort.bam.bai