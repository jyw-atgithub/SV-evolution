#!/bin/bash
## path
raw="/home/jenyuw/SV-project/raw"
qc_report="/home/jenyuw/SV-project/result/qc_report"
trimmed="/home/jenyuw/SV-project/result/trimmed"
assemble="/home/jenyuw/SV-project/result/assemble"
aligned_bam="/home/jenyuw/SV-project/result/aligned_bam"
polishing="/home/jenyuw/SV-project/result/polishing"
ref_genome="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"
## prep
source ~/.bashrc
nT=10


## The initial QC report SHORT reads
conda activate qc #this contain fastqc, porechop_abi, chopper
fastqc -t $nT -f fastq -o ${qc_report} ${raw}/nv1*_illumina_r?.fastq
#parallel "fastqc -f fastq -o ${qc_report} " ::: ${raw}/*.fastq
#fastq has built-in parallel option
wait
## The initial QC report LONG reads
#python longQC.py sampleqc -x ont-ligation -o ${qc_report} nv107_combined.fastq
#not working now


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


## Long-read assembly
conda activate assemble #this include flye, canu, bwa, fastp, trimmomatic, assembly-stats, mummer(4), bwa-mem2
# --nano-raw, The expected error rate is 10-15%
for i in $(ls ${trimmed}/*.trimmed.fastq)
do
name=$(basename ${i}|sed s/".trimmed.fastq"//g)
flye --threads $nT --genome-size 170m --nano-raw ${i} --out-dir ${assemble}/${name}
done

## Short-read mapping with sortingraw="/home/jenyuw/SV-project/raw"
qc_report="/home/jenyuw/SV-project/result/qc_report"
trimmed="/home/jenyuw/SV-project/result/trimmed"
assemble="/home/jenyuw/SV-project/result/assemble"
aligned_bam="/home/jenyuw/SV-project/result/aligned_bam"
polishing="/home/jenyuw/SV-project/result/polishing"
ref_genome="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"
#bwa index ${ref_genome}

: <<'SKIP' # this does NOT work.
SKIP
nptest="/home/jenyuw/SV-project/np-test"

for j in $(ls ${raw}/nv1*_illumina_r1.fastq)
do
name=$(basename ${j} |sed "s/_illumina_r.*.fastq//g")
r2=$(echo $j |sed 's/r1/r2/')
#The rules (of using * and ?) in sed is different.
mkfifo ${name}.read1 ${name}.read2
fastp -i ${j} -I ${r2} -o ${nptest}/${name}.read1 -O ${nptest}/${name}.read2 --thread ${nT} --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_tail --average_qual 3
bwa mem -M -t ${nT} ${ref_genome} ${nptest}/${name}.read1 ${nptest}/${name}.read2 | samtools view -bh -|\
samtools sort -@ ${nT} -o ${nptest}/${name}.sort.bam -
done




for j in $(ls ${raw}/nv1*_illumina_r1.fastq)
do
echo $j
name=$(basename ${j} |sed "s/_illumina_r.*.fastq//g")
#The rules (of using * and ?) in sed is different.
echo $name
r2=$(echo $j |sed 's/r1/r2/')
echo $r2
fastp -i ${j} -I ${r2} -o ${trimmed}/${name}.trimmed.r1.fastq -O ${trimmed}/${name}.trimmed.r2.fastq \
--thread ${nT} --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_tail --average_qual 3

bwa mem -M -t ${nT} ${ref_genome} ${trimmed}/${name}.trimmed.r1.fastq ${trimmed}/${name}.trimmed.r2.fastq |\
samtools sort -@ ${nT} - -o ${aligned_bam}/${name}.sort.bam
done

assembly-stats
#dnadiff is a part of mummer
dnadiff –p out assembly.fasta ${ref_genome}


for k in $(ls ${assemble}/nv???/assembly.fasta)
do
name=$(echo $k | sed "s@${assemble}\/@@g; s@\/assembly.fasta@@g")
r1="${trimmed}/${name}.trimmed.r1.fastq"
r2="${trimmed}/${name}.trimmed.r2.fastq"

#bwa index $k
#bwa mem -t ${nT} $k $r1 $r2 |samtools view -S -b -h |\
#samtools sort -@ ${nT} -o ${polishing}/${name}_ILL2ONT.sort.bam
#samtools index ${polishing}/${name}_ILL2ONT.sort.bam

java -XX:+AggressiveHeap -jar /home/jenyuw/Software/pilon-1.24.jar --diploid \
--genome ${assemble}/${name}/assembly.fasta \
--frags ${polishing}/${name}_ILL2ONT.sort.bam \
--output ${name} --outdir ${polishing} 
done
#--threads is not supported by Pilon anymore
#Do NOT use the pilon installed by Anaconda, it will crash because of memory limit.
#Just download the precompiled jar file from the latest release on Github.

java -Xmx128G -jar /home/jenyuw/Software/pilon-1.24.jar --diploid \
--genome ${assemble}/nv107/assembly.fasta \
--frags ${polishing}/nv107_ILL2ONT.sort.bam \
--output nv107 --outdir ${polishing} 

java -XX:+AggressiveHeap -jar /home/jenyuw/Software/pilon-1.24.jar --diploid \
--genome ${assemble}/nv107/assembly.fasta \
--frags ${polishing}/nv107_ILL2ONT.sort.bam \
--output nv107 --outdir ${polishing} 


trimmed="/home/jenyuw/SV-project/result/trimmed"
ref1="/home/jenyuw/SV-project/bwamem2-test/dmel-all-chromosome-r6.49.fasta"
ref2="/home/jenyuw/SV-project/bwamem2-test/dmel-all-chromosome-r6.49-2.fasta"
source ~/.bashrc
nT=10

time bwa index ${ref1}
echo "bwa mem began at"
date
time bwa mem -t ${nT} ${ref1} ${trimmed}/nv107.trimmed.r1.fastq ${trimmed}/nv107.trimmed.r2.fastq >aln1.sam
date

time bwa-mem2 index ${ref2}
echo "bwa mem 2 began at"
date
time bwa-mem2 mem -t ${nT} ${ref1} ${trimmed}/nv107.trimmed.r1.fastq ${trimmed}/nv107.trimmed.r2.fastq >aln2.sam
date

java -Xmx128G -jar /home/jenyuw/Software/pilon-1.24.jar