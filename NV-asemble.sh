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
nT=10


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
for i in $(ls ${trimmed}/*.trimmed.fastq)
do
name=$(basename ${i}|sed s/".trimmed.fastq"//g)
flye --threads $nT --genome-size 170m --nano-raw ${i} --out-dir ${assemble}/${name}_Flye
done

## Long-read assembly with Canu
# use Canu to correct, trim and assemble, all at once


for i in $(ls ${trimmed}/*.trimmed.fastq)
do
name=$(basename ${i}|sed s/".trimmed.fastq"//g)

canu \
-p ${name}.corrected -d ${assemble}/${name}_canu \
genomeSize=170m \
maxInputCoverage=90 \
minReadLength=500 \
-raw -nanopore ${i}
done



/data/home/jenyuw/anaconda3/envs/assemble/bin/sqStoreCreate \
  -o ./nv107.corrected.seqStore.BUILDING \
  -minlength 1000 \
  -genomesize 170000000 \
  -coverage   200 \
  -bias       0 \
  -raw -nanopore nv107_combined /data/home/jenyuw/NV_reads/nv107/combined_fastq/nv107_combined.fastq \


canu [-haplotype|-correct|-trim] \
   [-s <assembly-specifications-file>] \
   -p <assembly-prefix> \
   -d <assembly-directory> \
   genomeSize=<number>[g|m|k] \
   [other-options] \
   [-trimmed|-untrimmed|-raw|-corrected] \
   [-pacbio|-nanopore|-pacbio-hifi] *fastq



## Short-read mapping with sorting

#bwa index ${ref_genome}

: <<'SKIP'
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


for k in $(ls ${assemble}/nv???_Flye/assembly.fasta)
do
name=$(echo $k | sed "s@${assemble}\/@@g; s@_Flye@@g;; s@\/assembly.fasta@@g")
r1="${trimmed}/${name}.trimmed.r1.fastq"
r2="${trimmed}/${name}.trimmed.r2.fastq"
#bwa index $k
#bwa mem -t ${nT} $k $r1 $r2 |samtools view -S -b -h |\
#samtools sort -@ ${nT} -o ${polishing}/${name}_ILL2ONT.sort.bam
#samtools index ${polishing}/${name}_ILL2ONT.sort.bam

java -XX:+AggressiveHeap -jar /home/jenyuw/Software/pilon-1.24.jar --diploid \
--genome ${assemble}/${name}/assembly.fasta \
--frags ${polishing}/${name}_ILL2ONT.sort.bam \
--output ${name}.polished --outdir ${polishing} 
done
#--threads is not supported by Pilon anymore
#Do NOT use the pilon installed by Anaconda, it will crash because of memory limit.
#Just download the precompiled jar file from the latest release on Github.
#"java -Xmx128G -jar" means giving the program an allowance of 128Gb memory.
#java -XX:+AggressiveHeap -jar means letting the program use as much memory as needed.



## scaffolding.sh
conda activate post-proc
for i in $(ls ${polishing}/*.polished.fasta)
do
name=$(basename $i |sed s/.polished.fasta//g)
echo $name
ragtag.py scaffold -t ${nT} -o ${scaffold}/${name} $ref_genome ${i}
done

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

conda activate busco
for i in $(ls ${assemble}/nv*_Flye/assembly.fasta)
do
strain=$(echo $i | gawk -F "\/" '{print $7}' 2>>/dev/null| sed s/_Flye//g)
echo $strain
busco -i ${i} --out_path ${busco_out} -o ${strain} -m genome --cpu 20 -l diptera_odb10
done

busco -i nv107.polished.fasta -o nv107.p -m genome --cpu 20 -l diptera_odb10


minimap2 -t 20 -a -B 5 -x map-ont \
/home/jenyuw/SV-project-backup/reference_genome/dmel-all-chromosome-r6.49.fasta \
/home/jenyuw/SV-project-backup/result/assembly/nv107_Flye_assembly.fasta |samtools view -b -h -@ 20 -o nv107_mapped.bam
samtools sort -@ 20 -o nv107_mapped.sort.bam nv107_mapped.bam


cuteSV --threads 20 \
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
nv107_mapped.sort.bam /home/jenyuw/SV-project-backup/reference_genome/dmel-all-chromosome-r6.49.fasta \
nv107.vcf .