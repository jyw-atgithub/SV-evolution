## Fragments


srun -A jje_lab -c 4 --pty --x11 bash -i
cd /pub/jenyuw/

sbank balance statement -a jje_lab
sbank balance statement -u jenyuw


: <<'SKIP'
SKIP

## nameed pipe test
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
## Failed, because NOT all programs can use named pipe


#!/bin/bash
guppy6="/home/jenyuw/SV-project/raw/guppy6"

for i in $(echo /home/jenyuw/NV_reads/nv107/rapid_20180430/fast5)
do
name=$(echo $i | gawk -F "/" ' {print $6}' )
guppy_basecaller -i $i -s ${guppy6}/nv107-${name}/ -c dna_r9.4.1_450bps_hac.cfg \
--recursive --min_qscore 7 --device "cuda:0" \
--detect_adapter --detect_mid_strand_adapter \
--detect_primer --detect_barcodes --enable_trim_barcodes
done

#virtualenv medaka --python=python3 --prompt "(medaka)"
#after installing Medaka with python virtualenv
cd /home/jenyuw/SV-project
. medaka/bin/activate
deactivate

## Parsing the sam file generated by Murmmer4
## origin: https://github.com/mummer4/mummer/issues/24
cat dmel-dsim.sam | sed 's/HD\ /HD/' | sed 's/1.0\ /1.0/' | sed 's/\tSO:coordinate/SO:coordinate/' |\
sed s'/VN1/VN:1/' | sed 's/HD/HD\t/' | sed 's/SO:unsorted/\tSO:unsorted/' | sed 's/@PG /@PG\t/' |\
sed 's/ PN/\tPN/' | sed 's/ VN/\tVN/' | sed 's/ CL/\tCL/' > /dev/null

cat dmel-dsim.sam | sed 's/HD\ /HD/; s/1.0\ /1.0/; s/\tSO:coordinate/SO:coordinate/; s/VN1/VN:1/; s/HD/HD\t/; s/SO:unsorted/\tSO:unsorted/; s/@PG /@PG\t/; s/ PN/\tPN/; s/ VN/\tVN/; s/ CL/\tCL/' |\
samtools view -@ 10 -h --reference ${dmel_ref} - |samtools sort -@ 10 -O bam -o corrected.bam


#!/bin/bash
#on HYDRA
trimmed="/home/jenyuw/SV-project-backup/result/trimmed"
assemble="/home/jenyuw/SV-project-backup/result/assemble"

conda activate assemble

##For NANOPORE
for i in $(ls ${trimmed}/SRR*.trimmed.fastq)
do
name=$(basename ${i}|sed s/".trimmed.fastq"//g)

canu -p ${name} -d ${assemble}/${name}_canu \
genomeSize=135m \
maxInputCoverage=90 \
minReadLength=500 \
maxThreads=60 \
-raw -nanopore ${i}
done

##For PACBIO
for i in $(ls ${trimmed}/*.trimmed.rn.fastq.gz)
do
name=$(basename ${i}|sed s/".trimmed.rn.fastq.gz"//g)

canu -p ${name} -d ${assemble}/${name}_canu \
genomeSize=135m \
maxInputCoverage=90 \
minReadLength=500 \
maxThreads=60 \
-raw -pacbio ${i}
done


grep -v ^#  nv107syri.vcf | gawk '{print $3}' | cut -c 1-3 |sort |uniq -c
     14 CPG
     23 CPL
  76584 DEL
    151 DUP
   4509 HDR
  71926 INS
    352 INV
    231 NOT
 708407 SNP
   4749 SYN
      5 TDM
     93 TRA

bcftools query -f '%INFO/SVTYPE \n'  pav_nv107.vcf.gz|sort -h | uniq -c
 107641 DEL 
 103812 INS 
     15 INV 
 916851 SNV 


prefetch -pcv ${strain}

fastq-dump --split-spot --stdout  /dfs7/jje/jenyuw/SV-project-temp/raw/${strain}/*.sra |\
pigz -p 15 -v  >/dfs7/jje/jenyuw/SV-project-temp/raw/${strain}_${tech}.fastq.gz


bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\t[ %GT] \n' truvari_merge.vcf>extraction.tsv

bedtools intersect -a truvari_merge.sort.vcf  -b dsim-dmel.vcf
bcftools isec --collapse none -n =2 -O v -o ancestral.vcf truvari_merge.sort.vcf.gz dsim-dmel.vcf.gz

truvari consistency truvari_merge.sort.vcf.gz dsim-dmel.vcf.gz



## virtualenv of python3 is installed with pip.
## virtualenv is different from venv
#python -m pip install --user virtualenv
virtualenv medaka --python=python3 --prompt "(medaka)"
##activate the virtial env
. medaka/bin/activate
pip install --upgrade pip
pip install medaka

## deactivate the virtialenv
deactivate