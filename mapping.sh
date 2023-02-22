#!/bin/bash
cd /home/jenyuw/DSPR_snp/raw_fq
conda activate everything

forward="/home/jenyuw/DSPR_snp/raw_fq/*_1.fq"
#forward="/home/jenyuw/DSPR_snp/raw_fq/test/*_1.fq"
n_thread=60
ref=/home/jenyuw/Reference_genome/Dmel6/dmel-all-chromosome-r6.46.fasta
ref_fai=/home/jenyuw/Reference_genome/Dmel6/dmel-all-chromosome-r6.46.fasta.fai

#INDEXING THE GENOME for mapping. Only needed once.
bwa index $ref
#INDEXING THE GENOME for bam file. A .fai file is created. Only needed once.
samtools faidx $ref

for item in $forward
do
#echo $item
out_name1=$(basename $item)
#echo $out_name1
strain=${out_name1%_*_*.fq}
series=${out_name1%_*fq}
sample=${out_name1:3:1}
#echo "series name is $series"
#echo "strain name is $strain"
#echo "sample number is $sample"
item2=${item//"_1.fq"/"_2.fq"}
#echo $item2
out_name2=$(basename $item2)
#echo $out_name2

#"CLEANING UP"
fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread $n_thread --html ./trimmed/$series.html -i $item -I $item2 -o ./trimmed/trim$out_name1 -O ./trimmed/trim$out_name2

#"MAPPING!"
bwa aln -t $n_thread -f "$out_name1".sai $ref ./trimmed/trim$out_name1
bwa aln -t $n_thread -f "$out_name2".sai $ref ./trimmed/trim$out_name2
bwa sampe -r "@RG\tID:"$series"\tPL:illumina\tLB:"$sample"\tSM:"$strain"" $ref ""$out_name1".sai" ""$out_name2".sai" ./trimmed/trim$out_name1 ./trimmed/trim$out_name2 > ""$series"_pe_align.sam"

#"Sorting THE ALIGNMENT"
#"samtools fixmate expects name-sorted input files, which we can achieve with samtools sort -n."
samtools sort -n -O sam -n -@ $n_thread ""$series"_pe_align.sam"|samtools fixmate -m -O bam -@ $n_thread "-" "-"|\
samtools sort -O bam -n -@ $n_thread -|samtools sort -o "$series"_pe_sortedbyc.bam -O bam -@ $n_thread -
#"Cleaning THE ALIGNMENT"
samtools markdup -s -S -r -m s -@ $n_thread "$series"_pe_sortedbyc.bam -|\
samtools view -h -b -q 20 -@ 59 "-">"./well_mapped/"$series"_pe_rmdup_qc.bam"
bamtools index -in "./well_mapped/"$series"_pe_rmdup_qc.bam"
qualimap bamqc -bam "./well_mapped/"$series"_pe_rmdup_qc.bam"
done

echo "Mapping is done!"

mkdir combined_bam
mkdir snp
var="empty"

for i in ./well_mapped/*.bam
do
file=$(basename $i)
strain=${file:0:2}
echo "$strain is strain"
if [ "$strain" = "$var" ]
then
echo "$strain has been merged"
else
var=$strain
echo "$var is var"
#merge SAM files
samtools merge -@ $n_thread -o ./combined_bam/"$strain".bam ./well_mapped/"$strain"_1_pe_rmdup_qc.bam ./well_mapped/"$strain"_2_pe_rmdup_qc.bam ./well_mapped/"$strain"_3_pe_rmdup_qc.bam
fi
done