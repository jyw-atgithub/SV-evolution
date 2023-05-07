#!/bin/bash
# working on DSPR
## Remember to check the cuteSV setting, which needs to match the sequencing tech!!!
## path
ref_genome="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"
raw="/home/jenyuw/SV-project/raw"
#raw_prj="/home/jenyuw/SV-project/raw/PRJNA929424"

trimmed="/home/jenyuw/SV-project/result/trimmed"
assemble="/home/jenyuw/SV-project/result/assemble"
aligned_bam="/home/jenyuw/SV-project/result/aligned_bam"
polishing="/home/jenyuw/SV-project/result/polishing"
canu_proc="/home/jenyuw/SV-project/result/canu_processing"
scaffold="/home/jenyuw/SV-project/result/scaffold"
SVs="/home/jenyuw/SV-project/result/SVs"
merged_SVs="/home/jenyuw/SV-project/result/merged_SVs"

## prep
source ~/.bashrc
nT=24
#sample_prefix="{A{1..7},AB8,B{1,2,3,4,6},ORE}"
# on THOTH
conda activate sv-calling

## Mapping-based, so we have to map the "trimmed&filtered" ONT reads to the reference genome.

for i in $(ls ${trimmed}/*.trimmed.rn.fastq.gz)
do
name=$(basename ${i}|sed s/".trimmed.rn.fastq.gz"//g)
echo "mapping ${name} trimmed reads to reference genome"

minimap2 -t ${nT} -a -x map-ont \
${ref_genome} $i |\
samtools view -b -h -@ ${nT} -o - |\
samtools sort -@ ${nT} -o ${aligned_bam}/${name}.trimmed-ref.sort.bam
samtools index -@ ${nT} ${aligned_bam}/${name}.trimmed-ref.sort.bam
done

for j in $(ls ${aligned_bam}/{A{1..7},AB8,B{1,2,3,4,6},ORE}.trimmed-ref.sort.bam)
do
name=$(basename ${j}|sed s/".trimmed-ref.sort.bam"//g)
echo "calling SVs of ${name} wiht mapping based methods"
done
cd ${SVs}

#cuteSV #--max_size was not recognized. only -L worked
## Remember to check the cuteSV setting, which needs to match the sequencing tech!!!
cuteSV --threads ${nT} --genotype --sample ${name}-cute \
--min_support 10 \
--min_size 50 --min_mapq 20 --min_read_len 500 \
--merge_del_threshold 270 --merge_ins_threshold 270 \ 
-L 100000 \
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5 \
${j} ${ref_genome} ${name}-cutesv.vcf ${SVs}

#sniffles
sniffles --threads ${nT} --allow-overwrite --sample-id ${name}-snif \
--minsupport 10 \
--minsvlen 50 --mapq 20 --min-alignment-length 500 \
--cluster-merge-pos 270 \
--max-del-seq-len 100000 \
--reference ${ref_genome} \
--input $j --vcf "${name}-sniffles.vcf"

#SVIM
svim alignment --sample ${name}-svim \
--min_mapq 20 --min_sv_size 50 \
--max_sv_size 100000 \
--distance_normalizer 900 --cluster_max_distance 0.3 \
${SVs}/${name}-SVIM $j ${ref_genome}

cp ${SVs}/${name}-SVIM/variants.vcf ${SVs}/${name}-SVIM.vcf
done

printf "" > ${SVs}/sample.namelist.txt
for i in $(ls ${SVs}/{A{1..7},AB8,B{1,2,3,4,6},ORE}-*.vcf)
do
name=$(basename $i | gawk -F "-" '{print $1}')
prog=$(basename $i | gawk -F "-" '{print $2}'|sed 's/.vcf//g')
echo $name $prog
bgzip -f -@ ${nT} -k ${i}
# Only the vcf from SVIM is not sorted while others are sorted. We sort all because of convenience.
bcftools sort -O z -o ${SVs}/${name}-${prog}.sort.vcf.gz ${i}
tabix -f -p vcf ${SVs}/${name}-${prog}.sort.vcf.gz

bcftools view --threads 8 -r 2L,2R,3L,3R,4,X,Y \
-i 'QUAL >= 10 && FILTER = "PASS"'  -O v -o - ${SVs}/${name}-${prog}.sort.vcf.gz |\
sed 's/SVTYPE=DUP:INT/SVTYPE=DUP/g ; s/SVTYPE=DUP:TANDEM/SVTYPE=DUP/g ' |\
bcftools view --thread 8 -O z -o ${SVs}/${name}-${prog}.filtered.vcf.gz
echo ${name} >> ${SVs}/sample.namelist.txt
bgzip -f -dk ${SVs}/${name}-${prog}.filtered.vcf.gz
done

cat ${SVs}/sample.namelist.txt | sort | uniq >${SVs}/sample.namelist.u.txt

while read i
do
echo "we are merging SVs ${i} called by 3 programs"
#name=$(echo ${i})
ls ${SVs}/${i}-*.filtered.vcf >sample_files
SURVIVOR merge sample_files 0.05 3 1 0 1 50 ${merged_SVs}/${i}.consensus.vcf
done <${SVs}/sample.namelist.u.txt