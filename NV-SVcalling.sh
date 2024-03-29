#!/bin/bash
## path
ref_genome="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"

raw="/home/jenyuw/SV-project/raw"
raw_prj="/home/jenyuw/SV-project/raw/PRJNA929424"

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
nT=20

# on THOTH
conda activate sv-calling
: <<'SKIP'
SKIP

## Assembly based
#### We need to use the scaffold!!####
## suggestion from syri: minimap2 -ax asm5 --eqx, --eqx is required
## install syri locally. Conda failed. "python3 setup.py install --user"

for i in $(ls ${scaffold}/nv*/ragtag.scaffold.fasta)
do
echo $i
name=$(echo $i|awk -F "/" '{print $7}')
echo ${name}

minimap2 -t ${nT} -a -x asm5 --cs --eqx \
${ref_genome} ${i} \
|samtools view -b -h -@ ${nT} -o -|samtools sort -@ ${nT} -o ${aligned_bam}/${name}.scfd-ref.sort.bam
samtools index ${aligned_bam}/${name}.scfd-ref.sort.bam

#SVIM-asm
echo "calling SVs of ${name} wiyh SVIM-asm"
svim-asm haploid --sample ${name} \
${SVs}/${name}-svim-asm ${aligned_bam}/${name}.scfd-ref.sort.bam ${ref_genome}

done


for i in $(ls ${polishing}/nv*.polished.pilon.1.fasta)
do
cd ${SVs}
#syri
syri -F B -c ${aligned_bam}/${name}.Flye-ref.sort.bam -r ${ref_genome} \
-q ${i} --dir ${SVs}/syri-result --prefix ${name} --samplename ${name}

# smartie-sv
done

bam="/home/jenyuw/SV-project/result/aligned_bam/nv107.scfd-ref.sort.bam"
ref="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"
asm="/home/jenyuw/SV-project/result/scaffold/nv107/ragtag.scaffold.fasta"
syri -F B -c ${bam} -r ${ref} -q ${asm} --dir /home/jenyuw/SV-project/temp2 --prefix nv107 --samplename nv107

fixchr -F B -c ${bam} -r ${ref} -q ${asm}



sed 's/"_RagTag"//g' qry.filtered.fa > qry.filtered.rn.fa
sed 's/>/>chr/g' ref.filtered.fa > ref.filtered.rn.fa

minimap2 -t 10 -a -x asm5 --cs --eqx ref.filtered.rn.fa qry.filtered.rn.fa |\
samtools view -b -h -@ 10 -o -|samtools sort -@ 10 -o nv107.filtered.scfd-ref.sort.bam
samtools index nv107.filtered.scfd-ref.sort.bam

syri -F B -c nv107.filtered.scfd-ref.sort.bam -r ref.filtered.rn.fa -q qry.filtered.rn.fa \
--dir /home/jenyuw/SV-project/temp2 --prefix nv107 --samplename nv107

grep -v ^#  nv107syri.vcf | gawk '{print $3}' | cut -c 1-3 |sort |uniq -c

## PAV in singularity
## requires a clean folder and a config.json
singularity run --bind "$(pwd):$(pwd)" --writable-tmpfs library://becklab/pav/pav:latest -c 20

## non-default settings of of the config.json: "map_threads"[12], "inv_sig_cluster_svlen_min"[4]



## Mapping-based, so we have to map the "trimmed&filtered" ONT reads to the reference genome.
## 9mapping-based.sh

for i in $(ls ${trimmed}/nv*.trimmed.fastq)
do
name=$(basename ${i}|sed s/".trimmed.fastq"//g)
echo "mapping ${name} trimmed reads to reference genome"

minimap2 -t ${nT} -a -x map-ont \
${ref_genome} $i |\
samtools view -b -h -@ ${nT} -o - |\
samtools sort -@ ${nT} -o ${aligned_bam}/${name}.trimmed-ref.sort.bam
samtools index -@ ${nT} ${aligned_bam}/${name}.trimmed-ref.sort.bam
done


for j in $(ls ${aligned_bam}/*.trimmed-ref.sort.bam)
do
name=$(basename ${j}|sed s/".trimmed-ref.sort.bam"//g)
echo "calling SVs of ${name} wiht mapping based methods"
cd ${SVs}
#sniffles
sniffles --threads ${nT} --allow-overwrite --sample-id ${name}-snif \
--minsupport 10 \
--minsvlen 50 --mapq 20 --min-alignment-length 500 \
--cluster-merge-pos 270 \
--max-del-seq-len 100000 \
--reference ${ref_genome} \
--input $j --vcf "${name}-sniffles.vcf"

#cuteSV #--max_size was not recognized. only -L worked
cuteSV --threads ${nT} --genotype --sample ${name}-cute \
--min_support 10 \
--min_size 50 --min_mapq 20 --min_read_len 500 \
--merge_del_threshold 270 --merge_ins_threshold 270 \ 
-L 100000 \
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
$j ${ref_genome} "${name}-cutesv.vcf" .

#SVIM
svim alignment --sample ${name}-svim \
--min_mapq 20 --min_sv_size 50 \
--max_sv_size 100000 \
--distance_normalizer 900 --cluster_max_distance 0.3 \
${SVs}/${name}-SVIM $j ${ref_genome}

cp ${SVs}/${name}-SVIM/variants.vcf ${SVs}/${name}-SVIM.vcf
done

## the --sample-id option (of sniffles) behaves weirdly, which let vcftoolz generate errors
## vcftoolz is not suitable to compare VCF containing SVs
#conda activate vcf-kit
## change the sample name
#vk rename --subst=SAMPLE:${name} ${name}-sniffles.vcf  > ${name}-sniffles.2.vcf

: <<'SKIP'

SKIP


## filter the SVs!  only major chromosomes
## the vcftools command is not used. 
#vcftools --vcf nv107-cutesv.vcf --chr 2L --chr 2R --chr 3L --chr 3R --chr 4 --chr X --chr Y --out nv107-cutesv --recode --recode-INFO-all

## do all fitering with bcftools
## Remember that the length of DELETION is NEGATIVE!!!
## For BND, no SVLEN is given!!!
## SVIM report DUP:INT and DUP:TANDEM, while cuteSV & sniffles only report DUP
## SVIM does NOT report SVLEN for INV, but cuteSV and sniffles do.

#filter_merge.sh

printf "" > ${SVs}/sample.namelist.txt
for i in $(ls ${SVs}/nv*.vcf)
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
echo ${i}
#name=$(echo ${i})
ls ${SVs}/${i}-*.filtered.vcf >sample_files
SURVIVOR merge sample_files 0.05 3 1 0 1 50 ${merged_SVs}/${i}.consensus.vcf
done <${SVs}/sample.namelist.u.txt

for i in $(ls ${merged_SVs}/*.consensus.vcf)
do
echo $i 
name=$(basename ${i} | gawk -F "." '{print $1}')
bcftools view --threads ${nT} --samples ${name}-cute ${i} -O v -o ${name}.consensus.cut.vcf
done 

ls ${merged_SVs}/*.consensus.cut.vcf >sample.con.list.txt
SURVIVOR merge sample.con.list.txt 0.05 1 1 0 1 50 all.consensus.vcf



bgzip -@ 8 -k
tabix -p vcf -0
bcftools sort -O z -o nv107-SVIM.sort.vcf.gz nv107-SVIM.vcf.gz
#bcftools view --threads 8 -r 2L,2R,3L,3R,4,X,Y \
#-i 'QUAL >= 10  && ( SVLEN >= 50 || SVLEN <= -50 || SVTYPE = "BND")' -O z -o nv107-cutesv.filtered.vcf.gz nv107-cutesv.vcf.gz

.   
bcftools view --threads 8 -r 2L,2R,3L,3R,4,X,Y \
-i 'QUAL >= 10 && FILTER = "PASS"' -O z -o nv107-cutesv.filtered.vcf.gz nv107-cutesv.vcf.gz


bcftools view --threads 8 -r 2L,2R,3L,3R,4,X,Y \
-i 'QUAL >= 10 && FILTER = "PASS"' -O z -o nv107-sniffles.filtered.vcf.gz nv107-sniffles.vcf.gz


bcftools view --threads 8 -r 2L,2R,3L,3R,4,X,Y \
-i 'QUAL >= 10 && FILTER = "PASS"'  -O v -o - nv107-SVIM.sort.vcf.gz |\
 sed 's/SVTYPE=DUP:INT/SVTYPE=DUP/g ; s/SVTYPE=DUP:TANDEM/SVTYPE=DUP/g ' |\
bcftools view --thread 8 -O z -o nv107-SVIM.filtered.vcf.gz




bcftools query -f 'cuteSV %CHROM  %POS   %INFO/SVTYPE %INFO/SVLEN \n' nv107-cutesv.filtered.vcf.gz >stats-nv107.tsv
bcftools query -f 'sniffles %CHROM  %POS   %INFO/SVTYPE %INFO/SVLEN \n' nv107-sniffles.filtered.vcf.gz >>stats-nv107.tsv
bcftools query -f 'SVIM %CHROM  %POS   %INFO/SVTYPE %INFO/SVLEN \n' nv107-SVIM.filtered.vcf.gz >>stats-nv107.tsv

bgzip -dk -f nv107-cutesv.filtered.vcf.gz
bgzip -dk -f nv107-sniffles.filtered.vcf.gz
bgzip -dk -f nv107-SVIM.filtered.vcf.gz

pd=`pwd`
## SURVIVOR does NOT accept vcf.gz !!! and it does not report any warning
ls $pd/nv107-*.filtered.vcf >sample_files
SURVIVOR merge sample_files 0.05 1 1 1 1 50 sample_merged.005.vcf
SURVIVOR merge sample_files 0.1 2 1 1 0 30 sample_merged.01.vcf
SURVIVOR merge sample_files 1000 2 1 1 0 30 sample_merged.1000.vcf
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/'  sample_merged.005.vcf | sed -e 's/\(.\)/\1 /g' > sample_merged_overlapp.005.txt
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/'  sample_merged.01.vcf | sed -e 's/\(.\)/\1 /g' > sample_merged_overlapp.01.txt
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/'  sample_merged.1000.vcf | sed -e 's/\(.\)/\1 /g' > sample_merged_overlapp.1000.txt

SURVIVOR merge sample_files 0.05 3 1 0 1 50 sample_merged.consensus.vcf
#bedtools intersect -a sample_merged.consensus.vcf  -b nv107-sniffles.vcf.gz