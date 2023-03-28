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
## prep
source ~/.bashrc
nT=10

# on THOTH
conda activate sv-calling
: <<'SKIP'
SKIP

## Assembly based
for i in $(ls ${polishing}/nv*.polished.pilon.1.fasta)
do
name=$(basename $i | sed s/.polished.pilon.1.fasta//g)
echo "mapping ${name} polished assembly to reference genome"
minimap2 -t ${nT} -a -x asm5 --cs \
${ref_genome} ${i} \
|samtools view -b -h -@ ${nT} -o -|samtools sort -@ ${nT} -o ${aligned_bam}/${name}.Flye-ref.sort.bam
samtools index ${aligned_bam}/${name}.Flye-ref.sort.bam
#SVIM-asm
echo "calling SVs of ${name} wiht assembly based methods"
svim-asm haploid --sample ${name} \
${SVs}/${name}-svim-asm ${aligned_bam}/${name}.Flye-ref.sort.bam ${ref_genome}
done


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

conda activate sv-calling
for j in $(ls ${aligned_bam}/*.trimmed-ref.sort.bam)
do
name=$(basename ${j}|sed s/".trimmed-ref.sort.bam"//g)
echo "calling SVs of ${name} wiht mapping based methods"
cd ${SVs}

#sniffles
sniffles --threads 16 --allow-overwrite --sample-id ${name}\
--minsvlen 50 --mapq 20  \
--reference ${ref_genome} \
--input $j --vcf "${name}-sniffles.vcf"

#cuteSV
cuteSV --threads 20 --genotype \
-l 50 -L 500000 \
-r 500 -q 20 -s 3 \
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
--sample ${name} \
$j ${ref_genome} "${name}-cutesv.vcf" .

#SVIM
svim alignment --sample ${name} \
--min_sv_size 50 \
${SVs}/${name}-SVIM $j ${ref_genome}
done

## the --sample-id option behaves weirdly, which let vcftoolz generate errors
conda activate vcf-kit
## change the sample name
vk rename --subst=SAMPLE:${name} ${name}-sniffles.vcf  > ${name}-sniffles.2.vcf

: <<'SKIP'

SKIP


## filter the SVs!  only major chromosomes
## not used. 
#vcftools --vcf nv107-cutesv.vcf --chr 2L --chr 2R --chr 3L --chr 3R --chr 4 --chr X --chr Y --out nv107-cutesv --recode --recode-INFO-all

## do all fitering with bcftools
## Remember that the length of DELETION is NEGATIVE!!!
## For BND, no SVLEN is given!!!
bgzip -@ -k
tabix -p vcf -0 
bcftools view --threads 8 -r 2L,2R,3L,3R,4,X,Y \
-i 'QUAL >= 10  && ( SVLEN >= 50 || SVLEN <= -50 || SVTYPE = "BND")' -O z -o nv107-cutesv.filtered.vcf.gz nv107-cutesv.vcf.gz
bcftools view --threads 8 -r 2L,2R,3L,3R,4,X,Y \
-i 'QUAL >= 10  && ( SVLEN >= 50 || SVLEN <= -50 || SVTYPE = "BND")' -O z -o nv107-sniffles.filtered.vcf.gz nv107-sniffles.2.vcf.gz
bcftools view --threads 8 -r 2L,2R,3L,3R,4,X,Y \
-i 'QUAL >= 10  && ( SVLEN >= 50 || SVLEN <= -50 || SVTYPE = "BND")' -O z -o nv107-SVIM.filtered.vcf.gz nv107-SVIM.2.vcf.gz



bcftools query -f 'cuteSV %CHROM  %POS   %INFO/SVTYPE %INFO/SVLEN \n' nv107-cutesv.filtered.vcf.gz >stats-nv107.tsv
bcftools query -f 'sniffles %CHROM  %POS   %INFO/SVTYPE %INFO/SVLEN \n' nv107-sniffles.filtered.vcf.gz >>stats-nv107.tsv
bcftools query -f 'SVIM %CHROM  %POS   %INFO/SVTYPE %INFO/SVLEN \n' nv107-SVIM.filtered.vcf.gz >>stats-nv107.tsv

bgzip -dk -f nv107-cutesv.filtered.vcf.gz
bgzip -dk -f nv107-sniffles.filtered.vcf.gz
bgzip -dk -f nv107-SVIM.filtered.vcf.gz

pd=`pwd`
## SURVIVOR does NOT accept vcf.gz !!! and it does not report any warning
ls $pd/nv107-*.filtered.vcf >sample_files
SURVIVOR merge sample_files 0.05 1 1 1 1 50 sample_merged.vcf
SURVIVOR merge sample_files 0.1 2 1 1 0 30 sample_merged.vcf
SURVIVOR merge sample_files 1000 2 1 1 0 30 sample_merged.vcf
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/'  sample_merged.vcf | sed -e 's/\(.\)/\1 /g' > sample_merged_overlapp.txt

