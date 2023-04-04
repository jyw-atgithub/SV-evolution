#!/bin/bash
cd /home/jenyuw/DSPR_snp/raw_fq
conda activate everything

n_thread=60
ref=/home/jenyuw/Reference_genome/Dmel6/dmel-all-chromosome-r6.46.fasta
ref_fai=/home/jenyuw/Reference_genome/Dmel6/dmel-all-chromosome-r6.46.fasta.fai


for j in ./combined_bam/*.bam
do
file=$(basename $j)
strain=${file%.bam}
echo $strain
#call SNPs, single-thread
freebayes -f $ref -g 200 $j >./snp/"$strain"_freebayes.vcf

done
echo "SNPs are called!"


for k in ./snp/*_freebayes.vcf
do
file=$(basename $k)
strain=${file:0:2}
echo $strain
bgzip -k -@ $n_thread $k
# index
tabix -p vcf "$k".gz
rtg vcfstats "$k".gz
rtg vcffilter -q 30 -i "$k".gz -o "-"|vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1"|\
bgzip > ./snp/$strain.2filtered.vcf.gz
#ANNOTATING. If "-t $n_thread" is used, no summary file will be generated.
snpEff -v BDGP6.32.105 ./snp/$strain.2filtered.vcf.gz >./snp/$strain.annotated.vcf &
done

echo "Finally it is the end!"


#######
#parallelization the freebayes
#parallel --verbose --jobs $n_thread ./calling.sh ::: together.bam  ::: together2_freebayes.vcf
#this does NOT work, because the parallel does not break jobs into tiny chunks.
#######


##################
# single file workflow
# "together.bam" is the merged bam file from all strains and reads. 
# Actually, Freebayes suggested to call the SNPs from all individuals in a population all at once instead of calling by each individual.
n_thread=60
ref=/home/jenyuw/Reference_genome/Dmel6/dmel-all-chromosome-r6.46.fasta
ref_fai=/home/jenyuw/Reference_genome/Dmel6/dmel-all-chromosome-r6.46.fasta.fai

#file location: /home/jenyuw/DSPR_snp/raw_fq/together
# index of bam file is required for freebayes-parallel
samtools index -@ $n_thread -b together.bam

#freebayes-parallel <(fasta_generate_regions.py $ref_fai 100000) $n_thread -f $ref together.bam > parr_out.vcf
#freebayes-parallel <(fasta_generate_regions.py $ref 100000) $n_thread -f $ref together.bam > parr_out.vcf
freebayes-parallel <(fasta_generate_regions.py $ref 100000) $n_thread -f $ref -g 200 together.bam > parr_out.vcf

bgzip -k -@ $n_thread parr_out.vcf
tabix -p vcf parr_out.vcf.gz
rtg vcfstats parr_out.vcf.gz
#-d --min-read-depth=INT, -q --min-quality=FLOAT
rtg vcffilter --all-samples -d 10 -q 30 -i parr_out.vcf.gz -o "-"|\
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1"|\
bgzip > parr_out.2filtered.vcf.gz

snpEff -v BDGP6.32.105 parr_out.2filtered.vcf.gz >parr_out.annotated.vcf

##################
# filtering with new criteria
bcftools filter -i ' QUAL >= 30 && INFO/DP > 20 && INFO/DP < 1000 &&  QUAL / INFO/AO > 10 && SAF > 0 && SAR > 0 && RPR > 1 && RPL > 1' \
--SnpGap 10 --threads 8 -r 2L,2R,3L,3R,4,X,Y  \
-o - parr_out.vcf.gz |\
bcftools filter --threads 8 -i ' MQM >= 30 && MQMR >= 30 ' -o - - |\
bcftools filter --threads 8 -S . -e 'FMT/DP<3 | FMT/GQ<20' -o - - |\
bcftools filter --threads 8 -e 'AC==0 || AC==AN' --SnpGap 10 -o - - |\
bcftools view --threads 8 -m2 -M2 -v snps -O z -o all_variants_filtered.vcf.gz

jenyuw@hydra:~/DSPR_snp/raw_fq/together$ bcftools view -H all_variants_filtered.vcf.gz | wc -l
1476844
[jenyuw@login-i15:/pub/jenyuw/EE283/DNAseq/results/SNP] $bcftools view -H all_variants_filtered.vcf.gz |wc -l
50474
jenyuw@hydra:~/DSPR_snp/official$ bcftools view -H --threads 4 DSPR.r6.SNPs.vcf.gz|wc -l
1776091
