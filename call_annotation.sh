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
#call SNPs, PARALLEL. NOT working now
#freebayes-parallel <(fasta_generate_regions.py dmel-all-chromosome-r6.46.fasta.fai 100000) 36 \
#-f $ref -g 200 $j>"$strain"_freebayes.vcf
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

# index of bam file is required for freebayes-parallel
samtools index -@ $n_thread -b together.bam

#freebayes-parallel <(fasta_generate_regions.py $ref_fai 100000) $n_thread -f $ref together.bam > parr_out.vcf
#freebayes-parallel <(fasta_generate_regions.py $ref 100000) $n_thread -f $ref together.bam > parr_out.vcf
freebayes-parallel <(fasta_generate_regions.py $ref 100000) $n_thread -f $ref -g 200 together.bam > parr_out.vcf

bgzip -k -@ $n_thread parr_out.vcf
tabix -p vcf parr_out.vcf.gz
rtg vcfstats parr_out.vcf.gz
rtg vcffilter --all-samples -d 10 -q 30 -i parr_out.vcf.gz -o "-"|\
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1"|\
bgzip > parr_out.2filtered.vcf.gz

snpEff -v BDGP6.32.105 parr_out.2filtered.vcf.gz >parr_out.annotated.vcf

##################