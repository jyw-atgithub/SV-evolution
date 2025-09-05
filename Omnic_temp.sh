#! /bin/bash
#SBATCH --job-name=fastp
#SBATCH -A jje_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=20
raw="/dfs7/jje/jenyuw/Assembling_ISO1/raw"
trimmed="/dfs7/jje/jenyuw/Assembling_ISO1/results/trimmed"
filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
fastp --verbose --detect_adapter_for_pe --overrepresentation_analysis --correction \
--cut_right --average_qual 20 --length_required 100 --thread 20 \
--html ${trimmed}/S33_qc.html \
-i ${raw}/Iso1_6_S33_R1_001.fastq.gz -I ${raw}/Iso1_6_S33_R2_001.fastq.gz \
-o ${trimmed}/S33_R1.fastq.gz -O ${trimmed}/S33_R2.fastq.gz

fastp --verbose --detect_adapter_for_pe --overrepresentation_analysis --correction \
--cut_right --average_qual 20 --length_required 100 --thread 20 \
--html ${trimmed}/S34_qc.html \
-i ${raw}/Iso1_7_S34_R1_001.fastq.gz -I ${raw}/Iso1_7_S34_R2_001.fastq.gz \
-o ${trimmed}/S34_R1.fastq.gz -O ${trimmed}/S34_R2.fastq.gz
cat ${trimmed}/S{33,34}_R1.fastq.gz > ${filtered}/omnic_R1.fastq.gz
cat ${trimmed}/S{33,34}_R2.fastq.gz > ${filtered}/omnic_R2.fastq.gz

#! /bin/bash
#SBATCH --job-name=omnic
#SBATCH -A jje_lab
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=6G
#SBATCH --tmp=250G                ## requesting local scratch
#SBATCH --constraint=fastscratch  ## requesting nodes with fast scratch in /tmp
#SBATCH --output="omnic_scfd.out"
nT=$SLURM_CPUS_PER_TASK

filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
assembly="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly"
omnic="/dfs7/jje/jenyuw/Assembling_ISO1/results/omnic"

genome="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly/hifiasm_14/ONT_14.fasta.gz"
cp $genome $omnic
genome=${omnic}/`basename $genome`

samtools faidx $genome
cut -f1,2 $genome.fai > ${genome}.path
bwa index $genome

bwa mem -5SP -T0 -t $nT $genome ${filtered}/omnic_R1.fastq.gz ${filtered}/omnic_R2.fastq.gz | samtools view -bS -@ $nT -o ${omnic}/read_alignment.bam

samtools view ${omnic}/read_alignment.bam |\
pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $nT --nproc-out $nT --chroms-path ${genome}.path | \
pairtools sort --tmpdir=$TMPDIR --nproc $nT |\
pairtools dedup --nproc-in $nT --nproc-out $nT --mark-dups --output-stats ${omnic}/stats.txt |\
pairtools split --nproc-in $nT --nproc-out $nT --output-pairs ${omnic}/mapped.pairs --output-sam - |\
samtools view -bS -@ $nT | \
samtools sort -n -@ $nT -o ${omnic}/mapped.PT.bam

#real scaffolding
yahs $genome ${omnic}/mapped.PT.bam -o ${omnic}/scfd -v 1

##--> it still creat misjoining.