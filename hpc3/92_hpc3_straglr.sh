#!/bin/bash

#SBATCH --job-name=straglr    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p highmem       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=24   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU
###SBATCH --tmp=100G                ## requesting 100 GB local scratch
###SBATCH --constraint=fastscratch  ## requesting nodes with fast scratch in /tmp

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
trimmed="/dfs7/jje/jenyuw/SV-project-temp/result/trimmed"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
tandem_repeat="/dfs7/jje/jenyuw/SV-project-temp/result/tandem_repeat"
source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK

file=`head -n $SLURM_ARRAY_TASK_ID ${trimmed}/namelist.txt |tail -n 1`
name=$(basename ${file}|sed s/".trimmed.fastq.gz"//g)
read_type=`echo ${name} | gawk -F "_" '{print $2}'`
echo "the file is ${file}"
echo "the name is ${name}"

##!!!Importent!! remeber to declare the array##
declare -A mapping_option=(["CLR"]='map-pb' ["hifi"]='map-hifi' ["ONT"]='map-ont')
echo "The mapping option is ${mapping_option[$read_type]}"

##Please use the option -Y to enable soft-clipping so that read sequences can be assessed directly from the BAM file.
minimap2 -t ${nT} -a -x ${mapping_option[$read_type]} -Y \
${ref_genome} ${file} |\
samtools view -b -h -@ ${nT} -o - |\
samtools sort -@ ${nT} -o ${aligned_bam}/${name}.trimmed-ref.SOFT.bam
samtools index -@ ${nT} ${aligned_bam}/${name}.trimmed-ref.SOFT.bam

#module load python/3.10.2 ## It produced some exceptions
module load python/3.8.0
cd ${tandem_repeat}
python3 /pub/jenyuw/Software/straglr/straglr.py ${aligned_bam}/${name}.trimmed-ref.SOFT.bam ${ref_genome} ${name} \
--nprocs ${nT} --min_ins_size 50 --max_str_len 100
#module unload python/3.10.2
module unload python/3.8.0

#python /pub/jenyuw/Software/straglr/straglr.py ${aligned_bam}/SRR9951099_ONT.trimmed-ref.SOFT.bam ${ref_genome} SRR9951099_ONT \
#--nprocs 16 --min_ins_size 50 --max_str_len 100


#Remove SVs overlapping with  STRs
dmel_ref="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
dsim_ref="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing/GCF_016746395.2_Dsim_3.1.fasta"
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
merged_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/merged_SVs"
con_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/consensus_SVs"

bedtools intersect -header -v -a ${con_SVs}/${name}.asm-2.tru_con.sort.vcf.gz -b ${tandem_repeat}/${name}.bed |\
bgzip -@ ${nT} > ${con_SVs}/${name}.asm-2.noSTR.tru_con.vcf.gz
bcftools index -f -t ${con_SVs}/${name}.asm-2.noSTR.tru_con.vcf.gz

bcftools merge -m none --threads ${nT} ${con_SVs}/*.asm-2.noSTR.tru_con.vcf.gz -O z -o - |\
bcftools sort --max-mem 4G -O z -o - > ${merged_SVs}/merge.asm-2.noSTR.vcf.gz
bcftools index -f -t ${merged_SVs}/merge.asm-2.noSTR.vcf.gz

module load python/3.10.2
##Do not use "-f" in truvari collapse, because this is an old function and less accurate
truvari collapse --sizemax 5000000 -k common \
-i ${merged_SVs}/merge.asm-2.noSTR.vcf.gz \
-c ${merged_SVs}/truvari.asm-2.noSTR.collapse.vcf.gz |\
bcftools sort --max-mem 4G |\
bgzip -@ ${nT} > ${merged_SVs}/truvari.asm-2.noSTR.vcf.gz

bcftools index -f -t ${merged_SVs}/truvari.asm-2.noSTR.vcf.gz

##Polarizing Again
ls ${con_SVs}/*.asm-2.noSTR.tru_con.vcf.gz >${polarizing}/conSV_noSTR.txt
ls ${polarizing}/dsim.final.vcf.gz >>${polarizing}/conSV_files.txt
bcftools merge -m none --threads ${nT} `cat ${polarizing}/conSV_files.txt` |\
bcftools sort --max-mem 4G  |bgzip -@ ${nT} > ${polarizing}/allandsim.asm.noSTR.sort.vcf.gz
bcftools index --threads ${nT} -t -f ${polarizing}/allandsim.asm.noSTR.sort.vcf.gz


module load python/3.10.2
truvari collapse -k common --sizemax 5000000 --sizemin 50 --refdist 700 --pctseq 0.85 --pctsize 0.85 \
-i ${polarizing}/allandsim.asm.noSTR.sort.vcf.gz -o ${polarizing}/allandsim.asm.noSTR.truvari.vcf.gz \
-c ${polarizing}/all2sim.asm.noSTR.collapsed.vcf
module unload python/3.10.2
printf "dsim_mumco2">${polarizing}/filter.txt
bcftools filter --threads ${nT} -i 'GT[@filter.txt]="hom"' ${polarizing}/allandsim.asm.noSTR.truvari.vcf.gz |\
bcftools view --threads ${nT} -i 'NumCollapsed >= 1'|\
bcftools sort -O z >  ${polarizing}/all2sim.asm.noSTR.onlysim.vcf.gz
bcftools index -t -f ${polarizing}/all2sim.asm.noSTR.onlysim.vcf.gz
# --> zcat ${polarizing}/all2sim.asm.noSTR.onlysim.vcf.gz |grep -v "#"|wc -l
# --> ????

bedtools subtract -A -header -a ${merged_SVs}/truvari.asm-2.noSTR.vcf.gz -b ${polarizing}/all2sim.asm.noSTR.onlysim.vcf.gz |\
bcftools sort -O v -o ${polarizing}/nochange.noSTR.vcf

bedtools subtract -A -a ${merged_SVs}/truvari.asm-2.noSTR.vcf.gz -b ${polarizing}/nochange.noSTR.vcf |\
sed s@'\.\/\.'@'hahaha'@g |sed s@'1\/1'@'\.\/\.'@g|sed s@'hahaha'@'1\/1'@g >${polarizing}/reversed.noSTR.txt

##################################################
#####Important!!! Drrmove the duplicated SVs######
##################################################
cat ${polarizing}/nochange.noSTR.vcf ${polarizing}/reversed.noSTR.txt |bcftools sort --max-mem 2G -O v|uniq|\
bgzip -@ ${nT} -c > ${polarizing}/polarized.asm.noSTR.vcf.gz
bcftools index -f -t ${polarizing}/polarized.asm.noSTR.vcf.gz