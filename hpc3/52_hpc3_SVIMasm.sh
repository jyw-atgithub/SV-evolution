#!/bin/bash

#SBATCH --job-name=svim-asm    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p highmem       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=12   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU

##############################################################################
###This ia a crude version. We need to feed descent assemblies afterwards. ###
##############################################################################
source ~/.bashrc

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
assemble="/dfs7/jje/jenyuw/SV-project-temp/result/assemble"
scaffold="/dfs7/jje/jenyuw/SV-project-temp/result/scaffold"
nT=$SLURM_CPUS_PER_TASK

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
echo "Yes, ARRAY_TASK_ID=1"
ls ${assemble}/*_Flye/assembly.fasta > ${assemble}/Flyelist.txt
else
echo "No need to list the file again"
fi

file=`head -n $SLURM_ARRAY_TASK_ID ${assemble}/Flyelist.txt |tail -n 1`
name=$(echo ${file} | cut -d '/' -f 8 |cut -d '_' -f 1-2)

conda activate sv-calling

ragtag.py scaffold -u -r -g 50 -t ${nT} -o ${scaffold}/rr${SLURM_ARRAY_TASK_ID}/ \
${ref_genome} ${file}
#--aligner nucmer --nucmer-params "--maxmatch -l 100 -c 500 -t ${nT}"
##I don't know why, but nucmer canNOT be multithreaded. 

mv ${scaffold}/rr${SLURM_ARRAY_TASK_ID}/ragtag.scaffold.fasta ${scaffold}/${name}.scaffold.fasta


minimap2 -t ${nT} -a -x asm5 --cs --eqx \
${ref_genome} ${scaffold}/${name}.scaffold.fasta \
|samtools view -b -h -@ ${nT} -o -|samtools sort -@ ${nT} -o ${aligned_bam}/${name}.Flye-ref.sort.bam
samtools index ${aligned_bam}/${name}.Flye-ref.sort.bam

svim-asm haploid --sample ${name} --min_sv_size 50 \
${SVs}/${name}_svimASM ${aligned_bam}/${name}.Flye-ref.sort.bam ${ref_genome}

mv ${SVs}/${name}_svimASM/variants.vcf ${SVs}/${name}.svimASM.vcf 
conda deactivate