#!/bin/bash

#SBATCH --job-name=comb    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=3   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU
source ~/.bashrc

SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
con_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/consensus_SVs"
nT=$SLURM_CPUS_PER_TASK

file=`head -n $SLURM_ARRAY_TASK_ID ${aligned_bam}/alignedlist.txt |tail -n 1`
name=$(basename ${file}|sed s/".trimmed-ref.sort.bam"//g)

#We are using bcftools v18 and samtools v18!! 

module load perl/5.34.1
cd ${con_SVs}

perl /pub/jenyuw/Software/combiSV/combiSV2.2.pl \
-sniffles ${SVs}/${name}.sniffles.filtered.vcf \
-cutesv ${SVs}/${name}.cutesv.filtered.vcf \
-svim ${SVs}/${name}.SVIM.filtered.vcf \
-c 3 -o ${name}

rm ${con_SVs}/*_${name}.vcf

echo " This is the end!"
