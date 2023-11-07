#!/bin/bash

#SBATCH --job-name=svim    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=8   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
trimmed="/dfs7/jje/jenyuw/SV-project-temp/result/trimmed"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"

source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK

file=`head -n $SLURM_ARRAY_TASK_ID ${aligned_bam}/alignedlist.txt |tail -n 1`
name=$(basename ${file}|sed s/".trimmed-ref.sort.bam"//g)
read_type=`echo ${name} | gawk -F "_" '{print $2}'`

echo "name is $name "

module load python/3.10.2
#module load samtools/1.15.1
#"svim reads"requires samtools, but "svim alignment" does not
##some parameter names are changed. Check the svim alignment --help messages.
svim alignment --sample ${name}-svim \
--min_mapq 20 --min_sv_size 50 \
--max_sv_size 10000000 \
--position_distance_normalizer 900 --cluster_max_distance 0.3 \
${SVs}/${name}_SVIM ${file} ${ref_genome}

cp ${SVs}/${name}_SVIM/variants.vcf ${SVs}/${name}.SVIM.vcf

module unload python/3.10.2

echo "This is the end"