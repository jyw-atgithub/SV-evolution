#!/bin/bash

#SBATCH --job-name=cute    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=2-62%1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=16   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=3G     # requesting memory per CPU

###--array=2-62%1 This make the system only run 1 array job at a time
###Let everything goes one by one

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
trimmed="/dfs7/jje/jenyuw/SV-project-temp/result/trimmed"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"

source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
echo "Yes, ARRAY_TASK_ID=1"
ls ${aligned_bam}/*.trimmed-ref.sort.bam >${aligned_bam}/alignedlist.txt
else
echo "No need to list the file again"
fi

file=`head -n $SLURM_ARRAY_TASK_ID ${aligned_bam}/alignedlist.txt |tail -n 1`
name=$(basename ${file}|sed s/".trimmed-ref.sort.bam"//g)
read_type=`echo ${name} | gawk -F "_" '{print $2}'`

echo "name is $name "

module load python/3.10.2

#waitingtime=`expr $SLURM_ARRAY_TASK_ID \* 2 \* \( $RANDOM % 30 \)`
#sleep ${waitingtime}
#remember to put back slash for escape characters
#due to cuteSV's behavior. the same intermediate files cannot exist in the same directory at the same time. 
#This does NOT allow real parallel computation, so let's wait a little bit. 
##This method failed.

cd ${SVs} # for unknown reason, cuteSV does save the file at our desired directory
cuteSV --threads ${nT} --genotype --sample ${name}_cute \
--min_support 10 \
--min_size 50 --min_mapq 20 --min_read_len 500 \
-L '-1' \
--merge_del_threshold 270 --merge_ins_threshold 270 \
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
 "${file}" "${ref_genome}" "${name}_cutesv.vcf" "${SVs}"

module unload python/3.10.2