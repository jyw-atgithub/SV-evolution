#!/bin/bash

#SBATCH --job-name=comb    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=12   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU
source ~/.bashrc

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
con_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/consensus_SVs"
nT=$SLURM_CPUS_PER_TASK

file=`head -n $SLURM_ARRAY_TASK_ID ${aligned_bam}/alignedlist.txt |tail -n 1`
name=$(basename ${file}|sed s/".trimmed-ref.sort.bam"//g)

#We are using bcftools v18 and samtools v18!!
## combiSV is not a good tool here. It's output content is too simple!!

ls ${SVs}/${name}.*.filtered.vcf.gz
for i in `ls ${SVs}/${name}.*.filtered.vcf`
do
bgzip -@ ${nT} -f -k ${i}
bcftools index -f -t ${i}.gz
done

bcftools merge -m none ${SVs}/${name}.*.filtered.vcf.gz | bcftools sort -m 2G -O z -o ${con_SVs}/${name}.3.vcf.gz
bcftools index -f -t ${con_SVs}/${name}.3.vcf.gz

#Truvari requires a .tbi index
# truvari can output to stdout
#--intra is only provided later than v4.2 (experimental)
module load python/3.10.2

truvari collapse --intra -k maxqual --sizemax 200000000 \
-i ${con_SVs}/${name}.3.vcf.gz \
-c ${con_SVs}/${name}.tru_collapsed.vcf -f ${ref_genome} |\
bcftools sort -m 2G |bgzip -@ ${nT} > ${con_SVs}/${name}.tru_con.sort.vcf.gz

bcftools index -f -t ${con_SVs}/${name}.tru_con.sort.vcf.gz
module unload python/3.10.2

#bcftools +fixref ${con_SVs}/A1_CLR.tru_con.sort.vcf.gz -- -f ${ref_genome}
#bcftools norm --check-ref e -f ${ref_genome} ${con_SVs}/A1_CLR.tru_con.sort.vcf.gz -O v -o /dev/null
bcftools +fixref ${con_SVs}/${name}.tru_con.sort.vcf.gz -O z -o ${con_SVs}/${name}.tru_con.fix.vcf.gz -- -f ${ref_genome} -m top
bcftools index -f -t ${con_SVs}/${name}.tru_con.fix.vcf.gz
#truvari collapse -r 500 -p 0.95 -P 0.95 -s 50 -S 100000 #A draft human pangenome reference

echo " This is the end!"