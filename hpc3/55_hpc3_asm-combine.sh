#!/bin/bash

#SBATCH --job-name=asm-combine    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p highmem       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=4   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=10G     # requesting memory per CPU

source ~/.bashrc

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
con_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/consensus_SVs"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
assemble="/dfs7/jje/jenyuw/SV-project-temp/result/assemble"
scaffold="/dfs7/jje/jenyuw/SV-project-temp/result/scaffold"
nT=$SLURM_CPUS_PER_TASK

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
ls ${scaffold}/*/ragtag.scaffold.fasta > ${scaffold}/scfd_list.txt
fi

file=`head -n $SLURM_ARRAY_TASK_ID ${scaffold}/scfd_list.txt |tail -n 1`
name=$(echo ${file} | cut -d '/' -f 8)
echo "the name is ${name}"

bcftools index -f -t ${SVs}/${name}.mumco.filtered.vcf.gz
bcftools index -f -t ${SVs}/${name}.svimASM.filtered.vcf.gz

bcftools merge -m none ${SVs}/${name}.mumco.filtered.vcf.gz ${SVs}/${name}.svimASM.filtered.vcf.gz |\
bcftools sort -m 2G -O z -o ${con_SVs}/${name}.asm-2.vcf.gz
bcftools index -f -t ${con_SVs}/${name}.asm-2.vcf.gz

#Truvari requires a .tbi index
# truvari can output to stdout
#--intra is only provided later than v4.2 (experimental)
module load python/3.10.2

truvari collapse --intra -k maxqual --sizemax 200000000 \
-i ${con_SVs}/${name}.asm-2.vcf.gz \
-c ${con_SVs}/${name}.asm-2.tru_collapsed.vcf -f ${ref_genome} |\
bcftools sort -m 2G |bgzip -@ ${nT} > ${con_SVs}/${name}.asm-2.tru_con.sort.vcf.gz

bcftools index -f -t ${con_SVs}/${name}.asm-2.tru_con.sort.vcf.gz
module unload python/3.10.2

#bcftools +fixref ${con_SVs}/A1_CLR.tru_con.sort.vcf.gz -- -f ${ref_genome}
#bcftools norm --check-ref e -f ${ref_genome} ${con_SVs}/A1_CLR.tru_con.sort.vcf.gz -O v -o /dev/null
#bcftools +fixref ${con_SVs}/${name}.tru_con.sort.vcf.gz -O z -o ${con_SVs}/${name}.tru_con.fix.vcf.gz -- -f ${ref_genome} -m top
#bcftools index -f -t ${con_SVs}/${name}.tru_con.fix.vcf.gz
#truvari collapse -r 500 -p 0.95 -P 0.95 -s 50 -S 100000 #A draft human pangenome reference

echo " This is the end!"