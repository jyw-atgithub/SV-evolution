#!/bin/bash

#SBATCH --job-name="ASMst&fil"    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p highmem      ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=6   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU
source ~/.bashrc
scaffold="/dfs7/jje/jenyuw/SV-project-temp/result/scaffold"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
nT=$SLURM_CPUS_PER_TASK


if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
ls ${scaffold}/*/ragtag.scaffold.fasta > ${scaffold}/scfd_list.txt
fi

file=`head -n $SLURM_ARRAY_TASK_ID ${scaffold}/scfd_list.txt |tail -n 1`
name=$(echo ${file} | cut -d '/' -f 8)
echo "the name is ${name}"


for i in ${SVs}/${name}.svimASM.vcf 
do
bgzip -f --keep -@ ${nT} ${i}
bcftools sort --write-index --max-mem 4G -O z -o ${i}.sort.gz ${i}.gz
bcftools view --threads ${nT} -r 2L,2R,3L,3R,4,X,Y \
-i 'FILTER = "PASS"' -O v -o - ${i}.sort.gz |\
sed 's/DUP_TANDEM/DUP/g; s/DUP:TANDEM/DUP/g; s/DUP_INT/DUP/g; s/DUP:INT/DUP/g; s/BND/TRA/g' |\
bcftools view --threads ${nT} -O z -o ${SVs}/${name}.svimASM.filtered.vcf.gz
bgzip -f -dk ${SVs}/${name}.svimASM.filtered.vcf.gz
rm ${i}.sort.gz
rm ${i}.sort.gz.csi
done

for i in ${SVs}/${name}.mumco.good.sort.vcf.gz 
do
bcftools view --threads ${nT} -r 2L,2R,3L,3R,4,X,Y \
-i 'FILTER = "PASS"' -O v -o - ${i} |\
sed 's/DUP_TANDEM/DUP/g; s/DUP:TANDEM/DUP/g; s/DUP_INT/DUP/g; s/DUP:INT/DUP/g; s/BND/TRA/g' |\
bcftools view --threads ${nT} -O z -o ${SVs}/${name}.mumco.filtered.vcf.gz
bgzip -f -dk ${SVs}/${name}.mumco.filtered.vcf.gz
done

echo "This is the end!"
