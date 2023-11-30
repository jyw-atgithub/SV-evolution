#!/bin/bash

#SBATCH --job-name="ASMst&fil"    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=12   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU
source ~/.bashrc

SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
nT=$SLURM_CPUS_PER_TASK


if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
ls ${SVs}/*.svimASM.vcf >${SVs}/svimASMlist.txt
fi

file=`head -n $SLURM_ARRAY_TASK_ID ${SVs}/svimASMlist.txt |tail -n 1`
name=`basename ${file}|cut -d "." -f 1`
prog=`basename ${file}| cut -d "." -f 2`
echo "the file is ${file}"
echo "the name is ${name}"
echo "the program is ${prog}"


bgzip -f --keep -@ ${nT} ${file}
bcftools sort --write-index --max-mem 4G -O z -o ${SVs}/${name}.${prog}.sort.vcf.gz ${file}.gz

bcftools view --threads ${nT} -r 2L,2R,3L,3R,4,X,Y \
-i 'FILTER = "PASS"'  -O v -o - ${SVs}/${name}.${prog}.sort.vcf.gz |\
##SVIMasm do NOT give quality scores!
sed 's/DUP_TANDEM/DUP/g; s/DUP:TANDEM/DUP/g; s/DUP_INT/DUP/g; s/DUP:INT/DUP/g' |\
bcftools view --threads ${nT} -O v -o ${SVs}/${name}.${prog}.filtered.vcf
bgzip -f -dk ${SVs}/${name}.${prog}.filtered.vcf.gz

rm ${SVs}/${name}.${prog}.sort.vcf.gz
rm ${SVs}/${name}.${prog}.sort.vcf.gz.csi


echo "This is the end!"
