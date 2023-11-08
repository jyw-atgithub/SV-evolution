#!/bin/bash

#SBATCH --job-name="fil&com"    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=24   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=4G     # requesting memory per CPU
source ~/.bashrc

SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
nT=$SLURM_CPUS_PER_TASK

file=`head -n $SLURM_ARRAY_TASK_ID ${aligned_bam}/alignedlist.txt |tail -n 1`
name=$(basename ${file}|sed s/".trimmed-ref.sort.bam"//g)

#We are using bcftools v18 and samtools v18!! 

for i in `ls ${SVs}/${name}.*.vcf `
do
prog=`basename ${i} | gawk -F "." '{print $2}' `
bgzip -f --keep -@ ${nT} ${i}
# Only the vcf from SVIM is not sorted while others are sorted. We sort all because of convenience.
bcftools sort --write-index --max-mem 10000M -O z -o ${SVs}/${name}.${prog}.sort.vcf.gz ${i}.gz
#tabix -f -p vcf ${SVs}/${name}.${prog}.sort.vcf.gz #no need, because "bcftools sort --write-index" generate .csi index

bcftools view --threads ${nT} -r 2L,2R,3L,3R,4,X,Y \
-i 'QUAL >= 10 && FILTER = "PASS"'  -O v -o - ${SVs}/${name}.${prog}.sort.vcf.gz |\
sed 's/SVTYPE=DUP:INT/SVTYPE=DUP/g ; s/SVTYPE=DUP:TANDEM/SVTYPE=DUP/g ' |\
bcftools view --threads ${nT} -O v -o ${SVs}/${name}.${prog}.filtered.vcf
#bgzip -f -dk ${SVs}/${name}.${prog}.filtered.vcf.gz

rm ${SVs}/${name}.${prog}.sort.vcf.gz
rm ${SVs}/${name}.${prog}.sort.vcf.gz.csi
done


echo "This is the end!"