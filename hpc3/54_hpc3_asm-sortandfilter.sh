#!/bin/bash

#SBATCH --job-name="ASMst&fil"    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p highmem      ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=6   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU
source ~/.bashrc
trimmed="/dfs7/jje/jenyuw/SV-project-temp/result/trimmed"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
nT=$SLURM_CPUS_PER_TASK

file=`head -n $SLURM_ARRAY_TASK_ID ${trimmed}/namelist.txt |tail -n 1`
name=$(basename ${file}|sed s/".trimmed.fastq.gz"//g)
read_type=`echo ${name} | gawk -F "_" '{print $2}'`
echo "the file is ${file}"
echo "the name is ${name}"

#####The following sed commands contains literally tab, so copy and paste will not work. Edit it in nano manually.########
#####Truvari requires the quality score to be an integer and it can not be only "."########

for i in ${SVs}/${name}.mumco.good.sort.vcf.gz ${SVs}/${name}.svimASM.sort.gz
do
prog=`basename ${i} | gawk -F "." '{print $2}'`
bcftools view --threads ${nT} -r 2L,2R,3L,3R,4,X,Y \
-i 'FILTER = "PASS"' -O v -o - ${i} |\
sed 's/DUP_TANDEM/DUP/g; s/DUP:TANDEM/DUP/g; s/DUP_INT/DUP/g; s/DUP:INT/DUP/g; s/BND/TRA/g' |\
sed 's/ .       PASS/   30      PASS/g' |\
grep -v "NNNNNNNNNNNNNNNNNNNN" |bgzip -@ ${nT} -c > ${SVs}/${name}.${prog}.filtered.vcf.gz
done

##We only want the DUP and INS from mumco
##Change the quality score of mumco vcf. Better quality for duplication
mv ${SVs}/${name}.mumco.filtered.vcf.gz ${SVs}/${name}.mumco.filtered.ori.vcf.gz
zcat ${SVs}/${name}.mumco.filtered.ori.vcf.gz|grep  "#" |bgzip -@ ${nT} -c > ${SVs}/${name}.mumco.filtered.vcf.gz.header
zcat ${SVs}/${name}.mumco.filtered.ori.vcf.gz|grep -v "#" |gawk -v OFS='\t' '{
    if ($5=="\<DUP\>" || $5=="\<INS\>") 
    {
        print $1,$2,$3,$4,$5,"40",$7,$8,$9,$10
        }
    }'|\
grep -v "NNNNNNNNNNNNNNNNNNNN" |\
bgzip -@ ${nT} -c > ${SVs}/${name}.mumco.filtered.vcf.gz.body
cat ${SVs}/${name}.mumco.filtered.vcf.gz.header ${SVs}/${name}.mumco.filtered.vcf.gz.body > ${SVs}/${name}.mumco.filtered.vcf.gz
rm ${SVs}/${name}.mumco.filtered.vcf.gz.header ${SVs}/${name}.mumco.filtered.vcf.gz.body
#rm ${SVs}/${name}*.csi

echo "This is the end!"
