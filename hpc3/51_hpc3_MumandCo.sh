#!/bin/bash

#SBATCH --job-name=mumco    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=10   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU
#SBATCH --tmp=100G                ## requesting 100 GB local scratch
#SBATCH --constraint=fastscratch  ## requesting nodes with fast scratch in /tmp

source ~/.bashrc

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
assemble="/dfs7/jje/jenyuw/SV-project-temp/result/assemble"
scaffold="/dfs7/jje/jenyuw/SV-project-temp/result/scaffold"
nT=$SLURM_CPUS_PER_TASK

file=`head -n $SLURM_ARRAY_TASK_ID ${scaffold}/scfd_list.txt |tail -n 1`
name=$(echo ${file} | cut -d '/' -f 8)
echo "the file is ${file}"
echo "the name is ${name}"

cd ${SVs}

bash /pub/jenyuw/Software/MUMandCo-MUMandCov3.8/mumandco_v3.8.sh  \
-r ${ref_genome} -q ${scaffold}/${name}/ragtag.scaffold.fasta \
-g 135000000 -o ${name}_mumco -t ${nT} -ml 50
# Do not use -b

cd ${SVs}/${name}_mumco_output

cat ${SVs}/${name}_mumco_output/${name}_mumco.SVs_all.vcf|grep "#" |\
#tr '\n' '\f' |\
##'\f' is a special character that dpes not exist in the file, so it is used as a temporary separator
#sed -e 's/\f##contig/>\f##contig/g;s/\f##query=/>\f##query=/;s/.fasta>/.fasta/' |\
#sed -e 's/\f##query_contig/>\f##query_contig/g;s/\f##INFO=<ID=END/>\f##INFO=<ID=END/;s/.fasta>/.fasta/' |\
#tr '\f' '\n' |\
grep -v "##query_contig" |sed 's/ type=.*.;,/,/g' |tr -d " "|\
sed 's/qCHR/Q_CHROM/g; s/qSTART/Q_START/g; s/qEND/Q_END/g' |\
sed -E "s@##ALT=<ID=CONTR,Description=Contraction>\n@@g" |\
sed 's@ALT=<ID=DUP:TANDEM,Description=Tandem Duplication>@ALT=<ID=DUP,Description=Duplication>@g' >${SVs}/${name}_mumco_output/header.pre

## making VCF
grep -v "#" ${SVs}/${name}_mumco_output/${name}_mumco.SVs_all.vcf |\
sed 's/CONTR/DEL/g; s/qCHR/Q_CHROM/g; s/qSTART/Q_START/g; s/qEND/Q_END/g'|\
gawk '{print $1 "\t" $2 "\t" $1"_"$2 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 }' >${SVs}/${name}_mumco_output/container.pre
#contractions are treated as deletions!

cat ${SVs}/${name}_mumco_output/header.pre ${SVs}/${name}_mumco_output/container.pre >${SVs}/${name}_mumco_output/${name}_mumco.good.vcf
bgzip -@ ${nT} -f -k ${SVs}/${name}_mumco_output/${name}_mumco.good.vcf
bcftools sort -O z ${SVs}/${name}_mumco_output/${name}_mumco.good.vcf.gz >${SVs}/${name}_mumco.good.sort.vcf.gz
bcftools index -t -f ${SVs}/${name}_mumco.good.sort.vcf.gz

echo "This is the end!!"
