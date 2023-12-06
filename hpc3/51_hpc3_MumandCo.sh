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

cat ${name}_mumco.SVs_all.vcf|tr '\n' '\f' |\
sed -e 's/\f##contig/>\f##contig/g;s/\f##query=/>\f##query=/;s/.fasta>/.fasta/' |\
sed -e 's/\f##query_contig/>\f##query_contig/g;s/\f##INFO=<ID=END/>\f##INFO=<ID=END/;s/.fasta>/.fasta/' |\
tr '\f' '\n' |\
sed 's/ type=.*.;,/,/g' |grep "#" |grep -v "##ALT" >header.pre

## making VCF
printf "" >container.pre &&\
cat ${name}_mumco.SVs_all.withfragment.tsv |\
gawk ' (NR>1) && ($6=="insertion_mobile") {print $1 "\t" $3 "\t"$6"-"$5"-"$7 "\t" "." "\t" $9 "\t" "." "\t" "PASS" "\t" "SVTYPE=INS" ";SVLEN=" $5 ";qCHR=" $2 ";qSTART=" $7 ";qEND=" $8 "\t" "GT" "\t" "1/1"} ' >>container.pre &&\
cat ${name}_mumco.SVs_all.withfragment.tsv |\
gawk ' (NR>1) && ($6=="insertion_novel") {print $1 "\t" $3 "\t"$6"-"$5"-"$7 "\t" "." "\t" $9 "\t" "." "\t" "PASS" "\t" "SVTYPE=INS" ";SVLEN=" $5 ";qCHR=" $2 ";qSTART=" $7 ";qEND=" $8 "\t" "GT" "\t" "1/1"} '>>container.pre &&\
cat ${name}_mumco.SVs_all.withfragment.tsv |\
gawk ' (NR>1) && ($6=="deletion_mobile") {print $1 "\t" $3 "\t"$6"-"$5"-"$7 "\t" $9 "\t" "." "\t" "." "\t" "PASS" "\t" "SVTYPE=DEL" ";SVLEN=" $5 ";qCHR=" $2 ";qSTART=" $7 ";qEND=" $8 "\t" "GT" "\t" "1/1"} '>>container.pre &&\
cat ${name}_mumco.SVs_all.withfragment.tsv |\
gawk ' (NR>1) && ($6=="deletion_novel") {print $1 "\t" $3 "\t"$6"-"$5"-"$7 "\t" $9 "\t" "." "\t" "." "\t" "PASS" "\t" "SVTYPE=DEL" ";SVLEN=" $5 ";qCHR=" $2 ";qSTART=" $7 ";qEND=" $8 "\t" "GT" "\t" "1/1"} '>>container.pre &&\
cat ${name}_mumco.SVs_all.withfragment.tsv |\
gawk ' (NR>1) && ($6=="duplication") {print $1 "\t" $3 "\t"$6"-"$5"-"$7 "\t" $9 "\t" "." "\t" "." "\t" "PASS" "\t" "SVTYPE=DUP" ";SVLEN=" $5 ";qCHR=" $2 ";qSTART=" $7 ";qEND=" $8 "\t" "GT" "\t" "1/1"} '>>container.pre &&\
cat ${name}_mumco.SVs_all.withfragment.tsv |\
gawk ' (NR>1) && ($6=="inversion") {print $1 "\t" $3 "\t"$6"-"$5"-"$7 "\t" $9 "\t" "." "\t" "." "\t" "PASS" "\t" "SVTYPE=INV" ";SVLEN=" $5 ";qCHR=" $2 ";qSTART=" $7 ";qEND=" $8 "\t" "GT" "\t" "1/1"} '>>container.pre &&\
cat ${name}_mumco.SVs_all.withfragment.tsv |\
gawk ' (NR>1) && ($6=="transloc") {print $1 "\t" $3 "\t"$6"-"$5"-"$7 "\t" $9 "\t" "." "\t" "." "\t" "PASS" "\t" "SVTYPE=TRA" ";SVLEN=" $5 ";qCHR=" $2 ";qSTART=" $7 ";qEND=" $8 "\t" "GT" "\t" "1/1"} '>>container.pre &&\
cat ${name}_mumco.SVs_all.withfragment.tsv |\
gawk ' (NR>1) && ($6=="contraction") {print $1 "\t" $3 "\t"$6"-"$5"-"$7 "\t" $9 "\t" "." "\t" "." "\t" "PASS" "\t" "SVTYPE=DEL" ";SVLEN=" $5 ";qCHR=" $2 ";qSTART=" $7 ";qEND=" $8 "\t" "GT" "\t" "1/1"} '>>container.pre
#contractions are treated as deletions!


cat header.pre container.pre >${name}_mumco.good.vcf
bgzip -@ ${nT} -k ${name}_mumco.good.vcf
bcftools sort -O z ${name}_mumco.good.vcf.gz >${name}_mumco.good.sort.vcf.gz
bcftools index -t -f ${name}_mumco.good.sort.vcf.gz

#cat A4_CLR.SVs_all.withfragment.tsv |gawk ' NR>1 {print $6} '|sort -h|uniq -c|less -S
echo "We will have a good vcf!!"
echo "This is the end!!"
