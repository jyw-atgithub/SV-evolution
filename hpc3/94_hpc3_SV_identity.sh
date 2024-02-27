#!/bin/bash

#SBATCH --job-name=SV_identity    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p highmem      ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=20   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=10G     # requesting memory per CPU

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
trimmed="/dfs7/jje/jenyuw/SV-project-temp/result/trimmed"
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
TE="/dfs7/jje/jenyuw/SV-project-temp/result/TE_repeat"
source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK

: <<'SKIP'
########################This may cause redudant!!########################
cat ${trimmed}/namelist.txt|while read i
do
name=$(basename ${i} |sed s/".trimmed.fastq.gz"//g)
#TRA always gives error
################## bcftools consensus skips overlapping SVs, so we don't use it ##################
#bcftools consensus -i 'SVTYPE="DEL"' -s ${name}_mumco --fasta-ref ${ref_genome} --haplotype LR ${polarizing}/polarized.asm.vcf.gz >${TE}/${name}.SV.fa
#for DUP and DEL, just extract the REF allele
echo "Processing ${name}"
bcftools query -i 'SVTYPE="DUP" || SVTYPE="DEL"' -f '>%CHROM-%POS-%END-%ID\n%REF\n' \
-s ${name}_mumco ${polarizing}/polarized.asm.vcf.gz >${TE}/${name}.SV.fa
#for INS and INV, just extract the alt allele
bcftools query -i 'SVTYPE="INS" || SVTYPE="INV"' -f '>%CHROM-%POS-%END-%ID\n%ALT\n' \
-s ${name}_mumco ${polarizing}/polarized.asm.vcf.gz >>${TE}/${name}.SV.fa
done
cat ${TE}/*.SV.fa |gawk -F ";" '{print $1 }' >${TE}/SV.fa
#rm ${TE}/*.SV.fa
SKIP

#If use "--fasta-ref $ref_genome", "bcftools norm" will fail
###There are many repeated IDs of SVs representiting different alleles, so we have to rename the IDs
bcftools annotate --set-id '%CHROM-%POS-%INFO/SVTYPE-%INFO/CollapseId-%ID' ${polarizing}/polarized.asm.vcf.gz |\
bcftools sort --max-mem 4G -O v |\
##remove duplicated SVs and retain the first entry
bcftools norm --rm-dup all -O v |\
##shorten super long ID.
##[^\t]* meand any character except tab
sed 's/;svim[^\t]*\t/\t/g' |sed 's/svim_asm.//g'|bgzip -@ ${nT} -c> ${polarizing}/3corrected.polarized.asm.vcf.gz
bcftools index -f -t ${polarizing}/3corrected.polarized.asm.vcf.gz

#grep -v "#"|gawk '{print $3}'|uniq -d -c|less

bcftools query -i 'SVTYPE="DUP" || SVTYPE="DEL"' -f '>%ID\n%REF\n' \
${polarizing}/3corrected.polarized.asm.vcf.gz >${TE}/SV.DUP-DEL.fa
#for INS and INV, just extract the alt allele
bcftools query -i 'SVTYPE="INS" || SVTYPE="INV"' -f '>%ID\n%ALT\n' \
${polarizing}/3corrected.polarized.asm.vcf.gz >${TE}/SV.INS-INV.fa
cat ${TE}/SV.DUP-DEL.fa ${TE}/SV.INS-INV.fa >${TE}/SV.fa

module load singularity/3.11.3
singularity exec -B ${TE} \
-B /dfs7/jje/jenyuw/SV-project-temp/raw/Libraries:/opt/RepeatMasker/Libraries \
/pub/jenyuw/Software/dfam-tetools-latest.sif \
RepeatMasker -gff -s -xsmall \
-species "drosophila melanogaster" -dir ${TE}/extracted_SV3 \
${TE}/SV.fa

##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample
#2L-1055-DEL-1.1-DEL.3
cat ${TE}/extracted_SV3/SV.fa.out|gawk ' NR > 3 {print $0}' | while read line
do
chr_pos_id_type=`echo $line|gawk '{print $5}'|gawk -F "-" '{print $1 "\t" $2 "\t" $0 "\t" $3 "\t"}'`
order_supfam=`echo $line|gawk '{ if ( $11=="Simple_repeat") print $11 "\t" $10; else print $11}'|sed 's/\//\t/g;s/(//g;s/)n//g'`
echo -e "${chr_pos_id_type}${order_supfam}"
done |gawk ' { if( NF<6 ) print $0 "\t" $5 ; else print $0}'|sort -k 4,4 -k 5,5 -k 6,6 |uniq >${TE}/3SV_info.tsv
#cat SV_info.tsv|gawk ' { if( NF<5 ) print $0 "\t" $4 ; else print $0}'|sort -k 3,3 -k 4,4 -k 5,5|uniq >2SV_info.tsv

cat 3SV_info.tsv|sort -n -k 1,1 -k 2,2|gawk ' $1=="2L" {print $0}'  >SV_info.good.tsv
cat 3SV_info.tsv|sort -n -k 1,1 -k 2,2|gawk ' $1=="2R" {print $0}'  >>SV_info.good.tsv
cat 3SV_info.tsv|sort -n -k 1,1 -k 2,2|gawk ' $1=="3L" {print $0}'  >>SV_info.good.tsv
cat 3SV_info.tsv|sort -n -k 1,1 -k 2,2|gawk ' $1=="3R" {print $0}'  >>SV_info.good.tsv
cat 3SV_info.tsv|sort -n -k 1,1 -k 2,2|gawk ' $1=="4" {print $0}'  >>SV_info.good.tsv
cat 3SV_info.tsv|sort -n -k 1,1 -k 2,2|gawk ' $1=="X" {print $0}'  >>SV_info.good.tsv
cat 3SV_info.tsv|sort -n -k 1,1 -k 2,2|gawk ' $1=="Y" {print $0}'  >>SV_info.good.tsv
bgzip -k -f -@ ${nT} ${TE}/SV_info.good.tsv
tabix -f -s 1 -b 2 -e 2 ${TE}/SV_info.good.tsv.gz


#head -n 100 ${TE}/SV_info.tsv|sort -k 3,3 -k 4,4 -k 5,5|uniq |bgzip -c >${TE}/fake.info.tsv.gz
#tabix -f -s 1 -b 2 -e 2 ${TE}/fake.info.tsv.gz
#cat SV_info.tsv|gawk ' { if( NF<5 ) print $0 "\t" $4 ; else print $0}' >2SV_info.tsv
echo -e "##INFO=<ID=ORDER,Number=1,Type=String,Description=\"The order of the transposable element or a simple repeat\">
##INFO=<ID=SUPERFAMILY,Number=1,Type=String,Description=\"The SUPERFAMILY of the transposable element or the motif of a simple repeat\">" \
>${TE}/add_header.txt
zcat ${polarizing}/3corrected.polarized.asm.vcf.gz|\
bcftools annotate -a ${TE}/SV_info.good.tsv.gz -c 'CHROM,POS,ID,INFO/ORDER,INFO/SUPERFAMILY' --header-lines ${TE}/add_header.txt|less -S

#zcat ${polarizing}/corrected.polarized.asm.vcf.gz|head -n 10000| bcftools annotate -a ${TE}/fake.info.tsv.gz -c 'CHROM,POS,ID,INFO/ORDER,INFO/SUPERFAMILY' --header-lines ${TE}/add_header.txt|less -S

#1.CHROM   2.POS     3.ID 4.SVTYPE 5.ORDER 6.SUPERFAMILY 
#1.CHROM    2.POS 3.ID [ %GT]......
##split the VCF
bcftools query -f '%CHROM\t%POS\t%ID[ %GT]\n' ${polarizing}/3corrected.polarized.asm.vcf.gz >${TE}/SV_genotype.tsv
join -1 3 -2 3  <(sort -k 3,3 ${TE}/SV_info.good.tsv) <(sort -k 3,3 ${TE}/SV_genotype.tsv)|\
cut -d " " -f 2,3,1,4,5,8- |tr " " "\t" >${TE}/SV_genotype_info.tsv

module load R/4.2.2
Rscript -e '
library(dplyr)
library(tidyr)
t=read.table("SV_info.good.tsv",header=FALSE)
g=read.table("SV_genotype.tsv",header=FALSE)
colnames(t) <- c("CHROM","POS","ID","SVTYPE","ORDER","SUPERFAMILY")
colnames(g) <- c("CHROM","POS","ID", colnames(g)[4:ncol(g)])
t1=t %>% group_by(CHROM,POS,ID,SVTYPE)%>% summarise_at(vars(ORDER:SUPERFAMILY),paste, collapse=",")
all=g %>% left_join(t1,by=c("CHROM","POS","ID")) %>% replace_na(list(ORDER="not_repeat",SUPERFAMILY="not_repeat")) %>% 
relocate(SVTYPE, .after =ID) %>% relocate(ORDER, .after =SVTYPE)%>% relocate(SUPERFAMILY, .after =ORDER) 
head(all,25)
write.table(all,"repeat_type_genotype.tsv",quote=FALSE,sep="\t",row.names=TRUE)
'


: <<'SKIP'
cat ${TE}/extracted_SV2/SV.fa.out|gawk ' NR > 3 {print $0}' | while read line
do
chr_pos_id=`echo $line|gawk '{print $5}'|gawk -F "-" '{print $1 "\t" $2 "\t" $0 "\t"}'`
order_supfam=`echo $line|gawk '{ if ( $11=="Simple_repeat") print "ORDER="$11 ";" "SUPERFAMILY="$10; else print "ORDER="gensub("/",";SUPERFAMILY=","g",$11)}'|sed 's/(//g;s/)n//g'`
echo -e "${chr_pos_id}.\t.\t.\t.\t${order_supfam}\t."
done >${TE}/newbody.tsv

#-e The end column can be the same as the start column.

zcat ${polarizing}/corrected.polarized.asm.vcf.gz|grep "##" >${TE}/new_header.txt
echo -e "##INFO=<ID=ORDER,Number=1,Type=String,Description=\"The order of the transposable element or a simple repeat\">
##INFO=<ID=SUPERFAMILY,Number=1,Type=String,Description=\"The SUPERFAMILY of the transposable element or the motif of a simple repeat\">" \
>>${TE}/new_header.txt
zcat ${polarizing}/corrected.polarized.asm.vcf.gz|grep "#"|grep -v "##"|gawk 'BEGIN { OFS = "\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, "sample" }' >>${TE}/new_header.txt
##Keep failing parsing the new VCF
cat ${TE}/new_header.txt ${TE}/newbody.tsv |bcftools sort -O v >${TE}/annotation.vcf
SKIP