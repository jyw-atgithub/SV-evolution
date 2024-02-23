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
bcftools annotate --set-id '%CHROM-%POS-%INFO/CollapseId-%ID' ${polarizing}/polarized.asm.vcf.gz |\
bcftools sort --max-mem 4G -O v |\
##remove duplicated SVs and retain the first entry
bcftools norm --rm-dup all -O v |\
##shorten super long ID.
##[^\t]* meand any character except tab
sed 's/;svim[^\t]*\t/\t/g' |bgzip -@ ${nT} -c> ${polarizing}/corrected.polarized.asm.vcf.gz
bcftools index -f -t ${polarizing}/corrected.polarized.asm.vcf.gz

#grep -v "#"|gawk '{print $3}'|uniq -d -c|less

echo "Processing DUP & DEL"
bcftools query -i 'SVTYPE="DUP" || SVTYPE="DEL"' -f '>%ID\n%REF\n' \
${polarizing}/corrected.polarized.asm.vcf.gz >${TE}/SV.DUP-DEL.fa
#for INS and INV, just extract the alt allele
echo "Processing INS & INV"
bcftools query -i 'SVTYPE="INS" || SVTYPE="INV"' -f '>%ID\n%ALT\n' \
${polarizing}/corrected.polarized.asm.vcf.gz >${TE}/SV.INS-INV.fa
cat ${TE}/SV.DUP-DEL.fa ${TE}/SV.INS-INV.fa >${TE}/SV.fa

module load singularity/3.11.3
singularity exec -B ${TE} \
-B /dfs7/jje/jenyuw/SV-project-temp/raw/Libraries:/opt/RepeatMasker/Libraries \
/pub/jenyuw/Software/dfam-tetools-latest.sif \
RepeatMasker -gff -s -xsmall \
-species "drosophila melanogaster" -dir ${TE}/extracted_SV2 \
${TE}/SV.fa

##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample

head -n 100 ${TE}/extracted_SV2/SV.fa.out|gawk ' NR > 3 {print $0}' | while read line
do
chr_pos_id=`echo $line|gawk '{print $5}'|gawk -F "-" '{print $1 "\t" $2 "\t" $0 "\t"}'`
order_supfam=`echo $line|gawk '{ if ( $11=="Simple_repeat") print $11 "\t" $10; else print $11}'|sed 's/\//\t/g;s/(//g;s/)n//g'`
echo -e "${chr_pos_id}${order_supfam}"
done >${TE}/SV_info.tsv


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

echo -e "##INFO=<ID=ORDER,Number=1,Type=String,Description=\"The order of the transposable element or a simple repeat\">
##INFO=<ID=SUPERFAMILY,Number=1,Type=String,Description=\"The SUPERFAMILY of the transposable element or the motif of a simple repeat\">" \
>>${TE}/add_header.txt
zcat ${polarizing}/corrected.polarized.asm.vcf.gz| bcftools annotate -a ${TE}/SV_info.tsv.gz -c 'CHROM,POS,ID,INFO/ORDER,INFO/SUPERFAMILY' --header-lines ${TE}/add_header.txt|less -S

#1.CHROM   2.POS     3.ID 4.ORDER 5.SUPERFAMILY 
#1.CHROM    2.POS 3.ID [ %GT]
##split the VCF
bcftools query -f '%CHROM\t%POS\t%ID[ %GT]\n' ${polarizing}/corrected.polarized.asm.vcf.gz >${TE}/SV_genotype.tsv
join -1 3 -2 3  <(sort -k 3,3 ${TE}/SV_info.tsv) <(sort -k 3,3 ${TE}/SV_genotype.tsv)|cut -d " " -f 2,3,1,4,5,8- |tr " " "\t" >${TE}/SV_genotype_info.tsv