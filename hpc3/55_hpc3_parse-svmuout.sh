#!/bin/bash
#SBATCH --job-name=trs    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=4   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU


# transform svmu table to VCF4.3

##part1, create the header
sample_name="A4_CLR"
echo -e "
" >${sample_name}.header.txt

##part2, create the body
#1. REF_CHROM       2. REF_START       3. REF_END 4. SV_TYPE 5. Q_CHROM 6. Q_START 7. Q_END   8. ID      9. LEN     10. COV_REF COV_Q
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Hifi_60x_0.999_1_cute
input="sv.A4_CLR.txt"

printf "" >${sample_name}.body.txt
#($4 != "nCNV-R")&& ($4 != "CNV-R") Do not involve the nCNV-R and CNV-R (in the reference)
#&& (NR > 1) skip the header
#&& ($9 > 50) SVLEN longer than 50bp
#INS as INS
cat ${input}| gawk ' ($4 != "nCNV-R") && ($4 != "CNV-R") && (NR > 1) && ($9 > 50) {print $0}'|\
gawk ' $4=="INS" {print $1 "\t"  $2 "\t" "INS_" $8 "\t" "REF_seq" "\t" "ALT_seq" "\t" "30" "\t" "PASS" "\t" "SVTYPE=INS;SVLEN=" $9 ";END=" $3 ";Q_CHROM=" $5 ";Q_START=" $6 ";Q_END=" $7 "\t" "GT" "\t" "1/1"}' >>${sample_name}.body.txt
#DEL as DEL
cat ${input}| gawk ' ($4 != "nCNV-R") && ($4 != "CNV-R") && (NR > 1) && ($9 > 50) {print $0}'|\
gawk ' $4=="DEL" {print $1 "\t"  $2 "\t" "DEL_" $8 "\t" "REF_seq" "\t" "ALT_seq" "\t" "30" "\t" "PASS" "\t" "SVTYPE=DEL;SVLEN=" 0 - $9 ";END=" $3 ";Q_CHROM=" $5 ";Q_START=" $6 ";Q_END=" $7 "\t" "GT" "\t" "1/1"}' >>${sample_name}.body.txt
#CNV-Q as DUP
cat ${input}| gawk ' ($4 != "nCNV-R") && ($4 != "CNV-R") && (NR > 1) && ($9 > 50) {print $0}'|\
gawk ' ($4=="CNV-Q") && ($9 > 0) {print $1 "\t"  $2 "\t" "DUP_" $8 "\t" "REF_seq" "\t" "ALT_seq" "\t" "30" "\t" "PASS" "\t" "SVTYPE=DUP;SVLEN=" $9 ";END=" $3 ";Q_CHROM=" $5 ";Q_START=" $6 ";Q_END=" $7 "\t" "GT" "\t" "1/1"}' >>${sample_name}.body.txt
#nCNV-Q aslo as DUP
cat ${input}| gawk ' ($4 != "nCNV-R") && ($4 != "CNV-R") && (NR > 1) && ($9 > 50) {print $0}'|\
gawk ' ($4=="nCNV-Q") && ($9 > 0) {print $1 "\t"  $2 "\t" "DUP_" $8 "\t" "REF_seq" "\t" "ALT_seq" "\t" "30" "\t" "PASS" "\t" "SVTYPE=DUP;SVLEN=" $9 ";END=" $3 ";Q_CHROM=" $5 ";Q_START=" $6 ";Q_END=" $7 "\t" "GT" "\t" "1/1"}' >>${sample_name}.body.txt
##The length of INV were noted ad both Positive and NEgative !!
## only keep the INV with positive length
#INV as INV
cat ${input}| gawk ' ($4 != "nCNV-R") && ($4 != "CNV-R") && (NR > 1) && ($9 > 50) {print $0}'|\
gawk ' ($4=="INV") && ($9 > 0) {print $1 "\t"  $2 "\t" "INV_" $8 "\t" "REF_seq" "\t" "ALT_seq" "\t" "30" "\t" "PASS" "\t" "SVTYPE=INV;SVLEN=" $9 ";END=" $3 ";Q_START=" $6 ";Q_END=" $7 "\t" "GT" "\t" "1/1"}' >>${sample_name}.body.txt

ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
query_genome="/dfs7/jje/jenyuw/SV-project-temp/result/scaffold/A4_CLR/ragtag.scaffold.fasta"

printf "" >  ${sample_name}.good.body.txt

cat ${sample_name}.body.txt|while read line
do
    svtype=`echo $line|gawk '{print $4}'|gawk -F "_" '{print $1}'`
    ref_chrom=`echo $line|gawk '{print $1}'`
    ref_start=`echo $line|gawk '{print $2}'`
    ref_end=`echo $line|gawk '{print $8}'|gawk -F ";" '{print $3}'|gawk -F "=" '{print $2}'`
    q_chrom=`echo $line|gawk '{print $8}'|gawk -F ";" '{print $4}'|gawk -F "=" '{print $2}'`
    q_start=`echo $line|gawk '{print $8}'|gawk -F ";" '{print $5}'|gawk -F "=" '{print $2}'`
    q_end=`echo $line|gawk '{print $8}'|gawk -F ";" '{print $6}'|gawk -F "=" '{print $2}'`
    if [[ ${svtype} == "INV" ]];
    then
##
## Try named pipe for bedtools getfasta
## solve "bash: /usr/bin/sed: Argument list too long"
## solve "malformed BED entry at line 1. Start was greater than end. Exiting." -->> 從上方源頭的gawk修改，不要在這邊用if 判斷
##

    ##The coordinate of inversion is also "inversed"
    echo -e "${ref_chrom}\t${ref_start}\t$((ref_end + 1))">temp.1.bed
    real_r_seq=`bedtools getfasta -fi ${ref_genome} -bed temp.1.bed |gawk '{if (NR==2) print $1}'`
    line=`echo ${line}|sed s@"REF_seq"@${real_r_seq}@`
    ##So, switch the start and end
    echo -e "${q_chrom}\t${q_end}\t$((q_start + 1))">temp.2.bed
    real_q_seq=`bedtools getfasta -fi ${query_genome} -bed temp.2.bed |gawk '{if (NR==2) print $1}'`
    line=`echo ${line}|sed s@"ALT_seq"@${real_q_seq}@`
    echo ${line} >>${sample_name}.good.body.txt
    else
    echo -e "${ref_chrom}\t${ref_start}\t$((ref_end + 1))">temp.1.bed
    real_r_seq=`bedtools getfasta -fi ${ref_genome} -bed temp.1.bed |gawk '{if (NR==2) print $1}'`
    line=`echo ${line}|sed s@"REF_seq"@${real_r_seq}@`
    echo -e "${q_chrom}\t${q_start}\t$((q_end + 1))">temp.2.bed
    real_q_seq=`bedtools getfasta -fi ${query_genome} -bed temp.2.bed |gawk '{if (NR==2) print $1}'`
    line=`echo ${line}|sed s@"ALT_seq"@${real_q_seq}@`
    echo ${line} >>${sample_name}.good.body.txt
    fi
done


    echo ${line} >>${sample_name}.good.body.txt



init=`wc -l ${input}`
final`wc -l ${sample_name}.body.txt`
echo "There are" $init "lines in the input file"
echo "There are" $final "lines in the output file"
echo "Thus, there are" $((init-final)) "lines were filtered out"

cat ${input}| gawk ' ($4=="INV") && ($9 <= 0) {print $0}' > ${sample_name}.suspect_INV.txt
cat ${input}| gawk ' ($4 == "CNV-Q") && ($9 <= 0) {print $0}' > ${sample_name}.suspect_CNVQ.txt