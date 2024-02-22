#!/bin/bash

#SBATCH --job-name=calling
#SBATCH -A jje_lab
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=6G
#SBATCH --tmp=200G
#SBATCH --constraint=fastscratch

dmel_ref="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
dsim_ref="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing/GCF_016746395.2_Dsim_3.1.fasta"
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
SVs="/dfs7/jje/jenyuw/SV-project-temp/result/SVs"
merged_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/merged_SVs"
con_SVs="/dfs7/jje/jenyuw/SV-project-temp/result/consensus_SVs"
nT=$SLURM_CPUS_PER_TASK
source ~/.bashrc
## In interactive mode

: <<'SKIP'
nucmer -t ${nT} --maxmatch -l 100 -c 500 --sam-long=${polarizing}/mumm4.sam ${dmel_ref} ${dsim_ref}
cat ${polarizing}/mumm4.sam |\
sed 's/HD\ /HD/; s/1.0\ /1.0/; s/\tSO:coordinate/SO:coordinate/; s/VN1/VN:1/; s/HD/HD\t/; s/SO:unsorted/\tSO:unsorted/; s/@PG /@PG\t/; s/ PN/\tPN/; s/ VN/\tVN/; s/ CL/\tCL/' |\
sed 's@\t10\t@\t30\t@g'|\
samtools view -@ ${nT} -h --reference ${dmel_ref} - |samtools sort -@ ${nT} -O bam -o ${polarizing}/corrected.mumm4.bam
samtools index -@ ${nT} ${polarizing}/corrected.mumm4.bam

#nucmer -t ${nT} --maxmatch --delta=${polarizing}/mumm4.delta ${dmel_ref} ${dsim_ref}
#conda activate sv-calling
#lastz ${dmel_ref}[multiple] ${dsim_ref}[multiple] --chain --format=general:name1,strand1,start1,end1,name2,strand2,start2,end2 > lastz_out.txt
#conda deactivate
#svmu ${polarizing}/mumm4.delta ${dmel_ref} ${dsim_ref} l sam_lastz.txt svmu > svmu.tsv
#Does not work well
SKIP

##cuteSV
minimap2 -t ${nT} -a -x map-pb ${dmel_ref} ${polarizing}/dsim_pcbio.fastq.gz |\
samtools view -b -h -@ ${nT} -o - |\
samtools sort -@ ${nT} -o ${polarizing}/dsim_read-dmel.sort.bam
samtools index -@ ${nT} ${polarizing}/dsim_read-dmel.sort.bam
module load python/3.10.2
cuteSV --threads ${nT} --genotype --sample "dsim_cute" \
--min_support 10 \
--min_size 50 --min_mapq 20 --min_read_len 500 \
-L '-1' \
--merge_del_threshold 270 --merge_ins_threshold 270 \
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5 \
"${polarizing}/dsim_read-dmel.sort.bam" "${dmel_ref}" "dsim.cutesv.vcf" "${polarizing}"
module unload python/3.10.2

bgzip -k -@ ${nT} ${polarizing}/dsim.cutesv.vcf
bcftools index -t -f ${polarizing}/dsim.cutesv.vcf.gz
bcftools view --threads ${nT} -r 2L,2R,3L,3R,4,X,Y \
-i 'QUAL >= 20 && FILTER = "PASS"'  -O v -o - ${polarizing}/dsim.cutesv.vcf.gz |\
sed 's/DUP_TANDEM/DUP/g; s/DUP:TANDEM/DUP/g; s/DUP_INT/DUP/g; s/DUP:INT/DUP/g' |\
grep -v "NNNNNNNNNNNNNNNNNNNN" |bcftools sort --max-mem 4G -O z -o ${polarizing}/dsim.cutesv.filtered.vcf.gz

## minimap2-asm20
minimap2 -a -x asm20 -t ${nT} ${dmel_ref} ${dsim_ref} |\
samtools sort -@ ${nT} -O bam -o ${polarizing}/mm2-20.bam
samtools index -@ ${nT} ${polarizing}/mm2-20.bam

conda activate sv-calling
svim-asm haploid --sample "dsim_mm2_svimASM" --min_sv_size 50 \
${polarizing}/mm2-20_svimASM ${polarizing}/mm2-20.bam ${dmel_ref}
conda deactivate

bgzip -@ ${nT} -k -f ${polarizing}/mm2-20_svimASM/variants.vcf
bcftools index -t -f ${polarizing}/mm2-20_svimASM/variants.vcf.gz
bcftools view --threads ${nT} -r 2L,2R,3L,3R,4,X,Y \
-i 'FILTER = "PASS"' -O v -o - ${polarizing}/mm2-20_svimASM/variants.vcf.gz |\
sed 's/DUP_TANDEM/DUP/g; s/DUP:TANDEM/DUP/g; s/DUP_INT/DUP/g; s/DUP:INT/DUP/g; s/BND/TRA/g' |\
sed 's/ .   PASS/   30  PASS/g' |\
grep -v "NNNNNNNNNNNNNNNNNNNN" |bcftools sort --max-mem 4G -O z -o ${polarizing}/mm2-20.filtered.vcf.gz
bcftools index -f -t ${polarizing}/mm2-20.filtered.vcf.gz

## MUMandCo with "--maxmatch"
cd ${polarizing}

#Change the genome size to force MumandCo to use "nucmer --madmatch"
bash /pub/jenyuw/Software/MUMandCo-MUMandCov3.8/mumandco_v3.8.sh  \
-r ${dmel_ref} -q ${dsim_ref} \
-g 1000 -o dsim_mumco2 -t ${nT} -ml 50

cd ${polarizing}/dsim_mumco2_output
cat ${polarizing}/dsim_mumco2_output/dsim_mumco2.SVs_all.vcf|grep "#" |\
grep -v "##query_contig" |sed 's/ type=.*.;,/,/g' |tr -d " "|\
sed 's/qCHR/Q_CHROM/g; s/qSTART/Q_START/g; s/qEND/Q_END/g' |\
sed 's/##INFO=<ID=END,Number=1,Type=Integer,Description=Endpositioninthereferencegenome>/##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">/g' |\
sed 's/##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=DifferenceinlengthbetweenREFandALTalleles>/##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">/g' |\
sed 's/##INFO=<ID=SVTYPE,Number=1,Type=String,Description=Typeofstructuralvariant>/##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">/g' |\
sed 's/##INFO=<ID=Q_CHROM,Number=1,Type=String,Description=Chromosomeinquerygenome>/##INFO=<ID=Q_CHROM,Number=1,Type=String,Description=Chromosome in query genome>/g' |\
sed 's/##INFO=<ID=Q_START,Number=1,Type=Integer,Description=Startpositioninquerygenome>/##INFO=<ID=Q_START,Number=1,Type=Integer,Description="Start position in query genome">/g' |\
sed 's/##INFO=<ID=Q_END,Number=1,Type=Integer,Description=Endpositioninquerygenome>/##INFO=<ID=Q_END,Number=1,Type=Integer,Description="End position in query genome">/g' |\
sed -E "s/##ALT=<ID=CONTR,Description=Contraction>\n//g" |\
sed 's/##ALT=<ID=DUP:TANDEM,Description=TandemDuplication>/##ALT=<ID=DUP,Description="Duplication">/g' |\
sed 's/##ALT=<ID=DEL,Description=Deletion>/##ALT=<ID=DEL,Description="Deletion">/g' |\
sed 's/##ALT=<ID=INS,Description=Insertionofnovelsequence>/##ALT=<ID=INS,Description="Insertion of novel sequence">/g' |\
sed 's/##ALT=<ID=INV,Description=Inversion>/##ALT=<ID=INV,Description="Inversion">/g' |\
sed 's/##ALT=<ID=TRA,Description=Regioninvolvedintranslocationalternativetoabreakendposition>/##ALT=<ID=TRA,Description="Breakend">/g' \
>${polarizing}/dsim_mumco2_output/header.pre

## making VCF
grep -v "#" ${polarizing}/dsim_mumco2_output/dsim_mumco2.SVs_all.vcf |\
sed 's/CONTR/DEL/g; s/qCHR/Q_CHROM/g; s/qSTART/Q_START/g; s/qEND/Q_END/g'|\
gawk '{print $1 "\t" $2 "\t" $1"_"$2 "\t" $4 "\t" $5 "\t" 40 "\t" $7 "\t" $8 "\t" $9 "\t" $10 }' \
>${polarizing}/dsim_mumco2_output/container.pre


cat ${polarizing}/dsim_mumco2_output/{header,container}.pre |\
grep -v "NNNNNNNNNNNNNNNNNNNN" |\
bgzip -@ ${nT} -c > ${polarizing}/dsim_mumco2_output/good.vcf.gz
bcftools index -f -t  ${polarizing}/dsim_mumco2_output/good.vcf.gz

bcftools view --threads ${nT} -r 2L,2R,3L,3R,4,X,Y -i 'SVTYPE="DEL" || SVTYPE="INS"' \
-O z -o ${polarizing}/dsim_mumco2.filter.vcf.gz ${polarizing}/dsim_mumco2_output/good.vcf.gz
bcftools index -t -f ${polarizing}/dsim_mumco2.filter.vcf.gz

####################################################################
#Merge and collapse the SVs from two callers

bcftools merge -m none ${polarizing}/dsim_mumco2.filter.vcf.gz ${polarizing}/mm2-20.filtered.vcf.gz |\
bcftools sort -m 4G -O z -o ${polarizing}/dsim.2.vcf.gz
bcftools index -f -t ${polarizing}/dsim.2.vcf.gz

module load python/3.10.2

truvari collapse --intra -k maxqual --typeignore --sizemax 5000000 --sizemin 50 --refdist 600 --pctseq 0.8 --pctsize 0.8 \
-i ${polarizing}/dsim.2.vcf.gz \
-c ${polarizing}/dsim.collapsed.vcf.gz |\
bcftools sort -m 4G |bgzip -@ ${nT} > ${polarizing}/dsim.final.vcf.gz
bcftools index -f -t ${polarizing}/dsim.final.vcf.gz
module unload python/3.10.2
echo "there are " `zcat ${polarizing}/dsim.final.vcf.gz| grep -v "#" |wc -l` " SVs in the final VCF"

: <<'SKIP'
## minimap2-asm10
minimap2 -a -x asm10 -t ${nT} ${dmel_ref} ${dsim_ref} |\
samtools sort -@ ${nT} -O bam -o ${polarizing}/mm2-10.bam
samtools index -@ ${nT} ${polarizing}/mm2-10.bam
## minimap2-asm5
minimap2 -a -x asm10 -t ${nT} ${dmel_ref} ${dsim_ref} |\
samtools sort -@ ${nT} -O bam -o ${polarizing}/mm2-5.bam
samtools index -@ ${nT} ${polarizing}/mm2-5.bam
conda activate sv-calling
lra global -CONTIG ${dmel_ref}
lra align -CONTIG --refineBreakpoints ${dmel_ref} ${dsim_ref} -t 16 -p s |
samtools sort -@ ${nT} -O bam -o ${polarizing}/lra.bam
samtools index -@ ${nT} ${polarizing}/lra.bam
conda deactivate
svim-asm haploid --sample "dsim_lra_svimASM" --min_sv_size 50 \
${polarizing}/lra_svimASM ${polarizing}/lra.bam ${dmel_ref}

svim-asm haploid --sample "dsim_mumm4_svimASM" --min_sv_size 50 \
${polarizing}/mumm4_svimASM ${polarizing}/corrected.mumm4.bam ${dmel_ref}

svim-asm haploid --sample "dsim_mm2-10_svimASM" --min_sv_size 50 \
${polarizing}/mm2-10_svimASM ${polarizing}/mm2-10.bam ${dmel_ref}

svim-asm haploid --sample "dsim_mm2-5_svimASM" --min_sv_size 50 \
${polarizing}/mm2-5_svimASM ${polarizing}/mm2-5.bam ${dmel_ref}
SKIP

## We choose the SV called by minimap2-asm20,svim-asm and MumandCo


## create the collapsed (overlapping) SVs between the population (all) and D. simulans
#Avoid collapsing the SVs from samples TWICE

