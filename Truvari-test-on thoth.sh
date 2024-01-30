#!/bin/bash
ref_genome="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"
source /home/jenyuw/truvari_env/bin/activate

truvari collapse --sizemax 200000000 -k first \
-i merge.asm-2.sort.vcf.gz \
-c truvari.asm-2.noseq.collapse.vcf.gz -f ${ref_genome} \
--pctseq 0 |\
bcftools sort --max-mem 4G |\
bgzip -@ 4 > truvari.asm-2.noseq.vcf.gz

truvari collapse --sizemax 200000000 -k first \
-i merge.asm-2.sort.vcf.gz \
-c truvari.asm-2.noseq-2.collapse.vcf.gz -f ${ref_genome} \
--pctseq 0 --refdist 100 --minhaplen 30 --pctsize 0.98 |\
bcftools sort --max-mem 4G |\
bgzip -@ 4 > truvari.asm-2.noseq-2.vcf.gz


chr="2L"
truvari collapse --sizemax 200000000 -k first \
-i merge.asm-2.sort.vcf.gz \
-c truvari.asm-2.${chr}.collapse.vcf.gz -f ${ref_genome} --bed ${chr}.bed|\
bcftools sort --max-mem 4G |\
bgzip -@ 4 > truvari.asm-2.${chr}.vcf.gz

chr="2R"
truvari collapse --sizemax 200000000 -k first \
-i merge.asm-2.sort.vcf.gz \
-c truvari.asm-2.${chr}.collapse.vcf.gz -f ${ref_genome} --bed ${chr}.bed|\
bcftools sort --max-mem 4G |\
bgzip -@ 4 > truvari.asm-2.${chr}.vcf.gz

chr="3L"
truvari collapse --sizemax 200000000 -k first \
-i merge.asm-2.sort.vcf.gz \
-c truvari.asm-2.${chr}.collapse.vcf.gz -f ${ref_genome} --bed ${chr}.bed|\
bcftools sort --max-mem 4G |\
bgzip -@ 4 > truvari.asm-2.${chr}.vcf.gz

chr="3R"
truvari collapse --sizemax 200000000 -k first \
-i merge.asm-2.sort.vcf.gz \
-c truvari.asm-2.${chr}.collapse.vcf.gz -f ${ref_genome} --bed ${chr}.bed|\
bcftools sort --max-mem 4G |\
bgzip -@ 4 > truvari.asm-2.${chr}.vcf.gz

chr="4"
truvari collapse --sizemax 200000000 -k first \
-i merge.asm-2.sort.vcf.gz \
-c truvari.asm-2.${chr}.collapse.vcf.gz -f ${ref_genome} --bed ${chr}.bed|\
bcftools sort --max-mem 4G |\
bgzip -@ 4 > truvari.asm-2.${chr}.vcf.gz

chr="X"
truvari collapse --sizemax 200000000 -k first \
-i merge.asm-2.sort.vcf.gz \
-c truvari.asm-2.${chr}.collapse.vcf.gz -f ${ref_genome} --bed ${chr}.bed|\
bcftools sort --max-mem 4G |\
bgzip -@ 4 > truvari.asm-2.${chr}.vcf.gz

chr="Y"
truvari collapse --sizemax 200000000 -k first \
-i merge.asm-2.sort.vcf.gz \
-c truvari.asm-2.${chr}.collapse.vcf.gz -f ${ref_genome} --bed ${chr}.bed|\
bcftools sort --max-mem 4G |\
bgzip -@ 4 > truvari.asm-2.${chr}.vcf.gz