#!/bin/bash
#nohup ./Long-SNPcalling.sh >Long-SNPcalling.out & [1] 469525


## path
ref_genome="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"
aligned_bam="/home/jenyuw/SV-project/result/aligned_bam"
SVs="/home/jenyuw/SV-project/result/SVs"
merged_SVs="/home/jenyuw/SV-project/result/merged_SVs"
SNP="/home/jenyuw/SV-project/result/SNP"
merged_SNP="/home/jenyuw/SV-project/result/merged_SNP"
## prep
source ~/.bashrc
nT=30


## absolute output path prefix
## Use Clair3 for Nanoport data

for i in $(ls ${aligned_bam}/{nv*,SRR*}.trimmed-ref.sort.bam)
do
name=$(basename $i|sed 's/.trimmed-ref.sort.bam//g' )
echo ${name}
singularity exec /home/jenyuw/Software/clair3_latest.sif \
/opt/bin/run_clair3.sh \
--bam_fn=${aligned_bam}/${name}.trimmed-ref.sort.bam  --ref_fn=${ref_genome} \
--threads=${nT} --platform="ont" \
--model_path="/opt/models/r941_prom_sup_g5014" --ctg_name="2L,2R,3L,3R,4,X,Y" \
--remove_intermediate_dir \
--sample_name=${name} \
--output=${SNP}/${name}-clair3
cp ${SNP}/${name}-clair3/merge_output.vcf.gz ${SNP}/${name}-clair3.snps.vcf.gz
done


#conda activate long-snp
## Failed, status D
#for i in $(ls ${aligned_bam}/{A*,B*,ORE}.trimmed-ref.sort.bam)
#do
#name=$(basename $i|sed 's/.trimmed-ref.sort.bam//g' )
#NanoCaller --bam ${i} --ref ${ref_genome} --cpu ${nT} \
#--mode snps --regions 2L 2R 3L 3R 4 X Y \
#--preset clr --mincov 4 --maxcov 160 \
#--sample ${name} \
#--output ${SNP} --prefix ${name}-nanocaller
#done

## Use NanoCaller for Pacbio CLR data
for i in $(ls ${aligned_bam}/{A*,B*,ORE}.trimmed-ref.sort.bam)
do
name=$(basename $i|sed 's/.trimmed-ref.sort.bam//g' )
singularity exec -e /home/jenyuw/Software/nanocaller_latest.sif \
NanoCaller --bam ${i} --ref ${ref_genome} --cpu ${nT} \
--mode snps --regions 2L 2R 3L 3R 4 X Y \
--preset clr --mincov 4 --maxcov 160 \
--sample ${name} \
--output ${SNP} --prefix ${name}-nanocaller
done
## Here, the output file is called [prefix].snps.vcf.gz

#singularity exec -e --pwd /app nanocaller_latest.sif NanoCaller --help

for i in $(ls ${SNP}/*.snps.vcf.gz)
do
name=$(basename ${i}|sed 's/.snps.vcf.gz//g')
echo $name
bcftools filter --threads 8 -i 'DP > 10 && QUAL >= 30' ${i} | bcftools view -m 2 -M 2 -v snps -O z > ${name}.filtered.snps.vcf.gz
tabix -p vcf ${name}.filtered.snps.vcf.gz
done

#bcftools merge --threads 8 --missing-to-ref ${SNP}/*.filtered.snps.vcf.gz -O z > ${merged_SNP}/all.snps.vcf.gz
# Avoid using --missing-to-ref, because this does look like a good assumption
bcftools merge --threads 8 ${SNP}/*.filtered.snps.vcf.gz -O z > ${merged_SNP}/all.snps.vcf.gz
tabix -p vcf ${merged_SNP}/all.snps.vcf.gz

conda activate everything
snpEff -v BDGP6.32.105 ${merged_SNP}/all.snps.vcf.gz | bgzip -@ 8 -c  >all.snps.annotated.vcf.gz
