#!/bin/bash
#nohup ./Long-SNPcalling.sh >Long-SNPcalling.out & [1] 469525


## path
ref_genome="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"
aligned_bam="/home/jenyuw/SV-project/result/aligned_bam"
SVs="/home/jenyuw/SV-project/result/SVs"
merged_SVs="/home/jenyuw/SV-project/result/merged_SVs"
SNP="/home/jenyuw/SV-project/result/SNP"
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
done

conda activate long-snp

## Use NanoCaller for Pacbio CLR data
for i in $(ls ${aligned_bam}/{A*,B*,ORE}.trimmed-ref.sort.bam)
do
name=$(basename $i|sed 's/.trimmed-ref.sort.bam//g' )
NanoCaller --bam ${i} --ref ${ref_genome} --cpu ${nT} \
--mode snps --regions 2L 2R 3L 3R 4 X Y \
--preset clr --mincov 4 --maxcov 160 \
--sample ${name} \
--output ${SNP} --prefix ${name}-nanocaller
done