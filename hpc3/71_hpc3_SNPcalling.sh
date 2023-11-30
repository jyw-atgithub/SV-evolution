#!/bin/bash

#SBATCH --job-name=snp    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=2-62%1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=60   ## number of cores the job needs
#SBATCH --mem-per-cpu=6G     # requesting 6 GB memory per CPU, the max

#Because nanocaller produce intermediate files of the same name, we can only run one sample at a time

## path
ref_genome="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
SNPs="/dfs7/jje/jenyuw/SV-project-temp/result/SNPs"
clair3="/pub/jenyuw/Software/clair3_latest.sif"
nanocaller="/pub/jenyuw/Software/nanocaller_3.4.1.sif"

## prep
#source ~/.bashrc
nT=$SLURM_CPUS_PER_TASK

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
ls ${aligned_bam}/*.trimmed-ref.sort.bam >${aligned_bam}/trimmed-ref.list
fi

file=$(head -n $SLURM_ARRAY_TASK_ID ${aligned_bam}/trimmed-ref.list|tail -n 1)
name=$(basename ${file}|sed 's/.trimmed-ref.sort.bam//g' )
read_type=$(echo ${name}|cut -d"_" -f 2)
declare -A preset=(["CLR"]='clr' ["hifi"]='ccs' ["ONT"]='ont')
echo "The preset is ${preset[$read_type]}"

# Clair3 and NanoCaller are installed in singularity
module load singularity/3.11.3
# For simplicity, Use NanoCaller for both ONT and  Pacbio CLR data
## absolute output path prefix
singularity exec \
--bind ${aligned_bam},${SNPs} \
/pub/jenyuw/Software/nanocaller_3.4.1.sif \
NanoCaller --bam ${file} --ref ${ref_genome} --cpu ${nT} \
--mode snps --regions 2L 2R 3L 3R 4 X Y \
--preset ${preset[$read_type]} --mincov 4 --maxcov 120 \
--sample ${name} \
--output ${SNPs} --prefix ${name}_nanocaller
## Here, the output file is called [prefix].snps.vcf.gz
module unload singularity/3.11.3
wait


bcftools filter --threads ${nT} -i 'DP > 10 && QUAL >= 30' ${SNPs}/${name}_nanocaller.snps.vcf.gz |\
bcftools view --threads ${nT} -m 2 -M 2 -v snps -O z > ${SNPs}/${name}.filtered.snps.vcf.gz
sleep 10
wait
bcftools index -t -f ${SNPs}/${name}.filtered.snps.vcf.gz

echo "This is the end"