#! /bin/bash

#SBATCH --job-name="correct2"
#SBATCH -A jje_lab
#SBATCH -p hugemem
#SBATCH --cpus-per-task=40
#SBATCH --output="dorado_correct_2.out"
source ~/.bashrc

filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
corrected="/dfs7/jje/jenyuw/Assembling_ISO1/results/corrected"
nT=$SLURM_CPUS_PER_TASK
#dorado correct --verbose --to-paf ${filtered}/all_60k10k.fastq.gz > ${corrected}/all_60k10k.paf
dorado correct --verbose --to-paf ${filtered}/all_50k8k.renamed.fastq.gz > ${corrected}/all_50k8k.paf

###Then, use our desktop to perform the GPU-intense part.
#! /bin/bash
#in the same folder
export LD_LIBRARY_PATH="/home/jenyuwang/Software/dorado-1.0.2-linux-x64/lib"
#time dorado correct all_60k10k.fastq.gz --verbose -x cuda:all --from-paf all_60k10k.paf |bgzip -@ 4 -c >all_60k10k_correct.fasta.gz
dorado correct all_50k8k.renamed.fastq.gz --verbose -x cuda:all --from-paf all_50k8k.paf |bgzip -@ 4 -c >all_50k8k_correct.fasta.gz

########## ON THOTH BEGIN #######
#! /bin/bash
dorado correct  -x cuda:all --verbose all_60k10k.fastq.gz >all_60k10k_correct.fasta
# -x cuda:all ==> memory is not sufficient!! The GPU spects is improper.
########## ON THOTH END #######


#! /bin/bash
#SBATCH --job-name=verkko5
#SBATCH -A jje_lab
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --cpus-per-task=56
#SBATCH --mem-per-cpu=6G
#SBATCH --output="verkko_5.out"
source ~/.bashrc
raw="/dfs7/jje/jenyuw/Assembling_ISO1/raw"
filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
corrected="/dfs7/jje/jenyuw/Assembling_ISO1/results/corrected"
assembly="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly"
micromamba activate verkko
verkko -d ${assembly}/verkko_5 --hifi ${corrected}/all_60k10k_correct.fasta.gz --nano ${filtered}/all_50k8k.renamed.fastq.gz \
--local-memory 320 --unitig-abundance 4 --haploid
micromamba deactivate


#! /bin/bash
#SBATCH --job-name=hifiasm_10
#SBATCH -A jje_lab
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --cpus-per-task=60
#SBATCH --mem-per-cpu=6G
#SBATCH --output="hifiasm_10.out"
source ~/.bashrc
filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
corrected="/dfs7/jje/jenyuw/Assembling_ISO1/results/corrected"
assembly="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly"
nT=$SLURM_CPUS_PER_TASK
mkdir ${assembly}/hifiasm_10
cd ${assembly}/hifiasm_10
hifiasm -l 0 --primary -t ${nT} --ul-rate 0.02 \
-o ONT_10 --ul ${filtered}/all_50k8k.renamed.fastq.gz ${corrected}/all_60k10k_correct.fasta.gz


#! /bin/bash
#SBATCH --job-name=hifiasm11
#SBATCH -A jje_lab
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --cpus-per-task=60
#SBATCH --mem-per-cpu=6G
#SBATCH --output="hifiasm_11.out"
source ~/.bashrc
filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
corrected="/dfs7/jje/jenyuw/Assembling_ISO1/results/corrected"
assembly="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly"
nT=$SLURM_CPUS_PER_TASK
mkdir ${assembly}/hifiasm_11
cd ${assembly}/hifiasm_11
hifiasm -l 0 --primary -t ${nT} --ul-rate 0.02 -o ONT_11 \
--ont ${filtered}/all_60k10k.fastq.gz \
--ul ${filtered}/all_50k8k.renamed.fastq.gz ${corrected}/all_60k10k_correct.fasta.gz


#! /bin/bash
#SBATCH --job-name=verkko6
#SBATCH -A jje_lab
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem-per-cpu=6G
#SBATCH --time=7-00
#SBATCH --output="verkko_6.out"
#SBATCH -e verkko_6.err

source ~/.bashrc
raw="/dfs7/jje/jenyuw/Assembling_ISO1/raw"
filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
corrected="/dfs7/jje/jenyuw/Assembling_ISO1/results/corrected"
assembly="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly"
micromamba activate verkko2
verkko -d ${assembly}/verkko_6 --hifi ${filtered}/all_50k8k_corrected.filtered_40k.fasta  --nano ${filtered}/all_60k10k.fastq.gz \
--unitig-abundance 4 --haploid  \
--mbg-run 56 330 96 \
--local --local-memory 330 --local-cpus 56
#--snakeopts "--unlock --default-resources mem_gb=330 n_cpus=56"
#--lay_run <ncpus> <mem-in-gb> <time-in-h> DOES NOT WORK!
micromamba deactivate

#! /bin/bash
#SBATCH --job-name=hifiasm12
#SBATCH -A jje_lab 
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --cpus-per-task=60
#SBATCH --mem-per-cpu=6G
#SBATCH --output="hifiasm_12.out"
source ~/.bashrc
filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
corrected="/dfs7/jje/jenyuw/Assembling_ISO1/results/corrected"
assembly="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly"
nT=$SLURM_CPUS_PER_TASK
mkdir ${assembly}/hifiasm_12
cd ${assembly}/hifiasm_12
hifiasm -l 0 --primary -t ${nT} --ul-rate 0.02 -o ONT_12 \
--ul ${filtered}/all_60k10k.fastq.gz ${filtered}/all_50k8k_corrected.filtered_40k.fasta
#--> Terrible Result!!


#! /bin/bash
#SBATCH --job-name=hifiasm13
#SBATCH -A jje_lab
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --cpus-per-task=60
#SBATCH --mem-per-cpu=6G
#SBATCH --output="hifiasm_13.out"
source ~/.bashrc
filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
corrected="/dfs7/jje/jenyuw/Assembling_ISO1/results/corrected"
assembly="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly"
nT=$SLURM_CPUS_PER_TASK
mkdir ${assembly}/hifiasm_13
cd ${assembly}/hifiasm_13
hifiasm -l 0 --primary -t ${nT} --ul-rate 0.02 -o ONT_13 --ont \
--ul ${filtered}/all_60k10k.fastq.gz ${filtered}/regular_100k.fastq.gz


#! /bin/bash
#SBATCH --job-name=hifiasm14
#SBATCH -A jje_lab
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --cpus-per-task=60
#SBATCH --mem-per-cpu=6G
#SBATCH --output="hifiasm_14.out"
source ~/.bashrc
filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
corrected="/dfs7/jje/jenyuw/Assembling_ISO1/results/corrected"
assembly="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly"
nT=$SLURM_CPUS_PER_TASK
mkdir ${assembly}/hifiasm_14
cd ${assembly}/hifiasm_14
hifiasm -l 0 --primary -t ${nT} --ul-rate 0.02 -o ONT_14 \
--ul ${filtered}/all_60k10k.fastq.gz ${corrected}/all_60k10k_correct.fasta.gz

#! /bin/bash
#SBATCH --job-name=hifiasm15
#SBATCH -A jje_lab
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --cpus-per-task=60
#SBATCH --mem-per-cpu=6G
#SBATCH --output="hifiasm_15.out"
source ~/.bashrc
filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
corrected="/dfs7/jje/jenyuw/Assembling_ISO1/results/corrected"
assembly="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly"
nT=$SLURM_CPUS_PER_TASK
mkdir ${assembly}/hifiasm_15
cd ${assembly}/hifiasm_15
hifiasm -l 0 --primary -t ${nT}  -o ONT_15 \
--ul ${filtered}/all_50k8k.renamed.fastq.gz ${corrected}/all_50k8k_correct.fasta.gz

#################################################################################################################


#! /bin/bash

#SBATCH --job-name="correct3"
#SBATCH -A jje_lab
#SBATCH -p hugemem
#SBATCH --cpus-per-task=40
#SBATCH --output="dorado_correct_3.out"
source ~/.bashrc
raw="/dfs7/jje/jenyuw/Assembling_ISO1/raw"
filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
corrected="/dfs7/jje/jenyuw/Assembling_ISO1/results/corrected"
nT=$SLURM_CPUS_PER_TASK
#dorado correct --verbose --to-paf ${filtered}/all_60k10k.fastq.gz > ${corrected}/all_60k10k.paf
#dorado correct --verbose --to-paf ${filtered}/all_50k8k.renamed.fastq.gz > ${corrected}/all_50k8k.paf
dorado correct --verbose --to-paf ${raw}/sup_simplex_Q10.fastq.gz > ${corrected}/sup_simplex_Q10.paf


#! /bin/bash

#SBATCH --job-name=verkko2
#SBATCH -A jje_lab
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --cpus-per-task=50
#SBATCH --mem-per-cpu=5G
#SBATCH --output="verkko_2.out"

source ~/.bashrc

raw="/dfs7/jje/jenyuw/Assembling_ISO1/raw"
filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
assembly="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly"

micromamba activate verkko
verkko -d ${assembly}/verkko_2 --hifi ${filtered}/hifi_corrected-adaptive-3k.fastq --nano ${filtered}/9k_all.fastq.gz \
--local-memory 320 --unitig-abundance 4 --haploid
micromamba deactivate


#! /bin/bash

#SBATCH --job-name=verkko3
#SBATCH -A jje_lab
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --cpus-per-task=50
#SBATCH --mem-per-cpu=5G
#SBATCH --output="verkko_3.out"

source ~/.bashrc

filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
assembly="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly"

micromamba activate verkko
verkko -d ${assembly}/verkko_3 --hifi ${filtered}/hifi_corrected-30k5k-Q15.fastq.gz --nano ${filtered}/9k_all.fastq.gz \
--local-memory 320 --unitig-abundance 4 --haploid
micromamba deactivate

#! /bin/bash

#SBATCH --job-name=ont_only
#SBATCH -A jje_lab
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=4G
#SBATCH --output="hifiasm_ONT_correct.out"

source ~/.bashrc

filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
assembly="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly"
corrected="/dfs7/jje/jenyuw/Assembling_ISO1/results/corrected"

nT=$SLURM_CPUS_PER_TASK

mkdir ${assembly}/hifiasm_ONT
cd ${assembly}/hifiasm_ONT

hifiasm -l 0 --primary -N 150 -D 8 -u 1 -t ${nT} \
-o ONT_correct --ul ${filtered}/9k_all.fastq.gz ${corrected}/corrected_30k5k_Q15.fasta


#! /bin/bash

#SBATCH --job-name=flye
#SBATCH -A jje_lab
#SBATCH -p highmem
#SBATCH --array=1
#SBATCH --cpus-per-task=28
#SBATCH --mem-per-cpu=10G
source ~/.bashrc

filtered="/dfs7/jje/jenyuw/Assembling_ISO1/results/filtered"
assembly="/dfs7/jje/jenyuw/Assembling_ISO1/results/assembly"
nT=$SLURM_CPUS_PER_TASK


module load python/3.10.2
#flye 2.9.5
flye --nano-hq ${filtered}/all_30k5k_Q15.fastq.gz --out-dir ${assembly}/flye_2 \
--threads ${nT} --genome-size 180m --iterations 3 --read-error 0.015 \
--meta --keep-haplotypes --scaffold

module unload python/3.10.2