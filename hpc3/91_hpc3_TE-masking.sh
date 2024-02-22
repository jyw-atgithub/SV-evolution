#!/bin/bash

#SBATCH --job-name=TE    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=30   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU
#SBATCH --tmp=100G                ## requesting 100 GB local scratch
#SBATCH --constraint=fastscratch  ## requesting nodes with fast scratch in /tmp

dmel_ref="/dfs7/jje/jenyuw/SV-project-temp/reference/dmel-all-chromosome-r6.49.fasta"
dsim_ref="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing/GCF_016746395.2_Dsim_3.1.fasta"
TE="/dfs7/jje/jenyuw/SV-project-temp/result/TE_repeat"
nT=$SLURM_CPUS_PER_TASK
source ~/.bashrc

cp ${dmel_ref} ${TE}/dmel-r6.49.fasta
cd ${TE}

####### Avoid use the DROSOPHILA TRANSPOSON CANONICAL SEQUENCES #######
####### It does not involve Helitron #######
#module load singularity/3.11.3
#singularity exec -B ${TE} /pub/jenyuw/Software/dfam-tetools-latest.sif RepeatMasker -gff -s -xsmall \
#-lib "${TE}/D_mel_transposon_sequence_set.fa"  -dir ${TE}/dmel-r649 \
#${TE}/dmel-r6.49.fasta

module load singularity/3.11.3
singularity exec -B ${TE} \
-B /dfs7/jje/jenyuw/SV-project-temp/raw/Libraries:/opt/RepeatMasker/Libraries \
/pub/jenyuw/Software/dfam-tetools-latest.sif \
RepeatMasker -gff -s -xsmall \
-species "drosophila melanogaster" -dir ${TE}/dmel-r649-dfam \
${TE}/dmel-r6.49.fasta

##The process of building a custom library
##REF: https://github.com/Dfam-consortium/TETools
: <<'SKIP'
singularity run -B /dfs7/jje/jenyuw/SV-project-temp/raw /pub/jenyuw/Software/dfam-tetools-latest.sif
cp -r /opt/RepeatMasker/Libraries/ ./
exit
chown -R $USER ./Libraries/

cd /dfs7/jje/jenyuw/SV-project-temp/raw
wget https://www.dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.4.h5.gz
unpigz -k -p 6 dfam38_full.4.h5.gz
mv dfam38_full.4.h5 ./Libraries/famdb/

singularity run -B /dfs7/jje/jenyuw/SV-project-temp/raw/Libraries:/opt/RepeatMasker/Libraries \
/pub/jenyuw/Software/dfam-tetools-latest.sif
cd /opt/RepeatMasker/
rm ./Libraries/famdb/rmlib.config
./tetoolsDfamUpdate.pl
exit
SKIP