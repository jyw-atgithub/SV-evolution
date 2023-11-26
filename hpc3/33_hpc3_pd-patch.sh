#!/bin/bash

#SBATCH --job-name="pd&patch"    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=24   ## number of cores the job needs
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU

## input "polished Flye asm" --> purge_dups --> ragtag patch with "polished Nextdenovo asm" 
## --> purge_dups --> ragtag patch with "polished Canu asm" --> purge_dups --> "final purged asm"

## path

purge_dups="/dfs7/jje/jenyuw/SV-project-temp/result/purge_dups"
trimmed="/dfs7/jje/jenyuw/SV-project-temp/result/trimmed"
assemble="/dfs7/jje/jenyuw/SV-project-temp/result/assemble"
polishing="/dfs7/jje/jenyuw/SV-project-temp/result/polishing"
patched="/dfs7/jje/jenyuw/SV-project-temp/result/patched"
pd_scripts="/pub/jenyuw/Software/purge_dups/scripts"
pd_bin="/pub/jenyuw/Software/purge_dups/bin"

nT=$SLURM_CPUS_PER_TASK

source ~/.bashrc

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
echo "Yes, ARRAY_TASK_ID=1"
ls ${polishing}/*flye.nextpolish* > ${polishing}/flye_list.txt
else
echo "No need to list the file again"
fi

file=`head -n $SLURM_ARRAY_TASK_ID ${polishing}/flye_list.txt |tail -n 1`
name=$(basename ${file}|gawk -F "." '{print $1}')
read_type=`echo ${name} | gawk -F "_" '{print $2}'`
echo "the file is ${file}"
echo "the read type is ${read_type}"
declare -A mapping_option=(["CLR"]='map-pb' ["hifi"]='map-hifi' ["ONT"]='map-ont')
echo "The mapping option is ${mapping_option[$read_type]}"

function purge {
    # $1 is number of threads
    # $2 is mapping option, 'map-pb', 'map-hifi' or 'map-ont'
    # $3 is primary assembly
    # $4 is the trimmed reads
    # $5 is the prefix
    #prefix=`basename $3 | gawk -F "." '{print $1 "." $2 "." $3}'`
    prefix=$5
    echo "the prefix is ${prefix}"
    mkdir ${purge_dups}/${prefix}.polished.purged
    cd ${purge_dups}/${prefix}.polished.purged

    minimap2 -t $1 -x $2 $3 $4 | pigz -p 10 -c - > ${prefix}.paf.gz
    ${pd_bin}/pbcstat ${prefix}.paf.gz
    ${pd_bin}/calcuts PB.stat > cutoffs 2>calcults.log
    ${pd_bin}/split_fa $3 > ${prefix}.split
    minimap2 -t $1 -x asm5 -DP ${prefix}.split ${prefix}.split |pigz -p 10 -c - > ${prefix}.split.self.paf.gz
    ${pd_bin}/purge_dups -2 -T cutoffs -c PB.base.cov ${prefix}.split.self.paf.gz > dups.bed 2> purge_dups.log

    ${pd_bin}/get_seqs -e dups.bed $3
}



if [[ -f ${polishing}/${name}.canu.nextpolish.fasta ]]
then
echo "${polishing}/${name}.canu.nextpolish.fasta" "Exists"
else
echo "polished canu assembly of ${name} does not exist"
fi

# 1st purge_dups
purge ${nT} ${mapping_option[$read_type]} ${file} ${trimmed}/${name}.trimmed.fastq.gz ${name}.flye


#check the existance of polished assmblies and then patch with ragtag
if [[ -f ${polishing}/${name}.nextdenovo-45.nextpolish.fasta ]]
then
echo "${polishing}/${name}.nextdenovo-45.nextpolish.fasta" "Exists"
conda activate ragtag
ragtag.py patch -w -u -o ${patched}/${name}_1 --aligner 'nucmer' ${purge_dups}/${name}.flye.polished.purged/purged.fa ${polishing}/${name}.nextdenovo-45.nextpolish.fasta
conda deactivate
else
echo "polished nextdenovo assembly of ${name} does not exist"
exit
fi

# 2nd purge_dups
purge ${nT} ${mapping_option[$read_type]} ${patched}/${name}_1/ragtag.patch.fasta ${trimmed}/${name}.trimmed.fastq.gz ${name}.nd

echo "This is the end!!"