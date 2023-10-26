#!/bin/bash

#SBATCH --job-name=pd_ori    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name: standard or highmem
#SBATCH --array=1     ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=30    ## number of cores the job needs
#SBATCH --mem-per-cpu=6G

#"#SBATCH --mem=200G", "#SBATCH --nodes=1", "#SBATCH --ntasks=1"

## path
trimmed="/pub/jenyuw/SV-project-temp/result/trimmed"
assemble="/pub/jenyuw/SV-project-temp/result/assemble"
polishing="/pub/jenyuw/SV-project-temp/result/polishing"
purge_dups="/pub/jenyuw/SV-project-temp/result/purge_dups"
pd_scripts="/pub/jenyuw/Software/purge_dups/scripts"
pd_bin="/pub/jenyuw/Software/purge_dups/bin"

nT=$SLURM_CPUS_PER_TASK

source ~/.bashrc

##make name list
##Because of the lazy way of generating file list. Always submit --array=1 at the begining and then submit the rest

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
echo "Yes, ARRAY_TASK_ID=1"
ls ${polishing}/*.nextpolish.fasta| gawk -F "/" '{print $7}' | sed "s/.nextpolish.fasta//g" >${purge_dups}/namelist.txt
else
echo "No need to list the file again"
fi

n=`cat ${purge_dups}/namelist.txt|wc -l `
echo "there are ${n} names"
name=`head -n $SLURM_ARRAY_TASK_ID ${purge_dups}/namelist.txt |tail -n 1`
echo "name is ${name}"


function purge {
  # $1 is number of threads
  # $2 is mapping option
  # $3 is primary assembly
  # $4 is the trimmed reads
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


STR="A1 A2 A3 A4 A5 A6 A7 AB8 B1 B2 B3 B4 B6 ORE"
STR2="nv107 nv109 nv115"
strain=`echo ${name} |gawk -F "." '{print $1}'`
echo "strain is ${strain}"
if grep -q ${strain} <<< "$STR"
##for DSPR samples
then
  mapping_option=(["CLR"]="map-pb" ["HIFI"]="asm20" ["ONT"]="map-ont")
  read_type=CLR
  purge ${nT} ${mapping_option[$read_type]} \
  ${polishing}/${name}.nextpolish.fasta ${trimmed}/${strain}_${read_type}.trimmed.fastq.gz ${name}
elif grep -q ${strain} <<< "$STR2"
##for nv samples
then
  mapping_option=(["CLR"]="map-pb" ["HIFI"]="asm20" ["ONT"]="map-ont")
  read_type=ONT
  purge ${nT} ${mapping_option[$read_type]} \
  ${polishing}/${name}.nextpolish.fasta ${trimmed}/${strain}_${read_type}.trimmed.fastq.gz ${name}
elif [[ ${name} =~ "SRR232"* ]]
##for SRR samples
then
  mapping_option=(["CLR"]="map-pb" ["HIFI"]="asm20" ["ONT"]="map-ont")
  read_type=ONT
  purge ${nT} ${mapping_option[$read_type]} \
  ${polishing}/${name}.nextpolish.fasta ${trimmed}/${name}.trimmed.fastq.gz ${name}
else
echo "failed"
fi
