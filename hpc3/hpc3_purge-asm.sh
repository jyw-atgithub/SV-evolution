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
ls ${trimmed}/*.trimmed.fastq.gz| gawk -F "/" '{print $7}' | sed "s/.trimmed.fastq.gz//g ; s/.rn//g" >${purge_dups}/namelist.txt
else
echo "No need to list the file again"
fi

n=`cat ${purge_dups}/namelist.txt |wc -l `
echo "there are ${n}" "names"
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
  mkdir ${purge_dups}/${prefix}.purged
  cd ${purge_dups}/${prefix}.purged

  minimap2 -t $1 -x $2 $3 $4 | pigz -p 10 -c - > ${prefix}.paf.gz
  ${pd_bin}/pbcstat ${prefix}.paf.gz
  ${pd_bin}/calcuts PB.stat > cutoffs 2>calcults.log
  ${pd_bin}/split_fa $3 > ${prefix}.split
  minimap2 -t $1 -x asm5 -DP ${prefix}.split ${prefix}.split |pigz -p 10 -c - > ${prefix}.split.self.paf.gz
  ${pd_bin}/purge_dups -2 -T cutoffs -c PB.base.cov ${prefix}.split.self.paf.gz > dups.bed 2> purge_dups.log

  ${pd_bin}/get_seqs -e dups.bed $3
}



#### UGLY solution for poor nomenclature

if (( $SLURM_ARRAY_TASK_ID < 18 ))
then

##for nv and DSPR samples
assemblers="Flye canu nextdenovo"
for i in `echo ${assemblers}`
do
  read_type=`echo $name |gawk -F "_" '{print $2}' `
  echo "The read type is ${read_type}"
  name=`echo $name |gawk -F "_" '{print $1}' `
  echo "The new name is ${name}"

  mapping_option=(["CLR"]="map-pb" ["HIFI"]="asm20" ["ONT"]="map-ont")
  if [[ $i == "Flye" ]]
  then
  purge ${nT} ${mapping_option[$read_type]} \
  ${assemble}/${name}_${i}/assembly.fasta ${trimmed}/${name}_${read_type}.trimmed.fastq.gz ${name}_${read_type}_$i
  elif [[ $i == "canu" ]]
  then
  purge ${nT} ${mapping_option[$read_type]} \
  ${assemble}/${name}_${i}/${name}.contigs.fasta ${trimmed}/${name}_${read_type}.trimmed.fastq.gz ${name}_${read_type}_$i
  elif [[ $i == "nextdenovo" ]]
  then
  purge ${nT} ${mapping_option[$read_type]} \
  ${assemble}/${name}_${i}/03.ctg_graph/nd.asm.fasta ${trimmed}/${name}_${read_type}.trimmed.fastq.gz ${name}_${read_type}_$i
  else
  echo "NO such assembler ${i} was used"
  fi
echo "FOR nv and DSPR samples FINISHED"
name=`head -n $SLURM_ARRAY_TASK_ID ${purge_dups}/namelist.txt |tail -n 1`
done 

elif (( $SLURM_ARRAY_TASK_ID > 17 ))
then
##for SRR samples
assemblers="Flye canu nextdenovo"
for i in `echo ${assemblers}`
do
  read_type=`echo $name |gawk -F "_" '{print $2}' `
  echo "The read type is ${read_type}"
  mapping_option=(["CLR"]="map-pb" ["HIFI"]="asm20" ["ONT"]="map-ont")
  if [[ $i == "Flye" ]]
  then
  purge ${nT} ${mapping_option[$read_type]} \
  ${assemble}/${name}_${i}/assembly.fasta ${trimmed}/${name}.trimmed.fastq.gz ${name}_$i
  elif [[ $i == "canu" ]]
  then
  purge ${nT} ${mapping_option[$read_type]} \
  ${assemble}/${name}_${i}/${name}.contigs.fasta ${trimmed}/${name}.trimmed.fastq.gz ${name}_$i
  elif [[ $i == "nextdenovo" ]]
  then
  purge ${nT} ${mapping_option[$read_type]} \
  ${assemble}/${name}_${i}/03.ctg_graph/nd.asm.fasta ${trimmed}/${name}.trimmed.fastq.gz ${name}_$i
  else
  echo "NO such assembler ${i} was used"
  fi
echo "FOR SRR samples FINISHED"
done 

fi 

