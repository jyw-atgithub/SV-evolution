#!/bin/bash

#SBATCH --job-name=pd_ori    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name: standard or highmem
#SBATCH --array=1     ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=24    ## number of cores the job needs
#SBATCH --mem-per-cpu=6G

#"#SBATCH --mem=200G", "#SBATCH --nodes=1", "#SBATCH --ntasks=1"

## path
trimmed="/pub/jenyuw/SV-project-temp/result/trimmed"
assemble="/pub/jenyuw/SV-project-temp/result/assemble"
polishing="/pub/jenyuw/SV-project-temp/result/polishing"
purge_dups="/pub/jenyuw/SV-project-temp/result/purge_dups"
pd_scripts="/pub/jenyuw/Software/purge_dups/scripts"

nT=$SLURM_CPUS_PER_TASK

source ~/.bashrc

##make name list
##Because of the lazy way of generating file list. Always submit --array=1 at the begining and then submit the rest

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
then
echo "yes, array=1"
ls ${trimmed}/*.trimmed.fastq.gz| gawk -F "/" '{print $7}' | sed "s/.trimmed.fastq.gz//g ; s/.rn//g" >${purge_dups}/namelist.txt
else
echo "no need to list the file again"
fi

name=`head -n $SLURM_ARRAY_TASK_ID ${purge_dups}/namelist.txt |tail -n 1`
echo "name is ${name}"

## Purge_dups!!
#step 1
ls ${trimmed}/${name}*.trimmed.fastq.gz >${purge_dups}/${name}-read.fofn

echo "the content of ${name}-read.fofn is "
cat ${purge_dups}/${name}-read.fofn

#step 2
assemblers="Flye canu nextdenovo"
for i in `echo ${assemblers}`
do
  if [[ $i == "Flye" ]]
  then
    python3 ${pd_scripts}/pd_config.py \
    -l ${purge_dups}/${name}_${i} \
    -n ${purge_dups}/${name}_${i}-config.json \
    ${assemble}/${name}_${i}/assembly.fasta \
    ${purge_dups}/${name}-read.fofn
  elif [[ $i == "canu" ]]
  then
    python3 ${pd_scripts}/pd_config.py \
    -l ${purge_dups}/${name}_${i} \
    -n ${purge_dups}/${name}_${i}-config.json \
    ${assemble}/${name}_${i}/${name}.contigs.fasta \
    ${purge_dups}/${name}-read.fofn
  elif [[ $i == "nextdenovo" ]]
  then
    python3 ${pd_scripts}/pd_config.py \
    -l ${purge_dups}/${name}_${i} \
    -n ${purge_dups}/${name}_${i}-config.json \
    ${assemble}/${name}_${i}/03.ctg_graph/nd.asm.fasta \
    ${purge_dups}/${name}-read.fofn
  else
  echo "NO such assembler ${i} was used"
  fi
done

#step 3
echo -e "
{
  \"cc\": {
    \"fofn\": \"auto-read.fofn\",
    \"isdip\": 1,
    \"core\": ${nT},
    \"mem\": 60000,
    \"queue\": \"normal\",
    \"mnmp_opt\": \"\",
    \"bwa_opt\": \"\",
    \"ispb\": 1,
    \"skip\": 0
  },
  \"sa\": {
    \"core\": ${nT},
    \"mem\": 60000,
    \"queue\": \"normal\"
  },
  \"busco\": {
    \"core\": ${nT},
    \"mem\": 60000,
    \"queue\": \"long\",
    \"skip\": 0,
    \"lineage\": \"diptera_odb10\",
    \"prefix\": \"assembly_purged\",
    \"tmpdir\": \"busco_tmp\"
  },
  \"pd\": {
    \"mem\": 60000,
    \"queue\": \"normal\"
  },
  \"gs\": {
    \"mem\": 60000,
    \"oe\": 1
  },
  \"kcp\": {
    \"core\": ${nT},
    \"mem\": 60000,
    \"fofn\": \"\",
    \"prefix\": \"assembly_purged_kcm\",
    \"tmpdir\": \"kcp_tmp\",
    \"skip\": 1
  },
  \"ref\": \"path-to-assembly.fasta\",
  \"out_dir\": \"working-dir-set-in-pd_config.py-l\"
}" >${purge_dups}/config-template."$SLURM_ARRAY_TASK_ID".txt

assemblers="Flye canu nextdenovo"
for i in `echo ${assemblers}`
do
  if [[ $i == "Flye" ]]
  then
    cat ${purge_dups}/config-template."$SLURM_ARRAY_TASK_ID".txt |\
    sed "s@auto-read.fofn@${purge_dups}/${name}_${i}/pb.fofn@" |\
    sed "s@path-to-assembly.fasta@"${purge_dups}/${name}_${i}/assembly.fasta"@ " |\
    sed "s@working-dir-set-in-pd_config.py-l@"${purge_dups}/${name}_${i}"@ " \
    >${purge_dups}/${name}_${i}-config.json
  elif [[ $i == "canu" ]]
  then
    cat ${purge_dups}/config-template."$SLURM_ARRAY_TASK_ID".txt |\
    sed "s@auto-read.fofn@${purge_dups}/${name}_${i}/pb.fofn@" |\
    sed "s@path-to-assembly.fasta@"${purge_dups}/${name}_${i}/${name}.contigs.fasta"@ " |\
    sed "s@working-dir-set-in-pd_config.py-l@"${purge_dups}/${name}_${i}"@ " \
    >${purge_dups}/${name}_${i}-config.json
  elif [[ $i == "nextdenovo" ]]
  then
    cat ${purge_dups}/config-template."$SLURM_ARRAY_TASK_ID".txt |\
    sed "s@auto-read.fofn@${purge_dups}/${name}_${i}/pb.fofn@" |\
    sed "s@path-to-assembly.fasta@"${purge_dups}/${name}_${i}/nd.asm.fasta"@ " |\
    sed "s@working-dir-set-in-pd_config.py-l@"${purge_dups}/${name}_${i}"@ " \
    >${purge_dups}/${name}_${i}-config.json
  else
  echo "NO such assembler ${i} was used"
  fi
done


#step 4
assemblers="Flye canu nextdenovo"
for i in `echo ${assemblers}`
do
python3 ${pd_scripts}/run_purge_dups.py \
    --platform bash --wait 1 --retries 3 \
    ${purge_dups}/${name}_${i}-config.json /pub/jenyuw/Software/purge_dups/bin ${name}_${i}
done

