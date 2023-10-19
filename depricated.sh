## Unused

##polishing, without illumina reads or FAST5. --> Medaka or Racon.
##Medaka is designed to be used on Flye assembly directly!!
##Medaka was installed with "pip" because Anaconda kept failing.
conda activate post-proc


for i in $(ls ${assemble}/SRR*_ONT_Flye/assembly.fasta)
do
echo $i 
name=$(echo $i | gawk -F "\/" '{print $7}' 2>>/dev/null| sed s/_ONT_Flye//g)
echo $name

medaka_consensus -i ${raw}/PRJNA929424/${name}_ONT.fastq.gz -d ${i} -o ${polishing}/${name}-medaka-1 \
-t ${nT} -m r941_min_hac_g507
done


for i in /home/jenyuw/SV-project/result/assemble/SRR23269563_ONT_Flye/assembly.fasta
do
echo $i 
name=$(echo $i | gawk -F "\/" '{print $7}' 2>>/dev/null| sed s/_ONT_Flye//g)
echo $name

medaka_consensus -i /home/jenyuw/SV-project/raw/PRJNA929424/SRR23269563_ONT.fastq.gz \
-d /home/jenyuw/SV-project/result/assemble/SRR23269563_ONT_Flye/assembly.fasta \
-o ${polishing}/SRR23269563-medaka-1/ \
-t ${nT} -m r941_min_hac_g507
done


minimap2 -t ${nT} -B 5 -a -x map-ont \
$i ${trimmed}/${name}_ONT.trimmed.fastq |\
samtools view -b -h -@ ${nT} -o - |\
samtools sort -@ ${nT} -o ${aligned_bam}/${name}_ONT.trimmed_assembly.sort.bam
samtools index -@ ${nT} ${aligned_bam}/${name}_ONT.trimmed_assembly.sort.bam

medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${nT}\
-m r941_min_hac_g507

## purge_dups
purge_dups="/home/jenyuw/SV-project/result/purge_dups"
pd_scripts="/home/jenyuw/Software/purge_dups/scripts"

ls ${trimmed}/nv107.trimmed.fastq >${purge_dups}/read.fofn

python3 ${pd_scripts}/pd_config.py \
    -l ${purge_dups}/nv107.canu.nextpolish \
    -n ${purge_dups}/TESTconfig.json \
    ${polishing}/nv107.canu.nextpolish.fasta \
    ${purge_dups}/read.fofn
nT="2"
echo -e "
{
  \"cc\": {
    \"fofn\": \"${purge_dups}/read.fofn\",
    \"isdip\": 1,
    \"core\": ${nT},
    \"mem\": 60000,
    \"queue\": \"normal\",
    \"mnmp_opt\": \"-t ${nT} -x map-ont \",
    \"bwa_opt\": \"\",
    \"ispb\": 1,
    \"skip\": 0
  },
  \"sa\": {
    \"core\": ${nT},
    \"mem\": 40000,
    \"queue\": \"normal\"
  },
  \"busco\": {
    \"core\": ${nT},
    \"mem\": 60000,
    \"queue\": \"long\",
    \"skip\": 0,
    \"lineage\": \"diptera\",
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
  \"ref\": "\"/home/jenyuw/SV-project/result/purge_dups/nv107_Flye/assembly.fasta\"",
  \"out_dir\": "\"${purge_dups}/nv107_Flye\""
}" >${purge_dups}/TESTconfig.json

python3 ${pd_scripts}/run_purge_dups.py \
    --platform bash --wait 1 --retries 3 \
    ${purge_dups}/TESTconfig.json /home/jenyuw/Software/purge_dups/bin SPID