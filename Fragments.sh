## Fragments
: <<'SKIP'
SKIP

## nameed pipe test
nptest="/home/jenyuw/SV-project/np-test"
for j in $(ls ${raw}/nv1*_illumina_r1.fastq)
do
name=$(basename ${j} |sed "s/_illumina_r.*.fastq//g")
r2=$(echo $j |sed 's/r1/r2/')
#The rules (of using * and ?) in sed is different.
mkfifo ${name}.read1 ${name}.read2
fastp -i ${j} -I ${r2} -o ${nptest}/${name}.read1 -O ${nptest}/${name}.read2 --thread ${nT} --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_tail --average_qual 3
bwa mem -M -t ${nT} ${ref_genome} ${nptest}/${name}.read1 ${nptest}/${name}.read2 | samtools view -bh -|\
samtools sort -@ ${nT} -o ${nptest}/${name}.sort.bam -
done
## Failed, because NOT all programs can use named pipe
