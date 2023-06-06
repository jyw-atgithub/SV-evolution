#!/bin/bash
## path
ref_genome="/home/jenyuw/SV-project/reference_genome/dmel-all-chromosome-r6.49.fasta"
raw="/home/jenyuw/SV-project/raw"
qc_report="/home/jenyuw/SV-project/result/qc_report"
trimmed="/home/jenyuw/SV-project/result/trimmed"
assemble="/home/jenyuw/SV-project/result/assemble"
aligned_bam="/home/jenyuw/SV-project/result/aligned_bam"
polishing="/home/jenyuw/SV-project/result/polishing"
canu_proc="/home/jenyuw/SV-project/result/canu_processing"
scaffold="/home/jenyuw/SV-project/result/scaffold"
busco_out="/home/jenyuw/SV-project/result/busco_out"
## prep
source ~/.bashrc
nT=30



##Polishing with NextPolish
##nv samples
for j in $(ls ${assemble}/nv*_Flye/assembly.fasta)
do
name=`echo $j|gawk -F "/" '{print $7}'|sed s/_Flye//`
echo $name
#Set input and parameters
round=3
read=${trimmed}/${name}.trimmed.fastq
read_type=ont #{clr,hifi,ont}
mapping_option=(["clr"]="map-pb" ["hifi"]="asm20" ["ont"]="map-ont")
input=${j}

    for ((i=1; i<=${round};i++)); do
        minimap2 -ax ${mapping_option[$read_type]} -t ${nT} ${input} ${read} |\
        samtools sort - -m 2g --threads ${nT} -o ${aligned_bam}/${name}.trimmed-Flye.sort.bam;
        samtools index ${aligned_bam}/${name}.trimmed-Flye.sort.bam;
        ls ${aligned_bam}/${name}.trimmed-Flye.sort.bam > ${polishing}/lgs.sort.bam.fofn;
        python3 /home/jenyuw/Software/NextPolish/lib/nextpolish2.py -g ${input} -l ${polishing}/lgs.sort.bam.fofn \
        -r ${read_type} -p ${nT} -sp -o ${polishing}/${name}.nextpolish.fasta;
        # Finally polished genome file: ${name}.nextpolish.fasta
        if ((i!=${round}));then
            mv ${polishing}/${name}.nextpolish.fasta ${polishing}/${name}.nextpolishtmp.fasta;
            input=${polishing}/${name}.nextpolishtmp.fasta;
        fi;
    done;
done

conda activate assemble #this is for bwa-mem2

##Polishing with Pilon
##nv samples
for k in $(ls ${polishing}/nv*.nextpolish.fasta)
do
name=$(echo $k | sed "s@${polishing}\/@@g; s@.nextpolish.fasta@@g")
r1="${trimmed}/${name}.trimmed.r1.fastq"
r2="${trimmed}/${name}.trimmed.r2.fastq"
 
round=3
input=${k}

for ((i=1; i<=${round};i++)); do
bwa-mem2 index ${input}
bwa-mem2 mem -t ${nT} ${input} $r1 $r2 |samtools view -S -b -h |\
samtools sort -@ ${nT} -o ${aligned_bam}/${name}.ILL-Flye.sort.bam
samtools index ${aligned_bam}/${name}.ILL-Flye.sort.bam

java -XX:+AggressiveHeap -jar /home/jenyuw/Software/pilon-1.24.jar --diploid --minqual 7 \
--genome ${assemble}/${name}/assembly.fasta \
--frags ${aligned_bam}/${name}.ILL-Flye.sort.bam \
--output ${name}.pilon --outdir ${polishing}
    if ((i!=${round}));then
        mv ${polishing}/${name}.pilon.fasta ${polishing}/${name}.pilontmp.fasta;
        input=${polishing}/${name}.pilontmp.fasta;
    fi;
done
done
#--threads is not supported by Pilon anymore
#Do NOT use the pilon installed by Anaconda, it will crash because of memory limit.
#Just download the precompiled jar file from the latest release on Github.
#"java -Xmx128G -jar" means giving the program an allowance of 128Gb memory.
#java -XX:+AggressiveHeap -jar means letting the program use as much memory as needed.





## Short-read mapping with sorting
#bwa index ${ref_genome}
for j in $(ls ${raw}/nv1*_illumina_r1.fastq)
do
echo $j
name=$(basename ${j} |sed "s/_illumina_r.*.fastq//g")
#The rules (of using * and ?) in sed is different.
echo $name
r2=$(echo $j |sed 's/r1/r2/')
echo $r2
fastp -i ${j} -I ${r2} -o ${trimmed}/${name}.trimmed.r1.fastq -O ${trimmed}/${name}.trimmed.r2.fastq \
--thread ${nT} --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_tail --average_qual 3

bwa mem -M -t ${nT} ${ref_genome} ${trimmed}/${name}.trimmed.r1.fastq ${trimmed}/${name}.trimmed.r2.fastq |\
samtools sort -@ ${nT} - -o ${aligned_bam}/${name}.sort.bam
done





## scaffolding.sh
conda activate post-proc #this contain Racon, Ragtag
for i in $(ls ${polishing}/*.polished.fasta)
do
name=$(basename $i |sed s/.polished.fasta//g)
echo $name
ragtag.py scaffold -t ${nT} -o ${scaffold}/${name} $ref_genome ${i}
done


conda activate busco
for i in $(ls ${assemble}/nv*_Flye/assembly.fasta)
do
strain=$(echo $i | gawk -F "\/" '{print $7}' 2>>/dev/null| sed s/_Flye//g)
echo $strain
busco -i ${i} --out_path ${busco_out} -o ${strain} -m genome --cpu 8 -l diptera_odb10
done

busco -i nv107.polished.fasta -o nv107.p -m genome --cpu 20 -l diptera_odb10

busco -i nv107.polished.fasta -o nv107.p -m genome --cpu 20 -l diptera_odb10