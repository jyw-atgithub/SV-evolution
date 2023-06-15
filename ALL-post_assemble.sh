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


${assembler}

##Polishing with NextPolish
##nv samples
assembler="Flye"
for j in $(ls ${assemble}/nv*_${assembler}/assembly.fasta)
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
        samtools sort - -m 2g --threads ${nT} -o ${aligned_bam}/${name}.trimmed-${assembler}.sort.bam;
        samtools index ${aligned_bam}/${name}.trimmed-${assembler}.sort.bam;
        ls ${aligned_bam}/${name}.trimmed-${assembler}.sort.bam > ${polishing}/lgs.sort.bam.fofn;
        python3 /home/jenyuw/Software/NextPolish/lib/nextpolish2.py -g ${input} -l ${polishing}/lgs.sort.bam.fofn \
        -r ${read_type} -p ${nT} -sp -o ${polishing}/${name}.${assembler}.nextpolish.fasta;
        # Finally polished genome file: ${name}.nextpolish.fasta
        if ((i!=${round}));then
            mv ${polishing}/${name}.nextpolish.fasta ${polishing}/${name}.${assembler}.nextpolishtmp.fasta;
            input=${polishing}/${name}.${assembler}.nextpolishtmp.fasta;
        fi;
    done;

rm ${polishing}/${name}.${assembler}.nextpolishtmp.fasta
rm ${aligned_bam}/${name}.trimmed-${assembler}.sort.bam
rm ${aligned_bam}/${name}.trimmed-${assembler}.sort.bam.bai
done

##According to the author, making the loops manually is faster than using the package.
##The BWA contained in the original package is broken. We need to `git clone` the original bwa, so the installation (make) can be success. 
##because there both python2 and python3, change the shebang line as "#!/usr/bin/env python3"

conda activate assemble #this is for bwa-mem2

## trim the short-reads, only need once
: << 'SKIP'
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
done
SKIP

##Polishing with Pilon
##nv samples
for k in $(ls ${polishing}/nv*.${assembler}.nextpolish.fasta)
do
name=$(echo $k | sed "s@${polishing}\/@@g; s@${assembler}.nextpolish.fasta@@g")
r1="${trimmed}/${name}.trimmed.r1.fastq"
r2="${trimmed}/${name}.trimmed.r2.fastq"
 
round=3
input=${k}

    for ((i=1; i<=${round};i++)); do
    bwa-mem2 index ${input}
    bwa-mem2 mem -t ${nT} ${input} $r1 $r2 |samtools view -S -b -h |\
    samtools sort -@ ${nT} -o ${aligned_bam}/${name}.ILL-nextpolish.sort.bam
    samtools index ${aligned_bam}/${name}.ILL-nextpolish.sort.bam

    java -XX:+AggressiveHeap -jar /home/jenyuw/Software/pilon-1.24.jar --diploid --minqual 7 \
    --genome ${input} \
    --frags ${aligned_bam}/${name}.ILL-nextpolish.sort.bam \
    --output ${name}.${assembler}.pilon --outdir ${polishing}
        if ((i!=${round}));then
            mv ${polishing}/${name}.${assembler}.pilon.fasta ${polishing}/${name}.${assembler}.pilontmp.fasta;
            input=${polishing}/${name}.${assembler}.pilontmp.fasta;
        fi;
    done
rm ${aligned_bam}/${name}.ILL-nextpolish.sort.bam
rm ${aligned_bam}/${name}.ILL-nextpolish.sort.bam.bai
rm ${polishing}/${name}.${assembler}.pilontmp.fasta
done

#--threads is not supported by Pilon anymore
#Do NOT use the pilon installed by Anaconda, it will crash because of memory limit.
#Just download the precompiled jar file from the latest release on Github.
#"java -Xmx128G -jar" means giving the program an allowance of 128Gb memory.
#java -XX:+AggressiveHeap -jar means letting the program use as much memory as needed.


##NCBI-Nanopore samples
conda activate post-proc #this contain Racon, Ragtag
##Polishing with racon
assembler=Flye

for k in $(ls ${assemble}/*_ONT_${assembler}/assembly.fasta)
do
name=$(echo $k | gawk -F "/" '{print $7}' | sed "s/_${assembler}//g")
read=${trimmed}/${name}.trimmed.fastq

round=3
input=${k}

    for ((i=1; i<=${round};i++))
    do
    echo "round $i"
    minimap2 -x map-ont -t ${nT} -o ${aligned_bam}/${name}.trimmed-${assembler}.paf ${input} ${read}
    #samtools sort - -m 2g --threads ${nT} -o ${aligned_bam}/${name}.trimmed-${assembler}.sort.bam
    #samtools index ${aligned_bam}/${name}.trimmed-${assembler}.sort.bam
    racon -t ${nT} ${read} ${aligned_bam}/${name}.trimmed-${assembler}.paf ${input} >${polishing}/${name}.${assembler}.racon.fasta
        if ((i!=${round}));then
        mv ${polishing}/${name}.${assembler}.racon.fasta ${polishing}/${name}.${assembler}.racontmp.fasta;
        input=${polishing}/${name}.${assembler}.racontmp.fasta;
        fi;
    done
rm ${aligned_bam}/${name}.trimmed-${assembler}.paf
rm ${polishing}/${name}.${assembler}.racontmp.fasta
done
conda deactivate

##Then, polishing with Nextpolish
for j in $(ls ${polishing}/*_ONT.${assembler}.racon.fasta)
do
name=$(echo $j|gawk -F "/" '{print $7}'|sed "s/.${assembler}.racon.fasta//")
echo $name
#Set input and parameters
round=3
read=${trimmed}/${name}.trimmed.fastq
read_type=ont #{clr,hifi,ont}
mapping_option=(["clr"]="map-pb" ["hifi"]="asm20" ["ont"]="map-ont")
input=${j}

    for ((i=1; i<=${round};i++)); do
        minimap2 -ax ${mapping_option[$read_type]} -t ${nT} ${input} ${read} |\
        samtools sort - -m 2g --threads ${nT} -o ${aligned_bam}/${name}.trimmed-${assembler}.sort.bam;
        samtools index ${aligned_bam}/${name}.trimmed-${assembler}.sort.bam;
        ls ${aligned_bam}/${name}.trimmed-${assembler}.sort.bam > ${polishing}/lgs.sort.bam.fofn;
        python3 /home/jenyuw/Software/NextPolish/lib/nextpolish2.py -g ${input} -l ${polishing}/lgs.sort.bam.fofn \
        -r ${read_type} -p ${nT} -sp -o ${polishing}/${name}.${assembler}.nextpolish.fasta;
        # Finally polished genome file: ${name}.nextpolish.fasta
        if ((i!=${round}));then
            mv ${polishing}/${name}.${assembler}.nextpolish.fasta ${polishing}/${name}.${assembler}.nextpolishtmp.fasta;
            input=${polishing}/${name}.${assembler}.nextpolishtmp.fasta;
        fi;
    done;
rm ${polishing}/${name}.${assembler}.nextpolishtmp.fasta
rm ${aligned_bam}/${name}.trimmed-${assembler}.sort.bam
rm ${aligned_bam}/${name}.trimmed-${assembler}.sort.bam.bai
done


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




##Patching

for i in $(ls ${polishing}/)
do
ragtag.py patch <target.fa> <query.fa>
  -o PATH              output directory [./ragtag_output]
  -w                   overwrite intermediate files
done

## scaffolding.sh

for i in $(ls ${polishing}/*.polished.fasta)
do
name=$(basename $i |sed s/.polished.fasta//g)
echo $name
ragtag.py scaffold -t ${nT} -o ${scaffold}/${name} $ref_genome ${i}
done



conda activate busco
# BUSCO score of non-polished assembly
assembler="Flye"
for i in $(ls ${assemble}/*_${assembler}/assembly.fasta)
do
strain=$(echo $i | gawk -F "/" '{print $7}'| sed "s/_${assembler}//g")
echo $strain
busco -i ${i} --out_path ${busco_out} -o ${strain}_${assembler} -m genome --cpu 20 -l diptera_odb10
done

assembler="canu"
for i in $(ls ${assemble}/*_${assembler}/*.corrected.contigs.fasta)
do
strain=$(echo $i | gawk -F "/" '{print $7}'| sed "s/_${assembler}//g")
echo $strain
busco -i ${i} --out_path ${busco_out} -o ${strain}_${assembler} -m genome --cpu 20 -l diptera_odb10
done

assembler="nextdenovo nextdenovo-30"
for a in `echo $assembler`
do
    for i in $(ls ${assemble}/*_${a}/03.ctg_graph/nd.asm.fasta)
    do
    strain=$(echo $i | gawk -F "/" '{print $7}'| sed "s/_${a}//g")
    echo $strain
    busco -i ${i} --out_path ${busco_out} -o ${strain}_${a} -m genome --cpu 20 -l diptera_odb10
    done
done

assembler="canu Flye nextdenovo nextdenovo-30"
for a in $(echo ${assembler} )
do
    for i in $(ls ${busco_out}/*_${a}/short_summary.specific.diptera_odb10.*.txt)
    do
    name=$(echo $i| gawk -F "/" ' {print $7}' | sed "s/_${a}//g")
    grep "C:" $i|awk '{print substr($0, 4, 5)}'|tr -d "\n" && echo -e "\t${name}\t${a}"
    done
done 

#awk '{print substr(s, i, n)}' substr function accepts three arguments
#s: input string
#i: start index, 1-based!
#n: length