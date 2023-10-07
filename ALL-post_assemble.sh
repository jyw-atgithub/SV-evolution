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
patched="/home/jenyuw/SV-project/result/patched_contigs"
scaffold="/home/jenyuw/SV-project/result/scaffold"
busco_out="/home/jenyuw/SV-project/result/busco_out"
## prep
source ~/.bashrc
nT=20


## changing the polishing strategy. Simplify everything, only 1 kind workflows!!
### Racon*3  --> NextPolish*3
### if illumina reads exists, then ADDitional polishing with POLAC*3 at the FRONT
### Maybe conpare POLAC with Pilon if illumina reads exists,
## After polishing, purge_dups --> rag-tag patch --> rag-tag scaffold
## Medaka??? No, set aside temporarily because it's only suggested to works on Flye assembly. Also, it requires the experiment info.
## POLCA performs ok but it also takes Illumina or PacBio Hifi
: <<'SKIP'
assembler="Flye"
bash /home/jenyuw/Software/MaSuRCA-4.1.0/bin/polca.sh \
-a ${assemble}/nv107_${assembler}/assembly.fasta \
-r ${trimmed}/nv107.trimmed.r1.fastq \
-r ${trimmed}/nv107.trimmed.r2.fastq \
-t ${nT} -m 4G
SKIP

function polish_Ra {
# THREE input argumants:"path" tech rounds
# ${i} will be one of the assemblers
for k in $(ls $1 2> /dev/null)
do
#echo "first arg is " $1
#echo "second arg is " $2
echo "k is " $k
name=$(echo $k | gawk -F "/" '{print $7}' | sed "s/_${i}//g")
#echo $name
read=${trimmed}/${name}.trimmed.fastq
read_type=$2
mapping_option=(["clr"]="map-pb" ["hifi"]="asm20" ["ont"]="map-ont")
  if [[ $2 != "clr" && $2 != "hifi" && $2 != "ont" ]]
  then
    echo "The second argument can only be one of \"clr, hifi, ont\""
  fi
round=$3
input=${k}
  for ((count=1; count<=${round};count++))
  do
  echo "round $count"
  echo "input is" $input
  minimap2 -x ${mapping_option[$read_type]} -t ${nT} -o ${aligned_bam}/${name}.trimmed-${i}.paf ${input} ${read}
  racon -t ${nT} ${read} ${aligned_bam}/${name}.trimmed-${i}.paf ${input} >${polishing}/${name}.${i}.racon.fasta
    if ((${count}!=${round}))
    then
      mv ${polishing}/${name}.${i}.racon.fasta ${polishing}/${name}.${i}.racontmp.fasta
      input=${polishing}/${name}.${i}.racontmp.fasta
    fi
  done
rm ${aligned_bam}/${name}.trimmed-${i}.paf
rm ${polishing}/${name}.${i}.racontmp.fasta
done
}

function polish_Np {
# THREE input argumants:"path" tech rounds
# ${i} will be one of the assemblers
for j in $(ls $1 2> /dev/null)
do
echo "polish_Np starts"
name=`echo $j|gawk -F "/" '{print $7}'|gawk -F "." '{print $1}'`
echo "name is $name"
round=$3
read=${trimmed}/${name}.trimmed.fastq
read_type=$2 #{clr,hifi,ont}
mapping_option=(["clr"]="map-pb" ["hifi"]="asm20" ["ont"]="map-ont")
  if [[ $2 != "clr" && $2 != "hifi" && $2 != "ont" ]]
  then
    echo "The second argument can only be one of \"clr, hifi, ont\""
  fi
input=${j}
  for ((count=1; count<=${round};count++))
  do
  echo "round $count"
  echo "input is" $input
  minimap2 -ax ${mapping_option[$read_type]} -t ${nT} ${input} ${read} |\
  samtools sort - -m 2g --threads ${nT} -o ${aligned_bam}/${name}.trimmed-${i}.sort.bam
  samtools index ${aligned_bam}/${name}.trimmed-${i}.sort.bam
  ls ${aligned_bam}/${name}.trimmed-${i}.sort.bam > ${polishing}/lgs.sort.bam.fofn
  python3 /home/jenyuw/Software/NextPolish/lib/nextpolish2.py -g ${input} -l ${polishing}/lgs.sort.bam.fofn \
  -r ${read_type} -p ${nT} -sp -o ${polishing}/${name}.${i}.nextpolish.fasta
    if ((${count}!=${round}));then
        mv ${polishing}/${name}.${i}.nextpolish.fasta ${polishing}/${name}.${i}.nextpolishtmp.fasta;
        input=${polishing}/${name}.${i}.nextpolishtmp.fasta;
    fi;
  done
rm ${polishing}/${name}.${i}.nextpolishtmp.fasta
rm ${aligned_bam}/${name}.trimmed-${i}.sort.bam
rm ${aligned_bam}/${name}.trimmed-${i}.sort.bam.bai
done
}


##Actual Polishing with racon, NCBI ONT data
conda activate post-proc
assembler="nextdenovo"

for i in `echo $assembler`
do
  if [[ $i == "Flye" ]]
  then
  echo "racon $i now"
  #ls ${assemble}/*_ONT_${i}/assembly.fasta 2> /dev/null
  polish_Ra "${assemble}/*_ONT_${i}/assembly.fasta" "ont" "3"
  elif [[ $i == "canu" ]]
  then
  echo "racon $i now"
  #ls ${assemble}/*_ONT_${i}/*.corrected.contigs.fasta 2> /dev/null
  polish_Ra "${assemble}/*_ONT_${i}/*.contigs.fasta" "ont" "3"
  elif [[ $i == "nextdenovo" ]]
  then
  echo "racon $i now"
  #ls ${assemble}/*_ONT_${i}/03.ctg_graph/nd.asm.fasta 2> /dev/null
  polish_Ra "${assemble}/*_ONT_${i}/03.ctg_graph/nd.asm.fasta" "ont" "3"
  else
  echo "NO such assembler was used"
  fi
done
conda deactivate

## Actual Polishing with NextPolish, NCBI ONT data
assembler="nextdenovo"

for i in $(echo $assembler)
do
  if [[ $i == "Flye" ]]
  then
  echo "Nextpolish $i now"
  #ls ${assemble}/*_ONT_${i}/assembly.fasta 2> /dev/null
  polish_Np "${polishing}/*_ONT.${i}.racon.fasta" "ont" "3"
  elif [[ $i == "canu" ]]
  then
  echo "Nestpolish $i now"
  #ls ${assemble}/*_ONT_${i}/*.corrected.contigs.fasta 2> /dev/null
  polish_Np "${polishing}/*_ONT.${i}.racon.fasta" "ont" "3"
  elif [[ $i == "nextdenovo" ]]
  then
  echo "Nextpolish $i now"
  #ls ${assemble}/*_ONT_${i}/03.ctg_graph/nd.asm.fasta 2> /dev/null
  polish_Np "${polishing}/*_ONT.${i}.racon.fasta" "ont" "3"
  else
  echo "NO such assembler was used"
  fi
done


##nv samples
conda activate post-proc
assembler="nextdenovo"

for i in `echo $assembler`
do
  if [[ $i == "Flye" ]]
  then
  echo "racon $i now"
  polish_Ra "${assemble}/nv*_${i}/assembly.fasta" "ont" "3"
  elif [[ $i == "canu" ]]
  then
  echo "racon $i now"
  polish_Ra "${assemble}/nv*_${i}/*.contigs.fasta" "ont" "3"
  elif [[ $i == "nextdenovo" ]]
  then
  echo "racon $i now"
  polish_Ra "${assemble}/nv*_${i}/03.ctg_graph/nd.asm.fasta" "ont" "3"
  else
  echo "NO such assembler was used"
  fi
done
conda deactivate

assembler="nextdenovo"
for i in $(echo $assembler)
do
  if [[ $i == "Flye" ]]
  then
  echo "Nextpolish $i now"
  polish_Np "${polishing}/nv*.${i}.racon.fasta" "ont" "3"
  elif [[ $i == "canu" ]]
  then
  echo "Nestpolish $i now"
  polish_Np "${polishing}/nv*.${i}.racon.fasta" "ont" "3"
  elif [[ $i == "nextdenovo" ]]
  then
  echo "Nextpolish $i now"
  polish_Np "${polishing}/nv*.${i}.racon.fasta" "ont" "3"
  else
  echo "NO such assembler was used"
  fi
done

##DSPR samples

conda activate post-proc
assembler="nextdenovo"
prefix="{A{1..7},AB8,B{{1..4},6}}"
for i in `echo $assembler`
do
  if [[ $i == "Flye" ]]
  then
  echo "racon $i now"
  polish_Ra "${assemble}/${prefix}_${i}/assembly.fasta" "clr" "3"
  elif [[ $i == "canu" ]]
  then
  echo "racon $i now"
  polish_Ra "${assemble}/${prefix}_${i}/*.contigs.fasta" "clr" "3"
  elif [[ $i == "nextdenovo-30" ]]
  then
  echo "racon $i now"
  polish_Ra "${assemble}/${prefix}_${i}/03.ctg_graph/nd.asm.fasta" "clr" "3"
  else
  echo "NO such assembler was used"
  fi
done
conda deactivate

assembler="nextdenovo"
for i in $(echo $assembler)
do
  if [[ $i == "Flye" ]]
  then
  echo "Nextpolish $i now"
  polish_Np "${polishing}/${prefix}.${i}.racon.fasta" "clr" "3"
  elif [[ $i == "canu" ]]
  then
  echo "Nestpolish $i now"
  polish_Np "${polishing}/${prefix}.${i}.racon.fasta" "clr" "3"
  elif [[ $i == "nextdenovo-30" ]]
  then
  echo "Nextpolish $i now"
  polish_Np "${polishing}/${prefix}.${i}.racon.fasta" "clr" "3"
  else
  echo "NO such assembler was used"
  fi
done



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

##For NextPolish
##According to the author, making the loops manually is faster than using the package.
##The BWA contained in the original package is broken. We need to `git clone` the original bwa, so the installation (make) can be success. 
##because there both python2 and python3, change the shebang line as "#!/usr/bin/env python3"
## Ref:https://nextpolish.readthedocs.io/en/latest/TUTORIAL.html

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
purge_dups="/home/jenyuw/SV-project/result/purge_dups"
pd_scripts="/home/jenyuw/Software/purge_dups/scripts"

ls ${trimmed}/nv107.trimmed.fastq >${purge_dups}/read.fofn
python3 ${pd_scripts}/pd_config.py \
    -l ${purge_dups}/nv107_Flye \
    -n ${purge_dups}/TESTconfig.json \
    ${assemble}/nv107_Flye/assembly.fasta \
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
    \"mnmp_opt\": \"\",
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
  \"ref\": "\"/home/jenyuw/SV-project/result/purge_dups/nv107_Flye/assembly.fasta\"",
  \"out_dir\": "\"${purge_dups}/nv107_Flye\""
}" >${purge_dups}/TESTconfig.json

python3 ${pd_scripts}/run_purge_dups.py \
    --platform bash --wait 1 --retries 3 \
    ${purge_dups}/TESTconfig.json /home/jenyuw/Software/purge_dups/bin nv107

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
busco -i ${i} --out_path ${busco_out} -o ${strain}_${assembler} -m genome --cpu 8 -l diptera_odb10
done

assembler="canu"
for i in $(ls ${assemble}/*_${assembler}/*.contigs.fasta)
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


# BUSCO score of POLISHED assembly
for i in $(ls ${polishing}/*.fasta)
do
name=$(basename ${i}|sed s/.fasta// )
echo $name 
busco -i ${i} --out_path ${busco_out} -o ${name} -m genome --cpu 20 -l diptera_odb10
done

#collect the BUSCO scores
## original assembly
assembler="canu Flye nextdenovo nextdenovo-30"
for a in $(echo ${assembler} )
do
    for i in $(ls ${busco_out}/*_${a}/short_summary.specific.diptera_odb10.*.txt)
    do
    name=$(echo $i| gawk -F "/" ' {print $7}' | sed "s/_${a}//g")
    grep "C:" $i|awk '{print substr($0, 4, 5)}'|tr -d "\n"|tr -d "["|tr -d "%" && echo -e "\t${name}\t${a}\t"original"\t"busco""
    done
done
## Polished assembly
assembler="canu Flye nextdenovo nextdenovo-30"
for a in $(echo ${assembler} )
do
    for i in $(ls ${busco_out}/*.${a}.*/short_summary.specific.diptera_odb10.*.txt 2>/dev/null)
    do
    strain=$(echo $i| gawk -F "/" ' {print $7}' | gawk -F "." ' {print $1}')
    ass=$(echo $i| gawk -F "/" ' {print $7}' | gawk -F "." ' {print $2}')
    polisher=$(echo $i| gawk -F "/" ' {print $7}' | gawk -F "." ' {print $3}')
    grep "C:" $i|awk '{print substr($0, 4, 5)}'|tr -d "\n" |tr -d "[" |tr -d "%"&& echo -e "\t${strain}\t${ass}\t${polisher}\t"busco""
    done
done


#awk '{print substr(s, i, n)}' substr function accepts three arguments
#s: input string
#i: start index, 1-based!
#n: length

python3 /home/jenyuw/Software/compleasm_kit/compleasm.py list --local --library_path /home/jenyuw/Software/compleasm_kit/mb_downloads/

## # Compleasm score of original assembly
assembler="Flye"
for i in $(ls ${assemble}/*_${assembler}/assembly.fasta)
do
strain=$(echo $i | gawk -F "/" '{print $7}'| sed "s/_${assembler}//g")
echo $strain
python3 /home/jenyuw/Software/compleasm_kit/compleasm.py run  -a $i -o ${compleasm_out}/${strain}.${assembler}.original -t $nT \
-l diptera -L "/home/jenyuw/Software/compleasm_kit/mb_downloads/"
done

assembler="canu"
for i in $(ls ${assemble}/*_${assembler}/*.contigs.fasta)
do
strain=$(echo $i | gawk -F "/" '{print $7}'| sed "s/_${assembler}//g")
echo $strain
python3 /home/jenyuw/Software/compleasm_kit/compleasm.py run  -a $i -o ${compleasm_out}/${strain}.${assembler}.original -t $nT \
-l diptera -L "/home/jenyuw/Software/compleasm_kit/mb_downloads/"
done

assembler="nextdenovo nextdenovo-30"
for a in `echo $assembler`
do
    for i in $(ls ${assemble}/*_${a}/03.ctg_graph/nd.asm.fasta)
    do
    strain=$(echo $i | gawk -F "/" '{print $7}'| sed "s/_${a}//g")
    echo $strain
    python3 /home/jenyuw/Software/compleasm_kit/compleasm.py run \
    -a $i -o ${compleasm_out}/${strain}.${a}.original -t $nT \
    -l diptera -L "/home/jenyuw/Software/compleasm_kit/mb_downloads/"
    done
done

## Polished assembly
for i in $(ls ${polishing}/*.fasta)
do
name=$(basename ${i}|sed s/.fasta// )
python3 /home/jenyuw/Software/compleasm_kit/compleasm.py run  -a $i -o ${compleasm_out}/${name} -t 8 \
-l diptera -L "/home/jenyuw/Software/compleasm_kit/mb_downloads/"
done

## Collect the Compleasm scores
for i in $(ls ${compleasm_out}/*/summary.txt)
do
strain=$(echo $i|gawk -F "/" '{print $7}'|gawk -F "." '{print $1}')
ass=$(echo $i|gawk -F "/" '{print $7}'|gawk -F "." '{print $2}')
polisher=$(echo $i|gawk -F "/" '{print $7}'|gawk -F "." '{print $3}')
grep "S:" $i |gawk -F "%" '{print $1}'|tr -d "\n|sed s/"S:"//g"
echo -e "\t${strain}\t${ass}\t${polisher}\t"compleasm""
done