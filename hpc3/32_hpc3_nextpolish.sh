#!/bin/bash

#SBATCH --job-name=np    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p highmem        ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=24   ## number of cores the job needs
#SBATCH --mem-per-cpu=10G     # requesting memory per CPU

# racon --> nextpolisher +/- POLCA
### Racon*3  --> NextPolish*3

trimmed="/dfs7/jje/jenyuw/SV-project-temp/result/trimmed"
assemble="/dfs7/jje/jenyuw/SV-project-temp/result/assemble"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
polishing="/dfs7/jje/jenyuw/SV-project-temp/result/polishing"
nplib="/pub/jenyuw/Software/NextPolish/lib"
source ~/.bashrc

nT=$SLURM_CPUS_PER_TASK
echo "cpu number is $SLURM_CPUS_PER_TASK"



function polish_Np {
# THREE input argumants:"path" tech rounds
# ${i} will be one of the assemblers
for j in $(ls $1 2> /dev/null)
do
echo "polish_Np starts"
name=`echo $j|gawk -F "/" '{print $8}'|gawk -F "." '{print $1}'`
echo "name is $name"
round=$3
read=${trimmed}/${name}.trimmed.fastq.gz
read_type=$2 
echo "the second argument is $2, is read_type $read_type"
declare -A mapping_option=(["CLR"]='map-pb' ["hifi"]='asm20' ["ONT"]='map-ont')
echo "The mapping option is ${mapping_option[$read_type]}"
    if [[ $2 != "CLR" && $2 != "hifi" && $2 != "ONT" ]]
    then
    echo "The second argument can only be one of \"CLR, hifi, ONT\""
    fi
input=${j}
    for ((count=1; count<=${round};count++))
    do
    echo "round $count"
    echo "input is" $input
    minimap2 -a -x ${mapping_option[$read_type]} -t ${nT} ${input} ${read} |\
    samtools sort - -m 2g --threads ${nT} -o ${aligned_bam}/${name}.trimmed-${i}.sort.bam
    samtools index ${aligned_bam}/${name}.trimmed-${i}.sort.bam
    ls ${aligned_bam}/${name}.trimmed-${i}.sort.bam > ${polishing}/lgs.sort.bam.fofn
    python3 ${nplib}/nextpolish2.py -g ${input} -l ${polishing}/lgs.sort.bam.fofn \
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

file=`head -n $SLURM_ARRAY_TASK_ID ${trimmed}/namelist.txt |tail -n 1`
name=$(basename ${file}|sed s/".trimmed.fastq.gz"//g)
read_type=`echo ${name} | gawk -F "_" '{print $2}'`
echo "the file is ${file}"
echo "the read type is ${read_type}"

module load python/3.10.2

assembler="Flye"
for i in $(echo $assembler)
do
    if [[ $i == "Flye" ]]
    then
    echo "Nextpolish $i now"
    polish_Np "${polishing}/${name}.${i}.racon.fasta" "${read_type}" "3"
    elif [[ $i == "canu" ]]
    then
    echo "Nestpolish $i now"
    polish_Np "${polishing}/${name}.${i}.racon.fasta" "${read_type}" "3"
    elif [[ $i == "nextdenovo" ]]
    then
    echo "Nextpolish $i now"
    polish_Np "${polishing}/${name}.${i}.racon.fasta" "${read_type}" "3"
    else
    echo "NO such assembler was used"
    fi
done

module unload python/3.10.2
