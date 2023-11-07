#!/bin/bash

#SBATCH --job-name=racon    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p highmem       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=36   ## number of cores the job needs
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU

# racon --> nextpolisher +/- POLCA
### Racon*3  --> NextPolish*3

trimmed="/dfs7/jje/jenyuw/SV-project-temp/result/trimmed"
assemble="/dfs7/jje/jenyuw/SV-project-temp/result/assemble"
aligned_bam="/dfs7/jje/jenyuw/SV-project-temp/result/aligned_bam"
polishing="/dfs7/jje/jenyuw/SV-project-temp/result/polishing"

source ~/.bashrc

nT=$SLURM_CPUS_PER_TASK
echo "cpu number is $SLURM_CPUS_PER_TASK"

function polish_Ra {
# THREE input argumants:"path" tech rounds
# ${i} will be one of the assemblers
for k in $(ls $1 2> /dev/null)
do
#echo "first arg is " $1
#echo "second arg is " $2
echo "k is " $k
name=$(echo $k | gawk -F "/" '{print $8}' | sed "s/_${i}//g")
echo $name
read=${trimmed}/${name}.trimmed.fastq.gz
read_type=$2
mapping_option=(["CLR"]="map-pb" ["hifi"]="asm20" ["ONT"]="map-ont")
    if [[ $2 != "CLR" && $2 != "hifi" && $2 != "ONT" ]]
    then
    echo "The second argument can only be one of \"CLR, hifi, ONT\""
    fi
round=$3
input=${k}
    for ((count=1; count<=${round};count++))
    do
    echo "round $count"
    echo "input is $input"
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


file=`head -n $SLURM_ARRAY_TASK_ID ${trimmed}/namelist.txt |tail -n 1`
name=$(basename ${file}|sed s/".trimmed.fastq.gz"//g)
read_type=`echo ${name} | gawk -F "_" '{print $2}'`
echo "the file is ${file}"
echo "the read type is ${read_type}"


conda activate post-proc
assembler="Flye"
for i in `echo $assembler`
do
    if [[ $i == "Flye" ]]
    then
    echo "racon $i assembly now"
    polish_Ra "${assemble}/${name}_${i}/assembly.fasta" "${read_type}" "3"
    elif [[ $i == "canu" ]]
    then
    echo "racon $i assembly now"
    polish_Ra "${assemble}/${name}_${i}/*.contigs.fasta" "${read_type}" "3"
    elif [[ $i == "nextdenovo" ]]
    then
    echo "racon $i assembly now"
    polish_Ra "${assemble}/${name}_${i}/03.ctg_graph/nd.asm.fasta" "${read_type}" "3"
    else
    echo "NO such assembler was used"
    fi
done
conda deactivate
echo "It is the end!!"