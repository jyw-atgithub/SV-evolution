#!/bin/bash

#SBATCH --job-name=polish    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=36   ## number of cores the job needs
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU

# racon --> nextpolisher +/- POLCA

trimmed="/home/jenyuw/SV-project/result/trimmed"
assemble="/home/jenyuw/SV-project/result/assemble"
aligned_bam="/home/jenyuw/SV-project/result/aligned_bam"
polishing="/home/jenyuw/SV-project/result/polishing"

### Racon*3  --> NextPolish*3


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


module load python/3.10.2