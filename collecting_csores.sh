busco_out="/dfs7/jje/jenyuw/SV-project-temp/result/busco_out"
assembler="canu flye nextdenovo-45"
for a in $(echo ${assembler} )
do
    for i in $(ls ${busco_out}/*_${a}_busco_ori/short_summary.specific.diptera_odb10.*_${a}_busco_ori.txt)
    do
    name=$(echo $i| gawk -F "/" ' {print $8}' | sed "s/_${a}_busco_ori//g")
    grep "C:" $i|awk '{print substr($0, 4, 5)}'|tr -d "\n"|tr -d "["|tr -d "%" && echo -e "\t${name}\t${a}\t"original"\t"busco""
    done
done


for i in $(ls ${busco_out}/*_final_busco_final/short_summary.specific.diptera_odb10.*_final_busco_final.txt 2>/dev/null)
do
strain=$(echo $i| gawk -F "/" ' {print $8}' | gawk -F "_" ' {print $1"_"$2}')
ass=$(echo $i| gawk -F "/" ' {print $8}' | gawk -F "_" ' {print $3}')
polisher=$(echo $i| gawk -F "/" ' {print $8}' | gawk -F "_" ' {print $3}')
grep "C:" $i|awk '{print substr($0, 4, 5)}'|tr -d "\n" |tr -d "[" |tr -d "%"&& echo -e "\t${strain}\t${ass}\t${polisher}\t"busco""
done


busco_out="/dfs7/jje/jenyuw/SV-project-temp/result/busco_out"
assembler="canu flye nextdenovo-45"
for a in $(echo ${assembler} )
do
    for i in $(ls ${busco_out}/*_${a}_busco_ori/short_summary.specific.diptera_odb10.*_${a}_busco_ori.txt)
    do
    name=$(echo $i| gawk -F "/" ' {print $8}' | sed "s/_${a}_busco_ori//g")
    grep "Complete and duplicated" $i|gawk '{print $1}'|tr -d "\n" && echo -e "\t${name}\t${a}\t"original"\t"DUP""
    done
done

for i in $(ls ${busco_out}/*_final_busco_final/short_summary.specific.diptera_odb10.*_final_busco_final.txt 2>/dev/null)
do
strain=$(echo $i| gawk -F "/" ' {print $8}' | gawk -F "_" ' {print $1"_"$2}')
ass=$(echo $i| gawk -F "/" ' {print $8}' | gawk -F "_" ' {print $3}')
polisher=$(echo $i| gawk -F "/" ' {print $8}' | gawk -F "_" ' {print $3}')
grep "Complete and duplicated" $i|gawk '{print $1}'|tr -d "\n" &&  echo -e "\t${strain}\t${ass}\t${polisher}\t"DUP""
done
