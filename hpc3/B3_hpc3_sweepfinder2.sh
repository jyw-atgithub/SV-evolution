#!/bin/bash

#!/bin/bash

#in interactive mode
polarizing="/dfs7/jje/jenyuw/SV-project-temp/result/polarizing"
Merged_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/Merged_SNP"
SF2="/dfs7/jje/jenyuw/SV-project-temp/result/sweepfinder2"

#add the number of alternate alleles and number of samples (NS; is AC divided by 2) and then create a new copy
bcftools +fill-tags ${polarizing}/3corrected.polarized.asm.vcf.gz -O z -o ${SF2}/filled.SV.asm.vcf.gz -- -t AC,NS
bcftools index -f -t ${SF2}/filled.SV.asm.vcf.gz
#transform the SV vcf files into the required format
#Actually NS is what we want for "the allele count (x)"
#now, the sample number (n) is 62
#everything is unfolded, so all the values at folded colum is 0
#remove the rows with x=0
echo -e "position\tx\tn\tfolded" > ${SF2}/SV.asm.input.txt
bcftools query -r "2L,2R" -f '%POS\t%NS\t62\t0\n' ${SF2}/filled.SV.asm.vcf.gz|\
gawk '$2>0' >> ${SF2}/SV.asm.input.txt

#This is a super simple use
SweepFinder2 -sg 10000 ${SF2}/SV.asm.input.txt ${SF2}/SV.asm.output.txt