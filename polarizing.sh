#!/bin/bash

len=1000

bcftools view --thread 4 -i 'SVLEN<200 && SVLEN>90 && SVTYPE="INS" ' -r 2L:9000000-10000000 all.consensus-005.vcf.gz 2> /dev/null |\
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\n' |\
gawk -v len=$len ' $2-len > 0 {print $1 "\t" $2-len "\t" $3+len}' |sort|uniq|less -S