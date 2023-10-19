#!/bin/bash
guppy6="/home/jenyuw/SV-project/raw/guppy6"

for i in $(echo /home/jenyuw/NV_reads/nv115/{rapid_20180608,rapid_20181218,ligat_20180611}/fast5)
do
name=$(echo $i | gawk -F "/" ' {print $6}' )
guppy_basecaller -i $i -s ${guppy6}/nv115-${name}/ -c dna_r9.4.1_450bps_hac.cfg \
--recursive --min_qscore 7 --device "cuda:0" \
--detect_adapter --detect_mid_strand_adapter \
--detect_primer --detect_barcodes --enable_trim_barcodes
done

for i in $(echo /home/jenyuw/NV_reads/nv109/{ligat_20180605,ligat_20180606}/fast5)
do
name=$(echo $i | gawk -F "/" ' {print $6}' )
guppy_basecaller -i $i -s ${guppy6}/nv109-${name}/ -c dna_r9.4.1_450bps_hac.cfg \
--recursive --min_qscore 7 --device "cuda:0" \
--detect_adapter --detect_mid_strand_adapter \
--detect_primer --detect_barcodes --enable_trim_barcodes
done

for i in $(echo /home/jenyuw/NV_reads/nv107/{ligat_20180504,rapid_20180430,rapid_20180518}/fast5)
do
name=$(echo $i | gawk -F "/" ' {print $6}' )
guppy_basecaller -i $i -s ${guppy6}/nv107-${name}/ -c dna_r9.4.1_450bps_hac.cfg \
--recursive --min_qscore 7 --device "cuda:0" \
--detect_adapter --detect_mid_strand_adapter \
--detect_primer --detect_barcodes --enable_trim_barcodes
done


# Do NOT use CPU to call the bases, because it is too time consuming.
: <<'SKIP'
for i in $(echo /home/jenyuw/NV_reads/nv115/{rapid_20181218,ligat_20180611}/fast5)
do
name=$(echo $i | gawk -F "/" ' {print $6}' )
guppy_basecaller -i $i -s ${guppy6}/nv115-${name}/ -c dna_r9.4.1_450bps_hac.cfg \
--recursive --min_qscore 7 \
--num_callers 3 --cpu_threads_per_caller 8 \
--as_cpu_threads_per_scaler 8 \
--detect_adapter --detect_mid_strand_adapter \
--detect_primer --detect_barcodes --enable_trim_barcodes
done
SKIP