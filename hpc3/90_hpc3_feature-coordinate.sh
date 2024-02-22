#!/bin/bash
conda activate everything
cd /home/jenyuw/DSPR_snp/snp_processing

# ToDo:
# 1: Merge category (eg 5putr)
# 2: Subtract out all other genic categories (eg cds and 3putr).
# 00. Define files for use later.
basepath=/dfs7/jje/jenyuw/SV-project-temp/result/
#raw=$basepath/data/raw/
raw="/dfs7/jje/jenyuw/SV-project-temp/reference"
processed="${basepath}/processed_SNP"
dmel=dmel-all-r6.49
dmel_chr="dmel-all-chromosome-r6.49.fasta" #original
dmel_gff="dmel-all-r6.49.gff" #original from flybase



#counting features in gff
#grep -v ^## "dmel-all-filtered-r6.46.gff"| cut -f 3 | sort |uniq -c |sort -k 1|less

for featuretype in CDS intron exon five_prime_UTR three_prime_UTR gene ; do
  bioawk -c gff -v f=$featuretype \
    ' $feature == f { print $seqname "\t" $start-1 "\t" $end } ' \
    <(cat ${raw}/${dmel_gff}) \
  | bedtools merge \
    -i stdin \
  | sort -k1,1 -k2,2n -k3,3n \
  > $processed/$dmel.${featuretype}raw.bed
done


#modified som detailes: no cat for a single file.  For excluding complex loci. include more features??
bedtools subtract \
  -a $processed/$dmel.five_prime_UTRraw.bed \
  -b <(cat $processed/$dmel.{three_prime_UTR,CDS,intron}raw.bed) \
| bedtools merge -i stdin > $processed/$dmel.five_prime_UTRspans.bed

bedtools subtract \
  -a $processed/$dmel.three_prime_UTRraw.bed \
  -b <(cat $processed/$dmel.{five_prime_UTR,CDS,intron}raw.bed) \
| bedtools merge -i stdin >$processed/$dmel.three_prime_UTRspans.bed

bedtools subtract \
  -a $processed/$dmel.CDSraw.bed \
  -b <(cat $processed/$dmel.{five_prime_UTR,three_prime_UTR,intron}raw.bed) \
| bedtools merge -i stdin > $processed/$dmel.CDSspans.bed

bedtools subtract \
  -a $processed/$dmel.exonraw.bed \
  -b <(cat $processed/$dmel.{five_prime_UTR,three_prime_UTR,CDS,intron}raw.bed) \
| bedtools merge -i stdin > $processed/$dmel.exonspans.bed


# l=120: Smaller than or equal to 120bp is "not long introns".
bioawk -c gff -v l=120 \
  ' $feature == "intron" && $end - $start + 1 <= l { print $seqname "\t" $start - 1 "\t" $end } ' \
  <$raw/$dmel_gff \
| bedtools merge -i stdin \
> $processed/$dmel.notliraw.bed

# l=65: Greater than 65bp is "not short introns".
bioawk -c gff -v l=65 \
  ' $feature == "intron" && $end - $start + 1 > l { print $seqname "\t" $start - 1 "\t" $end } ' \
  <$raw/$dmel_gff \
| bedtools merge -i stdin \
> $processed/$dmel.notsiraw.bed

# Take raw short introns and remove any that overlap with introns that aren't short (bedtools subtract -A)
# Then subtract parts of the remaining short introns that overlap exons.
# Note that the second bioawk -c bed command is operating on bed coordinates (0 based 1/2 open intervals).
bioawk -c gff -v l=65 \
  ' $feature == "intron" && $end - $start + 1 <= l { print $seqname "\t" $start - 1 "\t" $end } ' \
  <$raw/$dmel_gff \
| bedtools subtract -A \
    -a stdin \
    -b $processed/$dmel.notsiraw.bed \
| bioawk -c bed -v s=7 -v e=30 \
    ' $end >= e { print $chrom "\t" $start+s "\t" $start+e } ' \
| bedtools subtract \
    -a stdin \
    -b $processed/$dmel.exonraw.bed \
| sort -k1,1 -k2,2n -k3,3n \
| bedtools merge -i stdin > $processed/$dmel.sispans.bed


bioawk -c gff -v l=120 -v s=7 \
  ' $feature == "intron" && $end - $start + 1 > l { print $seqname "\t" $start + s "\t" $end - s - 1 } ' \
<$raw/$dmel_gff \
| bedtools subtract \
    -a stdin \
    -b <(cat $processed/$dmel.{exon,notli}raw.bed) \
| sort -k1,1 -k2,2n -k3,3n \
| bedtools merge -i stdin > $processed/$dmel.lispans.bed


pl=4000
bedtools complement \
    -i $processed/$dmel.generaw.bed \
    -g <(faSize -detailed $raw/$dmel_chr | sort) \
| tee $processed/$dmel.igraw.bed \
| bioawk -c bed -v pl=$pl ' { s=$start+pl; e=$end-pl; if (s<=e) { $start=s; $end=e; print; } } ' \
| tee $processed/$dmel.digraw.bed |cat> $processed/$dmel.digspans.bed

bedtools subtract \
  -a $processed/$dmel.igraw.bed \
  -b $processed/$dmel.digraw.bed \
| tee $processed/$dmel.pigraw.bed |cat > $processed/$dmel.pigspans.bed