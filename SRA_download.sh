# for Pacbio, ideally, there is no paired read, so we don't have to use --split3. 
# (Actually, there are some paired reads for ??? reasones.)
# --split3 or --split-spot will ignore technical reads, so --skip-technical is not required.
for i in /home/jenyuw/SV-project/raw/SRR7*/*.sra
do
name=`basename $i | sed 's/.sra/.fastq.gz/'`
fastq-dump --split-spot --stdout  $i | pigz -p 8 -v  >/home/jenyuw/SV-project/raw/A8_pacbio/$name
#with --stdout, the -O is ignored
done


nlist=$(gawk '{printf $2 " "}' /home/jenyuw/SV-project/raw/DSPR_NCBI_SRX_list.tsv)
nline=$(wc -l </home/jenyuw/SV-project/raw/DSPR_NCBI_SRX_list.tsv)

for i in $(echo {1..$nline})
do
echo $i
line=$(head -n $i | tail -n 1)
echo $line
done
# echo {1..${nline}} # result is {1..14}, not an array

##NCBI_DSPR_download.sh
while read p
do
echo "$p"
strain=$(echo "$p" | gawk '{print $2}')
SRX=$(echo "$p" | gawk '{print $1}')
echo $strain
echo $SRX
mkdir /home/jenyuw/SV-project/raw/${strain}_pacbio
prefetch $SRX -O /home/jenyuw/SV-project/raw/${strain}_pacbio
wait
done </home/jenyuw/SV-project/raw/DSPR_NCBI_SRX_list.tsv

DSPR_NCBI_SRX_list.tsv
-------------------------first download batch, Finished
SRX4728152      AB8
SRX4599325      A2
SRX4713156      A4
SRX4699481      B1
SRX4692442      A3
SRX4668629      ORE
------------------------------second download batch
SRX4722917      A6
SRX4668167      B6
SRX4661415      B4
SRX4646796      A7
SRX4645948      B3
SRX4632943      B2
SRX4620074      A1
SRX4619206      A5

##NCBI_DSPR_conversion.sh
## file conversion (fastq-dump)
for i in /home/jenyuw/SV-project/raw/*_pacbio/SRR*/*.sra
do
echo $i
strain=`echo $i | gawk -F "\/" '{print $6}' 2>>/dev/null`
echo $strain
name=`basename $i | sed 's/.sra/.fastq.gz/'`
echo $name
nameSRR=`basename $i | sed 's/.sra//'`
echo $nameSRR
folder=/home/jenyuw/SV-project/raw/*_pacbio/${nameSRR}
echo $folder
fastq-dump --split-spot --stdout  $i | pigz -p 8 -v  >/home/jenyuw/SV-project/raw/${strain}/${name}
#with --stdout, the -O is ignored
rm -r $folder
done



##NCBI_download

# prefetch PRJNA929424 -O /home/jenyuw/SV-project/raw/PRJNA929424/
prefetch PRJNA574592 -O /home/jenyuw/SV-project/raw/PRJNA574592/
# this is fine, but it includes data other than DNA seq so we will waste on time to download unwanted items.


PRJNA929424 SRR23269563 ONT Denmark
PRJNA929424 SRR23269565 ONT Finland
PRJNA929424 SRR23269574 ONT Spain
PRJNA929424 SRR23269564 ONT Sweden
PRJNA929424 SRR23269572 ONT Turkey
PRJNA929424 SRR23269573 ONT Ukraine
PRJNA929424 SRR23269571 ONT  Zambia
PRJNA574592 SRR10199290 RS_II  RAL208
PRJNA574592 SRR10199524 RS_II  RAL379
PRJNA574592 SRR10199525 RS_II  RAL399
PRJNA574592 SRR10214519 RS_II  RAL799
PRJNA574592 SRR10200370 RS_II  RAL427
PRJNA574592 SRR10208187 RS_II  RAL517

for i in /home/jenyuw/SV-project/raw/PRJNA929424/SRR*/*.sra
do
echo $i
name=`basename $i | sed 's/.sra/_ONT.fastq.gz/'`
echo $name
fastq-dump --split-spot --stdout  $i | pigz -p 8 -v  >/home/jenyuw/SV-project/raw/PRJNA929424/${name}

#with --stdout, the -O is ignored
done

PRJ="PRJNA574592"
tech="pacbio"
prefetch ${PRJ} -O /home/jenyuw/SV-project/raw/${PRJ}/
for i in /home/jenyuw/SV-project/raw/${PRJ}/SRR*/*.sra
do
echo $i
# remember to use double quotes when a variable is inside!!
name=`basename $i | sed "s/.sra/_${tech}.fastq.gz/"`
echo $name
fastq-dump --split-spot --stdout  $i | pigz -p 8 -v  >/home/jenyuw/SV-project/raw/${PRJ}/${name}
done

for i in /home/jenyuw/SV-project/raw/*_pacbio/SRR*.fastq.gz
do
echo $i
strain=`echo $i | gawk -F "\/" '{print $6}' 2>>/dev/null|sed 's/_pacbio//g'`
echo $strain
cat $i >> /home/jenyuw/SV-project/raw/${strain}_pacbio/${strain}_combine_pacbio.fastq.gz
done

Long and short read sequencing and assembly of Natural Drosophila melanogaster genomes
Accession: PRJNA559813
Nanopore + PacBio +Illumina!!

SRR12187555 ONT 
SRR12187556 ONT 
SRR12187557 ONT 
SRR12187558 ONT 
SRR12187559 ONT 
SRR12187560 ONT 
SRR12187563 ONT 
SRR12187574 ONT 
SRR12187576 ONT 
SRR12187577 ONT 
SRR12187578 ONT 
SRR12187579 ONT 
SRR12187580 ONT 
SRR12187581 ONT 
SRR12187582 ONT 
SRR12187583 ONT 
SRR12187584 ONT 
SRR12187585 ONT 
SRR12187586 ONT 
SRR9951089 ONT 
SRR9951091 ONT 
SRR9951094 ONT 
SRR9951095 ONT 
SRR9951096 ONT 
SRR9951098 ONT 
SRR9951099 ONT 
SRR9951101 ONT 
SRR9951104 ONT 
SRR9951107 ONT 
SRR9951108 ONT 
SRR9951093 PacBio RS II
SRR9951102 PacBio RS II
