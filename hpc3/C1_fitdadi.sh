#!/bin/bash
#BASH preperation
Merged_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/Merged_SNP"
processed_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/processed_SNP"
fit_dadi="/dfs7/jje/jenyuw/SV-project-temp/result/fit_dadi"
nT=$SLURM_CPUS_PER_TASK

cd /dfs7/jje/jenyuw/SV-project-temp/result/fit_dadi
module load python/3.10.2

python3 -m pip install dadi --user ##The only successfull way to install dadi for now
#python3 -m pip install dadi --prefix "/dfs7/jje/jenyuw/SV-project-temp/result/fit_dadi/bin"
git clone https://bitbucket.org/gutenkunstlab/dadi.git

#Don't install from source, because it didn't work
#git clone https://github.com/LohmuellerLab/fitdadi.git
#cp /dfs7/jje/jenyuw/SV-project-temp/result/fit_dadi/fitdadi_original/Selection.py /dfs7/jje/jenyuw/SV-project-temp/result/fit_dadi/dadi
#cp /dfs7/jje/jenyuw/SV-project-temp/result/fit_dadi/fitdadi_original/__init__.py /dfs7/jje/jenyuw/SV-project-temp/result/fit_dadi/dadi
#cd /dfs7/jje/jenyuw/SV-project-temp/result/fit_dadi/dadi
#python3 setup.py install --user

bcftools annotate -x FORMAT/SUPP polarized.asm.vcf.gz -O v |\
sed 's@1\/1@1\|1@g;s@.\/.@0\|0@g' >good.vcf

bcftools view --threads ${nT} -m2 -M2 -v snps -S EU_snp_popfile.TXT ${processed_SNP}/synSNPs.vcf.gz |\
sed 's@1\/1@1\|1@g;s@.\/.@0\|0@g;s@0\/1@0\|1@g' >EU.syn-snp.vcf
python3 /pub/jenyuw/Software/easySFS/easySFS.py -i EU.syn-snp.vcf -p EU_snp_popfile.tsv -a --preview
python3 /pub/jenyuw/Software/easySFS/easySFS.py -i EU.syn-snp.vcf -p EU_snp_popfile.tsv -a -f --proj 66 -o EU_sfs


bcftools view --threads ${nT} -m2 -M2 -v snps -S AF-EU_snp_popfile.TXT ${processed_SNP}/synSNPs.vcf.gz |\
sed 's@1\/1@1\|1@g;s@.\/.@0\|0@g;s@0\/1@0\|1@g' >AF-EU.syn-snp.vcf
python3 /pub/jenyuw/Software/easySFS/easySFS.py -i AF-EU.syn-snp.vcf -p AF-EU_snp_popfile.tsv -a --preview
python3 /pub/jenyuw/Software/easySFS/easySFS.py -i AF-EU.syn-snp.vcf -p AF-EU_snp_popfile.tsv -a -f --proj "66,10" -o AF-EU_sfs

./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 --proj 5,21,3,33 -a --ploidy 1 --unfolded -o good_sfs1
#./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 --proj 5,21,3,33 --ploidy 1 --unfolded -o good_sfs2
./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 --proj 5,21,3,33 -a --unfolded -o good_sfs3
#./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 --proj 5,21,3,33 --unfolded -o good_sfs4
./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 --proj 10,42,6,66 -a --unfolded -o good_sfs5
#./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 --proj 2,2,2,2 -a --ploidy 1 --unfolded -o good_sfs6
#./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 --proj 2,2,2,2 -a  --unfolded -o good_sfs7
./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 -a --unfolded -o good_sfs8

## Let's only work on the Europe population first. 

###!/usr/bin/python3
python3
#Real work in Python interactive mode
import dadi
import numpy
#import Selection
fs = dadi.Spectrum.from_file("EU-66.sfs")
thetaW = fs.Watterson_theta()
pi = fs.pi()
D = fs.Tajima_D()
#folded = fs.fold() #imput spectrum is already folded
fs.mask[0] = True
fs.mask[66] = True
sample = fs.sample()
quit()

module unload python/3.10.2


module load anaconda/2022.05
conda activate dadi-cli
dadi-cli Model --names #this displays all the models
dadi-cli Model --names bottlegrowth_1d #this displays the parameters #params = (nuB,nuF,T)
dadi-cli Model --names growth #params = (nu,T)
dadi-cli Model --names two_epoch #params = (nu,T)
dadi-cli Model --names three_epoch_inbreeding #(nuB,nuF,TB,TF,F)
dadi-cli Model --names three_epoch #(nuB,nuF,TB,TF)

#snm_1d (Standard neutral model, populations never diverge.)
#three_epoch
#three_epoch_inbreeding
cd /dfs7/jje/jenyuw/SV-project-temp/result/fit_dadi/EU_sfs/dadi
dadi-cli InferDM --fs "EU-66.sfs" --model growth --lbounds 0.00001 0  --ubounds 100000 10000  --output gr --optimizations 300
dadi-cli BestFit --input-prefix "gr.InferDM" --lbounds 0.00001 0  --ubounds 100000 10000

dadi-cli InferDM --fs "EU-66.sfs" --model two_epoch --lbounds 0.001 0.001  --ubounds 100 100  --output 2epo --optimizations 300
dadi-cli BestFit --input-prefix "2epo.InferDM" --lbounds 0.001 0.001  --ubounds 100 100

dadi-cli InferDM --fs "EU-66.sfs" --model bottlegrowth_1d --lbounds 0.01 0.01 0.01  --ubounds 100 100 100  --output bg1d --optimizations 100
dadi-cli InferDM --fs "EU-66.sfs" --model bottlegrowth_1d --lbounds 0.01 0.01 0.0001  --ubounds 100 100 100  --output bg1d --optimizations 1000
dadi-cli InferDM --fs "EU-66.sfs" --model bottlegrowth_1d --lbounds 0.000001 0.000001 0.000001  --ubounds 100000 100000 100000  --output bg1d --optimizations 1000
dadi-cli BestFit --input-prefix "bg1d.InferDM" --lbounds 0 0 0  --ubounds 100000 100000 100000

dadi-cli InferDM --fs "EU-66.sfs" --model three_epoch_inbreeding --lbounds 0.1 0.1 0.1 0.1 0.0001 --ubounds 100 100 100 100 0.99 --output 3epoin --optimizations 300
dadi-cli BestFit --input-prefix "3epoin.InferDM" --lbounds 0.1 0.1 0.1 0.1 0.0001 --ubounds 100 100 100 100 0.99

dadi-cli InferDM --fs "EU-66.sfs" --model three_epoch --lbounds 0.1 0.1 0.1 0.01  --ubounds 1000 100 1000 100  --output 3epo --optimizations 300
dadi-cli BestFit --input-prefix "3epo.InferDM" --lbounds 0.1 0.1 0.1 0.01  --ubounds 1000 100 1000 100

#no good convergence in European population with 1D models
#so let's try AF-EU population

cd /dfs7/jje/jenyuw/SV-project-temp/result/fit_dadi/AF-EU_sfs/dadi
dadi-cli Model --names bottlegrowth_split_mig #(nuB,nuF,m,T,Ts)
dadi-cli Model --names split_asym_mig # (nu1,nu2,T,m12,m21)

dadi-cli InferDM --fs EU-AF.sfs --model bottlegrowth_split_mig --lbounds 0.001 0.001 0.001 0.001 0.001 --ubounds 100 100 100 100 100 --output bgsm --optimizations 300
dadi-cli BestFit --input-prefix "bgsm.InferDM" --lbounds 0.001 0.001 0.001 0.001 0.001 --ubounds 100 100 100 100 100
# --> No convergence

dadi-cli InferDM --fs EU-AF.sfs --model split_asym_mig --lbounds 0.01 0.01 0.01 0.01 0.01 --ubounds 100 100 100 100 100 --output sam --optimizations 100
dadi-cli InferDM --fs EU-AF.sfs --model split_asym_mig --lbounds 0.1 0.1 0.1 0.1 0.1 --ubounds 1000 1000 1000 1000 1000 --output sam --optimizations 1200
dadi-cli BestFit --input-prefix "sam.InferDM" --lbounds 0 0 0 0 0  --ubounds 1000 1000 1000 1000 1000 # --> No convergence