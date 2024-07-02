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
python3 /pub/jenyuw/Software/easySFS/easySFS.py -i EU.syn-snp.vcf -p EU_snp_popfile.tsv -a -f --proj 33 -o down_EU_sfs

bcftools view --threads ${nT} -m2 -M2 -v snps -S AF-EU_snp_popfile.TXT ${processed_SNP}/synSNPs.vcf.gz |\
sed 's@1\/1@1\|1@g;s@.\/.@0\|0@g;s@0\/1@0\|1@g' >AF-EU.syn-snp.vcf
python3 /pub/jenyuw/Software/easySFS/easySFS.py -i AF-EU.syn-snp.vcf -p AF-EU_snp_popfile.tsv -a --preview
python3 /pub/jenyuw/Software/easySFS/easySFS.py -i AF-EU.syn-snp.vcf -p AF-EU_snp_popfile.tsv -a -f --proj "66,10" -o AF-EU_sfs

python3 /pub/jenyuw/Software/easySFS/easySFS.py  -i good.vcf -p popfile.txt --proj 5,21,3,33 -a --ploidy 1 --unfolded -o good_sfs1

python3 /pub/jenyuw/Software/easySFS/easySFS.py  -i good.vcf -p EU_sv_popfile.tsv --proj 33 -a --ploidy 1 --unfolded -o SV_EU_sfs


## Let's only work on the Europe population first. 

# ###!/usr/bin/python3
# python3
# #Real work in Python interactive mode
# import pickle, random
# import numpy as np
# import nlopt
# import dadi
# import dadi.DFE as DFE
# import matplotlib.pyplot as plt
# from dadi.DFE import *
# #import Selection
# fs = dadi.Spectrum.from_file("EU-66.sfs")
# thetaW = fs.Watterson_theta()
# pi = fs.pi()
# D = fs.Tajima_D()
# #folded = fs.fold() #imput spectrum is already folded
# fs.mask[0] = True
# fs.mask[66] = True
# sample = fs.sample()
# quit()


# module unload python/3.10.2

# dataset = 'EU-66'
# data_fs = dadi.Spectrum.from_file('EU-66.sfs')
# demog_params = [0.08581,0.00738072]
# ns = data_fs.sample_sizes
# theta0 = 46237
# #ns = [66]
# theta_ns = theta0 * 2.31

# spectra = DFE.Cache1D(demog_params, ns, DFE.DemogSelModels.growth_sel, 
# gamma_bounds=(1e-5, 1000), gamma_pts=100, verbose=True, pts=[10])
# pickle.dump(spectra, open('gr.bpkl','wb'))

# data = dadi.Spectrum.from_file('EU-66.sfs')
# cache1d = pickle.load(open('gr.bpkl','rb'))
# cache1d = DFE.Cache1D(demog_params, ns, DFE.DemogSelModels.growth_sel, 
# gamma_bounds=(1e-5, 1000), gamma_pts=100, verbose=True, pts=[10])
# dfe_func = cache1d.integrate
# sele_dist1d = DFE.PDFs.gamma
# func_args = [sele_dist1d, theta_ns]
# params = [0.1, 15000, 0.01]
# lower_bounds = [1e-2, 1e-2, 1e-3]
# upper_bounds = [10, 10000, 1]
# p0 = dadi.Misc.perturb_params(params, fold=0, upper_bound=upper_bounds,
#                               lower_bound=lower_bounds)
# popt = dadi.Inference.opt(p0, data, cache1d.integrate, pts=None,
#                                     func_args=func_args, 
#                                     lower_bound=lower_bounds, 
#                                     upper_bound=upper_bounds,
#                                     maxeval=400, multinom=False, verbose=100)

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
dadi-cli InferDM --fs "EU-66.sfs" --model growth --lbounds 0.00001 0.00001  --ubounds 100000 10000  --output gr --optimizations 300
dadi-cli InferDM --fs "EU-66.sfs" --model growth --lbounds 0.00001 0.00001  --ubounds 10 10  --output gr --optimizations 300
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

#Let fit the growth model with EU population
#Manually edit the .bestfit file, remove the misid tag
dadi-cli GenerateCache --model growth_sel --demo-popt "gr.InferDM.bestfits" --sample-size 33 --output gr_bpkl \
--grids 1000 5000 10000 --gamma-pts 50 --gamma-bounds 0.0000001 2000

dadi-cli Pdf --names gamma #params = [alpha, beta] = [shape, scale]
input_fs="/dfs7/jje/jenyuw/SV-project-temp/result/fit_dadi/SV_EU_sfs/dadi/EU-33.sfs"
dadi-cli InferDFE --fs ${input_fs} --cache1d gr_bpkl --demo-popt "gr.InferDM.bestfits" \
--pdf1d gamma --p0 2 1 0.5 --lbounds 0.001 0.001 0.001 --ubounds 100 100 100 \
--ratio 2.31 \
--output EU-SV_gr-gamma --optimizations 500 --maxeval 400 --check-convergence 50

dadi-cli Plot --fs ${input_fs} --demo-popt EU-SV_gr-gamma.InferDFE.bestfits  --output snps.vs.SV.pdf --model growth

input_fs="/dfs7/jje/jenyuw/SV-project-temp/result/fit_dadi/down_EU_sfs/dadi/EU-33.sfs"
dadi-cli InferDFE --fs ${input_fs} --cache1d gr_bpkl --demo-popt "gr.InferDM.bestfits" \
--pdf1d gamma --p0 1 1 --lbounds 0.01 0.01 --ubounds 100 100 \
--ratio 2.31 \
--output EU-SV_gr-gamma --optimizations 50 --maxeval 400 --check-convergence 50

dadi-cli Plot --fs ${input_fs} --demo-popt gr.InferDM.bestfits  --output snps.pdf --model growth


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