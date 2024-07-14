#!/bin/bash

fastDFE="/dfs7/jje/jenyuw/SV-project-temp/result/fastDFE"
processed_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/processed_SNP"
fit_dadi="/dfs7/jje/jenyuw/SV-project-temp/result/fit_dadi"
#transfer file from "fit-dadi"
#cp EU.syn-snp.vcf ../fastDFE/
#cp EU.nonsyn-snp.vcf ../fastDFE/
#cp good.vcf $fastDFE

bcftools view --threads ${nT} -m2 -M2 -v snps -S ${fit_dadi}/EU_snp_popfile.TXT ${processed_SNP}/three_prime_UTRSNPs.vcf.gz |\
sed 's@1\/1@1\|1@g;s@.\/.@0\|0@g;s@0\/1@0\|1@g' >${fastDFE}/EU.3pri-snp.vcf

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  A1_CLR_mumco    A2_CLR_mumco ......
#fastDFE does not accept SV as input so we have to deit the ref and alt alleles namually. 
grep "#" ${fastDFE}/good.vcf >${fastDFE}/modified.SV.vcf
grep -v "#" ${fastDFE}/good.vcf|gawk '{s=""; for (i=6; i<= NF; i++) s = s $i "\t"; print $1 "\t" $2 "\t" $3 "\t" "C" "\t" "A" "\t" s}' >>${fastDFE}/modified.SV.vcf


module load anaconda/2022.05
#conda create -n "fastdfe" python=3.12 pip
#pip install fastdfe
conda activate fastdfe
cd ${fastDFE}
python3
#python codes below
import fastdfe as fd
p = fd.Parser(
    vcf="EU.syn-snp.vcf",
    n=66,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)
q= fd.Parser(
    vcf="EU.nonsyn-snp.vcf",
    n=66,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)
r= fd.Parser(
    vcf="modified.SV.vcf",
    n=33,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)
p_3pri=fd.Parser(
    vcf="EU.3pri-snp.vcf",
    n=66,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)

sfs = p.parse()
sfs2=q.parse()
sfs3=r.parse()
sfs_3pri=p_3pri.parse()

sfs.plot()
sfs2.plot()
sfs3.plot()
sfs_3pri.plot()

inf = fd.BaseInference(
    sfs_neut=sfs,
    sfs_sel=sfs2,
    n_runs=20,
#    n_bootstraps=10,
#    do_bootstrap=True,
    bounds = {'p_b': (0,5)},
    folded=True,
)
# run inference
inf.run();
inf.bootstrap(n_samples=30)
# plot discretized DFE
inf.plot_discretized()
inf.plot_discretized(title='DFE of Nonsynonymous SNPs',
intervals=[float('-inf'),-128,-64,-32,-16,-8,-6,-4,-2,-1,0,1,2,4,8,16,32,float('inf')])
# plot SFS comparison
inf.plot_sfs_comparison()


##Use small introns as a not-so neutral control
inf_3pri = fd.BaseInference(
    sfs_neut=sfs,
    sfs_sel=sfs_3pri,
    n_runs=20,
    bounds = {'p_b': (0,5), 'b':(0.01,20)},
    folded=True,
)
inf_3pri.run();
inf_3pri.plot_discretized(title="DFE of SNPs in 3\' UTRs",
intervals=[float('-inf'),-128,-64,-32,-16,-8,-6,-4,-2,-1,0,1,2,4,8,16,32,float('inf')])

##Let's work on the SVs. Remember to change the sample size. 33 haplotypes of SVs but 66 haplotypes of SNPs.
p_down = fd.Parser(
    vcf="EU.syn-snp.vcf",
    n=33,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)
sfs_down = p_down.parse()
sfs_down=sfs_down.fold()
inf2 = fd.BaseInference(
    sfs_neut=sfs_down,
    sfs_sel=sfs3,
    n_runs=20,
    folded=True,
    bounds = {'b':(0.01, 100),'eps':(0, 1)}
)

inf2.run();
inf2.bootstrap(n_samples=30)
inf2.plot_discretized()
inf2.plot_discretized(title='DFE of all SVs',
intervals=[float('-inf'),-128,-64,-32,-16,-8,-6,-4,-2,-1,0,1,2,4,8,16,32,float('inf')])
inf2.plot_sfs_comparison()