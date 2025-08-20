#!/bin/bash

######邏輯需更改！！Demography用Eu加Am，然後，SV也應該只用這些族群的！

fastDFE="/dfs7/jje/jenyuw/SV-project-temp/result/fastDFE"
processed_SNP="/dfs7/jje/jenyuw/SV-project-temp/result/processed_SNP"
fit_dadi="/dfs7/jje/jenyuw/SV-project-temp/result/fit_dadi"
nT=$SLURM_CPUS_PER_TASK
#transfer file from "fit-dadi"
#cp EU.syn-snp.vcf ../fastDFE/
#cp EU.nonsyn-snp.vcf ../fastDFE/
#cp good.vcf $fastDFE

bcftools view --threads ${nT} -m2 -M2 -v snps -S ${fit_dadi}/EU_snp_popfile.TXT ${processed_SNP}/three_prime_UTRSNPs.vcf.gz |\
sed 's@1\/1@1\|1@g;s@.\/.@0\|0@g;s@0\/1@0\|1@g' >${fastDFE}/EU.3pri-snp.vcf

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  A1_CLR_mumco    A2_CLR_mumco ......
#fastDFE does not accept SV as input so we have to edit the ref and alt alleles namually. 
grep "#" ${fastDFE}/good.vcf >${fastDFE}/modified.SV.vcf
grep -v "#" ${fastDFE}/good.vcf|gawk '{s=""; for (i=6; i<= NF; i++) s = s $i "\t"; print $1 "\t" $2 "\t" $3 "\t" "C" "\t" "A" "\t" s}' >>${fastDFE}/modified.SV.vcf

grep -v -e "##reference" -e "##query0" ${fastDFE}/modified.SV.vcf|bcftools view --threads ${nT} -i 'SVTYPE="DUP"'  >${fastDFE}/modified.DUP_only.vcf
grep -v -e "##reference" -e "##query0" ${fastDFE}/modified.SV.vcf|bcftools view --threads ${nT} -i 'SVTYPE="DEL"'  >${fastDFE}/modified.DEL_only.vcf
grep -v -e "##reference" -e "##query0" ${fastDFE}/modified.SV.vcf|bcftools view --threads ${nT} -i 'SVTYPE="INS"'  >${fastDFE}/modified.INS_only.vcf


###Fix fucking unexpected fidplay problem###
xhost +

module load anaconda/2022.05
#conda create -n "fastdfe" python=3.12 pip
#pip install fastdfe
conda activate fastdfe
cd ${fastDFE}

python3
#python codes below. INTERACTIVE mode
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

##Let's work on the subsets of SVs.

p_down = fd.Parser(
    vcf="EU.syn-snp.vcf",
    n=33,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)
sfs_down = p_down.parse()
sfs_down=sfs_down.fold()
#sfs_down.plot()

dup= fd.Parser(
    vcf="modified.DUP_only.vcf",
    n=33,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)
sfs_dup=dup.parse()
sfs_dup.plot()

inf_dup = fd.BaseInference(
    sfs_neut=sfs_down,
    sfs_sel=sfs_dup,
    n_runs=20,
    folded=True,
    bounds = {'b':(0.01, 100),'S_b':(0.01, 1000),'eps':(0, 1)}
)
inf_dup.run();
inf_dup.bootstrap(n_samples=30)
inf_dup.plot_discretized()
inf_dup.plot_discretized(title='DFE of duplications',
intervals=[float('-inf'),-128,-64,-32,-16,-8,-6,-4,-2,-1,0,1,2,4,8,16,32,float('inf')])
inf_dup.plot_sfs_comparison()


ins= fd.Parser(
    vcf="modified.INS_only.vcf",
    n=33,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)
sfs_ins=ins.parse()
sfs_ins.plot()
inf_ins = fd.BaseInference(
    sfs_neut=sfs_down,
    sfs_sel=sfs_ins,
    n_runs=20,
    folded=True,
    bounds = {'b':(0.01, 100),'S_b':(0.01, 1000),'eps':(0, 2)}
)
inf_ins.run();
inf_ins.bootstrap(n_samples=30)
#inf_ins.plot_discretized()
inf_ins.plot_discretized(title='DFE of insertions',
intervals=[float('-inf'),-128,-64,-32,-16,-8,-6,-4,-2,-1,0,1,2,4,8,16,32,float('inf')])
inf_ins.plot_sfs_comparison()


deletion= fd.Parser(
    vcf="modified.DEL_only.vcf",
    n=33,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)
sfs_deletion=deletion.parse()
sfs_deletion.plot()
inf_deletion = fd.BaseInference(
    sfs_neut=sfs_down,
    sfs_sel=sfs_deletion,
    n_runs=20,
    folded=True,
    bounds = {'b':(0.01, 100),'S_b':(0.01, 1000),'eps':(0, 2)}
)
inf_deletion.run();
inf_deletion.bootstrap(n_samples=30)
inf_deletion.plot_discretized()
inf_deletion.plot_discretized(title='DFE of deletion',
intervals=[float('-inf'),-128,-64,-32,-16,-8,-6,-4,-2,-1,0,1,2,4,8,16,32,float('inf')])
inf_deletion.plot_sfs_comparison()


##Let's work on the SVs associated with TEs
fastDFE="/dfs7/jje/jenyuw/SV-project-temp/result/fastDFE"
TE="/dfs7/jje/jenyuw/SV-project-temp/result/TE_repeat"
nT=$SLURM_CPUS_PER_TASK

cp ${TE}/SV_TE-annotated.vcf ${fastDFE}/

#make the VCF acceptable for fastDFE
grep "#" ${fastDFE}/SV_TE-annotated.vcf| grep -v -e "##reference" -e "##query0" >${fastDFE}/modified.SV_TE-annotated.vcf
grep -v "#" ${fastDFE}/SV_TE-annotated.vcf|\
gawk '{s=""; for (i=6; i<= NF; i++) s = s $i "\t"; print $1 "\t" $2 "\t" $3 "\t" "C" "\t" "A" "\t" s}' |\
sed 's@/@|@g'|sed 's/\.|\./0|0/g' >>${fastDFE}/modified.SV_TE-annotated.vcf

bcftools annotate -x 'FORMAT' ${fastDFE}/modified.SV_TE-annotated.vcf |\
tee >(bcftools view -i 'ORDER="LINE"' >${fastDFE}/modified.LINE.vcf) \
tee >(bcftools view -i 'ORDER="RC"' >${fastDFE}/modified.RC.vcf) \
tee >(bcftools view -i 'ORDER="Simple_repeat"' >${fastDFE}/modified.Simple_repeat.vcf) \
tee >(bcftools view -i 'ORDER="Satellite"' >${fastDFE}/modified.Satellite.vcf) \
>(bcftools view -i 'ORDER="DNA"' >${fastDFE}/modified.DNA.vcf) |\
bcftools view -i 'ORDER="LTR"' >${fastDFE}/modified.LTR.vcf


module load anaconda/2022.05
conda activate fastdfe
fastDFE="/dfs7/jje/jenyuw/SV-project-temp/result/fastDFE"
cd ${fastDFE}

python3
###########python codes below. INTERACTIVE mode########
import fastdfe as fd

p_down = fd.Parser(
    vcf="EU.syn-snp.vcf",
    n=33,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)
sfs_down = p_down.parse()
sfs_down=sfs_down.fold()
sfs_down.plot()

line= fd.Parser(
    vcf="modified.LINE.vcf",
    n=33,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)
sfs_line=line.parse()
sfs_line.plot()
inf_line = fd.BaseInference(
    sfs_neut=sfs_down,
    sfs_sel=sfs_line,
    n_runs=20,
    folded=True,
    bounds = {'b':(0.01, 100),'S_b':(0.1, 1000),'eps':(0, 1)}
)
inf_line.run();
inf_line.bootstrap(n_samples=30)
#inf_line.plot_discretized()
inf_line.plot_discretized(title='DFE of LINE',
intervals=[float('-inf'),-128,-64,-32,-16,-8,-6,-4,-2,-1,0,1,2,4,8,16,32,float('inf')])
inf_line.plot_sfs_comparison()



LTR= fd.Parser(
    vcf="modified.LTR.vcf",
    n=33,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)
sfs_LTR=LTR.parse()
sfs_LTR.plot()
inf_LTR = fd.BaseInference(
    sfs_neut=sfs_down,
    sfs_sel=sfs_LTR,
    n_runs=20,
    folded=True,
    bounds = {'b':(0.01, 100),'S_b':(0.1, 1000),'eps':(0, 1)}
)
inf_LTR.run();
inf_LTR.bootstrap(n_samples=30)
#inf_LTR.plot_discretized()
inf_LTR.plot_discretized(title='DFE of LTR',
intervals=[float('-inf'),-128,-64,-32,-16,-8,-6,-4,-2,-1,0,1,2,4,8,16,32,float('inf')])
inf_LTR.plot_sfs_comparison()



RC= fd.Parser(
    vcf="modified.RC.vcf",
    n=33,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)
sfs_RC=RC.parse()
sfs_RC.plot()

inf_RC = fd.BaseInference(
    sfs_neut=sfs_down,
    sfs_sel=sfs_RC,
    n_runs=20,
    folded=True,
    bounds = {'b':(0.01, 10),'S_b':(0.1, 1000),'eps':(0, 1)}
)
inf_RC.run();
inf_RC.bootstrap(n_samples=30)
#inf_RC.plot_discretized()
inf_RC.plot_discretized(title='DFE of RC',
intervals=[float('-inf'),-128,-64,-32,-16,-8,-6,-4,-2,-1,0,1,2,4,8,16,32,float('inf')])
inf_RC.plot_sfs_comparison()



DNA= fd.Parser(
    vcf="modified.DNA.vcf",
    n=33,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)
sfs_DNA=DNA.parse()
sfs_DNA.plot()

inf_DNA = fd.BaseInference(
    sfs_neut=sfs_down,
    sfs_sel=sfs_DNA,
    n_runs=20,
    folded=True,
    bounds = {'b':(0.01, 100),'S_b':(0.1, 1000),'eps':(0, 1)}
)
inf_DNA.run();
inf_DNA.bootstrap(n_samples=30)
#inf_DNA.plot_discretized()
inf_DNA.plot_discretized(title='DFE of DNA',
intervals=[float('-inf'),-128,-64,-32,-16,-8,-6,-4,-2,-1,0,1,2,4,8,16,32,float('inf')])
inf_DNA.plot_sfs_comparison()


Simple_repeat= fd.Parser(
    vcf="modified.Simple_repeat.vcf",
    n=33,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)
sfs_Simple_repeat=Simple_repeat.parse()
sfs_Simple_repeat.plot()

inf_Simple_repeat = fd.BaseInference(
    sfs_neut=sfs_down,
    sfs_sel=sfs_Simple_repeat,
    n_runs=20,
    folded=True,
    #bounds = {'b':(0.01, 100),'S_b':(0.1, 1000),'eps':(0, 1)}
)
inf_Simple_repeat.run();
inf_Simple_repeat.bootstrap(n_samples=30)
#inf_Simple_repeat.plot_discretized()
inf_Simple_repeat.plot_discretized(title='DFE of Simple_repeat',
intervals=[float('-inf'),-128,-64,-32,-16,-8,-6,-4,-2,-1,0,1,2,4,8,16,32,float('inf')])
inf_Simple_repeat.plot_sfs_comparison()



Satellite= fd.Parser(
    vcf="modified.Satellite.vcf",
    n=33,  # SFS sample size
    skip_non_polarized = False,
    subsample_mode ='probabilistic'
)
sfs_Satellite=Satellite.parse()
sfs_Satellite.plot()

inf_Satellite = fd.BaseInference(
    sfs_neut=sfs_down,
    sfs_sel=sfs_Satellite,
    n_runs=20,
    folded=True,
    bounds = {'b':(0.01, 100),'S_b':(0.1, 1000),'eps':(0, 1)}
)
inf_Satellite.run();
inf_Satellite.bootstrap(n_samples=30)
#inf_Satellite.plot_discretized()
inf_Satellite.plot_discretized(title='DFE of Satellite',
intervals=[float('-inf'),-128,-64,-32,-16,-8,-6,-4,-2,-1,0,1,2,4,8,16,32,float('inf')])
inf_Satellite.plot_sfs_comparison()