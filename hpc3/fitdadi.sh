#!/usr/bin/python3

#BASH preperation
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

./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 -a --preview --unfolded
./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 -a --ploidy 1 --preview --unfolded


./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 --proj 5,21,3,33 -a --ploidy 1 --unfolded -o good_sfs1
#./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 --proj 5,21,3,33 --ploidy 1 --unfolded -o good_sfs2
./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 --proj 5,21,3,33 -a --unfolded -o good_sfs3
#./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 --proj 5,21,3,33 --unfolded -o good_sfs4
./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 --proj 10,42,6,66 -a --unfolded -o good_sfs5
#./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 --proj 2,2,2,2 -a --ploidy 1 --unfolded -o good_sfs6
#./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 --proj 2,2,2,2 -a  --unfolded -o good_sfs7
./easySFS/easySFS.py  -i good.vcf -p popfile.txt2 -a --unfolded -o good_sfs8
python3
#Real work in Python interactive mode
import dadi
import numpy
import Selection

#/dfs7/jje/jenyuw/SV-project-temp/result/fit_dadi/fitdadi_original/manual_examples/example.py
#example.sfs
#nano /data/homezvol2/jenyuw/.local/lib/python3.10/site-packages/numpy/core/function_base.py

#bcftools query -l polarized.asm.vcf.gz

dd = dadi.Misc.make_data_dict_vcf("good.vcf", "popfile.txt2")
fs = dadi.Spectrum.from_data_dict(dd, ['AF', 'AM', 'AS', 'EU'], projections = [5, 21, 3, 33], polarized = True)