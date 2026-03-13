#covpock12_c_gcmc_pocket_analysis/pocket_conformations
ls -1 /lrlhps/users/c302900/ucn2/bfact/mxmd/apoucn_gcmc_pocket_analysis/pocket_conformations/*.maegz > confs.in
cwd=$(pwd)
while read F; do

Fsub=$(basename ${F} .maegz)
#echo ${Fsub}
mkdir ${Fsub}
cd ${Fsub}
#$SCHRODINGER/run split_structure.py -k -m pdb -many_files ${F} split.maegz
$SCHRODINGER/run split_structure.py -k -m molecule -many_files ${F} split.maegz

cp -p ../sitemap.in .
$SCHRODINGER/utilities/prepwizard -nometaltreat -disulfides -f S-OPLS -nopropka -noepik -WAIT -noprotassign -noimpref split_mol2.maegz -HOST localhost split_receptor1-prep.mae 
"${SCHRODINGER}/sitemap" sitemap.in -HOST sge_cpu -TMPLAUNCHDIR -j ${Fsub}

cd ${cwd}

#mv sitemap.log  ${Fsub}.log

done<confs.in


grep -A 1 Dscore apoucn_gcmc_*/apoucn*.log > Dscoreall.out
