#generate roation matrix to move mxmd output into docking reference frame
refdock=refboundxtal.maegz
mxmd=short_covpock9_gcmc_acetonitrile_frm152_cluster_11_pop_14.maegz
prefix=$(basename $mxmd .maegz)
$SCHRODINGER/utilities/structalign -asl 'protein AND (chain.name A)' $refdock $mxmd -matrix

#copy all probes from mxmd/pocketanalprobefolders (see copy_example.sh)
#$SCHRODINGER/run structcat.py probes/*.maegz probes_un/*.maegz -o cat_probe_top.maegz
$SCHRODINGER/run structcat.py probes/*.maegz -o cat_probe_top.maegz

$SCHRODINGER/run probe_act_scale.py

mkdir splits
$SCHRODINGER/run rot_poseview.py $prefix.rot cat_prob_top_scaled.maegz 
mv splits/cat_prob_top_scaled_ligand.maegz cat_probe_top_scaled_aligneddock.maegz

#$SCHRODINGER/run rot_poseview.py mxmd2dfe.rot cat_prob_top_scaled.maegz
#mv splits/cat_prob_top_scaled_ligand.maegz cat_probe_top_scaled_aligneddfe.maegz
