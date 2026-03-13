Preparation
1.) copy probes with activity metadata (out of mxmd analysis) into a probes/ folder
for example see copy_example.sh
2.) Align to model (docking /ref etc)
Copy the reference protein and one protein frame from the mxmd analysis (out of the pocket_analysis folder)
Edit refdock and mxmd entries in run_3dfield_prep.sh

Adjust the asl selection for alignment in run_3dfield_prep.sh
$SCHRODINGER/utilities/structalign -asl 'protein AND (chain.name R)' $refdock $mxmd -matrix

3) Run bash run_3dfield_prep.sh

this will generate all probes in one file for use in 3d-field qsar:
cat_probe_top_scaled_aligneddock.maegz
activity column will be: r_cluster_activity_linscale (this has been linearly scaled to give a large range of activities).

4) Run fqsar:
 bash runfqsargauss.sh
