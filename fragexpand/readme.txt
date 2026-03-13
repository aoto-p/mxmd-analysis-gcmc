PREPARATION:
=======
Before running site_batch.sh:

1.) Modify the first line of site_batch.sh to point to the pocket_conformation folder from mxmd analysis:
this line:
ls -1 /lrlhps/users/c302900/ucn2/bfact/mxmd/apoucn_gcmc_pocket_analysis/pocket_conformations/*.maegz > confs.in

2.) Check that structcat is pulling out the protein only - and that the maegz prepwizard is running on is the correct maegz (not ligand etc).
for example- modify which file is being run 'split_mol2.maegz' to be the receptor output of strutcat.
$SCHRODINGER/utilities/prepwizard -nometaltreat -disulfides -f S-OPLS -nopropka -noepik -WAIT -noprotassign -noimpref split_mol2.maegz -HOST localhost split_receptor1-prep.mae

3.) Create a ligandref.maegz -a ligand that superposes onto the pocket -the COM is used for the 3D location of where to run the sitemap analysis
Place ligandref.maegz in the same path as site_batch.sh

======
Before running dscore_anal.py

Prepare Dscoreall (output of bash script grep)- prior to running dscore_anal.py
1.) remove -- rows in vi Dscoreall with:
:g/--/d

2.) change only the first row to include name column and add a space at the colon (other rows will be ignored):
for example:
covpock12_c_gcmc_acetate_frm102_cluster_18_pop_39/covpock12_c_gcmc_acetate_frm102_cluster_18_pop_39.log:SiteScore size   Dscore  volume  exposure enclosure contact  phobic   philic   balance  don/acc  refdist  refmin   refavg  sitemin
change to:
Name SiteScore size   Dscore  volume  exposure enclosure contact  phobic   philic   balance  don/acc  refdist  refmin   refavg  sitemin
