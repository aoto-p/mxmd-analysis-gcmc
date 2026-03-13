#Check that structcat is pulling out the protein only - and that the maegz prepwizard is running on is the correct maegz (not ligand etc).
#create a refligand.maegz that superposes onto the pocket -this is used for the 3D location of where to run the sitemap analysis
#Prepare Dscoreall (output of bash script grep)- prior to running dscore_anal.py
#remove -- rows in vi Dscoreall with:
#:g/--/d
#change first row to include name column and add a space at the colon:
#for example:
#covpock12_c_gcmc_acetate_frm102_cluster_18_pop_39/covpock12_c_gcmc_acetate_frm102_cluster_18_pop_39.log:SiteScore size   Dscore  volume  exposure enclosure contact  phobic   philic   balance  don/acc  refdist  refmin   refavg  sitemin
#change to:
#Name SiteScore size   Dscore  volume  exposure enclosure contact  phobic   philic   balance  don/acc  refdist  refmin   refavg  sitemin
