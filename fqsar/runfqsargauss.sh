pre='dockedtopscalegauss'
phase_fqsar cat_probe_top_scaled_aligneddock.maegz ${pre}_fqsar_pred.mae "r_cluster_activity_linscale" -build -style gauss,gauss_r,gauss_ext -factors 6 -extend 3.0 -grid 1.0 -buff 2.0 -scut 30.0 -ecut 30.0 -sd 0.01 -pt 0.95 -osum ${pre}summary -omod ${pre}.fqsar

#allframe2-8_fqsar_pred.mae

