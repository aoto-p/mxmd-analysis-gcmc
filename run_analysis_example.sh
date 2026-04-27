module load schrodinger/2024.3
#copy the final trajectory tgz file for one replica in mxmd_3 into a separate folder and untar it
#for example: tar -xvf mxmd/mxmd_3/mxmd_mxmd_sulfamide-5_7-out.tgz
#open coordinates and find the pocket center (appoximately): for example -4.13, -2.25, 4.25 and use this as the --center values below
#adjust size in Angstroms of the box around the center to include for co-solvent probes
#edite probe names and PDB names -order between probes and solventsPDB must match

$SCHRODINGER/run pocket_analysis_gcmc_ligandpop.py --pname short_covpock9_gcmc --jname mxmd --size 9 --center -4.13 -2.25 4.25 --cluster --simulationStage 7 --probes imidazole pyrimidine acetonitrile isopropanol methylamine sulfamide --solventsPDB IMD P1R CCN IPA NME FUS
