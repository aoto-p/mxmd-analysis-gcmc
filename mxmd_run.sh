module load schrodinger/2024.3
host=sge_gpu_highpriority
#edit mxmd_npt.msj for restraints and mxmd analysis alignment asl
$SCHRODINGER/utilities/multisim 242_21_47_prime-outWT_apo.maegz -HOST sge_cpu:64 -SUBHOST $host -OPLSDIR custom_2024_3.opls -JOBNAME mxmd -maxjob 64 -m mxmd.msj
