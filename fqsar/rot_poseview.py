from schrodinger import structure
from schrodinger.application import ska
import subprocess as sp
import os
import numpy as np
import csv
import sys
import time

#reference = 'test_reference.mae'
#rote file
reference = sys.argv[1]
poseview = sys.argv[2]
#poseview = 'GIPR80nsnhmem2-20-Jan9_pv.maegz'

#cmd = ['$SCHRODINGER/utilities/structalign', '-asl \'protein AND (chain.name C)\'', reference, poseview, '-matrix' ]

#print(" ".join(cmd))
#p = sp.Popen(' '.join(cmd), stdout= sp.PIPE, shell = True)
#out,err = p.communicate()
#print(out,err)
poseview_prefix = os.path.splitext(os.path.basename(poseview))[0]
#poseview_prefix = os.path.splitext(os.path.basename(poseview))[0]
#pathpose = os.path.splitext(poseview)[0]
#wait for write
#print('prefix',poseview_prefix)
with open (reference) as trn_open:
#with open(pathpose+'.rot') as trn_open:
	next(trn_open)
#	trn_matrix = trn_open.read().splitlines()

	trn_reader = csv.reader(trn_open, delimiter=' ', quoting=csv.QUOTE_NONNUMERIC)
	trn_matrix = [arr for arr in [i for i in trn_reader]]
matrix = (trn_matrix[:3])
x = (trn_matrix[3:][0])


#matrix =[[0.897, 0.073, -0.435],[-0.440, 0.227, -0.869],[0.035, 0.971, 0.235]]
#x = [-18.214, 6.143, -2.991]

with structure.StructureReader(poseview) as reader:
#	recep = next(reader)
#	for st in reader:
#		transform.translate_structure(st, arr)
#		transform.rotate_structure(st,arr)
	with structure.StructureWriter('splits/'+str(poseview_prefix)+'_ligand.maegz') as writer:

#	with structure.StructureWriter(poseview_prefix+'_align.pdb') as writer:
		for st in reader:
			ska.transform_structure(st, matrix, x)
#			trns_st = transform.transform_structure(st, matrix)
			writer.append(st)
