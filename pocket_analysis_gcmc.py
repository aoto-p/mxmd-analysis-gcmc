import time
import csv
import queue
import os
import sys
from schrodinger.structutils import analyze
from schrodinger.utils import cmdline, fileutils
from schrodinger import structure
from schrodinger import project
from schrodinger.application.desmond.packages import traj
from schrodinger.application.desmond.packages import traj_util
from schrodinger.application.desmond.packages import topo
import os
import sys
import fileinput
import shutil
import subprocess as sp
import pandas as pd
import numpy as np
import argparse
import threading
import math
from schrodinger.application.desmond import cmj
from schrodinger.application.desmond import cms
from schrodinger.application.desmond import constants
from schrodinger.application.desmond import envir
from schrodinger.application.desmond import struc
from schrodinger.application.desmond import util
from schrodinger.application.desmond.cns_io import write_cns
from schrodinger.application.desmond.mxmd import mxmd_system_builder as msb
from schrodinger.structutils.analyze import evaluate_asl
from schrodinger.utils import subprocess
from schrodinger import structure
from schrodinger.structutils import measure
from schrodinger.application.desmond.packages import traj, topo
from schrodinger.application.desmond.packages import traj_util
from schrodinger.trajectory.trajectory_gui_dir import export_structures
from schrodinger.structure import StructureReader, StructureWriter
'''
MxMD Pocket analysis script 
author: lily steiner (contact for questions at lily.steiner@lilly.com)
features:
- finds centroid of the most frequently occuring population of the receptor conformation with solvent present in pocket of interest
- finds interactions with receptor made by solvent in pocket of interest
output: directory with up to 3 centroids representing the top populations of pocket & a csv file "interaction

workflow:
extract the hotspot frames fully -> save a version with one solvent per frame -> cluster on that version with 1 solvent per frame, full protein cluster -> also run interaction analysis on solvents 
##modified for gcmc, and outputs 10 cluster centroids and neareset neighbors
'''
#global:/
message_queue = queue.Queue()

def extract_trjs(probe,mxmd_results_dir,jname,probe_path,trj_stage_number,i):
    print('extracting trj ' + str(i) + ' files...')
    #0-9 for the probe runs. if more or less, adjust 
    #for i in range(0,total_reps):
        #first step is to get the files for trj2mae and put them in the probe dir(the cwd). 
    jobname = jname + '_mxmd_' + probe + '-' + str(i) + '_' + str(trj_stage_number)
    folder= mxmd_results_dir + '/' + jobname +'-out.tgz'
    outcms = jname + '_mxmd_'+probe+'-'+str(i)+'_' + str(trj_stage_number)+'/'+ jobname +'-out.cms'
    trj = jname + '_mxmd_'+probe+'-'+str(i)+'_' + str(trj_stage_number)+'/'+jobname +'.xtc'
    try:
       cmd = ['tar','-xzf',folder,'-C',mxmd_results_dir + '/']
       p = sp.Popen(' '.join(cmd), stdout = sp.PIPE, shell = True)
       out,err = p.communicate()    
    except:
       print("trajectory file missing")
    return 

def extract_frames_parallel(hotspot_frames_list,xyz_center,pocket_size,interaction,i,trj_stage_number,probe,jname,mxmd_dir,pocket_dir,solvent_pdb):  
    # outcms = str(i) + '-out.cms'
    # trj = str(i) + '.xtc'
    jobname = jname + '_mxmd_' + probe + '-' + str(i) + '_' + str(trj_stage_number)
    outcms = mxmd_dir + '/' + jname + '_mxmd_'+probe+'-'+str(i)+'_' + str(trj_stage_number)+'/'+ jobname +'-out.cms'
    trj = mxmd_dir + '/' + jname + '_mxmd_'+probe+'-'+str(i)+'_' + str(trj_stage_number)+'/'+jobname +'.xtc'
    #if its the first run, do this to get the same frames schrodinger has used for display of output hotspots
    #convert a subset of the trj frames to structures in a maegz file (cant do all, too much data/too slow). change hotspot_frames input to whatever way you want to slice 

    ##modified aug-2025: use parch to keep waters constant from gcmc runs -
    #set dew asl -sphere areound to keep waters-important for water analysis
    #uses hotspot frame as interval -assumes constant interval-should be correct unless frames are chosen differently
    interval = int(hotspot_frames_list[1])- int(hotspot_frames_list[0])
    slicedcmspre = mxmd_dir + '/' + jname + '_mxmd_'+probe+'-'+str(i)+'_' + str(trj_stage_number)+'/'+ jobname +'_hotslice'
    #run once per replica, modify below to extract consecutive frames from new slicde trajectory, use framelist only for bookkeeping
#    slicedcmspre
#$SCHRODINGER/run trj_parch.py -dew-asl 'res.num 100' -s :: -n 160 sliced_trjslice-out.cms sliced_trjslice_trj slice_dfe_parch 


#$SCHRODINGER/run trj2mae.py  -extract-asl "(all and NOT(res.ptype POPC) and NOT(water)) OR res.ptype SX1 OR res.ptype UNK" -align-asl "protein AND (chain.name C)" -out-format MAE slice_dfe_parch-out.cms slice_dfe_parch_trj/ imidazole_subset_slicedtrj_3_1041

    if not os.path.isfile(slicedcmspre+'-out.cms'):
        #print("1x extract hotspots to new cms for parching")    
        message_queue.put(f"1x extract hotspots to new cms for parching {slicedcmspre}")
        cmd = ['$SCHRODINGER/run','trj_parch.py','-dew-asl','\"res.num 100\"','-ref-mae',outcms,'-align-asl','"protein AND (chain.name C)"', '-s','::'+str(interval),'-n','160',outcms,trj,slicedcmspre]

#        cmd = ['$SCHRODINGER/run','trj_parch.py','-dew-asl','\"res.num 100\"','-ref-mae',outcms,'-align-asl','"protein AND (chain.name C)"', '-s','::'+str(interval),'-n','160',outcms,trj,slicedcmspre]
        p = sp.Popen(' '.join(cmd), stdout = sp.PIPE, shell = True)
        out,err = p.communicate()
        #print(out,err)
        #print(cmd)
        #message_queue.put(f"thread cmd {cmd}")
        p.wait()
        if p.returncode == 0:
            message_queue.put(f"trj_parch success!!! {cmd}")

#            done = True 
        else:
            print("trj_parch error") 
            message_queue.put(f"trj_parch ERROR {cmd}")

        message_queue.put(f"thread {out} Error:{err}")
    
    processes = []

#    handler_msg = threading.Thread(target=message_handler, daemon=True)
#    handler_msg.start()
    message_queue.put("starting extract frames parallel")
#    time.sleep(20)
    for j in range(len(hotspot_frames_list)):
#        thread = threading.Thread(target=extract_and_write, args=(hotspot_frames_list,j,xyz_center,pocket_size,i,interaction,outcms,trj,mxmd_dir,pocket_dir,probe,solvent_pdb))
        thread = threading.Thread(target=extract_and_write, args=(hotspot_frames_list,j,xyz_center,pocket_size,i,interaction,slicedcmspre+'-out.cms',slicedcmspre+'.xtc',mxmd_dir,pocket_dir,probe,solvent_pdb))

        processes.append(thread)
        thread.start()
    current_directory = os.getcwd()

    for process in processes:
        process.join()

#    message_queue.put("STOP")
#    handler_msg.join()
    return
'''
made so this can be multithreaded. "j" is the hotspot frame index, "i" is the trj number
'''
def extract_and_write(hotspot_frames_list,j,xyz_center,pocket_size,i,interaction,outcms,trj,mxmd_results_dir,pocket_dir,solv,solvent_pdb):
    check_subset = mxmd_results_dir + '/' + str(solv) + '_subset_' + str(i) +"_"+ str(hotspot_frames_list[j]) + '.maegz'
    os.chdir(mxmd_results_dir + "/")
    done = True 
    if not os.path.isfile(check_subset):
        done = False
        print("getting subset of" + str(i)+" frame:"+str(j)+"_"+str(hotspot_frames_list[j]))     
        cmd = ['$SCHRODINGER/run','trj2mae.py','-trj-frame-cutting',str(j),'-extract-asl','"(all and NOT(res.ptype POPC) and NOT(water)) OR res.ptype SX1 OR res.ptype UNK"','-align-asl','"protein AND (chain.name C)"','-out-format','MAE',outcms,trj,(mxmd_results_dir+'/'+str(solv)+'_subset_' + str(i) +"_"+ str(hotspot_frames_list[j]))]

#        cmd = ['$SCHRODINGER/run','trj2mae.py','-trj-frame-cutting',hotspot_frames_list[j],'-extract-asl','"(all and NOT(res.ptype POPC) and NOT(water)) OR res.ptype SX1 OR res.ptype UNK"','-align-asl','"protein AND (chain.name C)"','-out-format','MAE',outcms,trj,(str(solv)+'_subset_' + str(i) +"_"+ str(hotspot_frames_list[j]))]
        p = sp.Popen(' '.join(cmd), stdout = sp.PIPE, shell = True)
        out,err = p.communicate()
        print(out,err)
        print(cmd)
        message_queue.put(f"thread cmd {cmd}")
        p.wait()
        if p.returncode == 0:
            done = True 
        else:
            print("trj2mae error") 
            message_queue.put(f"trj2mae error")

        message_queue.put(f"thread {out} Error_trj2mae:{err}")

    #rewrite so its one probe per frame with probes in the region of interest
    os.chdir(pocket_dir)
    #print("slicing to get hotspot of subset")
    if done:    
        write_frames((mxmd_results_dir + '/'+str(solv) + '_subset_' +str(i) + "_" +  str(hotspot_frames_list[j])+'.maegz'),pocket_dir,xyz_center,pocket_size,hotspot_frames_list[j],i,solvent_pdb)

'''
runs the interaction analysis tool on each extracted frame, where ligand is the solvent in the pocket 
'''
def run_interaction_analysis(file,pdb):
    lig_asl = '"res.ptype ' + pdb +'"'
    #print(lig_asl)
    prot_asl = '"protein"'
    cmd = ['$SCHRODINGER/run','poseviewer_interactions_hub.py',file,'-lig_asl',lig_asl,'-rec_asl',prot_asl,'-csv']
    #print(' '.join(cmd))
    p = sp.Popen(' '.join(cmd), stdout = sp.PIPE, shell = True)
    out,err = p.communicate()

'''
given a coordinate and size, this function finds all atoms in the size range and finds which solvent molecule, if any, it belongs to. it keeps track of whole solvent molecules
in the return dictionary. it excludees waters and amino acids. output is dict with resnum as key, and value is pdb code and atom nums. example: {16: ['P1R', [5098, 5099, 5100, 5101, 5102, 5103, 5104, 5105, 5106, 5107]]} 
'''    
def get_probe_in_hotspot(st,xyz_center,pocket_size,solvent_pdb):
    #given xyz, get all probes within A angstroms 
    #get location of each of those probes 
    atoms_in_hotspot=measure.get_atoms_close_to_point(st,xyz_center,pocket_size) #get all atoms close to centroid of hotspot 
    hotspot_solvents = {}
    found = False
    if len(atoms_in_hotspot) > 0:
        found=True 
        for atom in atoms_in_hotspot:
            res_num = st.atom[atom].resnum #get the residue of an atom 
            pdb_code = st.atom[atom].pdbres
            pdb = pdb_code.strip()
            #print(pdb,solvent_pdb)
            if pdb == solvent_pdb:
            #if pdb not in ['UNK','SX1','ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS','MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                if res_num not in hotspot_solvents.keys():
                    hotspot_solvents[res_num] = [pdb]
        for residue_number,residue_name in hotspot_solvents.items():
            residue_atoms = [atom.index for atom in st.atom if atom.pdbres.strip() == residue_name[0] and atom.resnum == residue_number]
            hotspot_solvents[residue_number].append(residue_atoms)
    return hotspot_solvents, found

'''
for each solvent in the hotspot designated, generate a file that has the structure from that frame in the trj and the single extracted solvent in the pocket 
'''
def write_frames(input_file,pocket_dir,xyz_center,pocket_size,frame,run_num,solvent_pdb):
    with StructureReader(input_file) as reader:
            #print("in write_frames")
            #for each file in the sliced subset  - should only be 1 in this case 
            for structure in reader: 
                hotspot_solvents,found = get_probe_in_hotspot(structure,xyz_center,pocket_size,solvent_pdb) 
                #if there is any solvents even in that hotspot at this moment in time
                if len(hotspot_solvents) <= 0:
                    print('no hotspots found in', input_file)
                    message_queue.put(f"thread no hotposts {input_file}")

                if len(hotspot_solvents) > 0:
                    #need to create a new file for each solvent found in the hotspot
                    for residue in hotspot_solvents.keys():
                        output_file = pocket_dir + '/sliced_subset_' + str(run_num) + "_" + str(frame) + "_" + str(residue)+'.maegz'
                        with StructureWriter(output_file) as writer: 
                            new_st = structure.copy() 
                            protein_atoms = analyze.evaluate_asl(structure, 'protein') 
                            solvent_atoms = [] 
                            solvent_atoms += hotspot_solvents[residue][1]             
                            atoms2delete = []
                            ligand_atoms = []
                            ligands = analyze.find_ligands(structure)
                            if len(ligands)>0:
                                for ligand in ligands:
                                    lig_pdb= ligand.pdbres.strip()
                                    if lig_pdb == 'SX1' or lig_pdb == 'UNK':
                                        ligand_atoms = ligand.atom_indexes
                                        break
                            for atom in new_st.atom:
                                if atom.index not in protein_atoms+solvent_atoms+ligand_atoms:
                                    atoms2delete.append(atom.index)
                            
                            new_st.deleteAtoms(atoms2delete) 

                            writer.append(new_st)
'''
combine all the subset interaction files into a dataframe, then into a large csv file to be used for downstream analysis
'''
def pandas_intrxn_analysis(probe_path):
    directory = probe_path
    prefix = 'sliced_subset_'
    suffix = '.csv'
    files=[]
    # Loop through all files in the directory
    for filename in os.listdir(directory):
    # Check if the filename starts with the given prefix
        if filename.startswith(prefix) and filename.endswith(suffix):
            files.append(filename)
    #print(directory,files)
    df = pd.read_csv(files[0])
    probe_count = len(files)
    df['trj_Frame'] = files[0][14:21]
    for i in range(1,len(files)):
        df2 = pd.read_csv(files[i])
        df2['trj_Frame'] = files[i][14:21]
        df = df.append(df2)
    df.to_csv('interaction_analysis_combined.csv')

'''
run in the probe directory to combine each trj runs frames. 
'''
def concat_subsets(probe):
        #change this depending how many probe runs
        directory = '.'
        prefix = 'sliced_subset_'
        suffix = '.maegz'
        cmd = ['$SCHRODINGER/run','structcat.py','-o',(probe+'_poses.maegz')]
        subsets_found=False
        # Loop through all files in the directory
        for filename in os.listdir(directory):
            # Check if the filename starts with the given prefix
            if filename.startswith(prefix) and filename.endswith(suffix):
                cmd.append(filename)
                subsets_found=True
        if not subsets_found:
           print("No sliced files in directory, hotspot empty or pipeline error occured")
        else:
           p = sp.Popen(' '.join(cmd), stdout = sp.PIPE, shell = True)
           out,err = p.communicate()
        return subsets_found

'''
run schrodinger conformation clustering in probe folder. issue: this cant run with variable amounts of atoms, so i had to remove solvent. which means for final 
results we will have to go back and get that same frame but with the solvent. 
'''
def conformer_cluster(probe):
    print('submitted conformer_cluster.py for probe...')
    asl = "'protein'"
    cmd = ['$SCHRODINGER/run','conformer_cluster.py','-a',asl,'-HOST','sge_cpu','-n','0','-l','"Ward"','-rep','-o','clustered','-OVERWRITE',(probe +'_poses.maegz'),'-WAIT']
    p = sp.Popen(' '.join(cmd), stdout = sp.PIPE, stderr=sp.PIPE, shell = True)
    out,err = p.communicate()
    #print(out,err)



def extract_closest(df):
  framelst = []
  #find the 6 closest frames to centroid:
  df.sort_values(by=['DistanceToCentroid'], inplace=True)
  nfr = len(df.index)
  if nfr>6:
     nfr=6
  #frame = df.iloc[[0]]
  #print(frame['ItemLabel'].values)
  #print('cluster number:',frame.ClusterNumber.values)
  #print(f'cluster pop: {frame.ClusterItemCount.values}')
  #capture multiple centroids in a cluster
  dfcent = df.loc[(df['IsNearestToCentroid']==1)]
  cfr = len(dfcent.index)
  for i in range(0,cfr):
     frame = df.iloc[[i]]
     #print('centroid',frame['ItemLabel'].values)
     framelst.append([frame['ItemLabel'].values[0], frame.ClusterNumber.values[0]])
     #print('framelist-centroid',framelst)

  for i in range(cfr,nfr):
     frame = df.iloc[[i]]
     #print(frame['ItemLabel'].values)
     framelst.append([frame['ItemLabel'].values[0], frame.ClusterNumber.values[0]])
     #print('framelist',framelst)
  return(framelst)



'''
clustering analysis & finding centroids script from phillip aoto. run after schrodinger clustering script in the probe dir. 
'''
def extractN_clustering():
#    def extract_closest(df):
#        #find the closest frames to centroid:
#        df = df.sort_values(by=['DistanceToCentroid'])
#        frame = df.iloc[[0]]
#        frame_centroid = (frame['ItemLabel'].values)[0]
#        return(int(frame_centroid))

    datalsn = pd.read_csv('clustered_ligand1.grp', sep=',')
    print('analyzing clusters...')
    clustermax = np.max(datalsn['ClusterItemCount']) #find the max population amount
    # print(clustermax)
    cluster_order = datalsn.sort_values(by=['ClusterItemCount'],ascending=False) #sort by highest pop
    unique_labels = cluster_order['ClusterNumber'].drop_duplicates().to_list()[:10] #get the cluster labels for the 3 highest pop
   
    centroids = []
    cluster_pop = {}
    
    for cluster_num in unique_labels:
        maxidx = datalsn.index[datalsn['ClusterNumber']==cluster_num].tolist() #get the indices of the cluster
        datalsn0 = datalsn.iloc[maxidx] #get subset with just those of this cluster
#        population = datalsn0['ItemLabel'].to_list()
#        population = datalsn0['ItemLabel'].to_list()
#ClusterItemCount
        population = datalsn0['ClusterItemCount'].to_list()
        centroid = extract_closest(datalsn0)
#        centroids.append(centroid)
        centroids.extend(centroid)
        for centi in centroid:
            #cluster number and population
            cluster_pop[centi[1]]=population[0]
        #print(centroid,cluster_pop[centroid])
    with open("frames_rep.txt", "w", newline="") as fo:
        writer = csv.writer(fo)
        writer.writerows(centroids)

#    with open('frames_rep.txt', 'w') as fo:
#        fo.write('\n'.join(str(i) for i in centroids))
#        fo.write('\n')
    print('cluster_pop',cluster_pop)
    with open('frames_cluster.txt','w') as f:
        for key,value in cluster_pop.items():
            f.write(str(key) + ':' +str(value))

            #f.write((' '.join (str (i) for i in value)))
            f.write('\n')
    
    return centroids,cluster_pop

'''
this just extracts the found low conformation frame and puts it in a final directory so you dont have to dig for it yourself 
'''

def extract_centroid(centroid_frame,cluster_number,cluster_pop,probe,pocket_name,pocket_dir):
    subset = probe + '_poses.maegz'
    st = structure.StructureReader.read(subset, index=centroid_frame)
    out = (pocket_name + '_' + probe + '_frm' + str(centroid_frame)+"_cluster_"+str(cluster_number)+"_pop_"+str(cluster_pop)+".maegz") 
    with structure.StructureWriter(out) as writer:
        writer.append(st)

    shutil.copy((out),(os.path.join(pocket_dir,'pocket_conformations')))

'''
"main" -> puts together the other functions to create the basic workflow front to end 
'''
def per_probe_analysis(solv,pocket_dir,base,jname,pname,xyz_center,size,cluster,interaction,sim_stage,pdb,run_num=10):
        #global hotspot_frames_list 
        hotspot_frames_list = []
        print('running %s analysis...' %solv)
        os.chdir(pocket_dir)
        probe_path = (os.path.join(pocket_dir,str(solv)))
        try:
            os.mkdir(probe_path)
        except FileExistsError:
            pass
        os.chdir(probe_path)
        mxmd_results_path = base + '/'+jname + '_3'
        #extract the trajectory information . this line below checks if its already been extracted by looking at one of the stages. probably not the best but its good enough :') 
        #if you get a read_traj errorm, it could be because this file was already extracted so it didnt go through with the extraction 
        check_trj = mxmd_results_path  + '/' + jname + '_mxmd_'+probe+'-'+str(9)+'_' + str(sim_stage)+'/'+jname + '_mxmd_' + probe + '-' + str(9) + '_' + str(sim_stage)+'.xtc'

        if not os.path.isfile(check_trj):
            processes = []
            for i in range(run_num):  
                thread = threading.Thread(target=extract_trjs, args=(solv,mxmd_results_path,jname,probe_path,sim_stage,i))
                processes.append(thread)
                thread.start()
            for process in processes:
                process.join()
        else:
            print("trajectory files found")
        processes=[]
        #wont re-run  slicing if its already done
        files = os.listdir(os.getcwd())
        sliced_done = False
        out_files = []
        for file in files:
            if file.startswith('sliced') and file.endswith('.maegz'):
                out_files.append(file)
                sliced_done = True
        if not sliced_done:
            print("slicing trajectories")
            handler_msg = threading.Thread(target=message_handler, daemon=True)
            handler_msg.start()
            message_queue.put("starting slices parallel")


            for i in range(run_num):
                #generate hotspot_frames_list
                #NOTE: change hotspot_frames_list to whatever increment you want to slice up your trj with to analyze (cant analyze all of trj without crashing). this list is currently the list of frames schrodinger uses to generate mxmd maps.
                if i == 0:
                    if not os.path.isfile('hotspot_frames.txt'):
                        outcms = jname + '_mxmd_'+probe+'-'+str(i)+'_' + str(sim_stage)+'/'+jname + '_mxmd_' + probe + '-' + str(i) + '_' + str(sim_stage)+'-out.cms'
                        trj = mxmd_results_path  + '/'+jname + '_mxmd_'+probe+'-'+str(i)+'_' + str(sim_stage)+'/'+jname + '_mxmd_' + probe + '-' + str(i) + '_' + str(sim_stage)+'.xtc'
                        tmp_tr = traj.read_traj(trj)
                        gap = math.floor(len(tmp_tr)/20)
                        for k in range(0,len(tmp_tr),gap):
                            hotspot_frames_list.append(str(k))
                        hotspot_frames = ','.join(hotspot_frames_list)
                        with open('hotspot_frames.txt', 'w') as fo:
                            fo.write(hotspot_frames)   
                    else:             
                        with open('hotspot_frames.txt', 'r') as fi:
                            l = fi.readline()
                            hotspot_frames_list = [str(num) for num in l.split(',')]

                #now get the frames that have the solvent in hotspot
                print("slicing %s" %(str(i)))
                message_queue.put(f"slicing parallel {i}")

                thread = threading.Thread(target=extract_frames_parallel, args=(hotspot_frames_list,xyz_center,size,interaction,i,sim_stage,solv,jname,mxmd_results_path,probe_path,pdb))
                processes.append(thread)
                thread.start()   
            for process in processes:
                process.join()
            message_queue.put("STOP")
            handler_msg.join()
      
        os.chdir(probe_path)
        if interaction:
            files = os.listdir(os.getcwd())
            out_files = []
            for file in files:
                if file.startswith('sliced') and file.endswith('.maegz'):
                    out_files.append(file)
            source = os.path.join(pocket_dir,'poseviewer_interactions_hub.py')
            destination =os.path.join(probe_path,'poseviewer_interactions_hub.py')  # File path in Directory B
            shutil.copy(source, destination)
            for file in out_files:
                run_interaction_analysis(file,pdb)
            pandas_intrxn_analysis(probe_path)
        if cluster:
            #checks to continue workflow if something failed 
            if not os.path.isfile(str(solv)+'_poses.maegz'):

                print("combining files...")
                subsets_found = concat_subsets(solv)
                if subsets_found:
                  
            
                   if not os.path.isfile('clustered_ligand1.out'):
                        conformer_cluster(solv)

                   try:
                        centroids,cluster_pops = extractN_clustering()

                   except:
                        print("waiting for cluster run to be done")
                        while True:
                              if os.path.isfile('clustered_ligand1.out'):
                                  print("clustering done")
                                  break
                   centroids,cluster_pops = extractN_clustering()
        
                   for cent in centroids:
#                        print(cluster_pops[cent[1]],'clusterpops')
                        extract_centroid(int(cent[0]),int(cent[1]),int(cluster_pops[cent[1]]),solv,pname,pocket_dir)
    
        print('probe %s completed' %solv)

def message_handler():
    while True:
        message = message_queue.get()
        if message == "STOP":
            break
        print(message)
        

if __name__== '__main__':

    parser = argparse.ArgumentParser(description='run this script to analyze resiude-solvent interactions around a hotspot of interest. output are low energy mae snapshots of specified resiude conformations sampled during the 10 sims of a given probe',usage='move a copy of this script to the mxmd run directory. example of run command:  $SCHRODINGER/run pocket_analysis.py --pname cryptic_1 --jname mxmd --size 1.5 --center -0.8w1  9.15 -3.01 --cluster --interaction --simulationStage 7 --probes benzene',formatter_class=argparse.RawTextHelpFormatter)
    #argparse usage example: python md_workflow.py GIPR --analyze 30
    #all arguments with "--" are OPTIONAL, with exception to pname and jname. options: run only prep, run only MD, run only analysis, run md + analysis, run prep + md + analysis, 
    parser.add_argument('--pname',nargs=1,type=str,required=True,help='pocket name for the analysis folder')
    parser.add_argument('--jname',nargs=1,type=str,required=True,help='job name from the mxmd run')
    parser.add_argument('--probes',nargs='+',type=str,required=True,help='add the name of the probes you want to do the analysis for, separated by spaces. ex: --probes acetaldehyde pyrimidine')
    #parser.add_argument('--parallel',action="store_true",required=False,help='ignore - used automatically')
    parser.add_argument('--size',nargs=1,type=str,required=True,help='size in angstroms of the pocket of interest')
    parser.add_argument('--center',nargs='+',type=str,required=True,help='xyz coordinates of the center of the pocket of interest reference the out.cms. to get this value, open up the mxmd prj file and select the hotspots in the entry list that are in the region of interest. select a representative atom / probe close to the center, and to see its XYZ go to style->apply labels-> edit custom label-> xyz coordinate. example: --center -0.23 7.48 -6.65')
    parser.add_argument('--cluster',action="store_true",required=False,help='perform clustering on receptor conformation with hotspot probes in a specific pocket')
    parser.add_argument('--interaction',action="store_true",required=False,help='run a pose view interaction analysis on all probes in a hotspot in context of frame')
    parser.add_argument('--simulationStage',nargs=1,type=str,required=True,help='the stage that you can find the trjs for all the runs, usually 5,6,or 7.')
    parser.add_argument('--solventsPDB',nargs='+',type=str,required=True,help='the pdb code of each of the solvent probes youre running this one (i.e. --solventsPDB EOH CCN)')
    args = parser.parse_args()

    #extract user input
    probes = args.probes #list of the probes in the run (i.e [pyrimidine,isopropanol,acetonitrile])
    
    pname = args.pname[0] #pocket name
    jname = args.jname[0] #job name 
    #parallel = args.parallel
    cluster = args.cluster
    sim_stage = int(args.simulationStage[0])
    interaction = args.interaction
    pocket_size = float(args.size[0])
    pocket_center_raw = args.center
    pocket_center = [float(x) for x in pocket_center_raw]
    pdbs = [str(x) for x in args.solventsPDB]
    base = os.getcwd() #directory we're running this in 
    directory = pname + '_pocket_analysis'
    pocket_dir = (os.path.join(base,directory))
    
    # if parallel:
    #     per_probe_analysis(probes[0],pocket_dir,base,jname,pname,pocket_center,pocket_size,cluster,interaction,sim_stage)
    try: #create pocket directory
        os.mkdir(pocket_dir)
    except FileExistsError:
        pass
    try: #create the final set of conformation directory
        os.mkdir(os.path.join(pocket_dir,'pocket_conformations'))
    except FileExistsError:
        pass
    source = os.path.join(base,'poseviewer_interactions_hub.py')
    destination = os.path.join(pocket_dir,'poseviewer_interactions_hub.py')
    shutil.copy(source, destination)
    p = 0
    for probe in probes:
        print(probe)
        per_probe_analysis(probe,pocket_dir,base,jname,pname, pocket_center,pocket_size,cluster,interaction,sim_stage,pdbs[p])
        p+=1
