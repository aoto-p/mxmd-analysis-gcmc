
import os
import shutil 
import numpy as np
import subprocess as sp
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import argparse
from schrodinger.application.desmond.packages import analysis
from schrodinger.application.desmond.packages import topo
from schrodinger.application.desmond.packages import traj_util
from schrodinger.structure import create_new_structure
from schrodinger.application.desmond import cmj
from schrodinger.application.desmond import cms
from schrodinger.application.desmond import constants
from schrodinger.application.desmond import envir
from schrodinger.application.desmond import struc
from schrodinger.application.desmond import util
from schrodinger.application.desmond.cns_io import write_cns
from schrodinger.application.desmond.mxmd import mxmd_system_builder as msb
from schrodinger.structutils.analyze import evaluate_asl
from schrodinger import structure

def write_cns_mae_files(grid_data,normgrid_data,box_size,grid_spacing,center):
        """
        Write CNS for each probe type.

        :return: List of cns files for each probe
        """
        cns_files = []
        for probe, data in grid_data.items():
            
            cns_fn = str(probe) + '-water.cns'
            grid = normgrid_data[probe]
            write_cns(grid,
                      box_size,
                      grid_spacing,
                      cns_fn,
                      center=center)

            cns_files.append(cns_fn)
            
        return cns_files

def gen_normgrid_data(grid_data):
    """
    Generate normalized occupancy data
    """
    normgrid_data = {}
    for probe, data in grid_data.items():
        grid = np.zeros(data[0].shape)
        for g in data:
            grid += g
        grid /= len(data)
        normgrid_data[probe] = normalize_probe_occupancy_grid(grid)
    return normgrid_data
def normalize_probe_occupancy_grid(grid):
    '''
    Convert the grid counts to z-scores by normalizing the grid. Since the grid
    contains probe occupancies, this matrix is mostly sparse (with ~1% of non-
    zero values) we need to omit all the zero-containing values from calculating
    of mean and standard deviation.
    '''
    normgrid = np.zeros(grid.shape, dtype='float16')
    mask = grid != 0.0
    normgrid[mask] = grid[mask] - grid[mask].mean()
    normgrid /= grid[mask].std()
    return normgrid

def water_analysis(jobname, jname,cms_fname, mae_fname,probe,rep):
    st = structure.StructureReader.read(mae_fname)
    box_length = st.property.get('r_mxmd_box_length') 
    grid_spacing = st.property.get('r_mxmd_grid_spacing') #usually 0.5
    align_asl = '(protein and atom.ptype CA)'
    #read in cms file 
    msys_model, cms_model, tr = traj_util.read_cms_and_traj(cms_fname)
    #get cosolvent structure 
    solvent_ct = struc.component_structures(cms_model).solvent

    solvent_aids=cms_model.get_fullsys_ct_atom_index_range(cms_model.comp_ct.index(solvent_ct))

    solvent_aids_noh = [
        aid for aid in solvent_aids
        if cms_model.fsys_ct.atom[aid].atomic_number != 1
    ]
    solvent_mols = solvent_ct.mol_total

    solvent_probe='water'
    # extract reference structure and use it as a reference to align
    # trajectory and cms_model
    ref_mae = struc.get_reference_ct(cms_model.fsys_ct)
    ref_pos = ref_mae.extract(evaluate_asl(ref_mae, align_asl)).getXYZ()
    center_pos = np.mean(ref_pos, axis=0)
    ref_gids = topo.asl2gids(cms_model,
                            align_asl,
                            include_pseudoatoms=False)
    tr = topo.superimpose(msys_model, ref_gids, tr, ref_pos)
    cms_model = topo.superimpose_cms(msys_model, ref_gids, cms_model,
                                    ref_pos)



    grid_spacing = [grid_spacing] * 3
    box_length = [box_length] * 3

    vmw = analysis.VolumeMapper(cms_model,
                            aids=solvent_aids_noh,
                            spacing=grid_spacing,
                            length=box_length,
                            center=center_pos,
                            normalize=False)

    grid = analysis.analyze(tr,vmw)
    # reduce precision
    grid = grid.astype('uint16')
    normgrid=normalize_probe_occupancy_grid(grid)
    # Write .cns and .raw files. The later file will be used by the
    # cleanup stage to generate aggregate data for each probe type.
    _, cms_fname = os.path.split(cms_fname)
    out_cns_fname =  jname+'_3/'+jobname+'/' + probe+ '_'+str(i)+'-water-out.cns'
    out_raw_fname = jname+'_3/'+jobname+'/' + probe+ '_'+str(i)+'-water-out.raw'
    out_mae_fname = jname+'_3/'+jobname+'/' + probe+ '_'+str(i)+'-water-out.mae'
    out_probes_fname = jname+'_3/'+jobname+'/' + probe+ '_'+str(i)+'-water-probes.mae'
    cms_fname = os.path.join(os.getcwd(), cms_fname)
    write_cns(normgrid, box_length, grid_spacing, out_cns_fname,center=[center_pos[0],center_pos[1],center_pos[2]])
    # Save probe molecules from 20 frames into a CT for later use.
    solvent_probes_ct = create_new_structure()

    _pct = struc.delete_ffio_ff(solvent_ct)

    solvent_gids = topo.aids2gids(cms_model,
                                solvent_aids,
                                include_pseudoatoms=False)            
    nframes = len(tr)
    fr_interval = 1 if nframes < 20 else nframes // 20
    for fr in tr[::fr_interval]:
        _pct.setXYZ(fr.pos(solvent_gids))
        y = _pct.copy()
        for at in y.atom:
            at.property["s_user_ensemble_frame"] = "%s_%d" %(jobname, fr._index)
        solvent_probes_ct = solvent_probes_ct.merge(y)


    solvent_probes_ct.title = _pct.title

    solvent_probes_ct.property[constants.MXMD_COSOLVENT_PROBE] = solvent_probe
    
    solvent_probes_ct.write(out_probes_fname)


    ct = struc.delete_ffio_ff(
    struc.component_structures(cms_model).solute)

    # Set ct properties

    ct.property[constants.MXMD_COSOLVENT_PROBE] = solvent_probe
    ct.property[constants.MXMD_GRID_SPACING] = grid_spacing[0]
    ct.property[constants.MXMD_BOX_LENGTH] = box_length[0]
    ct.property[constants.MXMD_NUM_PROBES] = solvent_mols
    ct.property[constants.MXMD_CENTER_X] = center_pos[0]
    ct.property[constants.MXMD_CENTER_Y] = center_pos[1]
    ct.property[constants.MXMD_CENTER_Z] = center_pos[2]
    ct.write(out_mae_fname)
    


    with open(out_raw_fname, 'wb') as fh:
        np.save(fh, grid)
    
    #shutil.move(out_probes_fname, jname+'_3/'+jobname+'/')
    #shutil.move(out_cns_fname, jname+'_3/'+jobname+'/')
    #shutil.move(out_mae_fname, jname+'_3/'+jobname+'/')
    #shutil.move(out_raw_fname, jname+'_3/'+jobname+'/')
    center=[center_pos[0],center_pos[1],center_pos[2]]
    return normgrid, box_length, grid_spacing, center



if __name__== '__main__':

    parser = argparse.ArgumentParser(description='mxmd analysis that treats waters as solvent probes to analyze water hotspots from competition with solvent',usage='move this to mxmd project directory and run on completed mxmd run',formatter_class=argparse.RawTextHelpFormatter)
    #argparse usage example: python md_workflow.py GIPR --analyze 30
    #all arguments with "--" are OPTIONAL, with exception to pname and jname. options: run only prep, run only MD, run only analysis, run md + analysis, run prep + md + analysis, \
    parser.add_argument('--simulationStage',nargs=1,type=str,required=True,help='the stage that you can find the trjs for all the runs, usually 5,6,or 7.')
    parser.add_argument('--jname',nargs=1,type=str,required=True,help='job name from the mxmd run')
    parser.add_argument('--probes',nargs='+',type=str,required=True,help='add the name of the probes you want to do the analysis for, separated by spaces. ex: --probes acetaldehyde pyrimidine')
    #parser.add_argument('--parallel',action="store_true",required=False,help='ignore - used automatically')
    args = parser.parse_args()

    #extract user input
    probes = args.probes #list of the probes in the run (i.e [pyrimidine,isopropanol,acetonitrile])
    sim_stage = int(args.simulationStage[0])
    jname = args.jname[0] #job name 
    #parallel = args.parallel
    grid_data={}
    for probe in probes:
        grid_data[probe] = []
        for i in range(0,10):
            analysis_stage = sim_stage + 1
            jobname_out = jname + '_mxmd_' + probe + '-' + str(i) + '_' + str(analysis_stage)
            cms_fname = jname + '_3/' + jname + '_mxmd_'+probe+'-'+str(i)+'_' + str(sim_stage)+'/'+ jname + '_mxmd_' + probe + '-' + str(i) + '_' + str(sim_stage) +'-out.cms'
            mae_fname = jname + '_3/' + jname + '_mxmd_' + probe + '-' + str(i) + '-out.mae'
            #extract mxmd cleanup analysis folder for probes from run, to place the water hotspot output 
            cmd = ['tar','-xzf',jname + '_3/'+jobname_out + '-out.tgz','-C', jname + '_3/']
            try:
                p = sp.Popen(' '.join(cmd), stdout = sp.PIPE, shell = True)
                out,err = p.communicate() 
            except:
                print('trj already extracted')
            g, box_length, grid_spacing, center= water_analysis(jobname_out, jname,cms_fname, mae_fname,probe,i)
            grid_data[probe].append(g)
    normgrid_data = gen_normgrid_data(grid_data)
    cns_files = write_cns_mae_files(grid_data,normgrid_data,box_length,grid_spacing,center)
    #  for idx, cns_file in enumerate(cns_files)
    # set_isosurface(pt, cns_file, row_index, self.sigma, color, name,
    #                        comment)
