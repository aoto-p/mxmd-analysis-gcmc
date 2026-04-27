import matplotlib.pyplot as plt
import numpy as np
from schrodinger import structure
import pandas as pd
st = structure.StructureReader('cat_probe_top.maegz')
#out = (pocket_name + '_' + probe + '_frm' + str(centroid_frame)+"_cluster_"+str(cluster_number)+"_pop_"+str(cluster_pop)+".maegz") 
#for frame in structure.StructureReader('myfile.sdf'):
edict={}
for idx, frame in enumerate(st):
    edict[idx]=frame.property["r_cluster_energy"]

df = pd.DataFrame.from_dict(edict,orient='index', columns=['r_cluster_energy'])
#print(df)
# = (x - old_min) / (old_max - old_min) * (new_max - new_min) + new_min
df['r_cluster_energy_linscale'] = ((df.r_cluster_energy-df.r_cluster_energy.min()) / (df.r_cluster_energy.max()-df.r_cluster_energy.min()) )* (-2.5+12.5)-12.5
df['r_cluster_activity_linscale'] = (-1*np.log10(np.exp(df['r_cluster_energy_linscale']/(0.001987*300))))

print(df[['r_cluster_energy_linscale','r_cluster_energy']])
out='cat_prob_top_scaled.maegz'
st = structure.StructureReader('cat_probe_top.maegz')

with structure.StructureWriter(out) as writer:
    for idx, frame in enumerate(st):
        frame.property["r_cluster_activity_linscale"] = df.loc[idx,'r_cluster_activity_linscale']

#        print(frame.property["r_cluster_energy"])
#    st.property["r_cluster_num"] = float(cluster_number)
#    st.property["r_cluster_trjfr"] = float(centroid_frame)
#    st.property["r_cluster_pop"] = float(clpop)
#    st.property["r_cluster_energy"] = float(clene)
#    st.property["r_cluster_activity"] = float(clactivity)
        writer.append(frame)
#ax = df.plot.scatter(x='r_cluster_energy', y='r_cluster_energy_linscale', c='Red')
#ax = df['r_cluster_activity_linscale'].plot.hist(bins=12, alpha=0.5)

#plt.show()
