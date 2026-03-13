import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#Dscoreall from grep -A 1
#remove -- rows with:
#:g/--/d
#change first row to include name column
#skiprows=lambda x: x > 0 and x % 2 == 1 
#covpock12_c_gcmc_acetate_frm102_cluster_18_pop_39/covpock12_c_gcmc_acetate_frm102_cluster_18_pop_39.log:SiteScore size   Dscore  volume  exposure enclosure contact  phobic   philic   balance  don/acc  refdist  refmin   refavg  sitemin
#pareto front from chatgpt
def pareto_front(df, objectives, maximize=True):
    """
    Compute the Pareto front from a DataFrame.

    Parameters:
        df (pd.DataFrame): Input data.
        objectives (list): Column names representing objectives.
        maximize (bool or list): If True, higher is better for all objectives.
                                 If list, specify per-objective True/False.

    Returns:
        pd.DataFrame: Pareto front subset of df.
    """
    if isinstance(maximize, bool):
        maximize = [maximize] * len(objectives)
    if len(maximize) != len(objectives):
        raise ValueError("Length of maximize list must match objectives list.")

    # Convert to numpy for speed
    data = df[objectives].to_numpy()
    is_pareto = np.ones(data.shape[0], dtype=bool)

    for i, point in enumerate(data):
        if not is_pareto[i]:
            continue
        # Compare against all other points
        for j, other in enumerate(data):
            if i == j:
                continue
            better_or_equal = []
            strictly_better = []
            for k, maximize_k in enumerate(maximize):
                if maximize_k:
                    better_or_equal.append(other[k] >= point[k])
                    strictly_better.append(other[k] > point[k])
                else:
                    better_or_equal.append(other[k] <= point[k])
                    strictly_better.append(other[k] < point[k])
            # If other dominates point, mark as not Pareto
            if all(better_or_equal) and any(strictly_better):
                is_pareto[i] = False
                break

    return df[is_pareto].reset_index(drop=True)

df =pd.read_csv('Dscoreall.out',skiprows=lambda x:x>0 and x%2==0, sep=r'\s+')
#print(df.Dscore)
df = df.sort_values(by='Dscore', ascending=False)
#df0 = pd.read_csv('../../mxmd_charge/site_batch/Dscoresorted.csv')

#df = pd.concat([df,df0])
df = df.sort_values(by='Dscore', ascending=False)

df.to_csv('Dscoresorted.csv', index=False)

dfpar = pareto_front(df, objectives=['Dscore','volume'])
print(dfpar)
dfpar.to_csv('pareto_front_dscore_vol.csv', index=False)
ax = df.plot.scatter(x='Dscore', y='volume', c='DarkBlue')
dfpar.plot.scatter(x='Dscore', y='volume', c='Red', ax=ax)
plt.show()
