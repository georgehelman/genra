import numpy as np
import pandas as pd
from genrapred import *

def plot_worthy(pdobject):
    if isinstance(pdobject,pd.core.series.Series):
        pdobject=pdobject[pd.notnull(pdobject)]
        pdobject=pdobject[pdobject!=np.inf]
        return pdobject
    elif isinstance(pdobject,pd.core.frame.DataFrame):
        pdobject=pdobject[pdobject.notnull().all(axis='columns')]
        pdobject=pdobject[(pdobject!=np.inf).all(axis=1)]
        return pdobject
    
def wtavg(df,name,k,s):
    df=df[df['jaccard']>=s]
    df=df[df[name]!=np.inf]
    df=df[df[name].notnull()].iloc[0:k]
    if df.empty:
        return np.nan
    weights=list(df['jaccard'])
    values=list(df[name])
    return np.average(values,weights=weights)

def exact_k_wtavg(df,name,k,s):
    df=df[df['jaccard']>s]
    df=df[df[name]!=np.inf]
    df=df[df[name].notnull()].iloc[0:k]
    if len(df)<k:
        return np.nan
    weights=list(df['jaccard'])
    values=list(df[name])
    return np.average(values,weights=weights)

def wtvar(df,name,k):
    df=df[(df[name].notnull()) & (df[name]!=np.inf)].iloc[0:k]
    if df.empty:
        return np.nan
    weights=list(df['jaccard'])
    values=list(df[name])
    return sum([weights[i]**2*values[i] for i in range(len(values))])/sum(weights)**2

def find_neighbors(sids,DB):
    neighbors=[]
    for sid in sids:
        sid_neighbors=searchCollByFP(sid,s0=.05,SID=sids,DB=DB)
        if sid_neighbors:
            sid_neighbors = [{'target_sid':sid, \
                              'neighbor_sid':record['dsstox_sid'], \
                              'jaccard':record['jaccard']} \
                             for record in sid_neighbors]
            neighbors=neighbors+sid_neighbors
    neighbors=pd.DataFrame(neighbors)
    neighbors=neighbors[neighbors['target_sid']!=neighbors['neighbor_sid']]
    return neighbors

def make_aggregate_df(df,agg):
    assert agg=='min' or agg=='mean', 'agg must equal min or mean'
    if agg=='min':
        return df.pivot_table(index='dsstox_sid',columns='endpoint_category',values='pod_value_LM',aggfunc='min')
    else:
        return df.pivot_table(index='dsstox_sid',columns='endpoint_category',values='pod_value_LM',aggfunc='mean')
    
def make_neighbors_df(df,neighbors,agg):
    agg_df=make_aggregate_df(df,agg)
    neighbors=neighbors[neighbors['target_sid']!=neighbors['neighbor_sid']]
    neighbors=neighbors.merge(agg_df,left_on='neighbor_sid',right_index=True)
    neighbors_df=neighbors.sort_values('jaccard',ascending=False)
    return neighbors_df
    
