import numpy as np
import pandas as pd
import copy
import pymongo
import sys
import re

from db.mongo import *

def getFP(SID,fp='chm_mrgn',FP=None,DB=None,fill=None):

    col,ds = getColFPMap(fp)
    if not (ds and col): return col,ds

    Agg = [
            # Match chemicals in cluster
            {'$match': {
                     'dsstox_sid':{'$in':SID}}
            },
            # Include these fields
            {'$project':{'dsstox_sid':1,'name':1,'_id':0,
                        'fp':'$%s.ds'%ds},
            },
            # Unwind the fp 
            {'$unwind':"$fp"}
            ]

    if FP: Agg.append({'$match': {'fp':{'$in': FP}}})
    
    X = DB[col].aggregate(Agg,allowDiskUse=True)
    
    if not X: return
    try:
        R = pd.DataFrame(X['result'])
    except:
        R = pd.DataFrame(list(X))

    if R.shape[0]==0 or R.shape[1]==0: return pd.DataFrame()

    return pd.pivot_table(R,index=['dsstox_sid'],columns='fp',values='name',aggfunc=len,fill_value=fill)

def getChemBioSummary(SID,col=None,ds=None,fill=None,cls='fp'):
    if not cls: cls='fp'
    Agg = [
            # Match chemicals in cluster
            {'$match': {
                     'dsstox_sid':{'$in':SID}}
            },
            # Include these fields
            {'$project':{'dsstox_sid':1,'_id':0,
                        cls:'$'+ds},
            },
            {'$unwind':'$'+cls}

            ]
   
    X = col.aggregate(Agg,allowDiskUse=True)
    if not X: return
    try:
        R = X['result']
    except:
        R = list(X)

    return R

def getChemToxSummary(SID,col=None,ds='tox_fpp1.ds',fill=None,cls='tox_fp'):
    if not cls: cls='fp'
    Agg = [
            # Match chemicals in cluster
            {'$match': {
                     'dsstox_sid':{'$in':SID}}
            },
            # Include these fields
            {'$project':{'dsstox_sid':1,'_id':0,
                        cls:'$'+ds},
            },
            {'$unwind':'$'+cls}

            ]
   
    X = col.aggregate(Agg,allowDiskUse=True)
    if not X: return
    try:
        R = X['result']
    except:
        R = list(X)

    return R



def getChemSummary(SID,MDB=None,Colls=None):
    Agg_match= {'$match': {'dsstox_sid':{'$in':SID}}}
    Agg_proj = {'$project':{'dsstox_sid':1,'_id':0,'n':''}}
    
    Res = []
    for prop,db_coll in Colls.iteritems():
        Agg_proj['$project']['n'] = "$%(ds)s.%(n)s" % db_coll
        Agg = [Agg_match,Agg_proj]
        X = MDB[db_coll['coll']].aggregate(Agg,allowDiskUse=True)
        
        if not X: continue
        R = None
        try:
            R = X['result']
        except:
            R = list(X)
        
        R_df = pd.DataFrame(R)
        R_df['prop']=prop
        Res.append(R_df)
        
    X = pd.concat(Res)
    R = X.pivot_table(index='dsstox_sid',columns='prop',values='n',aggfunc=min)
    R.fillna(0,inplace=True)

    return R



def getChemToxCastNNSummary(sid,s0=0.01,k0=10,fp='chm_mrgn',MDB=None,col='bio_fp',ds='bio1',
                            sel_by=None,                            
                            dbg=False):
    
    Hits  = searchFP(sid,fp=fp,s0=s0,max_hits=k0,DB=MDB,sel_by=sel_by)
    if not Hits:
        return
    
    NN      = pd.DataFrame(Hits)
    k0      = NN.shape[0]

    SID = list(NN.dsstox_sid)

    AI1 = None
    R1  = None

    # Get assay information
    AI1= pd.DataFrame(list(MDB.bio_fp_info.find({'target_family':{'$regex':'^((?!background).)'}},dict(_id=0))))
    AI1.fillna('',inplace=True)

    # get bioactivity data for NN
    B1 = pd.DataFrame(getChemBioSummary(SID,col=MDB[col],ds='%s.ds'%ds,cls='bio_fp'))

    R1 = B1.merge(AI1,left_on='bio_fp',right_on='bio_fp').merge(NN,left_on='dsstox_sid',right_on='dsstox_sid')

    return R1

def getChemTox21NNSummary(sid,s0=0.01,k0=10,fp='chm_mrgn',MDB=None,col='tx21_fp',ds='t211',
                          sel_by=None,                            
                          dbg=False):
    
    Hits  = searchFP(sid,fp=fp,s0=s0,max_hits=k0,DB=MDB,sel_by=sel_by)
    if not Hits:
        return
    
    NN      = pd.DataFrame(Hits)
    k0      = NN.shape[0]

    SID = list(NN.dsstox_sid)

    AI1 = None
    R1  = None

    # Get assay information
    AI1= pd.DataFrame(list(MDB.bio_fp_info.find({'target_family':{'$regex':'^((?!background).)'}},dict(_id=0))))
    AI1.fillna('',inplace=True)

    # get bioactivity data for NN
    B1 = pd.DataFrame(getChemBioSummary(SID,col=MDB[col],ds='%s.ds'%ds,cls='bio_fp'))
    B1['bio_fp']=B1.bio_fp.str.lower()
    AI1['bio_fp']=AI1.bio_fp.str.lower()

    R1 = B1.merge(AI1,left_on='bio_fp',right_on='bio_fp').merge(NN,left_on='dsstox_sid',right_on='dsstox_sid')

    return R1

def getChemToxRefNNSummary(sid,s0=0.01,k0=10,fp='chm_mrgn',col='tox5_fp',ds='tox_fpp1',MDB=None,
                           sel_by=None,                            
                           dbg=False):
    
    Hits  = searchFP(sid,fp=fp,s0=s0,max_hits=k0,DB=MDB,sel_by=sel_by)
    if not Hits:
        return pd.DataFrame()
    
    NN      = pd.DataFrame(Hits)
    k0      = NN.shape[0]

    SID = list(NN.dsstox_sid)
    
    T1 = pd.DataFrame(getChemToxSummary(SID,col=MDB[col],ds='%s.ds'%ds,cls='tox_fp'))
    if T1.shape[0]==0: return pd.DataFrame() 
    T1['study']=T1.tox_fp.apply(lambda x: x.split(':')[0])
    T1['organ']=T1.tox_fp.apply(lambda x: x.split(':')[1].lower())
    R1=T1.merge(NN,left_on='dsstox_sid',right_on='dsstox_sid')
    
    return R1

def getChemToxRefNNInfo(sid,s0=0.01,k0=10,fp='chm_mrgn',
                        MDB=None,col='tox5_fp',ds_pos='tox_fpp1',ds_neg='tox_fpn1',
                        
                        dbg=False):
    
    Hits  = searchFP(sid,fp=fp,s0=s0,max_hits=k0,DB=MDB)
    if not Hits:
        return pd.DataFrame()
    
    NN      = pd.DataFrame(Hits)
    k0      = NN.shape[0]

    SID = list(NN.dsstox_sid)
    
    T1 = pd.DataFrame(getChemToxSummary(SID,col=MDB[col],ds='%s.ds'%ds,cls='tox_fp'))
    if T1.shape[0]==0: return pd.DataFrame() 
    T1['study']=T1.tox_fp.apply(lambda x: x.split(':')[0])
    T1['organ']=T1.tox_fp.apply(lambda x: x.split(':')[1].lower())
    R1=T1.merge(NN,left_on='dsstox_sid',right_on='dsstox_sid')
    
    return R1

def getChemToxRefNNData(sid,s0=0.01,k0=10,fp='chm_mrgn',tox=None,pos_min=0,neg_min=0,
                        MDB=None,col='tox5_fp',ds_pos='tox_fpp1',ds_neg='tox_fpn1',   
                        sel_by=None,                            
                        filt_by=None,
                        dbg=False):
    target_sid = sid
    # Get NN
    Hits  = searchFP(sid,fp=fp,s0=s0,max_hits=k0,DB=MDB,sel_by=sel_by)
    if not Hits:
        return pd.DataFrame()
    
    NN      = pd.DataFrame(Hits)
    k0      = NN.shape[0]
    NN['d'] = 1-NN.jaccard
    SID = list(NN.dsstox_sid)
    
    # Load Tox Data including doses
    TX_pos = []
    TX_neg = []
    for X in MDB[col].find({'dsstox_sid':{'$in':SID}},dict(_id=0,dsstox_sid=1,name=1,tox_q1=1,tox_fpn1=1)):
        sid = X.get('dsstox_sid')
        if not sid: continue
        name= X['name']
        for x in X['tox_q1']:
            x['name']=name
            x['dsstox_sid']=sid

        TX_pos += X['tox_q1']
        Y = []
        for y in X['tox_fpn1']['ds']:
            Y.append(dict(dsstox_sid=sid,name=name,neg=1,effect=y))
        TX_neg += Y 

    if len(TX_pos)==0 or len(TX_neg)==0: return
    
    # Prepare the data for pivot
    TX_pos = pd.DataFrame(TX_pos)
    TX_pos['dose_w_unit'] = TX_pos.apply(lambda x: "%(dose)10.3f %(dose_unit)s" % dict(x),axis=1)
    TX_neg = pd.DataFrame(TX_neg)

    # Pivot the positives and negatives from Tox
    Pos=TX_pos.pivot_table(index='dsstox_sid',columns='effect',values='dose_w_unit',aggfunc=min)
    Neg=TX_neg.pivot_table(index='dsstox_sid',columns='effect',values='neg',aggfunc=min)

    # Fill the values for the common tox effects
    J=Pos.columns.intersection(Neg.columns)
    PJ=Pos.columns.difference(Neg.columns)
    NJ=Neg.columns.difference(Pos.columns)

    # Combine the common endpoints first
    P1 = Pos[J]
    N1 = Neg[J]

    # Remaining Positive
    P2 = Pos[PJ]
    N2 = Neg[NJ]

    # Label negatives
    N2 = N2.where(N2!=1,'no_effect')
    P1 = P1.where(N1!=1,'no_effect')
    TX = pd.merge(P1,P2,left_index=True,right_index=True)
    TX = TX.merge(N2,left_index=True,right_index=True)

    # Merge with NN data
    TX = TX.reset_index()
    TX = TX.merge(NN,left_on='dsstox_sid',right_on='dsstox_sid',how='outer').set_index(list(NN.columns))
    # Filter out columns based on filt_by
    if filt_by!=None and len(filt_by)>1:
        I1 = TX.columns.str.contains(filt_by,case=False)
        if np.sum(I1)==0: return
        TX = TX[TX.columns[I1]]

    if TX.shape[1]==0: return

    # Skip all chemicals without any data but keep the target
    TX['tox_n']=TX.notnull().sum(axis=1)
    TX = TX.query("tox_n>0 or jaccard==1").drop('tox_n',axis=1)
    TX.sort_index(level=4,inplace=True)
    TX=TX.where(TX.notnull(),other='no_data')
    TX['cls']=['target']+['analog']*(TX.shape[0]-1)
    ind = list(TX.index.names)+['cls']
    TX = TX.reset_index().set_index(ind)
    
    # Summarise row labels
    RL = pd.DataFrame([dict(zip(TX.index.names,x)) for x in TX.index.values])
    RL.rename(columns=dict(name='label'),inplace=True)
    RL['svg_url']= RL.dsstox_sid.apply(lambda x: "/api/genra/v3/viewChemGlyph/%s.svg" % x)

    # Summarize column labels
    # count the number of chemicals with dose, no effect and no data
    TX2= TX.replace('\s?\d+\\.\d+.+',value='dose',regex=True)
    TX2= TX2.replace('^\s+',value='',regex=True)
    CL = pd.DataFrame(TX2.apply(pd.value_counts).T.to_dict('records'),dtype=np.int)
    CL.fillna(0,inplace=True)
    # Change the counts to integer
    for c in CL: CL[c] = CL[c].astype(np.int)

    CL['label']=TX2.columns
    CL['study']=CL.label.apply(lambda x: x.split(':')[0])
    CL['organ']=CL.label.apply(lambda x: x.split(':')[1])

    # If negative and positive mins are defined then filter by them
    q = None
    if pos_min>0 or neg_min>0:
        q = "dose>=%d" % pos_min
        q += "and no_effect>=%d" % pos_min
        CL = CL.query(q)
        TX = TX[CL.label]
    
    # Rows 
    RD = TX.reset_index(drop=True).T
    RD.columns = RL.dsstox_sid
    
    return dict(RL=CL,CL=RL,D=TX,
                rows=RD.to_dict('index'),
                row_labs=CL.to_dict('records'),
                col_labs=RL.to_dict('records'),
                nn_opts=dict(target=target_sid,s0=s0,k0=k0,fp=fp))


