import pandas as pd
import numpy as np
import pymongo
import sklearn.metrics as metrics
from scipy.spatial.distance import squareform,pdist

from util.misc import *


def openMongo(host="galaxy",user='ishah',passwd='ishah',db='genra_v1'):
    con2 = pymongo.MongoClient("mongodb://%s:%s@%s/%s" % (user,passwd,host,db))
    DB = con2[db]
    return DB

def getChemData(cid,MDB=None):
    if not MDB: return
    CID = [cid]
    X_chm = getChemFP(CID,col=MDB.chm_fp,ds='mrgn.ds',fill=0)
    X_ct  = getChemFP(CID,col=MDB.chm_fp,ds='chmtp1.ds',fill=0)
    X_bio = getChemFP(CID,col=MDB.bio_fp,ds='bio1.ds',fill=0)
    X_tox = getChemFP(CID,col=MDB.tox5_fp,ds='tox_fpp1.ds',fill=0)

    if X_bio.shape[0]>0 and X_chm.shape[0]>0:
        X_bc  = pd.merge(X_bio,X_chm,left_index=True,right_index=True)
    else:
        X_bc = pd.DataFrame()
    
    if X_bio.shape[0]>0 and X_ct.shape[0]>0:
        X_bct = pd.merge(X_bio,X_ct,left_index=True,right_index=True)
    else:
        X_bct = pd.DataFrame()

    R = {}
    if X_chm.shape[0]>0: R['chm']=X_chm
    if X_bio.shape[0]>0: R['bio']=X_bio
    if X_ct.shape[0]>0: R['ct']=X_ct
    if X_bc.shape[0]>0: R['bc']=X_bc
    if X_bct.shape[0]>0: R['bct']=X_bct
    if X_tox.shape[0]>0: R['tox']=X_tox

    return R

def getChemNNDataSet(cid,k0=20,s0=0.01,MDB=None,dbg=False):
    if not MDB: return

    NN = getSimChems(cid,DB=MDB,k0=k0,s0=s0,drop_self=False)
    n  = None
    try:
        n = NN.shape[0]
    except:
        print "  no hits"
        return 
    
    if n==0: return
    
    CID = list(NN.index)
    
    if dbg: print "Chemicals: ",len(CID)

    X_chm = getChemFP(CID,col=MDB.chm_fp,ds='mrgn.ds',fill=0)
    X_ct  = getChemFP(CID,col=MDB.chm_fp,ds='chmtp1.ds',fill=0)
    X_bio = getChemFP(CID,col=MDB.bio_fp,ds='bio1.ds',fill=0)

    if X_bio.shape[0]>0 and X_chm.shape[0]>0:
        X_bc  = pd.merge(X_bio,X_chm,left_index=True,right_index=True)
    else:
        X_bc = pd.DataFrame()
    
    if X_bio.shape[0]>0 and X_ct.shape[0]>0:
        X_bct = pd.merge(X_bio,X_ct,left_index=True,right_index=True)
    else:
        X_bct = pd.DataFrame()
    
    X_toxp = getChemFP(CID,col=MDB.tox5_fp,ds='tox_fpp1.ds')
    X_toxn = getChemFP(CID,col=MDB.tox5_fp,ds='tox_fpn1.ds')
    X_tox  = X_toxp.copy()
    X_tox[X_toxn==1]=0

    R = {}
    if X_chm.shape[0]>0: R['chm']=X_chm
    if X_bio.shape[0]>0: R['bio']=X_bio
    if X_ct.shape[0]>0: R['ct']=X_ct
    if X_bc.shape[0]>0: R['bc']=X_bc
    if X_bct.shape[0]>0: R['bct']=X_bct
    if X_tox.shape[0]>0: R['tox']=X_tox

    return R

def getDataSetSims(DS,metric='jaccard'):

    R = {}
    for k,X in DS.iteritems():
        if X.shape[0]<=1: continue
        S = 1-pd.DataFrame(squareform(pdist(X,metric)),columns=X.index,index=X.index)
        R[k]=S

    return R
    

    
def getToxDataSet(tox,MDB=None):
    if not MDB: return
    Pos = np.array([i['dsstox_cid'] for i in MDB.tox5_fp.find({'tox_fpp1.ds':tox,'dsstox_cid':{'$exists':1}},{'_id':0,'dsstox_cid':1})])
    Neg = np.array([i['dsstox_cid'] for i in MDB.tox5_fp.find({'tox_fpn1.ds':tox,'dsstox_cid':{'$exists':1}},{'_id':0,'dsstox_cid':1})])
    n_n = len(Neg)
    n_p = len(Pos)
    CID = list(np.concatenate((Pos,Neg)))
    X_tox = pd.DataFrame(np.concatenate((np.ones(n_p),np.zeros(n_n))),index=CID,columns=[tox])
    X_chm = getChemFP(CID,col=MDB.chm_fp,ds='mrgn.ds',fill=0)
    X_ct  = getChemFP(CID,col=MDB.chm_fp,ds='chmtp1.ds',fill=0)
    X_bio = getChemFP(CID,col=MDB.bio_fp,ds='bio1.ds',fill=0)
    X_bc  = pd.merge(X_bio,X_chm,left_index=True,right_index=True)
    X_bct = pd.merge(X_bio,X_ct,left_index=True,right_index=True)
    
    return dict(tox=X_tox,chm=X_chm,ct=X_ct,bio=X_bio,bc=X_bc,bct=X_bct)

def searchSpaceByFP(cid,target_cid=None,s0=0.5,col='chm_fp',fpn='mrgn',i1=0,i2=None,dbg=False,DB=None,col_comp=None):
    Q0 = DB[col].find_one({'dsstox_cid':cid})
    if not Q0: return
    C = col_comp.find_one({'dsstox_cid':Q0['dsstox_cid']})
    Q = Q0[fpn]
    qmin= int(s0*Q['n'])
    qmax= int(Q['n']/s0)
    if dbg: print '>Query: {} {} {}'.format(C['name'],C['dsstox_cid'],C['smiles'])

    Match1 = {'%s.n'%fpn:{'$gte':qmin, '$lte':qmax},
              '%s.ds'%fpn:{'$in':Q['ds']}
              }
    if target_cid:
        Match1['dsstox_cid']={'$in':target_cid}

    Agg = [
        {'$match': Match1},
        {'$project': 
             {'jaccard': 
                 {'$let':
                  {'vars': 
                   {'olap': {'$size':{'$setIntersection': ['$%s.ds'%fpn,Q['ds']] }}},
                   'in': {'$divide':['$$olap',
                                     {'$subtract': [{'$add':[Q['n'],'$%s.n'%fpn]},'$$olap'] }]}
                  }
                 },
              '_id':0,
              'dsstox_cid':1,
              'casrn':1,
              'name':1,
             }
        },
        {'$match':{'jaccard':{'$gte':s0}}},
        {'$sort': {'jaccard':-1}}
    ]

    #print qmin,qmax
    if pymongo.version[0]=='2':
        return DB[col].aggregate(Agg,allowDiskUse=True)['result']
    else:
        return DB[col].aggregate(Agg,allowDiskUse=True)


def getChemBioactivity(CID,bio=None):
    Agg = [
            # Match chemicals in cluster
            {'$match': {
                     'dsstox_cid':{'$in':CID}}
            },
            # Include these fields
            {'$project':{'dsstox_cid':1,'name':1,'_id':0,
                        'fp':'$bio1.ds'},
            },
            # Unwind the fp 
            {'$unwind':"$fp"}
            ]
    if bio: Agg.append({'$match': {'fp':{'$in': bio}}})
    
    X = DB.bio_fp.aggregate(Agg)

    R = pd.DataFrame(X['result'])
    if R.shape[0]==0: return R

    X = pd.pivot_table(R,index=['dsstox_cid'],columns='fp',values='name',aggfunc=len,fill_value=0)
    CID0 = list(set(CID).difference(X.index))
    X1 = np.zeros((len(CID0),X.shape[1]))
    X1[:,:]=None
    return pd.concat((X,pd.DataFrame(TX1,columns=TX2.columns,index=CID0)))



def getChemFP(CID,col=None,ds=None,fp=None,fill=None):
    Agg = [
            # Match chemicals in cluster
            {'$match': {
                     'dsstox_cid':{'$in':CID}}
            },
            # Include these fields
            {'$project':{'dsstox_cid':1,'name':1,'_id':0,
                        'fp':'$'+ds},
            },
            # Unwind the fp 
            {'$unwind':"$fp"}
            ]
    if fp: Agg.append({'$match': {'fp':{'$in': fp}}})
    
    X = col.aggregate(Agg,allowDiskUse=True)
    if not X: return
    try:
        R = pd.DataFrame(X['result'])
    except:
        R = pd.DataFrame(list(X))

    if R.shape[0]==0 or R.shape[1]==0: return pd.DataFrame()

    return pd.pivot_table(R,index=['dsstox_cid'],columns='fp',values='name',aggfunc=len,fill_value=fill)

def getChemTox(CID,col=None,ds=None,fp=None,fill=None):
    Agg = [
            # Match chemicals in cluster
            {'$match': {
                     'dsstox_cid':{'$in':CID}}
            },
            # Include these fields
            {'$project':{'dsstox_cid':1,'name':1,'_id':0,
                        'fp':'$'+ds},
            },
            # Unwind the fp 
            {'$unwind':"$fp"}
            ]
    if fp: Agg.append({'$match': {'fp':{'$in': fp}}})
    
    X = col.aggregate(Agg)
    R = pd.DataFrame(X['result'])
    
    if R.shape[0]==0: return R

    X = pd.pivot_table(R,index=['dsstox_cid'],columns='fp',values='name',aggfunc=len,fill_value=fill)
    TX2 = X.copy()
    ST = list(set([i.split(':')[0] for i in X.columns]))
    for st in ST:
        EFF = [i for i in X if i.startswith(st)]
        St_done = (X[EFF].notnull()).sum(axis=1)>0
 
        for E in EFF:
            Y = X[E]
            I = np.logical_and(Y.isnull(),St_done)
            Y[I] = 0
            TX2[E]=Y

    CID0 = list(set(CID).difference(TX2.index))
    TX1 = np.zeros((len(CID0),TX2.shape[1]))
    TX1[:,:]=None
    TX1 = pd.DataFrame(TX1,columns=TX2.columns,index=CID0)

    
    return pd.concat((TX2,TX1))

def getKNN(cid,Sim,k0=10,s0=0,sim=False,drop_self=True):
    S_i = Sim.ix[cid,:]
    if drop_self:
        S_i = S_i.drop(cid)
    NN  = None
    
    if k0 and s0:
        S_i = S_i[S_i>s0]
        S_i.sort()
        NN = S_i[-k0:]
    elif k0:
        S_i.sort()
        NN = S_i[-k0:]
    elif s0:
        NN = S_i[S_i>s0]
    
    if sim:
        return NN
    else:
        return NN.index
    

def getSimChems(cid,k0=10,s0=0.05,sim=False,drop_self=True,col='chm_fp',DB=None):
    if not DB: return
    H = searchSpaceByFP(cid,s0=s0,DB=DB,col_comp=DB['compound'])
    if not H: return
    H = list(H)
    if len(H)>k0:
        H = H[:k0]
    H = pd.DataFrame(H)
    H.set_index('dsstox_cid',inplace=True)

    if drop_self:
        H.drop(cid,axis=0,inplace=True)
    
    return H
        

def GenRAPerf(ba,Data=None,Sc=None,sim_type='chm',k0=5,s0=0.5,
              ret='perf',cl=None,
              wt=True,perm=0):
    """
    Systematically evaluate the automated read-across prediction of activities of chemical based on 
    nearest neighbours. Activity is inferred based on distance-weighted similarity and performance is evaluated
    using leave-one-out validation across the neighbours

    Inputs:
      ba  - the input descriptor/column of Data that is being predicted 
      Data- Full bioactivity + tox data matrix
      Sc  - Similarity matrix
      sim_type - 'chm','bio','bc' (not used in calc only returned in results)
      k0  - the maximum number of nearest neighbours to use
      s0  - the miminim similarity threshold

    Output: 
      A) ret == 'perf' returns dictionary with the performance summary:

      n_pos     - the number of positive chemicals in the entire neighbourhood
      n_neg     - the number of negative chemicals in the entire neightbourhood
      n_sim_pos - the number of positives for s0
      n_sim_neg - the number of negatives for s0 
      k0_e      - the actual number of nn used after s0 is applied (k0_e = n_sim_pos+n_sim_neg) 
      effect    - the same as ba input
      auc       - the  area under the curve based on a comparison of  
      f1_max    - the max f1 score across the entire ROC curve
      ba_max    - the ba_max corresponding to the f1_max
      sn_max    - the sn_max corresponding to the f1_max
      sp_max    - the sp_max corresponding to the f1_max
      sim_type  - the same as the input sim_type
      s_min     - the minimum value of neighbourhood similarity 
      s_max     - the maximum value of neighbourhood similarity 

      B) ret == 'all' return a list containing three items:-
      i)   The same dictionary containing the performance summary returned by ret=='perf'
      ii)  A dataframe containing the predictions for each chemical in the neighbourhood
      iii) A dataframe containing the ROC data 

    Usage:
      Find nearest neighbours within k0 and s0 using similarity metric
      iterate over each biactivity assay and infer activity


    """
    A = Data[ba]
    A = A[pd.notnull(A)]
    if len(A) < 2: return {}

    CID = [i for i in A.index]
    A_p = A.copy()
    A_t = A.copy()
    A_t[A>0]=1
    S_min = []    
    S_max = []
    S_mn = []
    N_pos= []
    N_neg= []
    n_sim_pos=0
    n_sim_neg=0
    CID1=[]
    for cid in A.index:
        S_i  = getKNN(cid,Sc.ix[CID,CID],k0=k0,s0=s0,sim=True)
        if S_i.shape[0]==0: continue
        CID1.append(cid)
        if A[cid]>0: 
            n_sim_pos += 1
        else:
            n_sim_neg += 1
        S_min.append(S_i.min())
        S_mn.append(S_i.mean())
        S_max.append(S_i.max())
        A_i  = A_t[S_i.index]
        if S_i.sum()==0:
            A_p[cid] = 0
        else:
            A_p[cid] = np.sum(A_i*S_i)/S_i.sum()

        # How many positives were there in the nearest neighbors
        N_pos.append(np.sum(A_i>0))
        N_neg.append(np.sum(A_i==0))
        #print cid,n_sim_pos,n_sim_neg
    
    if len(CID1)==0: return {}

    Res_chm = pd.DataFrame(zip(CID1,A_t,A_p,N_pos,N_neg,S_min,S_mn,S_max),
                           columns=['ID','a_t','a_p','n_pos','n_neg','s_min','s_ave','s_max']
                           )
    Res_chm['k0']=k0
    Res_chm['s0']=s0
    Res_chm['effect']=ba

    # Is this zero ? 
    # A_p.ix[pd.isnull(A_p)]=0.0

    # t0 is sometimes>1 but it makes no sense
    fpr,tpr,t0 = metrics.roc_curve(A_t,A_p,pos_label=1)
    tnr = 1-fpr
    Roc = pd.DataFrame(zip(fpr,tnr,tpr,t0),columns=['fpr','sp','sn','t0'])
    Roc.t0[Roc.t0>1]=1.0
    # This is what we tend to use but it's not that useful when sp or sn are 0
    Roc['ba']=0.5*(Roc.sn+Roc.sp)
    Roc['f1']= 2*(Roc.sn*Roc.sp)/(Roc.sn+Roc.sp)
    
    # Best sn,sp,f1
    #ba_max = Roc.ba.max()
    #P = Roc.ix[Roc.ba==ba_max,:]
    f1_max = Roc.f1.max()
    P = Roc.ix[Roc.f1==f1_max,:]
    
    P.rename(columns=dict(sn='sn_max',sp='sp_max',ba='ba_max',
                          f1='f1_max',
                          t0='t0_max',fpr='fpr_max'),inplace=True)
    # ROC AUC
    #print "new"
    if P.shape[0]>0:
        P_max = P.reset_index(drop=True).ix[0,:].to_dict()    
        auc = metrics.auc(fpr,tpr)
    else:
        P_max = {}
        auc   = None
        
    P_max['t_mdn']=np.median(A_p)
    P_max['t_mn']=np.min(A_p)
    P_max['t_mx']=np.max(A_p)

    Res = dict(auc=auc,effect=ba,k0=k0,s0=s0,n_pos=np.sum(A_t==1),n_neg=np.sum(A_t==0),
               n_sim_pos=n_sim_pos,n_sim_neg=n_sim_neg,k0_e=n_sim_pos+n_sim_neg,
               s_min=min(S_min),s_max=max(S_max),
               sim_type=sim_type,
               n_tot=len(A_t),n_cl=Data.shape[0]
               )
    
    # If there are no negatives in the vicinity i.e. all positives then 
    # predicted class should be positive so set auc = 1 
    if Res['n_sim_neg']==0 and Res['n_sim_pos']>0 and not auc:
        Res['auc'] = 1 
        perm= 0
    
    if cl!=None: Res['cl']=cl

    Res.update(P_max)

    if auc and perm: 
        Res['auc_pval']=permuteAUC(A_t,A_p,N=perm)

    if ret=='perf':
        return Res
    elif ret=='all':
        return Res,Res_chm,Roc


def permuteAUC(Y_t,Y_p,N=100,pos=1):
    AUC = []
    fpr,tpr,t0 = metrics.roc_curve(Y_t,Y_p,pos_label=pos)
    auc = metrics.auc(fpr,tpr)

    for i in range(N):
        Y_r = np.array(Y_t.copy())
        np.random.shuffle(Y_r)
        fpr,tpr,t0 = metrics.roc_curve(Y_r,Y_p,pos_label=pos)
        AUC.append(metrics.auc(fpr,tpr))

    p_val = 1.0*np.sum(np.array(AUC)>auc)/N

    return p_val

def GenRAPred(Q,ba,Data=None,Sc=None,sim_type='chm',wt=True,
              k0=5,s0=0.5,t0=1):

    """
    Predict activities of a chemical based on distance-weighted similarity.

    Inputs:
      qcid- Query cid for chemical
      ba  - the input descriptor/column of Data that is being predicted 
      Data- Full bioactivity + tox data matrix
      Sc  - Similarity matrix
      sim_type - 'chm','bio','bc' (not used in calc only returned in results)
      k0  - the maximum number of nearest neighbours to use
      s0  - the miminim similarity 
      t0  - the activity score threshold for 1 or 0 activity
    """

    A = Data[ba]
    A = A[pd.notnull(A)]
    CID = list(A.index)
    N_pos= []
    N_neg= []
    n_sim_pos=0
    n_sim_neg=0
    CID1=[]
    A_t = []
    A_s = []

    for cid in Q:
        a_s=None
        if cid not in Sc.index: 
            print cid," is not part of the similarity matrix"
            continue

        # If this is a chemical with known bioactivity ...
        if cid in CID: 
            A_t.append(A.ix[cid])
        else:
            A_t.append(None)
        
        # Search the neighbourhood of cid in the context of neighbours with known activities
        CIDi = [cid]+CID
        S_i  = getKNN(cid,Sc.ix[CIDi,CIDi],k0=k0,s0=s0,sim=True)
        if S_i.shape[0]==0: continue

        A_i  = A[S_i.index]
        if wt:
            if S_i.sum()==0:
                a_s = 0
            else:
                a_s = np.sum(A_i*S_i)/S_i.sum()
        else:
            a_s = np.sum(A_i)/len(A_i)
        # How many positives were there in the nearest neighbors
        N_pos.append(np.sum(A_i>0))
        N_neg.append(np.sum(A_i==0))
        A_s.append(a_s)
        
    #return A_s,CID1,A_t
    Res_chm = pd.DataFrame(zip(Q,A_t,A_s,N_pos,N_neg),
                           columns=['ID','a_t','a_s','n_pos','n_neg']
                           )
    Res_chm['k0']=k0
    Res_chm['s0']=s0
    Res_chm['t0']=t0
    Res_chm['effect']=ba
    Res_chm['a_p']=Res_chm.a_s
    Res_chm['sim_type']=sim_type
    Res_chm.a_p[Res_chm.a_s>=t0]=1
    Res_chm.a_p[Res_chm.a_s<t0]=0

    return Res_chm[['ID','effect','sim_type','k0','s0','t0','a_t','a_p','a_s','n_pos','n_neg']]
