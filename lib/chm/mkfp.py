from rdkit import Chem
from rdkit.DataStructs import *
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem
from rdkit.Chem.SaltRemover import SaltRemover
import numpy as np
import pandas as pd
import copy
import pymongo
import sys

def makeFPsForChems(CID,save=False):
    
    FP = []
    for C in CID: 
        FP.append(makeFP(C))

    return FP

def storeFPsForChems(CID,remove=True,DB_col=None):
    FP = makeFPsForChems(CID)
    r  = None
    if remove:
        #DB.chem_fp.remove({'dsstox_cid':{'$in':CID}})
        DB_col.insert([i for i in FP if i])
    else:
        for fp in FP:
            DB_col.update({'dsstox_cid':fp['dsstox_cid']},fp,dict(upsert=1))


    
def skipNulls(X):
    return {k:v for k,v in X.iteritems() if v and v==v}

def makeFP(dtxsid,save=False,replace=True,col_comp=None,col_chm_fp=None):
    if not (col_comp and col_chm_fp): return
    C = col_comp.find_one(dict(dsstox_sid=dtxsid),dict(_id=0))
    M = Chem.MolFromSmiles(str(C['smiles']))
    if not (dtxsid and M): return

    FPT = dict(httr=lambda i: AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(i),
               mrgn=lambda i: AllChem.GetMorganFingerprintAsBitVect(i,3,2048))


    FP = col_chm_fp.find_one({'dsstox_sid':dtxsid})
    if FP and not replace: return
    if FP: col_chm_fp.delete_one({'dsstox_sid':dtxsid})
        
    FP = {'dsstox_sid':dtxsid,
          'name':C['name'],
          'dsstox_cid':C['dsstox_cid'],
          'casrn':C['casrn']}
    
    for fpn,fp_func in FPT.iteritems():
        Y = {}
        V = fp_func(M)
        if V: 
            Y['ds']= list(np.core.defchararray.add(fpn+'_',np.where(V)[0].astype('string')))
            Y['n'] = len(Y['ds'])
            FP[fpn]=copy.copy(Y)
    
    if save:
        col_chm_fp.insert_one(FP)
    else:
        return FP
    
