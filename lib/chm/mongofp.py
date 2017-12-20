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

#con = pymongo.Connection("mongodb://ishah:ishah@localhost/txbrn_v1")
#DB = con['txbrn_v1']

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
    return dict(((k,v) for k,v in X.iteritems() if v!=0 or v!=None))

def makeFP(cid,save=False,MDB=None,col_chem_fp='chm_fp'):
    dsstox_cid = cid
    if not MDB: return
    C = MDB.compound.find_one(dict(dsstox_cid=cid),dict(_id=0))
    M = Chem.MolFromSmiles(str(C['smiles']))
    if not (dsstox_cid and M): return

    FPT = dict(httr=lambda i: AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(i),
               mrgn=lambda i: AllChem.GetMorganFingerprintAsBitVect(i,3,2048))


    FP = MDB.chm_fp.find_one({'dsstox_cid':dsstox_cid})
    if FP and replace:
        MDB.chm_fp.remove(FP)

    FP = {'dsstox_cid':dsstox_cid,
          'name':C['name'],
          'dsstox_sid':C['dsstox_sid'],
          'casrn':C['casrn']}
    
    for fpn,fp_func in FPT.iteritems():
        Y = {}
        V = fp_func(M)
        if V: 
            Y['ds']= list(np.core.defchararray.add(fpn+'_',np.where(V)[0].astype('string')))
            Y['n'] = len(Y['ds'])
            FP[fpn]=copy.copy(Y)
    
    if save:
        MDB.chem_fp.insert(FP)
    else:
        return FP
    
def getPerms(L=2048, n_p=100):
    return map(lambda _: np.random.permutation(L), range(n_p))

def FP2Vect(C,fp='mrgn'):
    V1 = C[fp].get('ds',None)
    if not V1: return
    V2 = [int(i.split('_')[1]) for i in V1]
    
    return V2

def getMinHash(V, perms):
    min_hash = []
    for p in perms:
        for idx, i in enumerate(p):
            if i in V:
                min_hash.append(idx)
                break            
    return min_hash


def hashToBuckets(min_hash, num_buckets=25, n_bits=2048):
    if len(min_hash) % num_buckets:
        raise Exception('number of buckets must be divisiable by the hash length')
    Buckets = []
    hash_per_bucket = int(len(min_hash) / num_buckets)
    num_bits = (n_bits-1).bit_length()
    if num_bits * hash_per_bucket > sys.maxint.bit_length():
        raise Exception('numbers are too large to produce valid buckets')
    for b in range(num_buckets):
        Buckets.append(reduce(lambda x,y: (x << num_bits) + y, 
                              min_hash[b:(b + hash_per_bucket)]))
    return Buckets


def searchFP(cid,s0=0.5,fpn='mrgn',i1=0,i2=None,dbg=False):
    if not i2: i2 = DB.chem_fp.count()
    Q0 = DB.chem_fp.find_one({'dsstox_cid':cid})
    if not Q0: return
    C = DB.compound.find_one({'dsstox_cid':Q0['dsstox_cid']})
    Q = Q0[fpn]
    qmin= int(s0*Q['n'])
    qmax= int(Q['n']/s0)
    if dbg: print '>Query: {} {} {}'.format(C['name'],C['dsstox_cid'],C['smiles'])

    Agg = [
        {'$match': {'mrgn.n':{'$gte':qmin, '$lte':qmax},
                    'mrgn.ds':{'$in':Q['ds']},
                    'i':{'$gte':i1, '$lte':i2}
                   }},
        {'$project': 
             {'tanimoto': 
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
        {'$match':{'tanimoto':{'$gte':s0}}},
        {'$sort': {'tanimoto':-1}}
    ]

    #print qmin,qmax
    return DB.chem_fp.aggregate(Agg)['result']

def searchCollByFP(cid,s0=0.5,col='chm_fp',fpn='mrgn',i1=0,i2=None,dbg=False,DB=None,
                   max_hits=100):
    Q0 = DB[col].find_one({'dsstox_cid':cid})
    if not Q0: return
    C = DB.compound.find_one({'dsstox_cid':Q0['dsstox_cid']})
    Q = Q0[fpn]
    qmin= int(s0*Q['n'])
    qmax= int(Q['n']/s0)
    if dbg: print '>Query: {} {} {}'.format(C['name'],C['dsstox_cid'],C['smiles'])

    Agg = [
        {'$match': {'%s.n'%fpn:{'$gte':qmin, '$lte':qmax},
                    '%s.ds'%fpn:{'$in':Q['ds']}
                   }},
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
        {'$sort': {'jaccard':-1}},
        {'$limit': max_hits}
    ]

    #print qmin,qmax
    return list(DB[col].aggregate(Agg))

def searchParallelFP(cid,s0=0.5,fpn='mrgn',n_par=40,dbg=False,par_mach=None):
    N = DB.chem_fp.find({'%s.n'%fpn:{'$gt':0}}).count()
    dn= int(N/n_par)
    I = [[i,i+dn] for i in range(0,N,dn)]
    
    X = par_mach.map(lambda (i,j): searchFP(cid,s0=s0,fpn=fpn,i1=i,i2=j), I)
    
    return sorted(reduce(lambda x,y: x+y, X),key=lambda x: x['jaccard'],reverse=True)
