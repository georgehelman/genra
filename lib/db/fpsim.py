import numpy as np
import pandas as pd
import copy
import pymongo
import sys
import re

from db.mongo import *

def searchCollByFP(sid,s0=0.5,col='chm_fp',fpn='mrgn',
                   SID=None,
                   i1=0,i2=None,dbg=False,DB=None,
                   max_hits=100):
    Q0 = DB[col].find_one({'dsstox_sid':sid})
    if not Q0: return
    C = DB.compound.find_one({'dsstox_sid':Q0['dsstox_sid']})
    Q = Q0[fpn]
    qmin= int(s0*Q['n'])
    qmax= int(Q['n']/s0)
    if dbg: print '>Query: {} {} {}'.format(C['name'],C['dsstox_sid'],C['smiles'])

    Match1 = {'%s.n'%fpn:{'$gte':qmin, '$lte':qmax},
                    '%s.ds'%fpn:{'$in':Q['ds']}
              }
    if SID: Match1['dsstox_sid']={'$in':SID}

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
              'dsstox_sid':1,
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


def searchFP(sid,fp='chm_mrgn',DB=None,sel_by=None,**kwargs):
    """
    sid: dsstox_sid
    fp: is the fingerprint name that matches COL
    sel_by:  Find similar chemicals that also have data in another data collection. 
             e.g. If sel_by is tox_txrf then search results are limited to the chemicals 
             with toxicity data 
    """

    col,ds = getColFPMap(fp)
    if not (ds and col): return

    SID_h= None
    if sel_by and COL.has_key(sel_by):
        col_sel =COL.get(sel_by)    
        SID_h = DB[col_sel].find({'dsstox_sid':{'$exists':1}}).distinct('dsstox_sid')
        if sid not in SID_h: SID_h.append(sid)

    # Find NN
    H  = searchCollByFP(sid=sid,col=col,fpn=ds,DB=DB,SID=SID_h,**kwargs)
    if not H: H=[] 
    
    Q = [i for i in H if i.get('dsstox_sid')==sid]
    H1= sorted([i for i in H if i.get('dsstox_sid')!=sid],key=lambda x: x['jaccard'],reverse=True)
    
    return Q+H1
