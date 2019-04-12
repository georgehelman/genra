from __future__ import division
import pymongo

con = pymongo.MongoClient("mongodb://ghelman:ghelman@pb.epa.gov/genra_dev_v4")
DB = con['genra_dev_v4']

def getColFPMap(fp):
    COL = dict(chm_mrgn='chm_fp',chm_httr='chm_fp',chm_ct='chemotypes',
               bio_txct='bio_fp',bio_tx21='tox21_fp',
               toxp_txrf='tox_fp',toxn_txrf='tox_fp')

    DS  = dict(chm_mrgn='mrgn',chm_httr='httr',chm_ct='chemotypes',
               bio_txct='bio1',bio_tx21='t211',
               toxp_txrf='tox_fpp1',toxn_txrf='tox_fpn1'
               )

    col=COL.get(fp)
    ds =DS.get(fp)

    return col,ds

def getfp(SID,fp='chm_mrgn'):
    col, ds = getColFPMap(fp)
    Agg = [
        # Match chemicals in cluster
        {'$match': {
            'dsstox_sid': {'$in': SID}}
        },
        # Include these fields
        {'$project': {'dsstox_sid': 1,
                      'fp': '$%s.ds' % ds},
         },
    ]

    X = DB[col].aggregate(Agg, allowDiskUse=True)

    return {record['dsstox_sid']:record['fp'] for record in X}

def jaccard(l1,l2):
    num = set(l1) & set(l2)
    denom = set(l1) | set(l2)
    return len(num)/len(denom)

def radial_sims(target,neighbors,fp='chm_mrgn'):
    #neighbors needs to already be ordered properly
    sids=neighbors[:]
    sids.append(target)
    fps=getfp(sids,fp=fp)
    sims=[jaccard(fps[target],fps[neighbor]) for neighbor in neighbors]
    return sims