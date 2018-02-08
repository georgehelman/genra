import pymongo
import os

def openMongo(host="pb.epa.gov",user=None,passwd=None,db=None):
    if not user or passwd:
        user,passwd = file(os.getenv('HOME')+'/.mngdb/passwd').read().strip().split(':')
        
    con2 = pymongo.MongoClient("mongodb://%s:%s@%s/%s" % (user,passwd,host,db))
    DB = con2[db]
    return DB

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