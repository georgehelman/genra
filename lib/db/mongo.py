import pymongo

def openMongo(host="pb.epa.gov",user='ishah',passwd='ishah',db='biospyder_v1'):
    con2 = pymongo.MongoClient("mongodb://%s:%s@%s/%s" % (user,passwd,host,db))
    DB = con2[db]
    return DB
