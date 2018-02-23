import pymongo

mongocon=pymongo.MongoClient("mongodb://ghelman:ghelman@pb.epa.gov/genra_dev_v4")
DB=mongocon['genra_dev_v4']
chemotypes_coll=DB['chemotypes']

for record in chemotypes_coll.find({}):
    chemotypes=record['chemotypes']
    record['chemotypes_temp']=record['chemotypes']
    del record['chemotypes']
    record['chemotypes']={'n':len(chemotypes),'ds':chemotypes}
    del record['chemotypes_temp']
    DB['schemotypes'].insert(record)
