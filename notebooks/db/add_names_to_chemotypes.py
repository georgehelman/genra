import pymongo

mongocon=pymongo.MongoClient("mongodb://ghelman:ghelman@pb.epa.gov/genra_dev_v4")
DB=mongocon['genra_dev_v4']
dsstox=DB['compounds']
chemotypes_coll=DB['chemotypes']

for record in list(chemotypes_coll.find({'dsstox_sid':{'$exists':True}})):
    sid=record['dsstox_sid']
    name=dsstox.find_one({'dsstox_sid':sid}).get('name',None)
    chemotypes_coll.update({'dsstox_sid':sid},{'$set':{'name':name}})
