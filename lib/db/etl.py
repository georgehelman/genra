from datetime import datetime
import mysql.connector
from collections import defaultdict
import pymongo
import jsonpatch

def query_mysql(connection,query):
    cursor = connection.cursor(dictionary=True)
    cursor.execute(query)
    ret = cursor.fetchall()
    cursor.close()
    if len(ret)==1:
        ret=ret[0]
    return ret

def query_mongo(db,query_dict,result_dict,limit=0):
    #Result dict specifies desired return fields
    # limit=0 equivalent to no limit
    return list(db.find(query_dict, result_dict).limit(limit))


def load_record(database,record,log_db=None):
    database.insert_one(record)
    if log_db:
        cid=record['dsstox_cid']
        load_record(log_db,{'dsstox_cid':cid})

def load_records(database,records): #records should be list of records
    if len(records)>0:
        database.insert_many(records)

def update_record(database,match,update,upsert=True):
    database.update_one(match,update,upsert=upsert)

def replace_record(database,match,record,upsert=True):
    database.replace_one(match,record,upsert)

def log_patch(db,patch,cid):
    time = datetime.utcnow().replace(microsecond=0)
    patch={'timestamp':time,
           'patch':patch}
    db.update_one({'dsstox_cid':cid},{'$push':{'patches':patch}})

#Helper functions
def keyify(str):
    return str.replace('.','-')

def default_dictify(query_results,field_name):
    default_dict=defaultdict(list)
    for row in query_results:
        default_dict[keyify(row['name'])].append(row['result_value'])
    return {field_name: default_dict}

def extract_field(list_of_dicts,field):
    if type(list_of_dicts) is dict:
        values=list_of_dicts.get(field,None)
    else:
        values = [c.get(field,None) for c in list_of_dicts]
    return values
