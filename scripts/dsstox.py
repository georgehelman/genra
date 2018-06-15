import os
import sys
import pymongo
TOP = '/'.join(os.getcwd().split('/')[:-1])+'/'
LIB = TOP+'lib'
sys.path.insert(0,LIB)
from utl.queries import get_new_chemicals,get_measured_properties,get_predicted_properties,get_synonyms
from db.etl import *

mysql_cnx = mysql.connector.connect(option_files='/share/home/ghelman/.my.cnf')

mongocon = pymongo.MongoClient("mongodb://ghelman:ghelman@pb.epa.gov/ghelman")
DB_dev = mongocon['ghelman']
compounds=DB_dev['compounds']
physprop=DB_dev['physprop']
logs=DB_dev['logs']


current_compounds=query_mongo(compounds,{},{'_id':0,'dsstox_sid':1})
current_sids=extract_field(current_compounds,'dsstox_sid')
new_compounds_query=get_new_chemicals(current_sids)
new_compounds=query_mysql(mysql_cnx,new_compounds_query)
new_physprop=[]
date=datetime.now()
for compound in new_compounds:
    cid=compound['dsstox_cid']
    sid=compound['dsstox_sid']
    compound.update({'created_at': date})
    compound.update({'updated_at': date})
    get_synonyms_query=get_synonyms(sid)
    synonyms=query_mysql(mysql_cnx,get_synonyms_query)
    compound.update({'synonyms': synonyms})
    get_measured_query=get_measured_properties(sid)
    get_predicted_query=get_predicted_properties(sid)
    measured_properties=query_mysql(mysql_cnx,get_measured_query)
    predicted_properties=query_mysql(mysql_cnx,get_predicted_query)
    new_physprop.append({'dsstox_sid':sid,'dsstox_cid':cid,'measured_properties':measured_properties,'predicted_properties':predicted_properties,'created_at':date,'updated_at':date})

load_records(compounds,new_compounds)
load_records(physprop,new_physprop)

with open('logfile','w') as f:
    f.writelines(new_compounds)
