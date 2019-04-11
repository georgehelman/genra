import os
import sys

TOP = '/'.join(os.getcwd().split('/')[:-1])+'/'
LIB = TOP+'lib'
if not LIB in sys.path: 
    sys.path.insert(0,LIB)

from utl.queries import get_invitrodb
import pymongo
from db.etl import *

mysql_cnx = mysql.connector.connect(option_files='/share/home/ghelman/.my.cnf')

mongocon = pymongo.MongoClient("mongodb://ghelman:ghelman@pb.epa.gov/genra_dev_v4")
DB = mongocon['genra_dev_v4']
invitrodb=DB['invitrodb']

current_invitro=query_mongo(invitrodb,{},{'_id':0,'assays.l5id':1})
l5ids=[str(l5id_field['l5id']) for l in current_invitro for l5id_field in l['assays']]
invitro_query=get_invitrodb(l5ids)
cursor=mysql_cnx.cursor(dictionary=True,buffered=True)
cursor.execute(invitro_query)

for assay in cursor:
	sid=assay['dsstox_sid']
	del assay['dsstox_sid']
	invitrodb.update_one({'dsstox_sid':sid},{'$push':{'assays':assay}},upsert=True)


