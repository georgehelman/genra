import pandas as pd
import datetime
import sys
import os
import pymongo
import mysql.connector

TOP = '/'.join(os.getcwd().split('/')[:-1])+'/'
LIB = TOP+'lib'
if not LIB in sys.path: 
    sys.path.insert(0,LIB)

from utl.queries import *

mysql_cnx = mysql.connector.connect(option_files='/share/home/ghelman/.my.cnf',database='ro_stg_dsstox')
curs=mysql_cnx.cursor(dictionary=True)
list_query=get_lists()
curs.execute(list_query)
results=curs.fetchall()
curs.close()

list_info=pd.DataFrame(results)

list_info['source_data_updated_at']=[datetime.datetime.combine(d,datetime.time.min) for d in list_info['source_data_updated_at']]
#Pymongo cannot encode date, needs datetime

metadata_columns=list(list_info.columns)
metadata_columns.remove('dsstox_sid')

chemical_lists=[]
for metadata, list_df in list_info.groupby(metadata_columns,as_index=True):
    list_dict={k:v for k,v in zip(metadata_columns,metadata)}
    list_dict['n']=len(list_df)
    list_dict['chemicals']=list(list_df['dsstox_sid'])
    chemical_lists.append(list_dict)

mongocon = pymongo.MongoClient("mongodb://ghelman:ghelman@pb.epa.gov/genra_dev_v4")
DB = mongocon['genra_dev_v4']
dsstox_lists=DB['chemical_lists']

dsstox_lists.drop()
dsstox_lists.insert_many(chemical_lists)
