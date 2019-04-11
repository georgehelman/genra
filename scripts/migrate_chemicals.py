import sys
sys.path.append('/share/home/ghelman/workspace/anaconda2/pkgs')
sys.path.append('/share/home/ghelman/dev/cronjobs/lib')
from datetime import datetime
import pymongo
import mysql.connector
import data_migration

mongocon = pymongo.MongoClient("mongodb://ghelman:ghelman@pb.epa.gov/ghelman")
DB_dev = mongocon['ghelman']
comps=DB_dev['compounds']
logs=DB_dev['logs']
tox=DB_dev['toxval']
meta=DB_dev['meta']
cnx = mysql.connector.connect(option_files='/share/home/ghelman/.my.cnf')

last_update=meta.find_one(sort=[("date",-1)])['date']
datetime.strftime(last_update,"%Y-%m-%d %H:%M:%S")

f=open('/share/home/ghelman/dev/cronjobs/import.log','w+',1)
data_migration.get_new_chemicals(cnx,comps,tox,logs,meta,last_update,f)
data_migration.get_new_predictions(cnx,comps,logs,meta,last_update)
data_migration.get_new_measurements(cnx,comps,logs,meta,last_update)
