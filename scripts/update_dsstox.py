import os
import sys
TOP = '/'.join(os.getcwd().split('/')[:-1])+'/'
LIB = TOP+'lib'
if not LIB in sys.path: 
    sys.path.insert(0,LIB)

from utl.queries import get_new_chemicals,get_measured_properties,get_predicted_properties,get_synonyms
import pymongo
from db.etl import *
from db.vis import smiles2svg
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

mysql_cnx = mysql.connector.connect(option_files='/share/home/ghelman/.my.cnf')

mongocon = pymongo.MongoClient("mongodb://ghelman:ghelman@pb.epa.gov/genra_dev_v4")
DB = mongocon['genra_dev_v4']
compounds=DB['compound']
physprop=DB['physprop']
fingerprints=DB['chm_fp']
chemotypes=DB['chemotypes']
logs=DB['logs']


current_compounds=query_mongo(compounds,{},{'_id':0,'dsstox_sid':1})
current_sids=extract_field(current_compounds,'dsstox_sid')
new_compounds_query=get_new_chemicals(current_sids)
new_compounds=query_mysql(mysql_cnx,new_compounds_query)
new_physprop=[]
new_fps=[]
date=datetime.now()
for compound in new_compounds:
    cid=compound['dsstox_cid']
    sid=compound['dsstox_sid']
    smiles=compound['smiles']
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
    viz=smiles2svg(smiles)
    compound.update({'viz':viz})
    #fps
    httr=AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(Chem.MolFromSmiles(smiles))
    httr=list(np.core.defchararray.add('httr_',np.where(httr)[0].astype('string')))
    mrgn=AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles),3,2048)
    mrgn=list(np.core.defchararray.add('mrgn_',np.where(mrgn)[0].astype('string')))
    new_fps.append({'dsstox_sid':sid,'mrgn':{'ds':mrgn,'n':len(mrgn)},'httr':{'ds':httr,'n':len(httr)}})

#Chemotypes
from subprocess import call

corina = '/opt/CORINA_Symphony/CORINA_Symphony_14698/bin/moses'
smile_file = '/share/home/ghelman/dev/read_across/applied/fingerprints/chemotypes/CORINA_Symphony/compounds.smi'
output = '/share/home/ghelman/dev/read_across/applied/fingerprints/chemotypes/CORINA_Symphony/results.txt'
descriptors = '/share/home/ghelman/dev/read_across/applied/fingerprints/chemotypes/toxprint_V2.0_r711.xml'

with open(smile_file, 'w') as f:
    f.write('\n'.join([record['smiles'] + '\t' + record['dsstox_sid'] for record in new_compounds]))

call([corina, '-N', 'symphony', 'batch', '-i', smile_file, '-o', output, 'descriptors', '-f', descriptors])

df = pd.read_csv(output, sep=';')
df = df.drop(['M_COMPOUND_HISTORY_[STRING]', 'M_CORINA_SYMPHONY_ERRORS_[STRING]'], axis=1)

fp_names = df.columns.values[1:df.shape[0]]
new_chemotypes = []
for (i, row) in df.iterrows():
    sid = row['M_NAME']
    fps_binary = row.drop('M_NAME')
    fps = [fp for (fp, b) in fps_binary.iteritems() if b]
    new_chemotypes.append({'dsstox_sid': sid, 'chemotypes': {'ds': fps, 'n': len(fps)}})

call(['rm',smile_file])
call(['rm',output])

load_records(compounds,new_compounds)
load_records(physprop,new_physprop)
load_records(fingerprints,new_fps)
load_records(chemotypes,new_chemotypes)
