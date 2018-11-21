from __future__ import print_function
from datetime import datetime
from collections import defaultdict

def get_new_chemicals(mysql_connection,compound_collection,tox_collection,log_collection,meta_collection,date,function_log): #not_list should be comma-separated list like in an IN sql query
    time_of_exec=datetime.utcnow().replace(microsecond=0)
    last_exec=datetime.strftime(date,"%Y-%m-%d %H:%M:%S") #date should be type datetime.datetime
    cursor=mysql_connection.cursor(dictionary=True)
    query="SELECT casrn, chemspider_id, dsstox_compound_id AS dsstox_cid, dsstox_substance_id AS dsstox_sid, generic_substances.id AS gsid, jchem_inchi_key AS inchi_key, acd_iupac_name AS iupac, mol_weight, preferred_name AS name, pubchem_cid, smiles FROM stg_dsstox.compounds INNER JOIN stg_dsstox.generic_substance_compounds ON compounds.id=generic_substance_compounds.fk_compound_id INNER JOIN stg_dsstox.generic_substances ON generic_substances.id=generic_substance_compounds.fk_generic_substance_id WHERE compounds.created_at>'" + last_exec + "'" 
    cursor.execute(query)
    new_compounds=cursor.fetchall()
    cursor.close()
    if len(new_compounds)>0:
	    for compound in new_compounds:
		dtxcid=compound['dsstox_cid']
		compound['measured_props']=get_measured_props(mysql_connection,dtxcid)
		compound['predicted_props']=get_predicted_props(mysql_connection,dtxcid)
		dtxsid=compound['dsstox_sid']
		compound['synonyms']=get_synonyms(mysql_connection,dtxsid)
		compound['created_at']=time_of_exec
		compound['updated_at']=time_of_exec
		#Compounds
		#prev_compound=compound_collection.find_one({'dsstox_cid':dtxcid}) #Need for logs later
		compound_collection.replace_one({'dsstox_cid':dtxcid},compound,upsert=True)
		#Toxvals
		toxvals=get_tox_data(mysql_connection,dtxcid)
		tox_collection.update_one({'dsstox_cid':dtxcid},{'$addToSet':{'toxvals': {'$each': toxvals}}},upsert=True)
		#Logs
		log={'date': time_of_exec,'change': 'Added to database'}
		log_collection.update_one({'dsstox_cid':dtxcid},{'$push': {'logs': log}},upsert=True)
		print(dtxcid, file=function_log)
	    meta_collection.insert_one({'date': time_of_exec, 'compounds_added': len(new_compounds)})

def get_measured_props(mysql_connection,dsstox_cid):
    cursor=mysql_connection.cursor(dictionary=True)
    query="SELECT name, dsstox_compound_id, result_value FROM stg_chemprop.measured_properties a INNER JOIN stg_dsstox.source_generic_substance_mappings b ON a.efk_dsstox_source_substance_id = b.fk_source_substance_id INNER JOIN stg_dsstox.generic_substances c on b.fk_generic_substance_id = c.id INNER JOIN stg_chemprop.endpoints ON endpoints.id= a.fk_endpoint_id INNER JOIN stg_dsstox.generic_substance_compounds e ON c.id=e.fk_generic_substance_id INNER JOIN stg_dsstox.compounds f ON e.fk_compound_id=f.id WHERE dsstox_compound_id = '" + dsstox_cid + "'"
    cursor.execute(query)
    ret=cursor.fetchall()
    cursor.close()
    measured_properties = defaultdict(list)
    for prop in ret:
        measured_properties[keyify(prop['name'])].append(prop['result_value'])
    return measured_properties

def get_predicted_props(mysql_connection,cid):
    cursor=mysql_connection.cursor(dictionary=True)
    query="SELECT name, dsstox_compound_id, result_value FROM stg_chemprop.qsar_predicted_properties INNER JOIN stg_dsstox.generic_substance_compounds ON qsar_predicted_properties.efk_dsstox_compound_id = generic_substance_compounds.fk_compound_id INNER JOIN stg_dsstox.generic_substances ON generic_substance_compounds.fk_generic_substance_id = generic_substances.id INNER JOIN stg_qsar.models ON qsar_predicted_properties.efk_qsar_model_id = models.id  INNER JOIN stg_dsstox.compounds ON qsar_predicted_properties.efk_dsstox_compound_id = compounds.id WHERE dsstox_compound_id = '" + cid + "'"
    cursor.execute(query)
    ret=cursor.fetchall()
    cursor.close()
    predicted_properties = defaultdict(list)
    for prop in ret:
        predicted_properties[keyify(prop['name'])].append(prop['result_value'])
    return predicted_properties

def get_tox_data(mysql_connection,cid):
    cursor=mysql_connection.cursor(dictionary=True)
    query="SELECT * FROM stg_toxval_v4.toxval INNER JOIN stg_toxval_v4.chemical ON stg_toxval_v4.chemical.chemical_id=stg_toxval_v4.toxval.chemical_id WHERE dsstox_compound_id = '" + cid + "'"
    cursor.execute(query)
    ret=cursor.fetchall()
    cursor.close()
    return ret

def get_synonyms(mysql_connection,sid):
    cursor=mysql_connection.cursor(dictionary=True)
    query="SELECT identifier FROM stg_dsstox.synonyms INNER JOIN stg_dsstox.generic_substances ON synonyms.fk_generic_substance_id = generic_substances.id WHERE dsstox_substance_id='" + sid + "'"
    cursor.execute(query)
    ret=cursor.fetchall()
    cursor.close()
    synonyms=[r['identifier'] for r in ret]
    return synonyms

def keyify(str):
    return str.replace('.','-')
    meta_collection.insert_one({'date': time_of_exec, 'compounds_added': len(new_compounds)})

def get_new_predictions(mysql_connection,compound_collection,log_collection,meta_collection,date):
    time_of_exec=datetime.utcnow().replace(microsecond=0)
    last_exec=datetime.strftime(date,"%Y-%m-%d %H:%M:%S") #date should be type datetime.datetime
    cursor=mysql_connection.cursor(dictionary=True)
    query="SELECT preferred_name, dsstox_compound_id AS dsstox_cid, name, result_value FROM stg_chemprop.qsar_predicted_properties INNER JOIN stg_dsstox.generic_substance_compounds ON qsar_predicted_properties.efk_dsstox_compound_id = generic_substance_compounds.fk_compound_id INNER JOIN stg_dsstox.generic_substances ON generic_substance_compounds.fk_generic_substance_id = generic_substances.id INNER JOIN stg_qsar.models ON qsar_predicted_properties.efk_qsar_model_id = models.id INNER JOIN stg_dsstox.compounds ON qsar_predicted_properties.efk_dsstox_compound_id = compounds.id WHERE qsar_predicted_properties.created_at >'" + last_exec + "'"
    cursor.execute(query)
    ret=cursor.fetchall()
    cursor.close()
    add_count=0
    for pred in ret:
        dtxcid=pred['dsstox_cid']
        prop_name=pred['name']
        value=pred['result_value']
        if compound_collection.find_one({'$and': [{'dsstox_cid':dtxcid},{('predicted_props.'+prop_name): {'$eq':value}}]}) is None:
            compound_collection.update_one({'dsstox_cid':dtxcid},{'$push':{('predicted_props.'+prop_name):value}})
            compound_collection.update_one({'dsstox_cid':dtxcid},{'$set': {'updated_at':time_of_exec}})                             
            log={"time":time_of_exec,"type":'add','collection':'compounds','compound':dtxcid,'field':'predicted_props','prop_name':prop_name,'value':value}
            log_collection.update_one({'dsstox_cid':dtxcid},{'$push':{'logs':log}})
            add_count+=1
    if add_count>0:
        meta_collection.insert_one({'date': time_of_exec, 'predictions_added': add_count})

def get_new_measurements(mysql_connection,compound_collection,log_collection,meta_collection,date):
    time_of_exec=datetime.utcnow().replace(microsecond=0)
    last_exec=datetime.strftime(date,"%Y-%m-%d %H:%M:%S") #date should be type datetime.datetime
    cursor=mysql_connection.cursor(dictionary=True)
    query="SELECT preferred_name, dsstox_compound_id AS dsstox_cid, name, result_value FROM stg_chemprop.measured_properties a INNER JOIN stg_dsstox.source_generic_substance_mappings b ON a.efk_dsstox_source_substance_id = b.fk_source_substance_id INNER JOIN stg_dsstox.generic_substances c on b.fk_generic_substance_id = c.id INNER JOIN stg_chemprop.endpoints ON endpoints.id= a.fk_endpoint_id INNER JOIN stg_dsstox.generic_substance_compounds d ON c.id=d.fk_generic_substance_id INNER JOIN stg_dsstox.compounds ON d.fk_compound_id=compounds.id WHERE a.created_at >'" + last_exec + "'"
    cursor.execute(query)
    ret=cursor.fetchall()
    cursor.close()
    add_count=0
    for meas in ret:
        dtxcid=meas['dsstox_cid']
        prop_name=meas['name']
        value=meas['result_value']
        if compound_collection.find_one({'$and': [{'dsstox_cid':dtxcid},{('measured_props.'+prop_name): {'$eq':value}}]}) is None:
            compound_collection.update_one({'dsstox_cid':dtxcid},{'$push':{('measured_props.'+prop_name):value}})
            compound_collection.update_one({'dsstox_cid':dtxcid},{'$set': {'updated_at':time_of_exec}})
            log={"time":time_of_exec,"type":'add','collection':'compounds','compound':dtxcid,'field':'measured_props','prop_name':prop_name,'value':value}
            log_collection.update_one({'dsstox_cid':dtxcid},{'$push':{'logs':log}})
            add_count+=1
    if add_count>0:
        meta_collection.insert_one({'date': time_of_exec, 'measurements_added': add_count})
