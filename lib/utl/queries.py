databases={
    'dsstox':'ro_stg_dsstox',
    'chemprop':'ro_stg_chemprop',
    'qsar':'ro_stg_qsar',
    'invitrodb':'dev_invitrodb'
}


def get_new_chemicals(l):
    in_list= '("'+'","'.join(l)+'")'
    return "SELECT casrn, chemspider_id, dsstox_compound_id AS dsstox_cid, dsstox_substance_id AS dsstox_sid, generic_substances.id AS gsid, " \
        "jchem_inchi_key AS inchi_key, acd_iupac_name AS iupac, mol_weight, preferred_name AS name, pubchem_cid, smiles " \
        "FROM {dsstox}.compounds INNER JOIN {dsstox}.generic_substance_compounds ON compounds.id=generic_substance_compounds.fk_compound_id " \
        "INNER JOIN {dsstox}.generic_substances ON generic_substances.id=generic_substance_compounds.fk_generic_substance_id " \
        "WHERE generic_substances.dsstox_substance_id NOT IN ".format(**databases) + in_list

def get_measured_properties(sid):
    return  "SELECT name, result_value FROM {chemprop}.measured_properties a " \
            "INNER JOIN {dsstox}.source_generic_substance_mappings b ON a.efk_dsstox_source_substance_id = b.fk_source_substance_id " \
            "INNER JOIN {dsstox}.generic_substances c on b.fk_generic_substance_id = c.id " \
            "INNER JOIN {chemprop}.endpoints ON endpoints.id= a.fk_endpoint_id " \
            "INNER JOIN {dsstox}.generic_substance_compounds e ON c.id=e.fk_generic_substance_id " \
            "INNER JOIN {dsstox}.compounds f ON e.fk_compound_id=f.id " \
            "WHERE c.dsstox_substance_id = '".format(**databases) + sid + "'"

def get_predicted_properties(sid):
    return  "SELECT name, result_value FROM {chemprop}.qsar_predicted_properties " \
            "INNER JOIN {dsstox}.generic_substance_compounds ON qsar_predicted_properties.efk_dsstox_compound_id = generic_substance_compounds.fk_compound_id " \
            "INNER JOIN {dsstox}.generic_substances ON generic_substance_compounds.fk_generic_substance_id = generic_substances.id " \
            "INNER JOIN {qsar}.models ON qsar_predicted_properties.efk_qsar_model_id = models.id " \
            "INNER JOIN {dsstox}.compounds ON qsar_predicted_properties.efk_dsstox_compound_id = compounds.id " \
            "WHERE generic_substances.dsstox_substance_id = '".format(**databases) + sid + "'"

def get_synonyms(sid):
    return "SELECT identifier FROM {dsstox}.synonyms " \
           "INNER JOIN {dsstox}.generic_substances ON synonyms.fk_generic_substance_id = generic_substances.id " \
           "WHERE dsstox_substance_id='".format(**databases) + sid + "'"

def get_invitrodb():
    return  "SELECT dsstox_substance_id AS dsstox_sid, c_casrn_id AS cas, preferred_name AS name, modl, hitc, fitc, coff, actp, modl_er, modl_tp, modl_ga, modl_gw, modl_la, modl_lw, modl_prob, modl_rmse, modl_acc, modl_acb, modl_ac10, bmad, resp_max, resp_min, max_mean, max_mean_conc, max_med, max_med_conc, logc_max, logc_min, cnst, hill, hcov, gnls, gcov, cnst_er, cnst_aic, cnst_rmse, cnst_prob, hill_tp, hill_tp_sd, hill_ga, hill_ga_sd, hill_gw, hill_gw_sd, hill_er, hill_er_sd, hill_aic, hill_rmse, hill_prob, gnls_tp, gnls_tp_sd, gnls_ga, gnls_ga_sd, gnls_gw, gnls_gw_sd, gnls_la, gnls_la_sd, gnls_lw, gnls_lw_sd, gnls_er,gnls_er_sd, gnls_aic, gnls_rmse, gnls_prob, nconc, npts, nrep, nmed_gtbl, tmpi"\
            " FROM {invitrodb}.level5"\
            " INNER JOIN {invitrodb}.level4 ON level4.l4id=level5.l4id"\
            " INNER JOIN {invitrodb}.sample ON level4.spid=sample.sa_sample_id"\
            " INNER JOIN {invitrodb}.casrn ON sample.sa_gsid=casrn.c_gsid_id"\
            " INNER JOIN {dsstox}.generic_substances ON casrn.c_casrn_id=generic_substances.casrn".format(**databases)

def get_sid_from_cid(cid):
    return "SELECT dsstox_substance_id AS dsstox_sid FROM {dsstox}.generic_substances " \
            "INNER JOIN {dsstox}.generic_substance_compounds ON generic_substances.id = generic_substance_compounds.fk_generic_substance_id " \
            "INNER JOIN {dsstox}.compounds ON compounds.id = generic_substance_compounds.fk_compound_id " \
            "WHERE dsstox_compound_id='".format(**databases) + cid + "'"

def get_cid_from_sid(sid):
    return  "SELECT dsstox_compound_id AS dsstox_cid FROM {dsstox}.generic_substances " \
            "INNER JOIN {dsstox}.generic_substance_compounds ON generic_substances.id = generic_substance_compounds.fk_generic_substance_id " \
            "INNER JOIN {dsstox}.compounds ON compounds.id = generic_substance_compounds.fk_compound_id " \
            "WHERE dsstox_substance_id='".format(**databases) + sid + "'"
            
def get_lists():
        return "SELECT a.name, a.label, a.long_description, ncct_contact, source_contact, source_contact_email,list_type, list_update_mechanism, list_accessibility, curation_complete, source_data_updated_at, a.created_by, a.updated_at, block_updates, a.short_description, d.name as source_level, dsstox_substance_id as dsstox_sid \
        from {dsstox}.chemical_lists as a \
        join {dsstox}.source_substances as b ON a.id=b.fk_chemical_list_id \
        join {dsstox}.source_generic_substance_mappings as bc on b.id=bc.fk_source_substance_id \
        join {dsstox}.generic_substances as c on bc.fk_generic_substance_id=c.id \
        join {dsstox}.list_types as d on a.fk_list_type_id=d.id".format(**databases)
