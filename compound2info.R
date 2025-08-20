
library(httr)
library(jsonlite)

compound_info = function(chem_name) {
  chem_search_url = paste0("https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name__icontains=", chem_name)
  chem_response = GET(chem_search_url)
  
  if (http_status(chem_response)$category == "Success") {
    chem_content = fromJSON(content(chem_response, "text", encoding = "UTF-8"))
  } else {
    stop("Failed to retrieve chemicals information from ChEMBL.")
  }
  
  chem_ids = chem_content$molecules$molecule_chembl_id
  chem_ids = Filter(function(x) !is.null(x) & !is.na(x), chem_ids)
  
  # find the drug indications associated with the parent compound
  numerical_ids =  as.numeric(sub("CHEMBL", "", chem_ids))
  parent_idx = which.min(numerical_ids) # parent chemical has the lowest CHEMBL ID
  parent_chem_id = chem_ids[parent_idx]
  indication_search_url = paste0("https://www.ebi.ac.uk/chembl/api/data/drug_indication.json?molecule_chembl_id__exact=", parent_chem_id)
  indication_response = GET(indication_search_url)
  if(http_status(indication_response)$category == "Success") {
    indication_content = fromJSON(content(indication_response, "text", encoding = "UTF-8"))
  } else {
    stop("Failed to get indication information from CHEMBL.")
  }
  indications = indication_content$drug_indications$efo_term
  indications = Filter(function(x) !is.null(x) & !is.na(x), indications)
  
  # ChEMBL API to get all mechanisms associated with compound family
  target_ids = c() # initialise
  mechanisms_of_action = c()
  for (id in chem_ids) {
    mech_search_url = paste0("https://www.ebi.ac.uk/chembl/api/data/mechanism.json?molecule_chembl_id=", id)
    mech_response = GET(mech_search_url)
    
    if(http_status(mech_response)$category == "Success") {
      mech_content = fromJSON(content(mech_response, "text", encoding = "UTF-8"))
      if (length(mech_content$mechanisms) > 0){ # if mechanism exists
        target_id = mech_content$mechanisms$target_chembl_id
        target_ids = c(target_ids, target_id) 
        mechanism = mech_content$mechanisms$mechanism_of_action
        mechanisms_of_action = c(mechanisms_of_action, mechanism)
      }
    } else {
      stop("Failed to retrieve mechanism information from ChEMBL.")
    }
  }
  target_ids = Filter(function(x) !is.null(x) & !is.na(x), target_ids)
  mechanisms_of_action = Filter(function(x) !is.null(x) & !is.na(x), mechanisms_of_action)
  
  
  # ChEMBL API to get gene accession IDS from target ChEMBL IDs
  accession_ids = c()
  for (id in target_ids) {
    target_search_url = paste0("https://www.ebi.ac.uk/chembl/api/data/target/", id, ".json")
    target_response = GET(target_search_url)
    
    if(http_status(target_response)$category == "Success") {
      target_content = fromJSON(content(target_response, "text", encoding = "UTF-8"))
      if (length(target_content$target_components) > 0) { # if target components exist
        accession_id = target_content$target_components$accession
        accession_ids = c(accession_ids, accession_id)
      }
    } else {
      stop("Failed to get target information from ChEMBL.")
    }
  }
  accession_ids = Filter(function(x) !is.null(x) & !is.na(x), accession_ids)
  
  # UniProt API to get gene names from accession IDs
  gene_names = c()
  for (id in accession_ids) {
    uniprot_query_url = paste0("https://rest.uniprot.org/uniprotkb/", id, ".json")
    uniprot_response = GET(uniprot_query_url)
    
    if(http_status(target_response)$category == "Success") {
      uniprot_content = fromJSON(content(uniprot_response, "text", encoding = "UTF-8"))
      gene_name = uniprot_content$genes$geneName$value
      gene_names = c(gene_names, gene_name)
    } else {
      stop("Failed to get information from UniProt.")
    }
  }
  gene_names = Filter(function(x) !is.null(x) & !is.na(x), gene_names)
  
  # tidy up, get rid of replicates
  indications = unique(indications)
  mechanisms_of_action = unique(mechanisms_of_action)
  gene_names = unique(gene_names)
  
  return(list(indications = indications, mechanisms = mechanisms_of_action, target_genes = gene_names))
}
