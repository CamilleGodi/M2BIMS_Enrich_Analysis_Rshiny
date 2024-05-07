################ GSEA - GO-terms ################################################

### GSEA - GO-terms

#' @description
#' Perform whole enrichment analysis on GO terms with GSEA
#' @param reactive_annotated_data : data.frame - data.frame from Rshiny reactive filter
#' @param organism_db : organism annotation db to yse (example : "org.Hs.eg.db")
#' @return enrichResults - Ready to plot enrichment analysis results
#'
#' @example do_gsea_go_terms(filtered_data,"org.Hs.eg.db")
#'
do_gsea_go_terms <- function(reactive_annotated_data, organism_db, reactive_abs) {
  gsea_ids <- prepare_gsea(reactive_annotated_data, abs = reactive_abs)
  gsea_go <- load_gsea_GO_enrichment(gsea_ids, organism_db)
  return(gsea_go)
}

#######################

### GSEA - KEGG

#' @description
#' Perform whole enrichment analysis on GO terms with GSEA
#' @param reactive_annotated_data : data.frame - data.frame from Rshiny reactive filter
#' @param organism_db : organism annotation db to yse (example : "org.Hs.eg.db")
#' @param universe : organism "universe" prepared with prepare_universe(reactive_annotated_data, organism_db, "ENSEMBL") 
#' @param kegg_organism_code : organism 3- nor 4-letters code (example : "hsa")sis results
#'
#' @example do_gsea_kegg(filtered_data, "org.Hs.eg.db", "hsa")
#'
do_gsea_kegg <- function(reactive_annotated_data, organism_db, kegg_organism_code, reactive_abs) {
  gsea_ids <- prepare_gsea(reactive_annotated_data, abs = reactive_abs)
  gsea_kegg <- load_gsea_kegg_enrichment(gene_list = gsea_ids,
                            organism_db = kegg_organism_code)
  return(gsea_kegg)
}

#######################

### GSEA - Reactome

#' @description
#' Perform whole enrichment analysis on GO terms with GSEA
#' @param reactive_annotated_data : data.frame - data.frame from Rshiny reactive filter
#' @param reactome_organism_name : organism name in reactome (example : "human")
#' @param universe : organism "universe" prepared with prepare_universe(reactive_annotated_data, organism_db, "ENSEMBL") 
#' @return enrichResults - Ready to plot enrichment analysis results
#'
#' @example do_gsea_reactome(filtered_data, "human")
#'
do_gsea_reactome <- function(reactive_annotated_data, reactome_organism_name, reactive_abs) {
  gsea_ids <- prepare_gsea(reactive_annotated_data, abs = reactive_abs)
  gsea_reactome <- load_gsea_reactome_enrichment(gene_list = gsea_ids,
                                    organism_db = reactome_organism_namee)
  return(gsea_reactome)
}