################ GSEA - GO-terms ################################################

### GSEA - GO-terms

#' @description
#' Perform whole enrichment analysis on GO terms with GSEA
#' @param reactive_annotated_data : data.frame - data.frame from Rshiny reactive filter
#' @param organism_db : organism annotation db to yse (example : "org.Hs.eg.db")
#' @param ontology : ontologies to consider for GO-term analysis. Should use "MF", "CC" and "BP" separately rather than using "ALL"
#' @param p_value_cutoff : p value cutoff
#' @param p_adj_cutoff :  adjusted p value cutoff
#' @param q_value_cutoff :  adjusted q value cutoff
#' @return enrichResults - Ready to plot enrichment analysis results
#'
#' @example do_gsea_go_terms(filtered_data,"org.Hs.eg.db", "MF", 0.01, 0.05, 0.05)
#'
do_gsea_go_terms <- function(reactive_annotated_data, organism_db, universe, ontology, p_value_cutoff, p_adj_cutoff, q_value_cutoff, min_GS_size, max_GS_size) {
  
  gsea_ids <- prepare_gsea(reactive_annotated_data,
                           metric = "log2FC",
                          abs = TRUE)
  
  gsea_go <- load_gsea_GO_enrichment(gsea_ids, 
                                     universe = universe, 
                                     organism_db = organism_db, 
                                     min_GS_size = 10, 
                                     max_GS_size = 50)
  
  gsea_go_after_filter <- filter_table_enrich_results(gsea_go, 
                                                     p_value_cutoff = p_value_cutoff, 
                                                     p_adj_cutoff   = p_adj_cutoff, 
                                                     q_value_cutoff = q_value_cutoff)
  
  gsea_go_after_filter_ontologies <- filter_go_enrich_results(gsea_go_after_filter, ontology = ontology)
  
  return(gsea_go_after_filter_ontologies)
}

#######################

### GSEA - KEGG
#' @description
#' Perform whole enrichment analysis on GO terms with GSEA
#' @param reactive_annotated_data : data.frame - data.frame from Rshiny reactive filter
#' @param organism_db : organism annotation db to yse (example : "org.Hs.eg.db")
#' @param kegg_organism_code : organism 3- nor 4-letters code (example : "hsa")
#' @param p_value_cutoff : p value cutoff
#' @param p_adj_cutoff :  adjusted p value cutoff
#' @param q_value_cutoff :  adjusted q value cutoff
#' @return enrichResults - Ready to plot enrichment analysis results
#'
#' @example do_gsea_kegg(filtered_data, "org.Hs.eg.db", "hsa", 0.01, 0.05, 0.05)
#'
do_gsea_kegg <- function(reactive_annotated_data, organism_db, universe, kegg_organism_code, p_value_cutoff, p_adj_cutoff, q_value_cutoff) {
  
  gsea_ids <- prepare_gsea(reactive_annotated_data,
                           metric = "log2FC",
                           abs = TRUE)
  
  gsea_kegg <- load_gsea_kegg_enrichment(gene_list = gsea_ids,
                                         organism_db = kegg_organism_code,
                                         universe = universe)
  
  gsea_kegg_after_filter <- filter_table_enrich_results(gsea_kegg, 
                                                       p_value_cutoff = p_value_cutoff, 
                                                       p_adj_cutoff   = p_adj_cutoff, 
                                                       q_value_cutoff = q_value_cutoff)
  
  return(gsea_kegg_after_filter)
}
