source("./utils/enrichment.R")

################ ORA - GO-terms ################################################

### ORA - GO-terms

#' @description
#' Perform whole enrichment analysis on GO terms with ORA
#' @param reactive_annotated_data : data.frame - data.frame from Rshiny reactive filter
#' @param organism_db : organism annotation db to yse (example : "org.Hs.eg.db")
#' @param ontology : ontologies to consider for GO-term analysis. Should use "MF", "CC" and "BP" separately rather than using "ALL"
#' @param p_value_cutoff : p value cutoff
#' @param p_adj_cutoff :  adjusted p value cutoff
#' @param q_value_cutoff :  adjusted q value cutoff
#' @return enrichResults - Ready to plot enrichment analysis results
#'
#' @example do_ora_go_terms(filtered_data,"org.Hs.eg.db", "MF", 0.01, 0.05, 0.05)
#'
do_ora_go_terms <- function(reactive_annotated_data, organism_db, ontology, p_value_cutoff, p_adj_cutoff, q_value_cutoff) {
  table_filtered_new_ids <- prepare_pipe(reactive_annotated_data, organism_db = organism_db, "ENSEMBL")
  
  ora_ids <- prepare_ora(table_filtered_new_ids)
  universe <- prepare_universe(reactive_annotated_data)
  
  ora_go <- load_ora_go(ora_ids, universe, organism_db)
  ora_go_after_filter <- filter_table_enrich_results(ora_go, 
                                                     p_value_cutoff = p_value_cutoff, 
                                                     p_adj_cutoff   = p_adj_cutoff, 
                                                     q_value_cutoff = q_value_cutoff)
  results <- filter_go_enrich_results(ora_go_after_filter, ontology = ontology)
  
  return(results)
}

#######################


# TMP : Tableau d'entré (non filtré etc), réactive "filtered_data"
table_filtered_tmp <-
  read.csv("./Example_files/filtered_data.csv", sep = ",")

# TMP : ext variables -> replace with the reactive variables
organism_library_tmp <- "org.Mm.eg.db"
p_value_cutoff_tmp <- 0.01
p_adj_cutoff_tmp   <- 0.05
q_value_cutoff_tmp <- 0.05
ontology_tmp <- c("MF", "CC")

res_tmp <- do_ora_go_terms(
    table_filtered_tmp,
    organism_library_tmp,
    ontology_tmp,
    p_value_cutoff_tmp,
    p_adj_cutoff_tmp,
    q_value_cutoff_tmp
  )

res_tmp %>% enrich_pagination

res_tmp %>% draw_dotplot(show_category = 30)

res_tmp %>% draw_cnetplot(
  category_label = 0.6,
  gene_list = res_tmp@result$geneID,
  category_color = "red",
  node_label = "category"
)

res_tmp %>% draw_treeplot(
  gradient_col = c("green", "black"),
  showCategory = 50,
  n_cluster = 10,
  label_words_n = 4,
  h_clust_method = "ward.D2"
)

res_tmp %>% draw_emapplot()
