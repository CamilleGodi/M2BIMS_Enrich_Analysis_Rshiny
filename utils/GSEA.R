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
                                    organism_db = reactome_organism_name)
  return(gsea_reactome)
}

#' @description
#' dessine le dotplot issu d'un enrichissement, pour une catégorie donnée trié selon ce qu'on choisi
#' @param ora_df: objet issus de l'enrichissement
#' @param order : paramètre selon lequel on veut ordonné, peut être pvalue, p.adjust, qvalue (dans ce cas ordre croissant), ou un autre paramètre numérique (dans ce cas décroissant)
#' @param category : catégorie qu'on veut montrer sur le dotplot (doit être numérique ou un ratio)
#' @param show_category : les n premieres categories a montrer
#' @param alpha : seuil de significativité pour accepter nos valeurs (important pour le trie)
#' @param title : titre
#' @param ylab : nom de l'axe y
#' @param xlab : nom de l'axe x
#' @param gradient_col : vecteur de couleur pour le gradient
#' @param y_text_size : taille du texte en ordonnée
#' 
#' @example draw_dotplot(ora.bp)
#' 
draw_gsea_dotplot = function(gsea_results_as_input,
                        order = "p.adjust",
                        category = "richFactor",
                        show_category = 10,
                        alpha = 0.05,
                        ylab = "Ontologies",
                        xlab = category,
                        title = paste("Dotplot of Gene Ontologies sorted by",order),
                        gradient_col = c("#f7ca64", "#46bac2", "#7e62a3"),
                        y_text_size = 10) {
  if(nrow(slot(gsea_results_as_input,"result")) < show_category){
    show_category = nrow(slot(gsea_results_as_input,"result"))
  }
  gsea_results = add_rich_factor_to_gsea(gsea_results_as_input)
  results_for_graph = slot(gsea_results,"result")[1:show_category,]
  df = data.frame("x" = results_for_graph[,category],
                  "y" = results_for_graph[,"Description"],
                  "order" = results_for_graph[,order],
                  "size" = results_for_graph[,"setSize"],
                  "color" = results_for_graph[,"p.adjust"],
                  "discriminant" = results_for_graph[,"NES"])
  df[,"x"][df[,"discriminant"] < 0] = -df[,"x"]
  ggplot2::ggplot(df, aes(x = x,y = fct_reorder(y,order))) +
    ggplot2::geom_point(aes(color = color, size = size)) +
    ggplot2::scale_y_discrete(label = function(y) stringr::str_trunc(y, 40))+
    ggplot2::scale_color_gradientn (
      colours = gradient_col,
      trans = "log10",
      guide = ggplot2::guide_colorbar(reverse =
                                        TRUE, order = 1)
    ) +
    ggplot2::theme(axis.text = element_text(size = 6))+
    ggplot2::labs(y = ylab, title = title, x = xlab,size = "Number of genes\nin set", color = "p.adjust") +
    ggplot2::geom_vline(xintercept = 0)+
    DOSE::theme_dose(12) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = y_text_size))
}

add_rich_factor_to_gsea = function(gsea_results_as_input){
  gsea_output = gsea_results_as_input
  core_enrichment = slot(gsea_output,"result")$core_enrichment
  slot(gsea_output,"result")$number_of_gene = sapply(core_enrichment,function(x) length(unlist(str_split(x,"\\/")))) %>% unname()
  slot(gsea_output,"result")$richFactor = slot(gsea_output,"result")$number_of_gene/slot(gsea_output,"result")$setSize
  return(gsea_output)
}

show_table_gsea = function(gsea_results_as_input){
  add_rich_factor_to_gsea(gsea_results_as_input) %>%
  slot("result") %>% 
    dplyr::select(Description,number_of_gene,setSize,richFactor,enrichmentScore,NES,pvalue,p.adjust,rank) %>%
    format(digits = 3) %>%
    # DT::datatable(options = list(scrollX = TRUE,pageLength = 25))  %>%
    # htmltools::tagList() %>%
    return()
}
