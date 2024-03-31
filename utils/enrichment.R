library(tidyverse)

#' @docType https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html


#' @description
#'
#' @param results_from_de_analysis dataframe from differential analysis
#' @example independent_filtering(analysis_of_cd4_vs_cd8)
#'
independent_filtering = function(results_from_de_analysis) {
  filtered_results_from_de_analysis = results_from_de_analysis %>% data.frame()   # on s'assure du type data.frame
  if ("padj" %in% colnames(filtered_results_from_de_analysis)) {
    filtered_results_from_de_analysis = filtered_results_from_de_analysis[!is.na(filtered_results_from_de_analysis[, "padj"]), ] #on retire les données "NA" qui poseront problème par la suite
  }
  return(filtered_results_from_de_analysis)
}

#' @description
#'
#' @param results_without_na data.frame - results of differential analysis without NA
#' @param from character - id known in the table to be converted, can be EMSEMBL | SYMBOL | ENTREZID
#' @param to character - type of id we want as output, can be EMSEMBL | SYMBOL | ENTREZID
#' @param organism_db character - organism in annotation_dbi, need for convertion, exemple : "Org.Mm.eg.db"
#'
conversion_table = function(vector_of_ids,
                            from,
                            organism_db) {
  if (from == "ENTREZID") {
    return(data.frame(vector_of_ids))
  } else{
    clusterProfiler::bitr(
      vector_of_ids,
      fromType = from,
      toType = "ENTREZID",
      OrgDb = organism_db
    ) %>% return()
  }
}

remove_duplicate = function(ids) {
  if (ncol(ids) > 1) {
    ids = ids[!duplicated(ids["ENTREZID"]), ]
    return(ids)
  }
  colnames(ids) = "ENTREZID"
  return(ids)
}

#' @description 
#' 
#' @param reactive_annotated_data data.frame - data.frame from Rshiny reactive filter
#' @param organism_db character - character of annotationDbi available database
#' @param from character - ID use as entry, might be ENTREZID|ENSEMBL|SYMBOL
#' 
prepare_pipe <- function(reactive_annotated_data,
                         organism_db,
                         from) {
  filtered_table <- reactive_annotated_data %>% independent_filtering()
  ids <- conversion_table(filtered_table[, "ID"],
                          organism_db = organism_db,
                          from = from)
  filtered_table %>% convert_results_ids(ids) %>% return()
}


convert_results_ids = function(results,
                               ids) {
  # on retire les doublons, lors de la conversion, certains id peuvent mener à un id unique et inversement
  
  if (ncol(ids) > 1) {
    ids = remove_duplicate(ids)
    results[, "ENTREZID"] = ids[match(results[, "ID"], ids[, 1]), "ENTREZID"]
  } else {
    results[, "ENTREZID"] = results[, "ID"]
  }
  results[!is.na(results[, "ENTREZID"]), ]  %>% return()
}

#### définir que diff_expressed n'existe pas
prepare_ora = function(results) {
  results[results[, "diff_expressed"] != "NO_DE", "ENTREZID"] %>% return()
}

prepare_gsea = function(results,
                        metric = "log2FC",
                        abs = FALSE) {
  if (abs) {
    setNames(abs(results[, metric]), results[, "ENTREZID"]) %>% sort(decreasing = TRUE) %>% return()
  } else {
    setNames(results[, metric], results[, "ENTREZID"]) %>% sort(decreasing = TRUE) %>% return()
  }
}

prepare_universe = function(initiale_data,
                            organism_db,
                            from) {
  initiale_data[, "ID"] %>%
    conversion_table(from = from, organism_db = organism_db) %>%
    remove_duplicate() %>%
    dplyr::select(ENTREZID) %>%
    unlist() %>%
    unname() %>%
    return()
}

#' @description
#' Launch enrichment analysis on GO terms with GSEA
#' @param gene_list : vecteur trié dans l'ordre décroissant des gènes selon la métrique voulu
#' @param organism_db : base de donnée à utiliser (exemple : "org.Hs.eg.db")
#' @param ont : ontologie à utiliser pour l'analyse des GO termes (de préférence utiliser séparément MF, CC et BP plutôt que ALL)
#' @param key_type : type d'identifiant utilisé
#' @param min_GS_size : taille minimal du gene set pris pour l'analyse
#' @param max_GS_size : taille maximal du gene set pris pour l'analyse
#' @param pvalue_cutoff : valeur de seuil minimale de la p value
#' @param p_adjust_method : méthode d'ajustement de la p-value
#' @param verbose voir défiler le text pendant le lancement de l'algo
#' @param n_perm_simple nombre de permutation, plus la valeur est élevé plus le process est long, mais plus les données sont exactes
#' @return vecteur trié dans l'ordre décroissant des gènes selon la métrique voulu
#'
#' @example load_gsea_GO_enrichment(results(dds),"org.Hs.eg.db")
#'
load_gsea_GO_enrichment = function(gene_list = list(),
                                   organism_db = character(),
                                   min_GS_size = 3,
                                   max_GS_size = 500) {
  output <- clusterProfiler::gseGO(
    geneList = gene_list,
    ont = "ALL",
    keyType = "ENTREZID",
    minGSSize = min_GS_size,
    maxGSSize = max_GS_size,
    pvalueCutoff = 1,
    verbose = FALSE,
    OrgDb = organism_db,
    nPermSimple = 10000,
    eps = 0,
    pAdjustMethod = "BH"
  )
  output %>% return()
}

#' @description
#' Launch enrichment analysis on KEGG with GSEA
#' @param gene_list : vecteur trié dans l'ordre décroissant des gènes selon la métrique voulu
#' @param organism_db : base de donnée KEGG utilisée (exemple : "hsa")
#' @param key_type : type d'identifiant utilisé
#' @param min_GS_size : taille minimal du gene set pris pour l'analyse
#' @param max_GS_size : taille maximal du gene set pris pour l'analyse
#' @param pvalue_cutoff : valeur de seuil maximale de la p value
#' @param p_adjust_method : méthode d'ajustement de la p-value
#' @param n_perm nombre de permutation effectuée pour le random walk
#' @return vecteur trié dans l'ordre décroissant des gènes selon la métrique voulu
#'
#' @example load_gsea_kegg_enrichment(results(dds),"hsa")

load_gsea_kegg_enrichment = function(gene_list = list(),
                                     organism_db = character(),
                                     min_GS_size = 3,
                                     max_GS_size = 500) {
  output <- clusterProfiler::gseKEGG(
    geneList = gene_list,
    organism = organism_db,
    minGSSize = min_GS_size,
    maxGSSize = max_GS_size,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    keyType = "ncbi-geneid"
  )
  output %>% return()
}

#'
#' @description
#' Launch enrichment analysis on GO terms with ORA
#' @param gene_list : vecteur de gènes différentielelment exprimés
#' @param universe : liste des gènes initiaux avant analyse DE (soit ça, soit OrgDb)
#' @param organism_db : base de donnée KEGG utilisée (exemple : "hsa")
#' @param key_type : type d'identifiant utilisé
#' @param readable : booléan, est-ce que le retour se fait sous gene symbols (TRUE) ou gene_id (FALSE)
#' @param ont : ontologie à utiliser pour l'analyse des GO termes (de préférence utiliser séparément MF, CC et BP plutôt que ALL)
#' @param min_GS_size : taille minimal du gene set pris pour l'analyse
#' @param max_GS_size : taille maximal du gene set pris pour l'analyse
#' @param p_adjust_method : méthode d'ajustement de la p-value
#' @param pvalue_cutoff : valeur de seuil maximale de la p value
#' @param qvalue_cutoff : valeur de seuil maximale de la q value
#' @return vecteur trié dans l'ordre décroissant des gènes selon la métrique voulu
#'
#' @example load_go_over_representation(genes_DE,"org.Hs.eg.db")

load_ora_go = function(gene_list = list(),
                       universe = vector(),
                       organism_db = character(),
                       min_GS_size = 10,
                       max_GS_size = 500) {
  output = clusterProfiler::enrichGO(
    gene = gene_list,
    universe = universe,
    OrgDb = organism_db,
    keyType = "ENTREZID",
    ont = "ALL",
    minGSSize = min_GS_size,
    maxGSSize = max_GS_size,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  slot(output, "method") = "ORA"
  output = clusterProfiler::mutate(output, Subset = as.numeric(sub("/\\d+", "", BgRatio)))
  output = clusterProfiler::mutate(output, RichFactor = Count / Subset)
  output %>% return()
}

#' Analyse la sur-représentation des GO terms dans le jeu de gènes différentiellement exprimés
#' @param gene_list : vecteur trié dans l'ordre décroissant des gènes selon la métrique voulu
#' @param universe : liste des gènes initiaux avant analyse DE (soit ça, soit OrgDb)
#' @param organism_db : base de donnée KEGG utilisée (exemple : "hsa") (soit ça, soit universe)
#' @param key_type : type d'identifiant utilisé
#' @param min_GS_size : taille minimal du gene set pris pour l'analyse
#' @param max_GS_size : taille maximal du gene set pris pour l'analyse
#' @param pvalue_cutoff : valeur de seuil maximale de la p value
#' @param qvalue_cutoff : valeur de seuil maximale de la q value
#' @param p_adjust_method : méthode d'ajustement de la p-value
#' @return vecteur trié dans l'ordre décroissant des gènes selon la métrique voulu
#'
#' @example loag_kegg_over_representation(gene_list,"org.Hs.eg.db")
#'
load_ora_kegg = function(gene_list = list(),
                         organism_db = character(),
                         universe = vector(),
                         min_GS_size = 10,
                         max_GS_size = 500) {
  output = clusterProfiler::enrichKEGG(
    gene = gene_list,
    organism = organism_db,
    keyType = "ncbi-geneid",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    universe = universe,
    minGSSize = min_GS_size,
    maxGSSize = max_GS_size,
    qvalueCutoff = 1
  )
  slot(output, "method") = "ORA"
  output = clusterProfiler::mutate(output, Subset = as.numeric(sub("/\\d+", "", BgRatio)))
  output = clusterProfiler::mutate(output, RichFactor = Count / Subset)
  
  output %>% return()
}

filter_table_enrich_results = function(enrich_results,
                                       p_value_cutoff = 0.05,
                                       p_adj_cutoff = 0.05,
                                       q_value_cutoff = 0.05) {
  new_results = slot(enrich_results, "result")
  slot(enrich_results, "result") = new_results[(new_results[, "pvalue"] < p_value_cutoff) &
                                                 (new_results[, "p.adjust"] < p_adj_cutoff) &
                                                 (new_results[, "qvalue"] < q_value_cutoff),]
  enrich_results %>% return()
}

#' @description
#'
#' @param go_enrich_results_filtered clusterprofiler objectf - results from a GO analysis after filter of p-value, p-adjust and q-value
#' @param ontology character - ontology to be used

filter_go_enrich_results = function(go_enrich_results_filtered,
                                    ontology = NULL) {
  # if no ontology was choose, show all ontology
  if (is.null(ontology)) {
    ontology = c("MF", "CC", "BP")
  }
  new_results = slot(go_enrich_results_filtered, "result")
  slot(go_enrich_results_filtered, "result") = new_results %>% dplyr::filter(ONTOLOGY %in% ontology)
  return(go_enrich_results_filtered)
}

#' @description
#' permet l'affichage du tableau de données résultant de "get_enrich" ou d'une fonction ORA/GSEA de manière plus lisible en rmd
#' @param get_enrich l'objet clusterprofiler résultante d'une des fonctions au dessus permettant de faire du ORA ou du GSEA
#' @param get_columns les noms des colonnes montrées après enrichissement
#' @example get_enrich(ora.bp,0.05) || get_enrich(ora.bp, ora.bp@pvalueCutoff)
enrich_pagination = function(get_enrich,
                             get_columns = c(
                               "ID",
                               "Description",
                               "GeneRatio",
                               "BgRatio",
                               "RichFactor",
                               "p.adjust",
                               "Count",
                               "Subset"
                             ),
                             to_scientific = c("p.adjust"),
                             alpha_cutoff = 0.05) {
  if (slot(get_enrich, "method") == "ORA") {
    res = slot(get_enrich, "result")[, get_columns]
    res = res[res[, "p.adjust"] < alpha_cutoff, ]
    res[to_scientific] = apply(
      X = res[to_scientific],
      MARGIN = 2,
      FUN = function(x)
        signif(x, 3)
    )
    res[to_scientific]  = format(res[to_scientific] , scientific = T)
    res = format(res, digits = 3)
  } else if (slot(get_enrich, "method") == "GSEA") {
    res = slot(get_enrich, "result")[, c("ID", "Description", "setSize", "enrichmentScore")]
  }
  res %>% return()
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

draw_dotplot = function(ora_df,
                        order = "p.adjust",
                        category = "GeneRatio",
                        show_category = 10,
                        alpha = 0.05,
                        ylab = "Ontologies",
                        xlab = category,
                        title = paste("Dotplot of Gene Ontologies sorted by", order),
                        gradient_col = c("#f7ca64", "#46bac2", "#7e62a3"),
                        y_text_size = 10) {
  if (slot(ora_df, "method") != "ORA") {
    return("error")
  } else {
    ora_df = slot(ora_df, "result")
    
    if (order %in% c("p.adjust", "pvalue", "qvalue")) {
      ora_df = ora_df[order(ora_df[, order]), ]
    } else {
      ora_df = ora_df[order(ora_df[, order], decreasing = TRUE), ]
    }
    ora_df = ora_df[1:show_category, ]
    if (category %in% c("GeneRatio", "BgRatio")) {
      ratios = strsplit(ora_df[, category], "/")
      res = lapply(ratios, function(calcul_of_ratio) {
        x = as.numeric(calcul_of_ratio[1])
        y = as.numeric(calcul_of_ratio[2])
        return(round(x / y, 4))
      })
      ora_df[, category] = res %>% unlist()
    }
    ora_df$category = ora_df[, category]
    ggplot2::ggplot(ora_df, aes(x = category, y = fct_reorder(Description, category))) +
      ggplot2::geom_point(aes(color = p.adjust, size = RichFactor)) +
      ggplot2::scale_y_discrete(
        label = function(y)
          stringr::str_trunc(y, 40)
      ) +
      ggplot2::scale_color_gradientn (
        colours = gradient_col,
        trans = "log10",
        guide = ggplot2::guide_colorbar(reverse =
                                          TRUE, order = 1)
      ) +
      ggplot2::theme(axis.text = element_text(size = 6)) +
      ggplot2::labs(y = ylab, title = title, x = xlab) +
      DOSE::theme_dose(12) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = y_text_size))
  }
}

#' @description
#' dessine le dotplot issu d'un enrichissement
#' @param enrich : objet clusterprofiler
#' @param show_category : les n premieres categories a montrer
#' @param title : titre
#' @param category_label : taille des noms des categories
#'
#' @example draw_emapplot(ora.bp)
#'
draw_emapplot = function(enrich,
                         show_category = 10,
                         title = "Carte d'enrichissement",
                         category_label = 0.7) {
  fig = enrichplot::emapplot(
    enrichplot::pairwise_termsim(enrich),
    showCategory = show_category,
    cex.params = list(category_label = category_label)
  ) +
    ggplot2::ggtitle(title)
  fig %>% return()
}

#' Permet de désiner un arbre en calculant les distances de cluster entre les différents groupes de fonctionnalités, selon le nombre de gène en commun
#' @param enrich resultats d'enrichissement (gsea ou enrichment) de GO ou KEGG
#' @param show_category nombre de catégories à montrer
#' @param hilight Booléan, est-ce qu'il faut surligner les branches en fonction de la couleur des clusters
#' @param n_cluster nombre de cluster que l'on souhaite représenter
#' @param label_words_n nombre de mots dans les clusters à afficher
#' @param gradient_name nom à donner à la légende e l'échelle de couleurs
#' @param h_clust_method Méthode de distance utilisée pour le clustering
#' @param gradient_col vector de longueur 2, couleurs du gradient des p-adjust
#' @param gradient_name nom de la légende du gradient
#' @param title titre
#'
#' @example draw_treeplot(gse_bp,gradient_col = c("green","black"), show_category = 50,n_cluster = 10,label_words_n = 4,h_clust_method = "ward.D2")
#'
draw_treeplot = function(enrich,
                         show_category = 30,
                         hilight = FALSE,
                         n_cluster = show_category / 5,
                         label_words_n = 5,
                         h_clust_method = "ward.D2",
                         gradient_col = c("red", "blue"),
                         gradient_name = "adjusted p-value",
                         title = "Carte d'enrichissement") {
  if (length(gradient_col) == 2) {
    if (gradient_col[1] == gradient_col[2]) {
      gradient_col = c("red", "blue")
    }
  } else {
    gradient_col = c("red", "blue")
  }
  fig = enrichplot::treeplot(
    enrichplot::pairwise_termsim(enrich),
    showCategory = show_category,
    hilight.params = list(hilight = hilight),
    cluster.params = list(
      method = h_clust_method,
      n = n_cluster,
      label_words_n = label_words_n
    )
  ) +
    ggplot2::scale_color_gradient(name = gradient_name,
                                  low = gradient_col[1],
                                  high = gradient_col[2]) +
    ggplot2::ggtitle(title)
  fig %>% return()
}

#' @param enrich resultats d'enrichissement (gsea ou enrichment) de GO ou KEGG
#' @param gene_list gene_list passée en argument pour obtenir le enrich (avec ou sans abs)
#' @param metrique métrique selon laquelle on veut colorer les points (seulement si results)
#' @param gradient_name nom à donner à la légende e l'échelle de couleurs (seulement si results)
#' @param showCategory nombre de catégories à montrer
#' @param node_label détermine le nom des points à afficher
#' @param title titre
#' @param category_label taille en pourcentage des écritures
#' @param category_node taille en pourcentage des points
#' @param size le nom de la légende des points
#' @param category_color couleurs des points des catégories
#'
#' @example > draw_cnetplot(gse_bp,category_label = 0.6,gene_list = gene_list, category_color = "red", node_label = "none")
draw_cnetplot = function(enrich,
                         gene_list = NULL,
                         metrique = "stat",
                         gradient_name = "Associated\ndata",
                         showCategory = 5,
                         node_label = c("category", "gene", "all", "none"),
                         title = "Netplot of category",
                         category_label = 0.6,
                         category_node = 0.7,
                         size_name = "number of edge",
                         category_color = "black") {
  if (length(node_label) > 1) {
    node_label = node_label[1]
  }
  
  if (!is.null(gene_list)) {
    foldChange = gene_list
  } else {
    foldChange = NULL
  }
  fig = enrichplot::cnetplot(
    enrich,
    showCategory = showCategory,
    color.params = list(
      foldChange = foldChange,
      edge = TRUE,
      category = category_color
    ),
    node_label = node_label,
    cex.params = list(category_label = category_label, category_node = category_node)
  ) +
    ggplot2::labs(title = title, size = size_name)
  
  if (!is.null(gene_list)) {
    fig = fig +
      ggplot2::scale_color_gradient2(name = gradient_name,
                                     low = "darkgreen",
                                     high = "darkred")
  }
  fig %>% return()
}

draw_ridgeplot = function(gse,
                          xlab = "Distribution des enrichissements",
                          ylab = "Nom des voies",
                          title = "Distribution de l'expression selon les résultats de GSEA",
                          y_text_size = 7,
                          ...) {
  enrichplot::ridgeplot(gse) +
    ggplot2::labs(x = xlab,
                  y = ylab,
                  title = title) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = y_text_size))
}

error_message = function(name,
                         alpha_enrichissement) {
  paste(
    "\n\nAucune ",
    name,
    "enrichie n'a été trouvé de façon significative au seuil alpha",
    alpha_enrichissement,
    "pour permettre l'affichage de ce graphique\n\n"
  ) %>% return()
}