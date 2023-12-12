library(tidyverse)

#' @docType https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html
#'
#' @description
#' converti les ID d'une base de donnée dans une autre
#' @param results : objet results de DESeq pour les conditions voulant être testées
#' @param metric : métrique à utiliser pour l'analyse, par défaut : valeur absolu de la statistique de test
#' @param abs : booléan, donne si TRUE : valeur absolu, si FALSE : valeur réel
#' @param convert : booléen, faut il convertir les ID dans les ID d'une autre base, par défaut : FALSE
#' @param from : le nom du type d'ID de départ ("SYMBOL","ENSEMBL","ENTREZID", "UNIPROT"...)
#' @param to : le nom du type d'ID de fin ("SYMBOL","ENSEMBL","ENTREZID", "UNIPROT"...)
#' @param organism_db : nom de la base de donnée utilisée pour faire les échanges d'ID
#' @return vecteur trié dans l'ordre décroissant des gènes selon la métrique voulu
#'
#' @example load_gsea_GO_enrichment(results(dds),organism_db = "org.Hs.eg.db")
#'
prepare_enrichment = function(results,
                                   metric = "stat",
                                   abs = TRUE,
                                   convert = FALSE,
                                   from,
                                   to = "ENTREZID",
                                   organism_db) {
  # retrait de toutes les valeurs n'ayant pas passé l'independent filtering et la distance de Cook
  if ("padj" %in% colnames(results)) {
    results = results[!is.na(results[, "padj"]),]
  }
  if (convert) {
    ids = clusterProfiler::bitr(
      rownames(results),
      fromType = from,
      toType = to,
      OrgDb = organism_db
    )
    ids = ids[!duplicated(ids[c(from)]),]
    results = results[rownames(results) %in% ids[, from],]
    rownames(results) = ids[, to]
  }
  if (abs) {
    output = setNames(abs(results[, metric]), rownames(results))
  } else {
    output = setNames(results[, metric], rownames(results))
  }
  return(sort(na.omit(output), decreasing = T))
}

#' @description
#' Lance l'enrichissement sur GO en utilisant GSEA
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
                                   ont = c("MF", "CC", "BP"),
                                   key_type = "ENTREZID",
                                   min_GS_size = 3,
                                   max_GS_size = 500,
                                   pvalue_cutoff = 0.05,
                                   p_adjust_method = "BH",
                                   verbose = FALSE,
                                   n_perm_simple = 20000) {
  if(length(genes_DE) > 0){
    if (length(key_type) > 1) {
      key_type = key_type[1]
    }
    if(length(ont) > 1){
      ont = ont[1]
    }
    output <- clusterProfiler::gseGO(
      geneList = gene_list,
      ont = ont,
      keyType = key_type,
      minGSSize = min_GS_size,
      maxGSSize = max_GS_size,
      pvalueCutoff = pvalue_cutoff,
      verbose = verbose,
      OrgDb = organism_db,
      nPermSimple = n_perm_simple,
      eps = 0,
      pAdjustMethod = p_adjust_method)
    if(!is.null(output)){
      slot(output, "method") = "GSEA"
    }
    else {
      setClass("data1", representation(method = "character",result = "data.frame"))
      output =new("data1", method = "GSEA",result = data.frame(matrix(nrow = 0, ncol = 3)))
    }
  } else {
    setClass("data1", representation(method = "character",result = "data.frame"))
    output =new("data1", method = "GSEA",result = data.frame(matrix(nrow = 0, ncol = 3)))
  }
  return(output)
}

#' @description
#' Lance l'enrichissement sur KEGG en utilisant GSEA
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
#'
load_gsea_kegg_enrichment = function(gene_list = list(),
                                     organism_db = character(),
                                     key_type = c("ncbi-geneid", "kegg", "ncbi-proteinid", "uniprot"),
                                     min_GS_size = 3,
                                     max_GS_size = 500,
                                     pvalue_cutoff = 0.05,
                                     # n_perm = 20000,
                                     p_adjust_method = "BH") {
  if (length(key_type) > 1) {
    key_type = key_type[1]
  }
  output <- clusterProfiler::gseKEGG(
    geneList = gene_list,
    organism = organism_db,
    minGSSize = minGS_size,
    maxGSSize = maxGS_size,
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = p_adjust_method,
    keyType = key_type
  )
  slot(output, "method") = "GSEA"
  return(output)
}

#'
#' @description
#' Analyse la sur-représentation des GO terms dans le jeu de gènes différentiellement exprimés
#' @param genes_DE : vecteur de gènes différentielelment exprimés
#' @param universe : liste des gènes initiaux avant analyse DE (soit ça, soit OrgDb)
#' @param organism_db : base de donnée KEGG utilisée (exemple : "hsa") (soit ça, soit universe)
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

load_go_over_representation = function(genes_DE = list(),
                                       universe = NA,
                                       organism_db = NA,
                                       key_type = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                                       readable = FALSE,
                                       ont = c("MF", "CC", "BP"),
                                       p_adjust_method = "BH",
                                       min_GS_size = 10,
                                       max_GS_size = 500,
                                       pvalue_cutoff = 0.05,
                                       qvalue_cutoff = 0.10) {
  if(length(genes_DE) > 0){
    if (length(key_type) > 1) {
      key_type = key_type[1]
    }
    if(length(ont) > 1){
      ont = ont[1]
    }
    output = clusterProfiler::enrichGO(
      gene = genes_DE,
      universe = universe,
      OrgDb = organism_db,
      keyType = key_type,
      readable = readable,
      ont = ont,
      minGSSize = min_GS_size,
      maxGSSize = max_GS_size,
      pAdjustMethod = p_adjust_method,
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff
    )
    if(!is.null(output)){
      slot(output, "method") = "ORA"
      output = clusterProfiler::mutate(output, Subset = as.numeric(sub("/\\d+", "", BgRatio)))
      output = clusterProfiler::mutate(output, RichFactor = Count / Subset)
    }
    else {
      setClass("data1", representation(method = "character",result = "data.frame"))
      output =new("data1", method = "ORA",result = data.frame(matrix(nrow = 0, ncol = 3)))
    }
  } else {
    setClass("data1", representation(method = "character",result = "data.frame"))
    output =new("data1", method = "ORA",result = data.frame(matrix(nrow = 0, ncol = 3)))
  }
  return(output)
}


#' @description
#' Convertir un vecteur de gènes d'une base de données dans une autre
#' @param gene_list : vecteur de gènes différentielelment exprimés
#' @param from : type d'identifiant de base de votre liste ("ENTREZID", "SYMBOL"...)
#' @param to : type d'identifiant de sortie ("ENTREZID","SYMBOL"...)
#' @param organism_db : base de donnée KEGG utilisée (exemple : "hsa") (soit ça, soit universe)
#' @return vecteur avec les ID convertis
#'
#' @note peut être faire une fonction pour convertir à partir du fichier gtfs
#' @example id_convert(genes_DE,from = "SYMBOL",to = "ENTREZID", organism_db = "org.Hs.eg.db")
#'
id_convert = function(gene_list,
                      from,
                      to = "ENTREZID",
                      organism_db) {
  ids = clusterProfiler::bitr(gene_list,
                              fromType = from,
                              toType = to,
                              OrgDb = organism_db)
  return(ids[, to])
}

#' Analyse la sur-représentation des GO terms dans le jeu de gènes différentiellement exprimés
#' @param genes_DE : vecteur trié dans l'ordre décroissant des gènes selon la métrique voulu
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
#' @example load_go_over_representation(genes_DE,"org.Hs.eg.db")
#'
loag_kegg_over_representation = function(genes_DE,
                                         organism_db = character(),
                                         key_type = c("ncbi-geneid", "kegg", "ncbi-proteinid", "uniprot"),
                                         pvalue_cutoff = 0.05,
                                         p_adjust_method = "BH",
                                         universe = NA,
                                         min_GS_size = 10,
                                         max_GS_size = 500,
                                         qvalue_cutoff = 0.2,
                                         use_internal_data = FALSE) {
  if (length(key_type) > 1) {
    key_type = key_type[1]
  }
  output = clusterProfiler::enrichKEGG(
    gene = genes_DE,
    organism = organism,
    keyType = key_type,
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = p_adjust_method,
    universe = universe,
    minGSSize = min_GS_size,
    maxGSSize = max_GS_size,
    qvalueCutoff = qvalue_cutoff,
    use_internal_data = use_internal_data
  )
  slot(output, "method") = "ORA"
  return(output)
}

#' @description
#' permet l'affichage du tableau de données résultant de "get_enrich" ou d'une fonction ORA/GSEA de manière plus lisible en rmd
#' @param get_enrich l'objet clusterprofiler résultante d'une des fonctions au dessus permettant de faire du ORA ou du GSEA
#'
#' @example get_enrich(ora.bp,0.05) || get_enrich(ora.bp, ora.bp@pvalueCutoff)
enrich_pagination = function(get_enrich,
                             get_columns = c("ID",
                                             "Description",
                                             "GeneRatio",
                                             "BgRatio",
                                             "RichFactor",
                                             "p.adjust",
                                             "Count",
                                             "Subset"),
                             to_scientific = c("p.adjust"),
                             alpha_cutoff = 0.05) {
  if (slot(get_enrich, "method") == "ORA") {
    res = slot(get_enrich, "result")[,get_columns]
    res = res[res[,"p.adjust"] < alpha_cutoff,]
    res[to_scientific] = apply(X = res[to_scientific], MARGIN = 2,FUN = function(x) signif(x,3))
    res[to_scientific]  = format(res[to_scientific] ,scientific = T)
    res = format(res,digits = 3)
  } else if (slot(get_enrich, "method") == "GSEA") {
    res = slot(get_enrich, "result")[, c("ID", "Description", "setSize", "enrichmentScore")]
  }
  return(res)
}


#' @description
#' mets des retours à la ligne dans les noms des descriptions trop longues (supérieur à 4 mots)
#' @param enrichGO : objet enrich issus de clusterprofiler
#' @note modifier pour rendre ajustable le nombre de mots par lignes voulus
#'
#' @example description_spacer(ora.bp)
#'
description_spacer = function(enrichGO) {
  list_desc = list_desc_iter = slot(enrichGO,"result")[, "Description"]
  for (i in 1:length(list_desc_iter)) {
    aux <- str_locate_all(list_desc_iter[[i]], " ")[[1]]
    if (str_count(list_desc_iter[[i]], " ") > 4) {
      pos <-
        aux[which(diff(aux[, 1] %/% ((aux[nrow(aux), 1] %/% 2) + 1)) == 1) + 1, 1]
      substr(list_desc[[i]], pos[1], pos[1]) <- "\n"
    }
  }
  slot(enrichGO,"result")[, "Description"] = list_desc
  return(enrichGO)
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
draw_dotplot = function(ora_df,
                        order = "p.adjust",
                        category = "GeneRatio",
                        show_category = 10,
                        alpha = 0.05,
                        ylab = "Ontologies",
                        xlab = category,
                        title = paste("Dotplot of Gene Ontologies sorted by",order),
                        gradient_col = c("#f7ca64", "#46bac2", "#7e62a3"),
                        y_text_size = 10) {
  if (slot(ora_df, "method") != "ORA") {
    return("error")
  } else {
    #slot(ora_df,"result") = name_cleaver(slot(ora_df,"result"),"Description",limit_of_letter = 60)
    ora_df = slot(ora_df,"result")[slot(ora_df,"result")[,"p.adjust"] < 0.05,]
    
    if(order %in% c("p.adjust","pvalue","qvalue")){
      ora_df = ora_df[order(ora_df[,order] ),]
    } else {
      ora_df = ora_df[order(ora_df[,order],decreasing = TRUE),]
    }
    ora_df = ora_df[1:show_category,]
    if(category %in% c("GeneRatio","BgRatio")){
      ratios = strsplit(ora_df[,category],"/")
      res = lapply(ratios, function(calcul_of_ratio){
        x = as.numeric(calcul_of_ratio[1])
        y = as.numeric(calcul_of_ratio[2])
        return(round(x/y,4))
      })
      ora_df[,category] = res %>% unlist()
    }
    ora_df$category = ora_df[,category]
    ggplot2::ggplot(ora_df, aes(x = category,y = fct_reorder(Description,category))) +
      #ggplot2::geom_bar(stat = "identity", aes(fill = p.adjust)) +
      ggplot2::geom_point(aes(color = p.adjust, size = RichFactor)) +
      ggplot2::scale_y_discrete(label = function(y) stringr::str_trunc(y, 40))+
      ggplot2::scale_color_gradientn (
        colours = gradient_col,
        trans = "log10",
        guide = ggplot2::guide_colorbar(reverse =
                                          TRUE, order = 1)
      ) +
      ggplot2::theme(axis.text = element_text(size = 6))+
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
    showCategory = showCategory,
    cex.params = list(category_label = category_label)
  ) +
    ggplot2::ggtitle(title)
  return(fig)
}

#' Permet de désiner un arbre en calculant les distances de cluster entre les différents groupes de fonctionnalités, selon le nombre de gène en commun
#' @param enrich resultats d'enrichissement (gsea ou enrichment) de GO ou KEGG
#' @param showCategory nombre de catégories à montrer
#' @param hilight Booléan, est-ce qu'il faut surligner les branches en fonction de la couleur des clusters
#' @param n_cluster nombre de cluster que l'on souhaite représenter
#' @param label_words_n nombre de mots dans les clusters à afficher
#' @param gradient_name nom à donner à la légende e l'échelle de couleurs
#' @param h_clust_method Méthode de distance utilisée pour le clustering
#' @param gradient_col vector de longueur 2, couleurs du gradient des p-adjust
#' @param gradient_name nom de la légende du gradient
#' @param title titre
#'
#' @example draw_treeplot(gse_bp,gradient_col = c("green","black"),showCategory = 50,n_cluster = 10,label_words_n = 4,h_clust_method = "ward.D2")
#'
draw_treeplot = function(enrich,
                         showCategory = 30,
                         hilight = FALSE,
                         n_cluster = showCategory / 5,
                         label_words_n = 5,
                         h_clust_method = "ward.D2",
                         gradient_col = c("red", "blue"),
                         gradient_name = "p-value ajustée",
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
    showCategory = showCategory,
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
  return(fig)
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
#' @example > draw_cnetplot(gse_bp,category_label = 0.6,results = res_t114a_wt,gene_list = gene_list, category_color = "red", node_label = "none")
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
  slot(enrich,"result") = name_cleaver(slot(enrich,"result"),"Description",limit_of_letter = 40)
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
  return(fig)
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

prepare_list_enrichment = function(res.DE, 
                            information_from_gtf,
                            method,
                            metric = "stat",
                            absolute_val = FALSE,
                            unique_id = params$type_of_key){
  df = merge_gtf_and_res(res.DE,information_from_gtf, trimmed = FALSE,unique_id = unique_id)
  df = df[df[,"Symbol"] %in% get_DE(res.DE),]
  if(unique_id == "ENTREZID"){
    res = setNames(df[,metric],df[,"Entrez-ID"])
  } else {
    res = setNames(df[,metric],df[,"ENSEMBL-ID"])
  }
  if(absolute_val == TRUE){
    res = abs(res)
  }
  if(method == "ORA"){
    return(res)
  }
  return(res[order(res, decreasing = TRUE)])
}

error_message = function(name,
                         alpha_enrichissement){
  return(paste("\n\nAucune ",name, " enrichie n'a été trouvé de façon significative au seuil alpha", alpha_enrichissement, "pour permettre l'affichage de ce graphique\n\n"))
}