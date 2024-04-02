library(org.Mm.eg.db)
source("./utils/enrichment.R")
#' Tableau d'entré (non filtré etc), réactive "filtered_data"
tableau = read.csv("./Example_files/filtered_data.csv", sep = ",")
tableau %>% head()

# appel du filtre indépendant
table_filtered <- tableau %>% independent_filtering()

table_filtered_new_ids = prepare_pipe(table_filtered,organism_db = "org.Mm.eg.db","ENSEMBL")

ora_ids = prepare_ora(table_filtered_new_ids)
gsea_ids = prepare_gsea(table_filtered_new_ids, metric = "log2FC",abs = TRUE)
universe = prepare_universe(tableau,"org.Mm.eg.db","ENSEMBL")

res = load_gsea_reactome_enrichment(gsea_ids, organism_db = "mouse")
res = load_ora_reactome(ora_ids,universe = universe, organism_db = "mouse")
gsea_go = load_gsea_GO_enrichment(gsea_ids, organism_db = "org.Mm.eg.db")
gsea_go_2 = load_gsea_GO_enrichment(
  gsea_ids,
  organism_db = "org.Mm.eg.db",
  min_GS_size = 10,
  max_GS_size = 50
)
ora_go = load_ora_go(ora_ids, universe, "org.Mm.eg.db")
gsea_kegg = load_gsea_kegg_enrichment(gsea_ids, "mmu")
ora_kegg = load_ora_kegg(gene_list = ora_ids,
                         organism_db = "mmu",
                         universe = universe)
ora_kegg_after_filter = filter_table_enrich_results(ora_kegg, p_value_cutoff = 0.2, 0.3, 0.3)
ora_kegg_after_filter %>% draw_cnetplot()
ora_go_after_filter = filter_table_enrich_results(ora_go, 0.2, 0.5, 0.5)
ora_go_after_filter@result$ONTOLOGY %>% unique()
a = filter_go_enrich_results(ora_go_after_filter, ontology = c("MF", "CC"))
ora_go_after_filter %>% draw_dotplot(show_category = 30)
ora_kegg@result %>% head()
gsea_kegg@result %>% head()
gsea_kegg %>% enrichplot::dotplot()
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
enrichplot::gseaplot2(gsea_go,1:5,pvalue_table = TRUE)
gsea_go@result %>% dplyr::select(NES) %>% unlist() %>% sort()
clusterProfiler::dotplot(gsea_go)
gsea_go@result[,1:6]


res %>% draw_cnetplot()
res %>% draw_dotplot()
res %>% draw_emapplot()
gsea_go_2 %>% draw_gsea_plot()
