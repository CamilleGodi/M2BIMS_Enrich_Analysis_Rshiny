library(org.Mm.eg.db)
source("./utils/enrichissement.R")
tableau = read.csv("./Example_files/exemple.csv",sep = ";")
table_filtered = tableau %>% independent_filtering()
ids = conversion_table(table_filtered,from = "ENSEMBL","ENTREZID",organism_db = "org.Mm.eg.db")
table_filtered_new_ids = convert_results_ids(results = table_filtered,ids)

ora_ids = prepare_ora(table_filtered_new_ids)
kegg_ids = prepare_gsea(table_filtered_new_ids)
c = prepare_universe(table_filtered_new_ids)

gsea_go = load_gsea_GO_enrichment(b,organism_db = "org.Mm.eg.db",ont = "ALL",key_type = "ENTREZID")
ora_go = load_ora_go(a,c,"org.Mm.eg.db","ENTREZID","MF")
ora_go@result
gsea_kegg = load_gsea_kegg_enrichment(b,"mmu","ncbi-geneid",pvalue_cutoff = 1)
ora_kegg = load_ora_kegg(gene_list = a,organism_db = "mmu",key_type = "ncbi-geneid",universe = c)
ora_kegg_after_filter = filter_table_enrich_results(ora_kegg,p_value_cutoff = 0.2,0.3,0.3)
ora_kegg_after_filter %>% draw_cnetplot()
ora_go_after_filter = filter_table_enrich_results(ora_go,0.2,0.5,0.5)
ora_go_after_filter %>% draw_dotplot()
ora_go_after_filter
