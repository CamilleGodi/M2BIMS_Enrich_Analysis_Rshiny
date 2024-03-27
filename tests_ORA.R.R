source("./utils/ORA.R")

# TMP : Tableau d'entré (non filtré etc), réactive "filtered_data"
table_filtered_tmp <- read.csv("./Example_files/filtered_data.csv", sep = ",")

# TMP : ext variables -> replace with the reactive variables
organism_library_go_tmp <- organism_conversion_table["Mus musculus", "annotation_db"]
p_value_cutoff_tmp <- 0.01
p_adj_cutoff_tmp   <- 0.05
q_value_cutoff_tmp <- 0.05
ontology_tmp <- "MF"


table_filtered_new_ids <- prepare_pipe(table_filtered_tmp, organism_db = organism_library_go_tmp, "ENSEMBL")

res_tmp <- do_ora_go_terms(
  table_filtered_new_ids,
  organism_library_go_tmp,
  ontology_tmp,
  p_value_cutoff_tmp,
  p_adj_cutoff_tmp,
  q_value_cutoff_tmp
)

res_tmp %>% enrich_pagination(alpha_cutoff = 0.05)

res_tmp %>% draw_dotplot(show_category = 30,
                         title = paste("ORA - GO termes -",
                                       ontology_tmp,
                                       "- Dot plot"))

res_tmp %>% draw_cnetplot(
  category_label = 0.6,
  gene_list = res_tmp@result$geneID,
  category_color = "red",
  node_label = "category",
  title = paste("ORA - GO termes -",
                ontology_tmp,
                "- CNET plot")
)

res_tmp %>% draw_treeplot(
  gradient_col = c("red", "blue"),
  show_category = 30,
  n_cluster = 10,
  label_words_n = 4,
  h_clust_method = "ward.D2",
  title = paste("ORA - GO termes -",
                ontology_tmp,
                "- Tree plot")
)

res_tmp %>% draw_emapplot(show_category = 10,
                          category_label = 0.5,
                          title = paste("ORA - GO termes -",
                                        ontology_tmp,
                                        "- EMAP plot"))