# M2BIMS_Enrich_Analysis_Rshiny

## Goal

The goal of this application is to facilitates functional enrichment analysis from differential expression results, and allow for quick and interactive visualization.

------------------------------------------------------------------------

## Dependancies

-   BiocManager

-   clusterProfiler

-   DT

-   shiny

-   shinydashboard

-   shinyalert

-   plotly

-   tidyverse

If not already existing on the computer, R will try to install those libraries itself.

------------------------------------------------------------------------

## Annotation libraries

Pre-installed annotation libraries from [Bioconductor](https://bioconductor.org/packages/3.18/data/annotation/) :

-   org.At.tair.db (*Arabidopsis thaliana*)

-   org.EcK12.eg.db (*Escherichia coli (K12)*)

-   org.Hs.eg.db (*Homo sapiens*)

-   org.Mm.eg.db" (*Mus musculus*)

-   org.Sc.sgd.db (*Saccharomyces cerevisiae*)

### How to add a new organism

In `ui.R` :

-   Add organism name in the list in line 36

In `global.R` :

-   Add the corresponding organism annotation library from [Bioconductor](https://bioconductor.org/packages/3.18/data/annotation/) in `load_libs_biocmanager`

-   Add the corresponding line to the conversion table : `organism_conversion_table <- add_organism_in_conversion_table(organism_conversion_table, "Organism name", "org.Xx.xxx.db", "kegg_code")`, with *kegg_code* the 3- or 4-letters KEGG organism code ([here](https://www.genome.jp/kegg/catalog/org_list.html)) and *org.Xx.xxx.db* the organism annotation library name

------------------------------------------------------------------------

## Authors

-   Victor BAILLEUL ( [victor.bailleul\@univ-rouen.fr](mailto:victor.bailleul@univ-rouen.fr){.email} )

-   Camille GODI ( [camille.godi\@univ-rouen.fr](mailto:camille.godi@univ-rouen.fr){.email} )

-   Benjamin MARSAC ( [benjamin.marsac\@univ-rouen.fr](mailto:benjamin.marsac@univ-rouen.fr){.email} )

-   Komlan Dieu-Donné TOTO ( [komlan-dieu-donne.toto\@univ-rouen.fr](mailto:komlan-dieu-donne.toto@univ-rouen.fr){.email} )

This app is the result of a group work in second year of Bioinformatics Master's Degree, 'BIMS', year 2023-2024s, Université de Rouen Normandie ( URN ).

## Acknowledgments :

We thank Solène Pety and Hélène Dauchel for their guidance and advices.
