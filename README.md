# M2BIMS_Enrich_Analysis_Rshiny

## Goal

The goal of this application is to facilitates functional enrichment analysis from differential expression results, and allow for quick and interactive visualization.

------------------------------------------------------------------------

## Installation from GitHub

In a terminal :

```bash
git clone git@github.com:CamilleGodi/M2BIMS_Enrich_Analysis_Rshiny.git
```

------------------------------------------------------------------------

## Dependancies

-   BiocManager

-   clusterProfiler

-   DT

-   shiny

-   shinydashboard

-   shinyalert

-   plotly

-   ReactomePA

-   tidyverse

If not already existing on the computer, the program will try to install those libraries itself.

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

-   Add the corresponding line to the conversion table : `organism_conversion_table <- add_organism_in_conversion_table(organism_conversion_table, "Organism name", "org.Xx.xxx.db", "kegg_code", "reactome_name")`, with *kegg_code* the 3- or 4-letters KEGG organism code ( see [here](https://www.genome.jp/kegg/catalog/org_list.html) ), *org.Xx.xxx.db* the organism annotation library name, and *reactome_name* the name of the organism in the ReactomePA package (NA if not available. Available : 'celegans', 'fly', 'human', 'mouse', 'rat', 'yeast' and 'zebrafish').

------------------------------------------------------------------------
          
## Tutorial

1/ Select a CSV (or CSV2) file. It must have the following columns and none others : 'GeneName', 'ID', 'baseMean', 'log2FC', 'pval', 'padj'.

2/ Select the scientific name of the organism from which the data originates. Tip : you can type to search in the box.

3/ (Optional) Explore your data through the 'Whole data inspection' tab.

4/ Perform desired analysis through the appropriate tab.

Note : any plot can be downloaded by doing right-click > 'save image as'.
          
------------------------------------------------------------------------
          
## Authors

-   Victor BAILLEUL ( [victor.bailleul\@univ-rouen.fr](mailto:victor.bailleul@univ-rouen.fr) )

-   Camille GODI ( [camille.godi\@univ-rouen.fr](mailto:camille.godi@univ-rouen.fr) )

-   Benjamin MARSAC ( [benjamin.marsac\@univ-rouen.fr](mailto:benjamin.marsac@univ-rouen.fr) )

-   Komlan Dieu-Donné TOTO ( [komlan-dieu-donne.toto\@univ-rouen.fr](mailto:komlan-dieu-donne.toto@univ-rouen.fr) )

This app is the result of a group work in second year of Bioinformatics Master's Degree, 'BIMS', year 2023-2024s, Université de Rouen Normandie ( URN ).

## Acknowledgments :

We thank Solène Pety and Hélène Dauchel for their guidance and advices.
