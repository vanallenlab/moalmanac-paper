# Comparison of PubMed IDs

PubMed IDs (PMIDs) for journal articles cited by the Molecular Oncology Almanac (MOAlmanac), CIViC, and OncoKB were obtained and compared using an Upset plot to look at overlap between sets. This analysis is limited as MOAlmanac also cites clinical guidelines and FDA approvals directly; however, the intent of the analysis is to illustrate that the literature space is vast and that no one knowledge base subsumes another. 

PMIDs were gathered and compared in the jupyter notebook `pmid-comparison.ipynb`. PMIDs from CIVIC were fetched via their API and were obtained for OncoKB relationships by listing the [300 biomarker-drug associations](https://www.oncokb.org/actionableGenes#levels=1,2,3,4,R1,R2) (65 genes, 52 tumor types, 6 levels of evidence; levels 1, 2, 3, 4, R1, and R2), downloading associations, and annotated for PMIDs through manual review. Summary statistics are then reported and visualized. A list of all PMIDs with a boolean for each database is present in `pmid-comparison.txt`. An UpSet plot is used to illustrate the intersection size of sets, `pmid-comparison.png` generated with `pmid-upset.R`. 


## References
1. [Griffith, M. et al. CIViC is a community knowledgebase for expert crowdsourcing the clinical interpretation of variants in cancer. Nat. Genet. 49, 170–174 (2017).](https://www.nature.com/articles/ng.3774)
2. [Chakravarty, D. et al. OncoKB: A Precision Oncology Knowledge Base. JCO Precis Oncol 2017, (2017).](https://ascopubs.org/doi/full/10.1200/po.17.00011)
3. [Jake R Conway, Alexander Lex, Nils Gehlenborg UpSetR: An R Package for the Visualization of Intersecting Sets and their Properties](https://doi.org/10.1093/bioinformatics/btx364)
4. [Pallarz, S. et al. Comparative Analysis of Public Knowledge Bases for Precision Oncology. JCO Precision Oncology 1–8 (2019).](https://ascopubs.org/doi/10.1200/PO.18.00371)
