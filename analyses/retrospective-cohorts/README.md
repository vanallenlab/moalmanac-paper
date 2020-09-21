# Retrospective cohorts
We evaluated 110 patients with metastatic melanoma and 150 patients with metastatic castration restistant prostate cancer for clinical actionability based on the Molecular Oncology Almanac and compared findings relative to PHIAL and TARGET. Tumor and normal whole-exome sequencing (WES) and RNA-seq were processed on the Broad Institute and Verily's [Terra platform](https://app.terra.bio/) that sits atop Google Cloud. 

Outputs, without germline, can be downloaded using the `download-from-terra.ipynb` jupyter notebook, located in the above directory. Files will be deposited into the `data/` subfolder of `2015-Robinson` and `2015-VanAllen` for the appropriate cohort. 

Figure 2 as well as a supplementary figure showing the number of clinically relevant features highlighted by both MOAlmanac and PHIAL&TARGET are generated in `1.illustrate-retrospective-cohorts`. All calculations appearing in the paper's section, "Evaluating expanded molecular profiling and actionability in two retrospective cohorts" are performed in `2.calculations-retrospective-cohorts`. 


## References
1. [Van Allen, E. M. et al. Genomic correlates of response to CTLA-4 blockade in metastatic melanoma. Science 350, 207–211 (2015).](http://science.sciencemag.org/content/350/6257/207.long)
2. [Robinson, D. et al. Integrative Clinical Genomics of Advanced Prostate Cancer. Cell 162, 454 (2015).](https://www.sciencedirect.com/science/article/pii/S0092867415005486?via%3Dihub)
3. [Van Allen, E. M. et al. Whole-exome sequencing and clinical interpretation of formalin-fixed, paraffin-embedded tumor samples to guide precision cancer medicine. Nat. Med. 20, 682–688 (2014).](https://www.nature.com/articles/nm.3559)
