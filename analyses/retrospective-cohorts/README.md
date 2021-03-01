# Retrospective cohorts
We evaluated 110 patients with metastatic melanoma, 150 patients with metastatic castration restistant prostate cancer, 100 patients with papillary renal cell carcinoma, and 59 pediatric patients with osteosarcoma for clinical actionability based on the Molecular Oncology Almanac and compared findings relative to PHIAL and TARGET. Tumor and normal whole-exome sequencing (WES) and RNA-seq were processed on the Broad Institute and Verily's [Terra platform](https://app.terra.bio/) that sits atop Google Cloud. 

Figure 2 and two supplementary figures showing the number of clinically relevant faetures highlighted by both MOAlmanac and PHIAL&TARGET are generated in `2.illustrate-retrospective-cohorts`. All calculations appearing in the paper's section, "Evaluating expanded molecular profiling and actionability in two retrospective cohorts" are performed in `3.calculations-retrospective-cohorts`.

Germline data has been redacted from all code and figures herein, so it is not entirely representative of calculations in the paper.

Figures produced and an excel sheet of the data used to produce them can be found in the `figures/figure-2/`, `figures/supplementary-figure-2/`, and `figures/supplementary-figure-3/` folder of this repository.

## References
1. [Van Allen, E. M. et al. Genomic correlates of response to CTLA-4 blockade in metastatic melanoma. Science 350, 207–211 (2015).](http://science.sciencemag.org/content/350/6257/207.long)
2. [Robinson, D. et al. Integrative Clinical Genomics of Advanced Prostate Cancer. Cell 162, 454 (2015).](https://www.sciencedirect.com/science/article/pii/S0092867415005486)
3. [Comprehensive Molecular Characterization of Papillary Renal-Cell Carcinoma. N Engl J Med 374, 135–145 (2016).](https://www.nejm.org/doi/full/10.1056/NEJMoa1505917)
4. [Perry, J. A. et al. Complementary genomic approaches highlight the PI3K/mTOR pathway as a common vulnerability in osteosarcoma. Proc Natl Acad Sci USA 111, E5564–E5573 (2014).](https://www.sciencedirect.com/science/article/pii/S0092867415005486)
5. [Van Allen, E. M. et al. Whole-exome sequencing and clinical interpretation of formalin-fixed, paraffin-embedded tumor samples to guide precision cancer medicine. Nat. Med. 20, 682–688 (2014).](https://www.nature.com/articles/nm.3559)
