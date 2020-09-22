# 2019 Sicklick
We compared clinical actions administered based on molecular profiles to patients in the [I-PREDICT](https://www.nature.com/articles/s41591-019-0407-5) study to those highlighted by the Molecular Oncology Almanac. This repository reproduces data processing from the original study, any calculations done on our part, as well as visualizations to display our findings. Jason Sicklick and authors did a _phenomenal_ job with preparing their paper supplement in a transparent manner and, because of their great efforts, we are able to compare our clinical interpretation methodology to real world evidence.

## Data processing
All variants considered were present in the supplementary text of the study. We extracted molecular features, therapies administered, and citations from the supplementary text. Disease ontologies were mapped from the supplementary text to [Oncotree](http://oncotree.mskcc.org/#/home) terms and codes, stored in the file `ontology.txt`. Variants were formatted using the notebook `0.format-variants-for-annotation.ipynb` to be annotated and evaluated with by the Molecular Oncology Almanac, `1.annotate-with-molecular-oncology-almanac.ipynb`. Formatted variants are found per patient in the `formatted_variants` folder and outputs from annotation with almanac are under `almanac_outputs`. 

## Review and annotation
After variants were annotated with the Molecular Oncology Almanac, three files were created for maunal review with the notebook `2.prepare-for-review.ipynb`. 

1. Citations were extracted from the supplementary text to the file `citations.not_annotated.txt`. They were then obtained, read, commented on, and categorized by [predictive implication level](https://moalmanac.org/about) and documented in the file `citations.annotated.txt`.

2. Molecular features considered by the study were merged with annotations made by the Molecular Oncology Almanac in `merged-features.not_annotated.txt` and then, using the author notes from the supplementary text, we annotated if the study targeted the feature in `merged-features.annotated.txt`

3. Therapies administered in the study and those highlighted by Molecular Oncology Almanac for therapeutic sensitivity were listed on a per patient basis in `therapies.not_annotated`. Predictive implication levels were annotated for each therapy per patient. For therapies highlighted by the Molecular Oncology Almanac, predictive implication levels were pulled from the outputs in `2.prepare-for-review`. For therapies administered in the study, the citations cited per patient were referenced again for the specific relationship between the therapeutic strategy or therapy and molecular features. Comments are written for all therapies per patient based on the citations cited. Each therapy administered was binned based on the evidence level or annotated as `no citation`, if the therapy was administered not on the basis of molecular features, or `citation listed not applicable`, if the citation(s) listed did not mention the therapy, strategy, or target. In some cases which would have resulted in the latter, we write in the comments that perhaps the authors forgot to cite a source and mention that in the comments. Therapies were then annotated for therapeutic strategy on a per case basis using expert view. Therapies were tagged with a boolean value if they were involved in a shared therapeutic strategy between what was administered in I-PREDICT and highlighted by the Molecular Oncology Almanac for the given patient.

## Results and visualization
Findings from the review and annotation were visualized and explored in `3.illustrate-sicklick-comparisons.ipynb`. In short, 

- (Supplementary figure) We visualize the number of therapies per patient that were either administered in I-PREDICT or highlighted by Molecular Oncology Almanac and color based on if the therapy is involved in a shared therapeutic strategy or not. 

- (A) Stacked bar charts are made of therapies administered to patients in I-PREDICT and highlighted by the Molecular Oncology Almanac by evidence level, and likewise annotated if they are involved in a shared therapeutic strategy. 

- (B) We create scatter plots to show the overlap of therapeutic strategies, specific therapies, and molecular features targeted. 

Any calculations used in the text of the paper pretaining to this topic are also located in this notebook.

## References
1. [Sicklick, J. K. et al. Molecular profiling of cancer patients enables personalized combination therapy: the I-PREDICT study. Nat. Med. 25, 744–750 (2019).](https://www.nature.com/articles/s41591-019-0407-5)
