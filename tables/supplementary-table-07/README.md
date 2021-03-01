# Supplementary Table 7
Profile-to-cell line matchmaking model performance.

Comparison of various models that were tested for matchmaking. Mean average precision, average precision at rank 1-5, and model description for all models used in study to assess genomic similarity between cancer cell lines for matchmaking.

Evaluation of significance of all models. Pairwise significance testing of all models evaluated. Our best performing model SNF: CGC & FDA was within the noise range of two other models, a multi-pass sort of first using agreement based measure of molecular features associated with an FDA approved therapy followed by agreement based sort of CGC genes mutated by any feature type (Multi-pass sort: FDA & CGC, p=0.415) and sorting cell lines by their mutant and wild type status of variants in order based on the somatic heuristic in MOAlmanac (Somatic tree, p=0.4397); however, SNF: CGC & FDA observed a stronger AP @ k = 1 in both cases, 0.191 versus 0.175 and 0.119, respectively. 
