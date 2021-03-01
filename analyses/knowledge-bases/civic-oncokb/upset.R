devtools::install_github("hms-dbmi/UpSetR")
library('UpSetR')
library('grid')

setwd('~/Github/moalmanac-paper/analyses/knowledge-bases/civic-oncokb/')

drugs <- read.csv('drug-comparison.txt', sep='\t', header=T)
genes <- read.csv('gene-comparison.txt', sep='\t', header=T)
pmids <- read.csv('pmid-comparison.txt', sep='\t', header=T)

png('~/Github/moalmanac-paper/figures/supplementary-figure-13/supplementary-figure-13.png', width=800, height=600)
upset(pmids, sets = c("CIViC", "MOAlmanac", "OncoKB"), sets.bar.color = "#56B4E9", order.by = "freq", text.scale=2)
grid.text("Comparison of PubMed IDs catalogued by\nCIViC, MOAlmanac, and OncoKB",x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()

png('~/Github/moalmanac-paper/figures/supplementary-figure-14/supplementary-figure-14.png', width=800, height=600)
upset(drugs, sets = c("CIViC", "MOAlmanac", "OncoKB"), sets.bar.color = "#56B4E9", order.by = "freq", text.scale=2)
grid.text("Comparison of therapies catalogued by\nCIViC, MOAlmanac, and OncoKB",x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()

png('~/Github/moalmanac-paper/figures/supplementary-figure-15/supplementary-figure-15.png', width=800, height=600)
upset(genes, sets = c("CIViC", "MOAlmanac", "OncoKB"), sets.bar.color = "#56B4E9", order.by = "freq", text.scale=2)
grid.text("Comparison of genes catalogued by\nCIViC, MOAlmanac, and OncoKB",x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()
