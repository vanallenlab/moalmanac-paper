devtools::install_github("hms-dbmi/UpSetR")
library('UpSetR')
library('grid')

setwd('~/Github/moalmanac-paper/analyses/knowledge-bases/civic-oncokb/')

drugs <- read.csv('drug-comparison.txt', sep='\t', header=T)
genes <- read.csv('gene-comparison.txt', sep='\t', header=T)
pmids <- read.csv('pmid-comparison.txt', sep='\t', header=T)

## PDF
pdf("~/Github/moalmanac-paper/extended-data-figures/extended-data-fig-4/extended-data-fig-4a.pdf", paper='a4')
upset(pmids, sets = c("CIViC", "MOAlmanac", "OncoKB"), sets.bar.color = "#56B4E9", order.by = "freq", text.scale=1.5,
      mainbar.y.label="PubMed IDs\nintersection size")
dev.off()

pdf("~/Github/moalmanac-paper/extended-data-figures/extended-data-fig-4/extended-data-fig-4b.pdf", paper='a4')
upset(drugs, sets = c("CIViC", "MOAlmanac", "OncoKB"), sets.bar.color = "#56B4E9", order.by = "freq", text.scale=1.5,
      mainbar.y.label="Therapies\nintersection size")
dev.off()

pdf("~/Github/moalmanac-paper/extended-data-figures/extended-data-fig-4/extended-data-fig-4c.pdf", paper='a4')
upset(genes, sets = c("CIViC", "MOAlmanac", "OncoKB"), sets.bar.color = "#56B4E9", order.by = "freq", text.scale=1.5,
      mainbar.y.label="Genes\nintersection size")
dev.off()

## PNG
png("~/Github/moalmanac-paper/extended-data-figures/extended-data-fig-4/extended-data-fig-4a.png", res=300, width=1000, height=750)
upset(pmids, sets = c("CIViC", "MOAlmanac", "OncoKB"), sets.bar.color = "#56B4E9", order.by = "freq", text.scale=1,
      mainbar.y.label="PubMed IDs\nintersection size")
dev.off()

png("~/Github/moalmanac-paper/extended-data-figures/extended-data-fig-4/extended-data-fig-4b.png", res=300, width=1000, height=750)
upset(drugs, sets = c("CIViC", "MOAlmanac", "OncoKB"), sets.bar.color = "#56B4E9", order.by = "freq", text.scale=1,
      mainbar.y.label="Therapies\nintersection size")
dev.off()

png("~/Github/moalmanac-paper/extended-data-figures/extended-data-fig-4/extended-data-fig-4c.png", res=300, width=1000, height=750)
upset(genes, sets = c("CIViC", "MOAlmanac", "OncoKB"), sets.bar.color = "#56B4E9", order.by = "freq", text.scale=1,
      mainbar.y.label="Genes\nintersection size")
dev.off()

