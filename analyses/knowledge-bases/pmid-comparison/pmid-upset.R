devtools::install_github("hms-dbmi/UpSetR")
library('UpSetR')

handle = 'pmid-comparison.txt'
pmids <- read.csv(handle, sep='\t', header=T)

png('../../../figures/pmid-comparison.png', width=800, height=600)
upset(pmids, sets = c("CIViC", "MOAlmanac", "OncoKB"), sets.bar.color = "#56B4E9",
      order.by = "freq", text.scale=2)
dev.off()
