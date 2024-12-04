rnaseq.dirs <- list.files('data/rnaseq/jovem')
print(rnaseq.dirs)

rnaseq.dirs[1]

raw_counts <- read.delim(paste0('./data/rnaseq/jovem/',  rnaseq.dirs[1], '/', list.files(paste0('data/rnaseq/jovem/', rnaseq.dirs[1]), pattern = '.tsv')), skip = 1)
raw_counts <- raw_counts[, c(1,4)]
head(raw_counts)

colnames(raw_counts)[2] <- rnaseq.dirs[1]
head(raw_counts)


maf.dirs <- list.files('data/maf/jovem')
print(maf.dirs)


maf <- read.delim(paste0('./data/maf/jovem/',  maf.dirs[1], '/', list.files(paste0('data/maf/jovem/', maf.dirs[1]), pattern = 'maf.gz')), skip = 7)
head(maf)
colnames(maf)
maf <- maf[, c(1,9,16)]
head(maf)