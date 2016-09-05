corum <- read.table("allComplexes.csv", sep="\t", h=T, quote="", comment.char="")
entrez <- strsplit(as.character(tab[,4]), ",")
corum[,5] <- unlist(lapply(entrez, length))
corum_formatted <- data.frame(V1 = rep(corum[,1], sapply(entrez, length)), V2 = unlist(entrez), V3=rep(corum[,5], sapply(entrez, length)))
