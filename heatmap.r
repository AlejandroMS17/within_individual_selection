#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

inputf = args[1]
outputf = args[2]

library(reshape2)

#inputf = "~/Documents/github/within_individual_selection/k21/distances.tab"

d = read.delim(inputf, header=F)
d = d[,c(1,2,3)]
names(d) = c("sample1", "sample2", "distance")

dst = acast(d, sample1 ~ sample2)
dst = data.matrix(dst)
diag(dst) = NA # set diagonals to NA to avoid colour washout, HT: https://twitter.com/BEDecato/status/847646772453285890
dim = ncol(dst)

samples = rownames(dst)
splits = unlist(strsplit(samples, ".fastq.gz"))

rownames(dst) = splits
colnames(dst) = splits

pdf(file=outputf)

image(1:dim, 1:dim, dst, axes = FALSE, xlab="", ylab="")
axis(1, 1:dim, rownames(dst), cex.axis = 0.5, las=3)
axis(2, 1:dim, rownames(dst), cex.axis = 0.5, las=1)
text(expand.grid(1:dim, 1:dim), sprintf("%.4f", dst), cex=0.4)

dev.off()