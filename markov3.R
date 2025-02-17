#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(feather)
library(markovchain)

df <-arrow::read_feather(args[1]) # "/Users/mathieuouellet/Desktop/AMP/AMP/src/polygon/test.feather")
singleMc<-markovchainFit(data=df, possibleStates=c('0','1','2'), byrow = FALSE )
tm <- singleMc$estimate@transitionMatrix
std <- singleMc$standardError


print(tm)

#print(paste(tm[1,1],tm[1,2],tm[2,1],tm[2,2],std[1,1],std[1,2],std[2,1],std[2,2], sep = ","))