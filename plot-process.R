#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1) {
  stop("What file should I plot? (By default the name is [source datafile]-posterior-[HyperTraPS parameters].txt.process")
}

message("Loading libraries...")

library(ggplot2)

df = read.csv(args[1])

ggplot(df, aes(x=Time,y=OriginalIndex,size=Probability)) + geom_point()
