## DESeq2 Analysis
## requires R package 'DESeq2', and a counts table
## in this case, I used HT-seq to generate counts based on 'gene' features

setwd("~/path/to/RNA_seq/files")
library('DESeq2')
all_counts <- read.delim("all_conditions_raw_counts.txt", sep='\t', header=T)
all_counts <- all_counts[!row.names(all_counts) == 'dps',] #remove dps
samples <- read.csv('samples.csv', header=T)

## in design formula, the "control" goes first, and the "variable of interest" second

dds <- DESeqDataSetFromMatrix(countData=all_counts, colData = all_counts[,c(1:8)], design = ~condition)
dds$condition <- factor(dds$mutant, levels=c("WT","dps"))
samples$reps <- c('A',"A", 'B','B', 'C','C','D','D')
collapseReplicates(dds, groupby=samples$reps)
dds1 <- DESeq(dds)
res <- results(dds1)
