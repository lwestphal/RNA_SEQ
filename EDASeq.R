##EDASeq read qualityfiles <- 
setwd('~/Desktop/Lacey/DATA/RNA_seq/tophat_bams')
library('EDASeq')
library(ShortRead)
#pull in un-aligned fastq reads
files <- dir(pattern="*\\.fastq$")
names(files) <- gsub("\\.fastq.*", "", basename(files))
met <- DataFrame(conditions=c(rep("dps",2), rep("wt",2),rep("dps",2), rep("wt",2)),
row.names=names(files))
fastq <- FastqFileList(files)
elementMetadata(fastq) <- met
fastq
#FastqFileList of length 8
#names(8): 24dps1_ALL 24dps2_ALL 24WT1_ALL 24WT2_ALL 4dps1_ALL 4dps2_ALL 4WT1_ALL 4WT2_ALL

# pull in aligned reads
setwd('~/Desktop/Lacey/DATA/RNA_seq/tophat_bams')
files <- dir(path='./tophat_bams/', pattern="*\\.bam$")
names(files) <- gsub("\\.bam", "", basename(files))
gt <- gsub(".*/", "", files)
gt <- gsub("_.*", "", gt)
lane <- gsub(".*(.)$", "\\1", gt)
geno <- gsub(".$", "", gt)
pd <- DataFrame(geno=geno, lane=lane,
row.names=paste(geno,lane,sep="."))
bfs <- BamFileList(files)
elementMetadata(bfs) <- pd
bfs
#BamFileList of length 8
#names(8): 24dps1_accepted_hits 24dps2_accepted_hits 24WT1_accepted_hits ... 4dps2_accepted_hits #4WT1_accepted_hits 4WT2_accepted_hits

colors <- c(rep(rgb(1,0,0,alpha=0.7),2),
rep(rgb(0,0,1,alpha=0.7),2),
rep(rgb(0,1,0,alpha=0.7),2),
rep(rgb(0,1,1,alpha=0.7),2))
barplot(bfs,las=2,col=colors, names.arg = c( "24dps1","24dps2","24WT1","24WT2","4dps1","4dps2","4WT1", "4WT2"))

dev.new()
plotQuality(bfs,col=colors,lty=1)
legend("topright",unique(elementMetadata(bfs)[,1]), fill=unique(colors))


counts <- read.delim('all_conditions_raw_counts.txt', sep='\t', stringsAsFactors=F, header=T)

