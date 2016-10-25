#edgeR dps RNA_seq data analysis

setwd('~/path/to/RNA_seq')
library('edgeR')

counts <- read.delim('all_conditions_counts.txt', sep='\t', header=T, stringsAsFactors=F)
group <- c(1,1,2,2,3,3,4,4)
y <- DGEList(counts=counts, group = group)
y$samples
       group lib.size norm.factors
WT41       1 22112058            1
WT42       1 63861448            1
WT241      2 20198850            1
WT242      2 18840134            1
d41      3 17323345            1
d42      3 17323345            1
d241     4 21237032            1
d242     4 20819735            1

#lib.size is a sum of all counts from one sample

#filter out lowly expressed genes (if lowly expressed in all conditions)

keep <- rowSums(cpm(y)>1) >=2
y <- y[keep, ,keep.lib.sizes=F]

#pairwise comparisons

#test for DE between 4hr WT and 4hr d:

et <- exactTest(y, pair=c("1","3"))
topTags(et)
Comparison of groups:  3-1 
         logFC    logCPM        PValue           FDR
ribB  3.893274 11.334951  0.000000e+00  0.000000e+00
dps  -3.080016 12.905690  0.000000e+00  0.000000e+00
astD -1.522062  9.440466 3.549882e-158 4.562782e-155
phoA -2.059049  5.372041 2.265466e-151 2.183909e-148
garK  1.455788  7.602721 3.866551e-148 2.981884e-145
ompF -1.751824  9.360731 5.628714e-148 3.617387e-145
astB -1.399226  8.776147 4.461848e-146 2.457841e-143
astA -1.524927  9.430567 9.564936e-145 4.610299e-142
yaiB -2.660329  6.288505 7.238038e-137 3.101097e-134
cydB  1.220191  9.602891 4.235693e-125 1.633283e-122

et4 <- exactTest(y, pair=c("1","3"))
comp4 <- et4$table
comp4_sig <- comp4[comp4$PValue < 0.05,]
length(comp4_sig$logFC)
#[1] 1574


#test for DE between 24hr WT and 24hr d:
et24 <- exactTest(y, pair=c("2","4"))
topTags(et)
Comparison of groups:  4-2 
          logFC    logCPM        PValue           FDR
ydcI -2.1762281  5.558690 1.911477e-212 7.370655e-209
dps  -1.9530969 12.905690 1.291264e-181 2.489557e-178
putA  2.0519035 10.634846 1.026807e-165 1.319789e-162
ynfM  1.4675942  8.137322 2.019695e-105 1.946986e-102
ompT -1.4435562  7.657235  4.193298e-93  3.233871e-90
malT -1.2424444 10.306173  4.342388e-90  2.790708e-87
hmp  -1.6359327  7.279399  1.527638e-87  8.415101e-85
cstA -1.0709137  8.754008  7.736593e-83  3.729038e-80
hdeA  1.4371081 12.517463  1.438751e-78  6.164248e-76
hepA  0.9462336  7.228994  2.145530e-65  8.273162e-63

comp24 <- et24$table
comp24_sig <- comp24[comp24$PValue < 0.05,]
length(comp24_sig$logFC)
#[1] 1889

plot(comp24$logFC, pch=20, xlab='genes with p-value < 0.05', ylab= 'log2 fold change', main='DE between 24hr WT and 24hr d')
abline(h = 2, col='red', lty=2)
abline(h = -2, col='red', lty=2)
abline(h = -1, col='red', lty=2)
abline(h = 1, col='red', lty=2)


comp24_sig <- comp24[comp24$PValue < 0.05,]
comp24_sig_neg <- comp24_sig[comp24_sig$logFC < -1,]
comp24_sig_pos <- comp24_sig[comp24_sig$logFC > 1,]
comp24_all <- rbind(comp24_sig_pos, comp24_sig_neg)
comp24_all

#test for DE between 24hr WT and 4hr WT:

etwt <- exactTest(y, pair=c("1","2"))
topTags(etwt)



compwt <- etwt$table
compwt_sig <- compwt[compwt$PValue < 0.05,]
length(compwt_sig$logFC)

comp4$index <- (1:length(comp4$logFC))
comp24$index <- (1:length(comp24$logFC))
compwt$index <- (1:length(compwt$logFC))
#[1] 1574

plot(compwt_all$index, compwt_all$logFC, col='black', xlab='gene index', ylab='log2FC', main='all DE genes')
>  points(comp24_all$index, comp24_all$logFC, col='blue')
>  points(comp4_all$index, comp4_all$logFC, col='red')
> legend('bottom', horiz=T, legend=c('d 4 v wt 4', 'd 24 v wt 24', 'wt 4 vs wt 24'), col=c('red', 'blue', 'black'), pch=1)

comp4_all_nonint <- comp4_all[!row.names(comp4_all) %in% wt4_inter,]

## plot along chromosome counts

counts4d1 <- read.delim('~/Desktop/4d1_counts.txt', header=F, stringsAsFactors=F, sep=',', col.names = c('index', 'gene', 'counts', 'ref'))

counts4d2 <- read.delim('~/Desktop/4d2_counts.txt', header=F, stringsAsFactors=F, sep=',', col.names = c('index', 'gene', 'counts', 'ref'))

counts24d1 <- read.delim('~/Desktop/24d1_counts.txt', header=F, stringsAsFactors=F, sep=',', col.names = c('index', 'gene', 'counts', 'ref'))

counts24d2 <- read.delim('~/Desktop/24d2_counts.txt', header=F, stringsAsFactors=F, sep=',', col.names = c('index', 'gene', 'counts', 'ref'))

counts4WT1 <- read.delim('~/Desktop/4WT1_counts.txt', header=F, stringsAsFactors=F, sep=',', col.names = c('index', 'gene', 'counts', 'ref'))

counts4WT2 <- read.delim('~/Desktop/4WT2_counts.txt', header=F, stringsAsFactors=F, sep=',', col.names = c('index', 'gene', 'counts', 'ref'))

counts24WT1 <- read.delim('~/Desktop/24WT1_counts.txt', header=F, stringsAsFactors=F, sep=',', col.names = c('index', 'gene', 'counts', 'ref'))

counts24WT2 <- read.delim('~/Desktop/24WT2_counts.txt', header=F, stringsAsFactors=F, sep=',', col.names = c('index', 'gene', 'counts', 'ref'))

layout(matrix(1:4))

plot(counts4d1$ref, counts4d1$counts, pch=20, main='4 d 1', xlab='reference position', ylab='counts', ylim=c(0,500000), type='h')
plot(counts4d2$ref, counts4d2$counts, pch=20, main='4 d 2', xlab='reference position', ylab='counts', ylim=c(0,500000), type='h')
plot(counts24d1$ref, counts24d1$counts, pch=20, main='24 d 1', xlab='reference position', ylab='counts', ylim=c(0,500000), type='h')
plot(counts24d2$ref, counts24d2$counts, pch=20, main='24 d 2', xlab='reference position', ylab='counts', ylim=c(0,500000), type='h')

plot(counts4WT1$ref, counts4WT1$counts, pch=20, main='4 WT 1', xlab='reference position', ylab='counts')
plot(counts4WT2$ref, counts4WT2$counts, pch=20, main='4 WT 2', xlab='reference position', ylab='counts' )
plot(counts24WT1$ref, counts24d1$counts, pch=20, main='24 WT 1', xlab='reference position', ylab='counts')
plot(counts24WT2$ref, counts24d2$counts, pch=20, main='24 WT 2', xlab='reference position', ylab='counts')


## # reads per sample ###
sum4d1 <- sum(counts4d1$counts)
sum4d2 <- sum(counts4d2$counts)
sum24d1 <- sum(counts24d1$counts)
sum24d2 <- sum(counts24d2$counts)
sum4WT1 <- sum(counts4WT1$counts)
sum4WT2 <- sum(counts4WT2$counts)
sum24WT1 <- sum(counts24WT1$counts)
sum24WT2 <- sum(counts24WT2$counts)


sample <- c('4d1', '4dp2', '24d1', '24d2', '4WT1', '4WT2', '24WT1', '24WT2')
sums <- c(sum4d1, sum4d2, sum24d1, sum24d2, sum4WT1, sum4WT2, sum24WT1, sum24WT2)
barplot(sums, names.arg=sample, main='# of reads per sample', col='salmon')

(848830, 849333
## d changes between d and WT at 4 hours ##



d <- c(counts4d1$counts[counts4d1$gene == 'd']/sum4d1,
counts4d2$counts[counts4d2$gene == 'd']/sum4d2,
counts24d2$counts[counts24d2$gene == 'd']/sum24d2,
counts24d1$counts[counts24d1$gene == 'd']/sum24d1,
counts24WT1$counts[counts24WT1$gene == 'd']/sum24WT1,
counts24WT2$counts[counts24WT2$gene == 'd']/sum24WT2,
counts4WT2$counts[counts4WT2$gene == 'd']/sum4WT2,
counts4WT1$counts[counts4WT1$gene == 'd']/sum4WT1)
barplot(d, names.arg=c('4d1', '4d2', '24d2', '24d1', '24WT1','24WT2', '4WT2','4WT1'), main='d counts / normalized', col='salmon')

## allD changes between d and WT at 4 hours ##

allD <- c(counts4d1$counts[counts4d1$gene == 'allD']/sum4d1,
counts4d2$counts[counts4d2$gene == 'allD']/sum4d2,
counts24d2$counts[counts24d2$gene == 'allD']/sum24d2,
counts24d1$counts[counts24d1$gene == 'allD']/sum24d1,
counts24WT1$counts[counts24WT1$gene == 'allD']/sum24WT1,
counts24WT2$counts[counts24WT2$gene == 'allD']/sum24WT2,
counts4WT2$counts[counts4WT2$gene == 'allD']/sum4WT2,
counts4WT1$counts[counts4WT1$gene == 'allD']/sum4WT1)
barplot(d, names.arg=c('4d1', '4d2', '24d2', '24d1', '24WT1','24WT2', '4WT2','4WT1'), main='allD counts / normalized', col='salmon')

barplot(counts4d1$counts, main='4 d 1', xlab='reference position', ylab='counts', xlim=c(848000, 850000), ylim=c(0,500000))

barplot(counts4d2$ref, counts4d2$counts, pch=20, main='4 d 2', xlab='reference position', ylab='counts', ylim=c(0,500000), xlim=c(848000, 850000),type='h')

barplot(counts4WT1$ref, counts4WT1$counts, xlim=c(848000, 850000),pch=20, main='4 WT 1', xlab='reference position', ylab='counts')
barplot(counts4WT2$ref, counts4WT2$counts, xlim=c(848000, 850000),pch=20, main='4 WT 2', xlab='reference position', ylab='counts' )

sapply(WT_d4$gene, function(x) {barplot(as.matrix(counts[row.names(counts) == x,][9:16]), main='normalized counts of dps', ylab='normalized counts', col='salmon')})


#operons for 4hour d/WT to check out

astABCDE, allDB, iscAS, fliTA, puuABCDEP

#operons for 24hr d/WT 

cysHIJ, gadABC, hdeABD, tnaAC, uxuAB, tdcAR, ompFT, fimIZ,yddAB, ydeMQ, ydeST,ydjNJ,
yheBD, ygcKW, yjfIK, yqeCJ



d_424 <- read.csv("~/Desktop/Lacey/DATA/RNA_seq/d_4_24_exactTest_withref.txt", stringsAsFactors=F, header=T)
WT_424 <- read.csv("~/Desktop/Lacey/DATA/RNA_seq/WT_4_24_exactTest_withref.txt", stringsAsFactors=F, header=T)


> both_intersect (4 hours with 4/24 d and WT)
 [1] "astA" "astB" "astC" "astD" "cspI" "entC" "fadB" "hisL" "hscB" "iscA" "iscS" "iscU" "puuA"
[14] "puuD" "ribB" "ycdM" "ydcI" "yigF" "yjiC" "ymfA" "allB" "cadA" "cydB" "garK" "gudD" "hypA"
[27] "phoA" "puuB" "puuC" "puuE" "tnaC" "yaiB" "ybbV" "ybgT"
> both_intersect <- intersect(time_intersect, WT_d24$gene)
> both_intersect (24 hours with 4/24 d and WT)
 [1] "cysI" "cysJ" "gntP" "proV" "uxuA" "uxuB" "ydcI" "ydeT" "yebV" "ygaR" "ygeF" "ylaC" "ynfM"
[14] "aphA" "blr"  "cbrC" "gadA" "gadB" "gadC" "hdeA" "hdeB" "hdeD" "hmp"  "malE" "malT" "mcrB"
[27] "ompT" "oppB" "paaB" "putA" "rbsB" "spr"  "srlA" "tdcA" "tnaA" "tnaC" "yaiB" "ybaJ" "ybbW"
[40] "ydgT" "yehD" "ygcK" "ygeW" "yqeC"

pdf('all_logFC_conditions.pdf')

df <- sapply(steve, function(x) {
	barplot(as.matrix(logFC_all[logFC_all$gene == x,][1:4]), ylim=c(-5,5), main=x, ylab='logFC',col='steelblue2') 
	abline(v=2.5, lwd=2, col='violetred4')
	abline(h=-1, lty=2)
	abline(h=1, lty=2)
	abline(h=0, lwd=2)
	text(1.25, 4, 'condition', col='violetred4')
	text(3.75, 4, 'time', col='violetred4')
     })
dev.off()	 
steve <- c("rpoN" ,"rpoD" ,"rpoS", "rpoH", "rpoN" ,"rpoE" ,"bfr",  "ftn" , "oxyR" ,"fecA" ,"fis" , "hns",  "hupA", "hupB", 
	  "ihfA" ,"ihfB" ,"fepE", 'fur', 'd')

