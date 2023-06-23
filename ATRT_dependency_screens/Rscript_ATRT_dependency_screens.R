rm(list = ls())
library(openxlsx)
library(ggplot2)
library(pheatmap)
library(CRISPRcleanR)
library(MAGeCKFlute)
library(mixOmics)
setwd('directory')

r1 = read.delim('scores-FP_GPP655_Merk_CP0043_20190409_sgRNA.txt')
r2a = read.delim('scores-FP_GPP843_Merk_20190906_CP0043_Plate1.txt')
r2b = read.delim('scores-FP_GPP843_Merk_20190906_CP0043_Plate2.txt')
r2c = read.delim('scores-FP_GPP843_Merk_20190906_CP0043_Plate3.txt')
r2d = read.delim('scores-FP_GPP843_Merk_20190906_CP0043_Plate4.txt')
r2e = read.delim('scores-FP_GPP843_Merk_20190906_CP0043_Plate5.txt')
r2f = read.delim('scores-FP_GPP843_Merk_20190906_CP0043_Plate6.txt')

## Sanity Checks if all tables are structured equally
SanityCheck1 = all.equal.character(r1$Construct_Barcode, r2a$Construct_Barcode, r2b$Construct_Barcode,
                                   r2c$Construct_Barcode, r2d$Construct_Barcode, r2e$Construct_Barcode,
                                   r2f$Construct_Barcode)
SanityCheck2 = all.equal(nrow(r1), nrow(r2a), nrow(r2b), nrow(r2c), nrow(r2d), nrow(r2e),nrow(r2f))
SanityCheck3 = (nrow(r1) == 77441)

## Delete barcodes and IDs, short tables
s2a = r2a[3:ncol(r2a)]
s2b = r2b[3:ncol(r2b)]
s2c = r2c[3:ncol(r2c)]
s2d = r2d[3:ncol(r2d)]
s2e = r2e[3:ncol(r2e)]
s2f = r2f[3:ncol(r2f)]

## Check if one sample is on two plates
if (colnames(s2a)[ncol(s2a)] == colnames(s2b)[1]){ 
  s2ab = s2a[ncol(s2a)]+s2b[1]
  colnames(s2ab) = colnames(s2b)[1]
  s2a = cbind(s2a[1:(ncol(s2a)-1)],s2ab) 
  s2b = s2b[2:ncol(s2b)] 
}

if (colnames(s2b)[ncol(s2b)] == colnames(s2c)[1]){ 
  s2bc = s2b[ncol(s2b)]+s2c[1]
  colnames(s2bc) = colnames(s2c)[1]
  s2b = cbind(s2b[1:(ncol(s2b)-1)],s2bc) 
  s2c = s2c[2:ncol(s2c)] 
}

if (colnames(s2c)[ncol(s2c)] == colnames(s2d)[1]){ 
  s2cd = s2c[ncol(s2c)]+s2d[1]
  colnames(s2cd) = colnames(s2d)[1]
  s2c = cbind(s2c[1:(ncol(s2c)-1)],s2cd) 
  s2d = s2d[2:ncol(s2d)] 
}

if (colnames(s2d)[ncol(s2d)] == colnames(s2e)[1]){ 
  s2de = s2d[ncol(s2d)]+s2e[1]
  colnames(s2de) = colnames(s2e)[1]
  s2d = cbind(s2d[1:(ncol(s2d)-1)],s2de) 
  s2e = s2e[2:ncol(s2e)] 
}

if (colnames(s2e)[ncol(s2e)] == colnames(s2f)[1]){ 
  s2ef = s2e[ncol(s2e)]+s2f[1]
  colnames(s2ef) = colnames(s2f)[1]
  s2e = cbind(s2e[1:(ncol(s2e)-1)],s2ef) 
  s2f = s2f[2:ncol(s2f)] 
}

## Merge all
raw_read_counts = cbind(r1,s2a,s2b,s2c,s2d,s2e,s2f)

#get gRNA Brunello code
code = read.delim('gRNA_to_Brunello_code.txt')
raw_read_counts = cbind(code, raw_read_counts)

#read in quality data from PoolQ
q1 = read.delim('quality-FP_GPP655_Merk_CP0043_20190409_sgRNA.txt')
q2a = read.delim('quality-FP_GPP843_Merk_20190906_CP0043_Plate1.txt')
q2b = read.delim('quality-FP_GPP843_Merk_20190906_CP0043_Plate2.txt')
q2c = read.delim('quality-FP_GPP843_Merk_20190906_CP0043_Plate3.txt')
q2d = read.delim('quality-FP_GPP843_Merk_20190906_CP0043_Plate4.txt')
q2e = read.delim('quality-FP_GPP843_Merk_20190906_CP0043_Plate5.txt')
q2f = read.delim('quality-FP_GPP843_Merk_20190906_CP0043_Plate6.txt')

## Merge lists to one complete data.frame
quality.merge = rbind.data.frame(q1,q2a,q2b,q2c,q2d,q2e,q2f)
quality.merge = cbind.data.frame(quality.merge$Condition, quality.merge$Matched.Construct, quality.merge$Matched.Sample)
colnames(quality.merge) = c('Condition', 'Construct_Matched', 'Sample_Matched')
quality.sum = aggregate(. ~ Condition, quality.merge, sum)

## Add column with proportion matched to library
quality.sum = cbind.data.frame(quality.sum,(as.numeric(quality.sum$Construct_Matched)/as.numeric(quality.sum$Sample_Matched)))
quality.sum = cbind.data.frame(quality.sum,(as.numeric(quality.sum$Sample_Matched)-as.numeric(quality.sum$Construct_Matched)))
colnames(quality.sum) = c('Condition', 'Construct_Matched', 'Sample_Matched', 'Fraction_Matched', 'reads_unmatched')

write.xlsx(raw_read_counts, 'directory/raw_read_counts_all.xlsx')
write.xlsx(quality.sum, 'directory/Read_counts_overview.xlsx')

#graph for matched/unmatched reads
library(reshape2)
quality.sum.long = melt(quality.sum[,c(1,2,5)], id.vars = c("Condition"))
quality.sum.long[order(quality.sum.long$Condition),]
Line = c(rep("BT12", times = 3), rep("BT16", times = 3), rep("CHLA02", times = 3),rep("CHLA04", times = 3), rep("CHLA05", times = 3),rep("CHLA06", times = 3),rep("CHLA266", times = 3),
         rep("BT12", times = 3), rep("BT16", times = 3), rep("CHLA02", times = 3),rep("CHLA04", times = 3), rep("CHLA05", times = 3),rep("CHLA06", times = 3),rep("CHLA266", times = 3))
quality.sum.long = cbind(quality.sum.long, Line)
quality.sum.long[,"Condition"] = as.factor(quality.sum.long[,"Condition"])
str(quality.sum.long)

setwd('directory')
pdf('Matched_and_Unmatched_Read_Counts_of_all_Unique_Replicates.pdf', height = 10, width = 15)
ggplot(quality.sum.long, aes(x=Condition, y=value, fill=variable, group=Line))+
  geom_bar(stat="identity")+
  facet_wrap(~Line, scales = "free_x", nrow = 1)+
  scale_fill_manual(values=c('#3B638E', '#B29F65'))+
  geom_hline(yintercept=38220500, color = "red")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
dev.off()


####check replicate correlation
#get LFC_uncorrected data for sgRNA level
gRNA_LFCs = read.xlsx('LFCs_all.xlsx')

#get gene level LFCs
library(CRISPRcleanR)
data("Brunello_Library")
sgFC_CHLA266_C = gRNA_LFCs$CHLA266_REP_C
names(sgFC_CHLA266_C) = gRNA_LFCs$sgRNA
geneFC_CHLA266_C = ccr.geneMeanFCs(sgFC_CHLA266_C,Brunello_Library)

gene_LFCs = cbind(gene_LFCs,geneFC_CHLA266_C)

#calculate variance across all replicates
variance_geneLFCs = cbind(gene_LFCs, variance = apply(gene_LFCs[,-1], 1, var))
#sort according to variance and take top 3% high variance genes
sort_variance_geneLFCs = variance_geneLFCs[order(variance_geneLFCs$variance),]
high_variance_geneLFCs = sort_variance_geneLFCs[18542:19114,-22]

###calculate correlation, make heatmap and boxplots
library("Hmisc")
library(pheatmap)
cor_gene_LFCs <- rcorr(as.matrix(high_variance_geneLFCs),type = c("pearson"))
cor_gene_LFCs


flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
mat_cor_gene_LFCs = flattenCorrMatrix(cor_gene_LFCs$r,cor_gene_LFCs$P)
write.xlsx(mat_cor_gene_LFCs, 'directory/Replicate_correlation_matrix.xlsx')

#heatmap
setwd('directory/output')
pal = colorRampPalette(c("#D6D6D6", "#27408B"))

pdf('Heatmap_Pearson_correlation_replicates.pdf', height = 10, width = 11)
pheatmap(cor_gene_LFCs$r, border_color = NA, color = pal(1000))
dev.off()

#boxplot
Pearson.intra = read.delim('Pearson_r_intra.txt')
Pearson.inter = read.delim('Pearson_r_inter.txt')
data.Pearson = as.data.frame(qpcR:::cbind.na(Pearson.intra$cor, Pearson.inter$cor))
colnames(data.Pearson) = c("cor.intra", "cor.inter")
data.Pearson.long = melt(data.Pearson)
data.Pearson.long = na.omit(data.Pearson.long)

mean(Pearson.intra$cor) #0.7794298
pal(1000)[(mean(Pearson.intra$cor)-min(data.Pearson.long$value))/((max(data.Pearson.long$value)-min(data.Pearson.long$value))/1000)] #6676A6
mean(Pearson.inter$cor) #0.5235352
pal(1000)[(mean(Pearson.inter$cor)-min(data.Pearson.long$value))/((max(data.Pearson.long$value)-min(data.Pearson.long$value))/1000)] #B1B6C6

pdf('Boxplot_Pearson_r.pdf', height = 10, width = 6)
ggplot(data.Pearson.long, aes(x=variable, y= value, fill=variable)) +
  stat_boxplot(geom ='errorbar')+
  geom_boxplot()+
  scale_fill_manual(values=c("#6676A6", "#B1B6C6"))+
  theme_classic()
dev.off()

###CRISPRcleanR analyses
#load data
getwd()
setwd('directory')
sgraw_counts_all = read.xlsx('raw_read_counts_all.xlsx')
sgraw_counts_all = cbind(sgraw_counts_all[,"Code"], sgraw_counts_all[,"Gene"], sgraw_counts_all[,7:28])
colnames(sgraw_counts_all)[1] = "sgRNA"
colnames(sgraw_counts_all)[2] = "Gene"
essential_sgRNAs = read.delim("essential_sgRNAs.txt", header=FALSE)
essential_genes = read.delim("essential_genes.txt", header=FALSE)
nonessential_sgRNAs = read.delim("nonessential_sgRNAs.txt", header=FALSE)
nonessential_genes = read.delim("nonessential_genes.txt", header=FALSE)



#do normalization, LFC, and correction of LFCs and read counts individually for each replicate  
#make dataframe for individual replicates and plasmid first
counts_BT12_REP_A = cbind(sgraw_counts_all[,"sgRNA"], sgraw_counts_all[,"Gene"],sgraw_counts_all[,"plasmid"], sgraw_counts_all[,"BT12_REP_A"] )
colnames(counts_BT12_REP_A) = c("sgRNA", "Gene", "plasmid", "BT12_REP_A")
normANDfcs<-ccr.NormfoldChanges(Dframe=counts_BT12_REP_A,
                                saveToFig = TRUE,
                                libraryAnnotation=Brunello_Library,
                                EXPname='BT12_REP_A')
write.csv(normANDfcs$norm_counts, file="BT12_REP_A_norm_counts.csv")
write.csv(normANDfcs$logFCs, file="BT12_REP_A_logFCs.csv")

#correct LFCs
gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs$logFCs,Brunello_Library)
correctedFCs <- ccr.GWclean(gwSortedFCs)

write.csv(correctedFCs$corrected_logFCs, file="BT12_REP_A_corrected_LFCs.csv")

#correct read counts
correctedCounts<-ccr.correctCounts("BT12_REP_A",
                                   normANDfcs$norm_counts,
                                   correctedFCs,
                                   Brunello_Library)
write.csv(correctedCounts, file="BT12_REP_A_corrected_counts.csv")

##make dataframe with all norm_reads, LFCs, norm_reads_corrected und LFCs_corrected
setwd('directory')
sgnorm_counts_uncorrected = read.xlsx('norm_counts_all.xlsx')
sgLFCs_uncorrected = read.xlsx('LFCs_all.xlsx')
sgnorm_counts_corrected = read.xlsx('norm_counts_corrected.xlsx')
sgLFCs_corrected = read.xlsx('LFCs_corrected.xlsx')

####generate gene-level LFC uncorrected and corrected
gLFCs_uncorrected = aggregate(sgLFCs_uncorrected[,-(1:2)], list(sgLFCs_uncorrected$Gene),mean)
colnames(gLFCs_uncorrected)[1] = 'gene'
gLFCs_corrected = aggregate(sgLFCs_corrected[,-(1:2)], list(sgLFCs_corrected$Gene),mean)
colnames(gLFCs_corrected)[1] = 'gene'

####generate average across cell line replicates for gene-level LCFs
gLFCs_uncorrected_avg = gLFCs_uncorrected
colnames(gLFCs_uncorrected_avg) = c(rep("gene", times=1), rep("BT12", times = 3), rep("BT16", times = 3), rep("CHLA02", times = 3),
                                    rep("CHLA04", times = 3), rep("CHLA05", times = 3),rep("CHLA06", times = 3),rep("CHLA266", times = 3))
gLFCs_uncorrected_avg = as.data.frame(sapply(unique(names(gLFCs_uncorrected_avg[,-1])), 
                                             function(col) rowMeans(gLFCs_uncorrected_avg[names(gLFCs_uncorrected_avg) == col])))
gLFCs_uncorrected_avg = gLFCs_uncorrected_avg[colSums(!is.na(gLFCs_uncorrected_avg)) > 0]
gLFCs_uncorrected_avg = cbind(gLFCs_uncorrected$gene, gLFCs_uncorrected_avg)
colnames(gLFCs_uncorrected_avg)[1] = "gene"

gLFCs_corrected_avg = gLFCs_corrected
colnames(gLFCs_corrected_avg) = c(rep("gene", times=1), rep("BT12", times = 3), rep("BT16", times = 3), rep("CHLA02", times = 3),
                                  rep("CHLA04", times = 3), rep("CHLA05", times = 3),rep("CHLA06", times = 3),rep("CHLA266", times = 3))
gLFCs_corrected_avg = as.data.frame(sapply(unique(names(gLFCs_corrected_avg[,-1])), 
                                           function(col) rowMeans(gLFCs_corrected_avg[names(gLFCs_corrected_avg) == col])))
gLFCs_corrected_avg = gLFCs_corrected_avg[colSums(!is.na(gLFCs_corrected_avg)) > 0]
gLFCs_corrected_avg = cbind(gLFCs_corrected$gene, gLFCs_corrected_avg)
colnames(gLFCs_corrected_avg)[1] = "gene"
write.csv(gLFCs_corrected_avg, file="All_ATRTs_gene_level_LFC_avg.csv")

#Make subsets for essentials and nonessentials for each line average, uncorrected and corrected
essentials_uncor_BT12 = subset(gLFCs_uncorrected_avg[,c("gene", "BT12")], gene %in% essential_genes$V1)
nonessentials_uncor_BT12 =subset(gLFCs_uncorrected_avg[,c("gene", "BT12")], gene %in% nonessential_genes$V1)
essentials_cor_BT12 = subset(gLFCs_corrected_avg[,c("gene", "BT12")], gene %in% essential_genes$V1)
nonessentials_cor_BT12 =subset(gLFCs_corrected_avg[,c("gene", "BT12")], gene %in% nonessential_genes$V1)

essentials_uncor_BT16 = subset(gLFCs_uncorrected_avg[,c("gene", "BT16")], gene %in% essential_genes$V1)
nonessentials_uncor_BT16 =subset(gLFCs_uncorrected_avg[,c("gene", "BT16")], gene %in% nonessential_genes$V1)
essentials_cor_BT16 = subset(gLFCs_corrected_avg[,c("gene", "BT16")], gene %in% essential_genes$V1)
nonessentials_cor_BT16 =subset(gLFCs_corrected_avg[,c("gene", "BT16")], gene %in% nonessential_genes$V1)

essentials_uncor_CHLA02 = subset(gLFCs_uncorrected_avg[,c("gene", "CHLA02")], gene %in% essential_genes$V1)
nonessentials_uncor_CHLA02 =subset(gLFCs_uncorrected_avg[,c("gene", "CHLA02")], gene %in% nonessential_genes$V1)
essentials_cor_CHLA02 = subset(gLFCs_corrected_avg[,c("gene", "CHLA02")], gene %in% essential_genes$V1)
nonessentials_cor_CHLA02 =subset(gLFCs_corrected_avg[,c("gene", "CHLA02")], gene %in% nonessential_genes$V1)

essentials_uncor_CHLA04 = subset(gLFCs_uncorrected_avg[,c("gene", "CHLA04")], gene %in% essential_genes$V1)
nonessentials_uncor_CHLA04 =subset(gLFCs_uncorrected_avg[,c("gene", "CHLA04")], gene %in% nonessential_genes$V1)
essentials_cor_CHLA04 = subset(gLFCs_corrected_avg[,c("gene", "CHLA04")], gene %in% essential_genes$V1)
nonessentials_cor_CHLA04 =subset(gLFCs_corrected_avg[,c("gene", "CHLA04")], gene %in% nonessential_genes$V1)

essentials_uncor_CHLA05 = subset(gLFCs_uncorrected_avg[,c("gene", "CHLA05")], gene %in% essential_genes$V1)
nonessentials_uncor_CHLA05 =subset(gLFCs_uncorrected_avg[,c("gene", "CHLA05")], gene %in% nonessential_genes$V1)
essentials_cor_CHLA05 = subset(gLFCs_corrected_avg[,c("gene", "CHLA05")], gene %in% essential_genes$V1)
nonessentials_cor_CHLA05 =subset(gLFCs_corrected_avg[,c("gene", "CHLA05")], gene %in% nonessential_genes$V1)

essentials_uncor_CHLA06 = subset(gLFCs_uncorrected_avg[,c("gene", "CHLA06")], gene %in% essential_genes$V1)
nonessentials_uncor_CHLA06 =subset(gLFCs_uncorrected_avg[,c("gene", "CHLA06")], gene %in% nonessential_genes$V1)
essentials_cor_CHLA06 = subset(gLFCs_corrected_avg[,c("gene", "CHLA06")], gene %in% essential_genes$V1)
nonessentials_cor_CHLA06 =subset(gLFCs_corrected_avg[,c("gene", "CHLA06")], gene %in% nonessential_genes$V1)

essentials_uncor_CHLA266 = subset(gLFCs_uncorrected_avg[,c("gene", "CHLA266")], gene %in% essential_genes$V1)
nonessentials_uncor_CHLA266 =subset(gLFCs_uncorrected_avg[,c("gene", "CHLA266")], gene %in% nonessential_genes$V1)
essentials_cor_CHLA266 = subset(gLFCs_corrected_avg[,c("gene", "CHLA266")], gene %in% essential_genes$V1)
nonessentials_cor_CHLA266 =subset(gLFCs_corrected_avg[,c("gene", "CHLA266")], gene %in% nonessential_genes$V1)

## Plots of distributions of LFC
setwd('directory')
pdf('Distributions_uncorrected_LFCs_BT12.pdf', width = 10, height=10)
plot(density(essentials_uncor_BT12$BT12), xlim = c(-6,1.2), ylim = c(0, 1.1), xlab = 'logfold change', col = 'transparent',
     main = 'Distribution uncorrected LFCs Controls', bty = 'n')
polygon(density(essentials_uncor_BT12$BT12), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(nonessentials_uncor_BT12$BT12), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
clip(x1 = -6, x2 = 1.2, y1 = -0.05, y2 = 1)
abline(v = median(essentials_uncor_BT12$BT12), col = rgb(0.8, 0, 0), lty = 2, lw = 3)
abline(v = median(nonessentials_uncor_BT12$BT12), col = rgb(0.2117647, 0.3921569, 0.5450980), lty = 2, lw = 3)
legend(-6, 1, c('Ess', 'Noness'), lty = 1, lw = 5,
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980)), alpha = .6), box.lty = 0, cex = 2)
dev.off()

pdf('Distributions_corrected_LFCs_BT12.pdf', width = 10, height=10)
plot(density(essentials_cor_BT12$BT12), xlim = c(-6,1.2), ylim = c(0, 1.1), xlab = 'logfold change', col = 'transparent',
     main = 'Distribution corrected LFCs Controls', bty = 'n')
polygon(density(essentials_cor_BT12$BT12), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(nonessentials_cor_BT12$BT12), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
clip(x1 = -6, x2 = 1.2, y1 = -0.05, y2 = 1)
abline(v = median(essentials_cor_BT12$BT12), col = rgb(0.8, 0, 0), lty = 2, lw = 3)
abline(v = median(nonessentials_cor_BT12$BT12), col = rgb(0.2117647, 0.3921569, 0.5450980), lty = 2, lw = 3)
legend(-6, 1, c('Ess', 'Noness'), lty = 1, lw = 5,
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980)), alpha = .6), box.lty = 0, cex = 2)
dev.off()

setwd('directory')
pdf('Distributions_uncorrected_LFCs_BT16.pdf', width = 10, height=10)
plot(density(essentials_uncor_BT16$BT16), xlim = c(-6,1.2), ylim = c(0, 1.1), xlab = 'logfold change', col = 'transparent',
     main = 'Distribution uncorrected LFCs Controls', bty = 'n')
polygon(density(essentials_uncor_BT16$BT16), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(nonessentials_uncor_BT16$BT16), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
clip(x1 = -6, x2 = 1.2, y1 = -0.05, y2 = 1)
abline(v = median(essentials_uncor_BT16$BT16), col = rgb(0.8, 0, 0), lty = 2, lw = 3)
abline(v = median(nonessentials_uncor_BT16$BT16), col = rgb(0.2117647, 0.3921569, 0.5450980), lty = 2, lw = 3)
legend(-6, 1, c('Ess', 'Noness'), lty = 1, lw = 5,
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980)), alpha = .6), box.lty = 0, cex = 2)
dev.off()

pdf('Distributions_corrected_LFCs_BT16.pdf', width = 10, height=10)
plot(density(essentials_cor_BT16$BT16), xlim = c(-6,1.2), ylim = c(0, 1.1), xlab = 'logfold change', col = 'transparent',
     main = 'Distribution corrected LFCs Controls', bty = 'n')
polygon(density(essentials_cor_BT16$BT16), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(nonessentials_cor_BT16$BT16), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
clip(x1 = -6, x2 = 1.2, y1 = -0.05, y2 = 1)
abline(v = median(essentials_cor_BT16$BT16), col = rgb(0.8, 0, 0), lty = 2, lw = 3)
abline(v = median(nonessentials_cor_BT16$BT16), col = rgb(0.2117647, 0.3921569, 0.5450980), lty = 2, lw = 3)
legend(-6, 1, c('Ess', 'Noness'), lty = 1, lw = 5,
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980)), alpha = .6), box.lty = 0, cex = 2)
dev.off()

pdf('Distributions_uncorrected_LFCs_CHLA06.pdf', width = 10, height=10)
plot(density(essentials_uncor_CHLA06$CHLA06), xlim = c(-6,1.2), ylim = c(0, 1.1), xlab = 'logfold change', col = 'transparent',
     main = 'Distribution uncorrected LFCs Controls', bty = 'n')
polygon(density(essentials_uncor_CHLA06$CHLA06), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(nonessentials_uncor_CHLA06$CHLA06), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
clip(x1 = -6, x2 = 1.2, y1 = -0.05, y2 = 1)
abline(v = median(essentials_uncor_CHLA06$CHLA06), col = rgb(0.8, 0, 0), lty = 2, lw = 3)
abline(v = median(nonessentials_uncor_CHLA06$CHLA06), col = rgb(0.2117647, 0.3921569, 0.5450980), lty = 2, lw = 3)
legend(-6, 1, c('Ess', 'Noness'), lty = 1, lw = 5,
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980)), alpha = .6), box.lty = 0, cex = 2)
dev.off()

pdf('Distributions_corrected_LFCs_CHLA06.pdf', width = 10, height=10)
plot(density(essentials_cor_CHLA06$CHLA06), xlim = c(-6,1.2), ylim = c(0, 1.1), xlab = 'logfold change', col = 'transparent',
     main = 'Distribution corrected LFCs Controls', bty = 'n')
polygon(density(essentials_cor_CHLA06$CHLA06), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(nonessentials_cor_CHLA06$CHLA06), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
clip(x1 = -6, x2 = 1.2, y1 = -0.05, y2 = 1)
abline(v = median(essentials_cor_CHLA06$CHLA06), col = rgb(0.8, 0, 0), lty = 2, lw = 3)
abline(v = median(nonessentials_cor_CHLA06$CHLA06), col = rgb(0.2117647, 0.3921569, 0.5450980), lty = 2, lw = 3)
legend(-6, 1, c('Ess', 'Noness'), lty = 1, lw = 5,
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980)), alpha = .6), box.lty = 0, cex = 2)
dev.off()

## Calculate NNMD as measure for distinction between positive and negative controls per replicate
#make subsets of only essential and nonessential genes for corrected gene-level LFCs on replicate level
essentials_cor_replicates = subset(gLFCs_corrected, gene %in% essential_genes$V1)
nonessentials_cor_replicates = subset(gLFCs_corrected, gene %in% nonessential_genes$V1)

NNMDs = ((apply(essentials_cor_replicates[,-1], 2, mean) - apply(nonessentials_cor_replicates[,-1], 2, mean))
         /apply(nonessentials_cor_replicates[,-1], 2, sd))
NNMDs

## Calculate Cohen's D per replicate
## Transform LFCs in FCs
CohenDs = (apply((2^nonessentials_cor_replicates[,-1]), 2, mean) - apply((2^essentials_cor_replicates[,-1]), 2, mean)) / 
  sqrt((((nrow(nonessentials_cor_replicates)-1)*apply((2^nonessentials_cor_replicates[,-1]), 2, var))+
          (nrow(essentials_cor_replicates)-1)*apply((2^essentials_cor_replicates[,-1]), 2, var)) /
         (nrow(nonessentials_cor_replicates)+nrow(essentials_cor_replicates)-2))
CohenDs

## F-measure on cell line level (harmonic mean of precision and recall at BF 5)
getwd()
setwd("directory/pr_files_BAGEL2")
files = list.files() 
F_mes = data.frame()
for (ii in 1:length(files)) { 
  if (grepl('.txt', files[ii], fixed = TRUE)) { 
    temp = read.table(files[ii], header = TRUE) 
    tempBF = subset(temp, BF >= 5) 
    if (tempBF$Recall[nrow(tempBF)] == 0 && tempBF$Precision[nrow(tempBF)] == 0){ 
      F_mes[nrow(F_mes)+1,1] = 0
    }
    else {
      F_mes[nrow(F_mes)+1,1] = 2*(tempBF$Recall[nrow(tempBF)]*tempBF$Precision[nrow(tempBF)])/
        (tempBF$Recall[nrow(tempBF)]+tempBF$Precision[nrow(tempBF)])
      
    }
    rownames(F_mes)[length(rownames(F_mes))] = sub('_pr\\.txt', '', files[ii])
  }
}

## Build data frame with measures of separation of gene controls, and plot data
setwd("directory")
rep_separation = cbind(NNMDs, CohenDs)
line_separation = cbind.data.frame(rep_separation, sub('_REP_.', '', rownames(rep_separation)))
colnames(line_separation) = c('NNMDs_AVG','CohenDs_AVG','Cell_Line') 
line_separation = aggregate(line_separation[,-3], list(line_separation$Cell_Line), mean)
colnames(line_separation)[1] = 'Cell_Line'
line_separation = cbind(line_separation, F_mes)
colnames(line_separation)[4] = 'Fmeasure'
write.xlsx(as.data.frame(rep_separation), 'Separation_replicates.xlsx', rowNames = TRUE)
write.xlsx(line_separation, "Separation_lines.xlsx", rowNames = TRUE)

rep_separation_df = as.data.frame(rep_separation)
Cell_line1 = c(rep("BT12", times = 3), rep("BT16", times = 3), rep("CHLA02", times = 3),rep("CHLA04", times = 3), rep("CHLA05", times = 3),rep("CHLA06", times = 3),rep("CHLA266", times = 3))
Cell_line1 = as.factor(Cell_line1)
rep_separation_df =cbind(rep_separation_df, Cell_line1)
line_separation=line_separation[,-1]
Cell_line2 = c(rep("BT12", times = 1), rep("BT16", times = 1), rep("CHLA02", times = 1),rep("CHLA04", times = 1), rep("CHLA05", times = 1),rep("CHLA06", times = 1),rep("CHLA266", times = 1))
Cell_line2 = as.factor(Cell_line2)
line_separation = cbind(line_separation, Cell_line2)

########adjust color code here to match figure 1,, adjust axes
setwd('directory')
pdf('NNMD_to_Cohen.pdf', width = 11, height=10)
ggplot(rep_separation_df, aes(x=NNMDs, y=CohenDs, group=Cell_line1, colour=Cell_line1))+
  geom_point(cex=6)+
  theme_classic()+
  scale_color_manual(values = c("#881899",
                                "#268c41",
                                "#080808",
                                "#2eb8b3",
                                "#2c6cb0",
                                "#de8110",
                                "#8057d4"))
dev.off()

pdf('NNMD_to_Fmeasure.pdf', width = 11, height=10)
ggplot(line_separation, aes(x=NNMDs_AVG, y=Fmeasure, group=Cell_line2, colour=Cell_line2))+
  geom_point(cex=6)+
  theme_classic()+
  scale_color_manual(values = c("#881899",
                                "#268c41",
                                "#080808",
                                "#2eb8b3",
                                "#2c6cb0",
                                "#de8110",
                                "#8057d4"))
dev.off()

pdf('Cohen_to_Fmeasure.pdf', width = 11, height=10)
ggplot(line_separation, aes(x=CohenDs_AVG, y=Fmeasure, group=Cell_line2, colour=Cell_line2))+
  geom_point(cex=6)+
  theme_classic()+
  scale_color_manual(values = c("#881899",
                                "#268c41",
                                "#080808",
                                "#2eb8b3",
                                "#2c6cb0",
                                "#de8110",
                                "#8057d4"))
dev.off()

###work on averaged data on 6 good ATRT lines for phenotypical validation
#make avg of corrected sgRNA-level LFCs across replicates, then make avg across all 6 good ATRT lines

sgLFCs_corrected_avg = sgLFCs_corrected
colnames(sgLFCs_corrected_avg) = c(rep("sgRNA", times=1),rep("gene", times=1), rep("BT12", times = 3), rep("BT16", times = 3), rep("CHLA02", times = 3),
                                   rep("CHLA04", times = 3), rep("CHLA05", times = 3),rep("CHLA06", times = 3),rep("CHLA266", times = 3))
sgLFCs_corrected_avg = as.data.frame(sapply(unique(names(sgLFCs_corrected_avg[,-(1:2)])), 
                                            function(col) rowMeans(sgLFCs_corrected_avg[names(sgLFCs_corrected_avg) == col])))
sgLFCs_corrected_avg = sgLFCs_corrected_avg[colSums(!is.na(sgLFCs_corrected_avg)) > 0]
sgLFCs_corrected_avg = cbind(sgLFCs_corrected$sgRNA, sgLFCs_corrected_avg)
colnames(sgLFCs_corrected_avg)[1] = "sgRNA"

sgLFCs_corrected_6lineavg = sgLFCs_corrected_avg[,-4]
sgLFCs_corrected_6lineavg = data.frame(sgRNA = sgLFCs_corrected_6lineavg$sgRNA, mean_LFC = apply(sgLFCs_corrected_6lineavg[,-1],1, mean))
sgLFCs_corrected_6lineavg.string = as.vector(sgLFCs_corrected_6lineavg$mean_LFC)
names(sgLFCs_corrected_6lineavg.string) = sgLFCs_corrected_6lineavg$sgRNA

#make gene level string across all good ATRT cell lines
geneFC.good.lines<-ccr.geneMeanFCs(sgLFCs_corrected_6lineavg.string,Brunello_Library)

#use CRISPRcleanR to check depletion of known pan-essential signatures
data(EssGenes.ribosomalProteins)
data(EssGenes.DNA_REPLICATION_cons)
data(EssGenes.KEGG_rna_polymerase)
data(EssGenes.PROTEASOME_cons)
data(EssGenes.SPLICEOSOME_cons)

SIGNATURES<-list(Ribosomal_Proteins=EssGenes.ribosomalProteins,
                 DNA_Replication = EssGenes.DNA_REPLICATION_cons,
                 RNA_polymerase = EssGenes.KEGG_rna_polymerase,
                 Proteasome = EssGenes.PROTEASOME_cons,
                 Spliceosome = EssGenes.SPLICEOSOME_cons,
                 CFG = essential_genes$V1,
                 non_essential = nonessential_genes$V1)

pdf('Depletion_sigantures_in_good_ATRT_lines.pdf', width = 13, height=9)
Recall_scores<-ccr.VisDepAndSig(FCsprofile = geneFC.good.lines,
                                SIGNATURES = SIGNATURES,
                                TITLE = 'AVG_good_ATRT_lines',
                                pIs = 6,
                                nIs = 7)
dev.off()


#perform functional analysis on MAGeCK MLE of 6 ATRT cell lines
#read data
getwd()
setwd("directory")
#gene level
gdata.all.ATRT = ReadBeta("All_ATRTs.gene_summary.txt")
genelist.all.ATRT = gdata.all.ATRT$ATRT_common
names(genelist.all.ATRT) = gdata.all.ATRT$Gene
genelist.all.ATRT = sort(genelist.all.ATRT, decreasing = TRUE)
head(genelist.all.ATRT)

GSEA.KEGG = EnrichAnalyzer(genelist.all.ATRT, method = "GSEA", type = "KEGG", organism = "hsa")
gseaplot2(GSEA.KEGG, geneSetID = which(GSEA.KEGG$NES<0)[1:5], base_size = 12)


####generate scaled LFCs (Dependency scores)
#get dataframe with corrected LFCs of essential and non-essential genes, also name LFC dataframe
essentials_cor_all = cbind(essentials_cor_BT12, essentials_cor_BT16$BT16, essentials_cor_CHLA04$CHLA04,
                           essentials_cor_CHLA05$CHLA05, essentials_cor_CHLA06$CHLA06, essentials_cor_CHLA266$CHLA266)
rownames(essentials_cor_all) = essentials_cor_all$gene
essentials_cor_all= essentials_cor_all[,-1]
colnames(essentials_cor_all)[6] ="CHLA266"
nonessentials_cor_all = cbind(nonessentials_cor_BT12, nonessentials_cor_BT16$BT16, nonessentials_cor_CHLA04$CHLA04,
                              nonessentials_cor_CHLA05$CHLA05, nonessentials_cor_CHLA06$CHLA06, 
                              nonessentials_cor_CHLA266$CHLA266)
rownames(nonessentials_cor_all) = nonessentials_cor_all$gene
nonessentials_cor_all= nonessentials_cor_all[,-1]
colnames(nonessentials_cor_all)[6] ="CHLA266"

gLFCs_final = gLFCs_corrected_avg
rownames(gLFCs_final) = gLFCs_final$gene
gLFCs_final = gLFCs_final[,-1]
gLFCs_final = gLFCs_final[,-3]

#scale so median of nonessentials is 0
nonessentials_median = apply(nonessentials_cor_all, 2, median)
nonessentials_median_df = as.data.frame(matrix(rep(nonessentials_median, 
                                                   each = nrow(gLFCs_final)), nrow = nrow(gLFCs_final)))
colnames(nonessentials_median_df) = colnames(gLFCs_final)
LFC_scale1 = gLFCs_final - nonessentials_median_df #step 1 of scaling

#scale so median of essential genes is -1
essentials_median = apply(essentials_cor_all, 2, median)
essentials_median_df = as.data.frame(matrix(rep(essentials_median, 
                                                each = nrow(LFC_scale1)),
                                            nrow = nrow(LFC_scale1)))
colnames(essentials_median_df) = colnames(LFC_scale1)
LFC_scale2 = LFC_scale1/(-essentials_median_df) # step 2 of scaling
Dep_score = LFC_scale2
setwd("directory")
write.csv(Dep_score, file="Dependency_scores_ATRTs.csv")


### perform integration of dependency, mRNA, and methylation data using mixOMICS package
###get data for ATRT cells
getwd()
setwd("directory")
#RNAseq data
data.RNAseq = read.csv("RNAseq_1.5STD.csv", header = T, sep= ",")
rownames(data.RNAseq) = data.RNAseq$X
data.RNAseq = data.RNAseq[,-1]
data.RNAseq.log = log(data.RNAseq, 2)
data.RNAseq.log = t(data.RNAseq.log)
dim(data.RNAseq.log)


#methylation data
data.methyl = read.csv("Methylation_1.5STD.csv", header = T, row.names = 1)
data.methyl = t(data.methyl)
dim(data.methyl)

#generate random methylation data

data.methyl.random = matrix(runif(68424),nrow=6)

#dependency data
data.DepScore = read.csv("DepScore_1.5STD.csv", header = T, row.names = 1)
data.DepScore = t(data.DepScore)
dim(data.DepScore)

#compile into one list
x <- list(RNASeq = data.RNAseq.log, Methylation = data.methyl, Dependency = data.DepScore)
x.random <- list(RNASeq = data.RNAseq.log, Methylation = data.methyl.random, Dependency = data.DepScore)
lapply(x, dim)
str(x)

# select arbitrary values of features to keep
list.keepX = c(100, 100)
list.keepY = c(100, 100)

# generate three pairwise PLS models
pls.RNAseq_methyl <- spls(x[["RNASeq"]], x[["Methylation"]], keepX = list.keepX, keepY = list.keepY)

pls.RNAseq_methyl_random <- spls(x.random[["RNASeq"]], x.random[["Methylation"]], keepX = list.keepX, keepY = list.keepY)

pls.RNAseq_DepScore <- spls(x[["RNASeq"]], x[["Dependency"]], keepX = list.keepX, keepY = list.keepY)

pls.Methylation_DepScore <- spls(x[["Methylation"]], x[["Dependency"]], keepX = list.keepX, keepY = list.keepY)
pls.Methylation_random_DepScore <- spls(x.random[["Methylation"]], x.random[["Dependency"]], keepX = list.keepX, keepY = list.keepY)

# plot features of first PLS
pair.Var.RNAseq_methyl = plotVar(pls.RNAseq_methyl, cutoff = 0.7, title = "RNASeq vs Methylation",
                                 legend = c("RNA", "Methylation"),
                                 var.names = F, style = 'graphics',
                                 pch = c(15, 16), cex = c(2.2,2),
                                 col = c("#0072B2", "#009E73"))
pair.Var.RNAseq_methyl = plotVar(pls.RNAseq_methyl, cutoff = 0.7, title = NA,
                                 var.names = F, style = 'graphics',
                                 pch = c(15, 16), cex = c(2.2,2),
                                 col = c("#0072B2", "#009E73"))
pair.Var.RNAseq_methyl.random = plotVar(pls.RNAseq_methyl_random, cutoff = 0.7, title = "RNASeq vs Methylation_random",
                                        legend = c("RNA", "Methylation_random"),
                                        var.names = F, style = 'graphics',
                                        pch = c(15, 16), cex = c(2.2,2),
                                        col = c("#0072B2", "#009E73"))

pair.Var.RNAseq_DepScore = plotVar(pls.RNAseq_DepScore, cutoff = 0.7, title = "RNASeq vs Dependency",
                                   legend = c("RNA", "Dependency"),
                                   var.names = F, style = 'graphics',
                                   pch = c(15, 17), cex = c(2.2,2),
                                   col = c("#0072B2", "#CC79A7"))
pair.Var.RNAseq_DepScore = plotVar(pls.RNAseq_DepScore, cutoff = 0.7, title = NA,
                                   var.names = F, style = 'graphics',
                                   pch = c(15, 17), cex = c(2.2,2),
                                   col = c("#0072B2", "#CC79A7"))

pair.Var.Methyl_DepScore = plotVar(pls.Methylation_DepScore, cutoff = 0.7, title = "Methylation vs Dependency",
                                   legend = c("Methylation", "Dependency"),
                                   var.names = F, style = 'graphics',
                                   pch = c(16, 17), cex = c(2.2,2),
                                   col = c("#009E73", "#CC79A7"))
pair.Var.Methyl_DepScore = plotVar(pls.Methylation_DepScore, cutoff = 0.7, title = NA,
                                   var.names = F, style = 'graphics',
                                   pch = c(16, 17), cex = c(2.2,2),
                                   col = c("#009E73", "#CC79A7"))
pair.Var.Methyl_random_DepScore = plotVar(pls.Methylation_random_DepScore, cutoff = 0.7, title = "Methylation_random vs Dependency",
                                          legend = c("Methylation_random", "Dependency"),
                                          var.names = F, style = 'graphics',
                                          pch = c(16, 17), cex = c(2.2,2),
                                          col = c("#009E73", "#CC79A7"))

setwd("directory")
write.table(pair.Var.RNAseq_methyl, file = "Pairwise_PLS_RNAseq-Methylation_100featuresperComponent.csv", sep = ",")
write.table(pair.Var.RNAseq_DepScore, file = "Pairwise_PLS_RNAseq-DepScore_100featuresperComponent.csv", sep = ",")
write.table(pair.Var.Methyl_DepScore, file = "Pairwise_PLS_Methylation-DepScore_100featuresperComponent.csv", sep = ",")

#multiblock analysis

design = matrix(1, ncol = length(x), nrow = length(x), # for square matrix filled with 0.5s
                dimnames = list(names(x), names(x)))
diag(design) = 0 # set diagonal to 0s

x.2 = list(RNASeq = x$RNASeq, 
           Methylation = x$Methylation)

list.keepX = list(RNASeq = rep(50, 2), Methylation = rep(50,2))
list.keepY = c(rep(50, 2))

final.mbspls.model = block.spls(X = x.2, Y = x$Dependency,
                                ncomp = 2, keepX = list.keepX,
                                keepY = list.keepY, design = 'full')

line = factor(final.mbspls.model$names$sample, levels = 1:6)

Var.plot.MB = plotVar(final.mbspls.model, var.names = FALSE,
                      legend = TRUE, cutoff = 0.7, style = "graphics",
                      pch = c(15,16,18), cex = c(2.2,2,1.8), col = c("#0072B2", "#009E73", "#CC79A7"))


Var.plot.MB.legend = plotVar(final.mbspls.model, var.names = FALSE,title = NA,
                             legend = F, cutoff = 0.7, style = "graphics",
                             pch = c(15,16,18), cex = c(2.2,2,1.8), col = c("#0072B2", "#009E73", "#CC79A7"))

setwd("directory")
write.table(Var.plot.MB, file = "Multiblock_PLS_RNAseq-Methylation-DepScore_100featuresperComponent.csv", sep = ",")


circosPlot.Multiblock = circosPlot(final.mbspls.model,
                                   cutoff = 0.9, group=line,
                                   Y.name = 'Y', blocks.link="Y", color.blocks=c("#0072B2", "#009E73", "#CC79A7"), size.variables=0.01, var.adj=1)
circosPlot(final.mbspls.model,cutoff = 0.9, group=line,
           Y.name = 'Y', blocks.link="Y", color.blocks=c("#0072B2", "#009E73", "#CC79A7"), 
           size.variables=0.01, var.adj=1, size.legend=0, size.labels=0)



###perform Pearson R correlation analysis for context-specific essentials RNASeq vs DepScore
#get expression and dependency data for context specific ATRT essentials
getwd()
setwd("directory")
data.exp.norm.Pearson = read.csv("Context_specific_essentials_norm_expression.csv", header = T, row.names = 1)
data.exp.norm.Pearson = data.exp.norm.Pearson + 1
data.exp.norm.Pearson = log(data.exp.norm.Pearson, 2)
data.exp.norm.Pearson = t(data.exp.norm.Pearson)
data.exp.norm.Pearson = as.data.frame(data.exp.norm.Pearson)

data.DepScore.Pearson = read.csv("Context_specific_essentials_DepScore.csv", header = T, row.names = 1)
data.DepScore.Pearson = t(data.DepScore.Pearson)
data.DepScore.Pearson = as.data.frame(data.DepScore.Pearson)

#linear correlation
cor.dep.exp <- sapply(1:ncol(data.exp.norm.Pearson), function(i) cor(data.exp.norm.Pearson[,i],data.DepScore.Pearson[,i]))
setwd("directory")
write.csv(cor.dep.exp,file="Pearson_r_Dep_Exp.csv")

# linear correlation of null distribution by bootstrap resampling the values in the dependency dataset
data.DepScore.Pearson.perm <- sapply(data.DepScore.Pearson,function(x){ sample(x) }) 
data.DepScore.Pearson.perm <- as.data.frame(data.DepScore.Pearson.perm)
cor.dep.exp.perm <- sapply(1:ncol(data.exp.norm.Pearson), function(i) cor(data.exp.norm.Pearson[,i],data.DepScore.Pearson.perm[,i]))
setwd("directory")
write.csv(cor.dep.exp.perm,file="Pearson_r_Dep_Exp.NULL.csv")

#make graph
Pearson_R = read.csv("Pearson_r_Dep_Exp.csv", header = T, row.names = 1)
Pearson_R = na.omit(Pearson_R)
Pearson_R_null = read.csv("Pearson_r_Dep_Exp.NULL.csv", header = T, row.names = 1)
Pearson_R_null = na.omit(Pearson_R_null)

median(Pearson_R$r)

p_null = density(Pearson_R_null$r, bw=0.15)
h = hist(Pearson_R$r,20, freq =TRUE, xlim=c(-1.5,1.5), yaxt="n", xaxt="n", ylim = c(0,150), col = "#0072B2", xlab = "Pearson r", main ='Correlation Expression Dependency')
clip(-0.09, -0.07,0,150)
abline(v=-0.08653528, col=c("black"), lty = c(2), lwd=c(3))
axis(2, pos=-1.5, lwd = 1)
axis(1, at=c(-1, 0, 1), lwd = 1)
w = (h$counts/ h$density)[1]
lines(p_null$x, p_null$y*w, col = "red", lty=2, lw=3)
legend('topright', 'null distribution\n(random permutation)', col = "red", lw=3,lty=2,bty='n', cex = 0.8)

#Wilcoxon rank sum (non-parametric test without normality assumption) test for RNASeq vs DepScore
#check for normality using Shapiro-Wilk test
hist(Pearson_R$r)
shapiro.test(Pearson_R$r)
hist(Pearson_R_null$r)
shapiro.test(Pearson_R_null$r)
wilcox.test(Pearson_R$r, Pearson_R_null$r)
#W = 1291952, p-value = 9.45e-05
#alternative hypothesis: true location shift is not equal to 0


###perform Pearson R correlation analysis for context-specific essentials for Methylation vs DepScore
#get methylation data for context specific ATRT essentials
getwd()
setwd("directory")
data.Methylation.Pearson = read.csv("Context_specific_essentials_Methylation.csv", header = T, row.names = 1)
data.Methylation.Pearson = t(data.Methylation.Pearson)
data.Methylation.Pearson = as.data.frame(data.Methylation.Pearson)

data.DepScore.Pearson.2 = read.csv("Context_specific_essentials_DepScore_for_Methylation.csv", header = T, row.names = 1)
data.DepScore.Pearson.2 = t(data.DepScore.Pearson.2)
data.DepScore.Pearson.2 = as.data.frame(data.DepScore.Pearson.2)

#linear correlation
cor.dep.methyl <- sapply(1:ncol(data.DepScore.Pearson.2), function(i) cor(data.Methylation.Pearson[,i],data.DepScore.Pearson.2[,i]))
setwd("directory")
write.csv(cor.dep.methyl,file="Pearson_r_Dep_Methyl.csv")

# linear correlation of null distribution by bootstrap resampling the values in the dependency dataset
data.DepScore.Pearson.perm.2 <- sapply(data.DepScore.Pearson.2,function(x){ sample(x) }) 
data.DepScore.Pearson.perm.2 <- as.data.frame(data.DepScore.Pearson.perm.2)
cor.dep.methyl.perm <- sapply(1:ncol(data.DepScore.Pearson.2), function(i) cor(data.Methylation.Pearson[,i],data.DepScore.Pearson.perm.2[,i]))
setwd("directory")
write.csv(cor.dep.methyl.perm,file="Pearson_r_Dep_Methyl.NULL.csv")

#make graph
Pearson_R_Dep_Methyl = read.csv("Pearson_r_Dep_Methyl.csv", header = T, row.names = 1)
Pearson_R_Dep_Methyl = na.omit(Pearson_R_Dep_Methyl)
Pearson_R_null_Dep_Methyl = read.csv("Pearson_r_Dep_Methyl.NULL.csv", header = T, row.names = 1)
Pearson_R_null_Dep_Methyl = na.omit(Pearson_R_null_Dep_Methyl)

mean(Pearson_R_Dep_Methyl$r)

p_null = density(Pearson_R_null_Dep_Methyl$r, bw=0.15)
h=hist(Pearson_R_Dep_Methyl$r,20, freq =TRUE, xlim=c(-1.5,1.5), yaxt="n", xaxt="n", ylim = c(0,70), col = "#009E73", xlab = "Pearson r", main ='Correlation Methylation Dependency')
clip(0.05, 0.07,0,70)
abline(v=0.06072267, col=c("black"), lty = c(2), lwd=c(3))
axis(2, pos=-1.5, lwd = 1)
axis(1, at=c(-1, 0, 1), lwd = 1)
w = (h$counts/ h$density)[1]
lines(p_null$x, p_null$y*w, col = "red", lty=2, lw=3)
legend('topright', 'null distribution\n(random permutation)', col = "red", lw=3,lty=2,bty='n', cex = 0.8)

#Wilcox test for Methyl vs DepScore
hist(Pearson_R_Dep_Methyl$r)
shapiro.test(Pearson_R_Dep_Methyl$r)
hist(Pearson_R_null_Dep_Methyl$r)
shapiro.test(Pearson_R_null_Dep_Methyl$r)
wilcox.test(Pearson_R_Dep_Methyl$r, Pearson_R_null_Dep_Methyl$r, alternative = c("greater"))
#W = 383356, p-value = 0.04026
#alternative hypothesis: true location shift is greater than 0
