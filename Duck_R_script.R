####################################### ~ Library loading ~ ####
library("DESeq2")
library(dplyr)
library("RColorBrewer")
library("gplots")
library( "genefilter" )
library("pheatmap")
library(tools)
library("ggplot2")
library(magrittr)
library("biomaRt")
library(apeglm)
library("genefilter")
library(tidyr)
library(goseq)
####################################### ~ DECLARE YOUR VARIABLES HERE ~ ####
myspecies <- "duck_sp"
mygenes <- c("DUSP10", "ANP32A", "TLR7")
directory <- "D:/05.OneDrive/OneDrive - The University of Queensland/Marcus_htseq/chicken_htseq/"
setwd("D:/05.OneDrive/OneDrive - The University of Queensland/Marcus_htseq/")
resultsdir <- "D:/05.OneDrive/OneDrive - The University of Queensland/Marcus_htseq/results/"
dir.create(resultsdir)

resultsdir
testgroup <- "duck"
dir.create(paste0(resultsdir, testgroup))
resultsdir <- paste0(resultsdir, testgroup, "/")
resultsdir
sample_table <- read.table(paste0("marcus_meta.txt"), sep = "\t", header=F)
length(sample_table)
dim(sample_table)

sample_table$V3 <- sub('txt', 'tsv', sample_table$V3) 
sample_table
sample_table <- sample_table[which(sample_table$V1 == 'du'),]
rownames(sample_table) <- sample_table$V3
rownames(sample_table)
setwd('D:/05.OneDrive/OneDrive - The University of Queensland/Marcus_htseq/duck_htseq/')
files <- sample_table
files
sampleFiles <- paste0(files$V3)
sampleFiles
all(file.exists(sampleFiles))

sampleTable <- data.frame(sampleName = files$V3, 
                          fileName = sampleFiles,
                          condition = as.factor(files$V2))
head(sampleTable)


filedir <- getwd()

sampleTable
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = filedir,
                                  design= ~ condition)

dimnames(dds)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]

ddsKEEP

ddsKEEP <- DESeq(ddsKEEP)
resultsNames(ddsKEEP)
vst <- vst(ddsKEEP)

nonnormalized_counts <- counts(dds)
write.table(nonnormalized_counts, file=paste(resultsdir, testgroup, "_non_normalized_counts.txt", sep = ""), sep="\t", quote=F, col.names=NA)
normalized_counts <- counts(ddsKEEP, normalized=TRUE)
write.table(normalized_counts, file=paste(resultsdir, testgroup, "_normalized_counts.txt", sep = ""), sep="\t", quote=F, col.names=NA)

####################################### ~ Duck vs Chicken ~ ################################################################################################# 
resname <- 'condition_infected_vs_control'
res <- results(ddsKEEP, name = resname)
table(is.na(res$padj))
sum(res$padj<0.01, na.rm = T)
res[which(res$padj<0.01),]
sum(res$log2FoldChange>2, na.rm = T)
res[which(res$log2FoldChange>2),]
resSig <- res[which(res$padj < 0.01),]
resSig
write.csv( as.data.frame(resSig), file=paste(resultsdir, testgroup, "_", myspecies,"_", resname, "_DEGS.csv", sep = "") )

####################################### ~ topgene start ####
sum(res$padj < 0.01, na.rm=TRUE)
summary(res)
sum(!is.na(res$pvalue))
resSig <- subset(res, padj < 0.01)
head(resSig[ order( resSig$log2FoldChange ), ])

sum(res$padj < 0.01, na.rm=TRUE)

assayed.genes <- rownames(ddsKEEP)
ensembl<- useMart("ensembl")
ensembl = useMart("ensembl", dataset = "applatyrhynchos_gene_ensembl")
sp <- useDataset("applatyrhynchos_gene_ensembl",mart = ensembl)
sp
#listDatasets(sp)
assayed.genes
genenames<- getBM(mart = sp, filters = "ensembl_gene_id", values = assayed.genes, 
                  attributes = c("ensembl_gene_id","external_gene_name"))
genenames$ensembl_gene_id
genenames$external_gene_name
genenames$external_gene_name <- ifelse(genenames$external_gene_name == "", genenames$ensembl_gene_id,
                                       genenames$external_gene_name)
genenames$external_gene_name
resultsNames(ddsKEEP)

res1 <- res
rownames(res) <- genenames$external_gene_name
topGene <- rownames(res)[which.min(res$padj)]
topGene
ddsKEEP1 <- ddsKEEP
rownames(ddsKEEP) <- genenames$external_gene_name
ddsKEEP <- ddsKEEP[!duplicated(rownames(ddsKEEP)),]

listAttributes(sp)


plotCounts(ddsKEEP, gene=topGene, intgroup=c("condition"))
mygene <- topGene
for (i in c('DUSP1', 'DUSP10', 'OASL', 'DUSP1')){
  plotCounts(ddsKEEP, gene=i, intgroup=c("condition"))
}
plotCounts(ddsKEEP, gene="OASL", intgroup=c("condition"))
d <- plotCounts(ddsKEEP,gene = mygene, intgroup = "condition", main = paste( myspecies, "gene -", mygene), returnData=TRUE)
p <- ggplot(d, aes(x=condition, y=count)) + 
  #geom_violin()  +
  geom_boxplot(colour = "red", fill = "orange", alpha = 0.2, 
               outlier.colour="black", outlier.shape=8, outlier.size=2, notch=F, notchwidth = 0.5) + 
  #geom_dotplot(binwidth = 50, binaxis='y', stackdir='center', dotsize=1)
  geom_point(position=position_jitter(w=0.1,h=0), colour = 'purple', size = 1) + 
  scale_y_log10(breaks=c(25,100,400)) + 
  theme(
    #panel background elements
    panel.background = element_rect(
      fill = "grey90",
      colour = "black",
      size = 1,
    ),
    legend.position= "bottom",
    plot.title = element_text(color="black", size=12, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=11, face="bold"),
    axis.title.y = element_text(color="#993333", size=11, face="bold")
  ) + 
  ggtitle(paste0("Condition: ", myspecies, " gene - ", mygene)) + xlab(testgroup) + 
  ylab("Noramlized gene count") +
  labs(fill = d$condition) +
  stat_summary(fun=mean, geom="point", shape=23, size=4) + scale_color_grey() +
  scale_fill_manual(values=c("#999999", "#E69F00"))
print(p)
ggsave(p, file=paste0(resultsdir, testgroup,"_", resname, mygene,".png", sep = ""), width = 14, height = 10, units = "cm")

####################################### ~ top-gene end ####
resSig= res[which(res$padj<0.01),]
head(resSig[which(resSig$log2FoldChange >= 2),])
tail(resSig[order(resSig$log2FoldChange),])
write.csv( as.data.frame(resSig), file=paste(resultsdir, testgroup, "_", myspecies, resname, "_DEGS.csv", sep = "") )
myspecies
(resSig_up= resSig[which(resSig$log2FoldChange >= 2),])
resSig_up=resSig_up[order(resSig_up$log2FoldChange),]
head(resSig_up)
write.csv( as.data.frame(resSig_up), file=paste(resultsdir, testgroup, "_", myspecies, "_", resname, "_UP_DEGS.csv", sep = "" ))
#file_up1 <- paste(resultsdir, testgroup, "_", myspecies, "Du_vs_Ch_UP_DEGS.csv", sep = "" )
#file_up1


resSig_down= resSig[which(resSig$log2FoldChange <= -2),]
resSig_down=resSig_down[order(resSig_down$log2FoldChange),]
head(resSig_down)
write.csv( as.data.frame(resSig_down), file=paste(resultsdir, testgroup, "_", myspecies, "_", resname, "_DOWN_DEGS.csv", sep = ""))
#file_down1 <- paste(resultsdir, testgroup, "_", myspecies, "Du_vs_Ch_DOWN_DEGS.csv", sep = "")
#file_down1

resD <- res
resDSort <- resD[order(resD$log2FoldChange, na.last = NA, decreasing = T),]
resDSort
topDESeq2 <- resDSort[1:100,]
topDESeq2
write.csv(topDESeq2, file=paste(resultsdir, testgroup, myspecies, "_", resname, "_topDESeq2_100_.csv", sep = ""))
topDESeq2

rld <- rlog(ddsKEEP)

topgenes <- head(rownames(resDSort),40)
mat <- assay(rld)[topgenes,]
(mat <- mat -rowMeans(mat))

col.pan <- colorpanel(100, "blue","white","red")
#Non scaled heatmap for topgenes

heatmap.2(mat, col=col.pan, Rowv=TRUE, scale="none", trace="none")

#Scaled heatmap for topgenes
scaled.mat<-t(scale(t(mat)))


pheatmap(scaled.mat, fontsize_row = 6, labels_col = d$condition,
         #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         color = col.pan, 
         scale = "row",
         main = "Heatmap for the top 40 differentially expressed genes",treeheight_col = 25, cluster_cols = T, 
         cluster_rows = F, clustering_distance_cols = "euclidean" )
heatmap.2(scaled.mat, col=col.pan, Rowv=TRUE, scale="row",
          trace="none" ,margins = c(10,8),cexRow=0.5, cexCol=1, keysize=1,labCol = d$condition, 
          main =  "Heatmap for the top 40 DEGs", dendrogram = "both", cex.main = 10)


pdf(paste0(resultsdir, testgroup, myspecies, "_", resname, "_VSD_scaled_topgenes_heatmap.pdf"))
#heatmap.2(scaled.mat, col=col.pan, Rowv=TRUE, scale="none",
#          trace="none", labRow= "",margins = c(10,8),cexRow=0.5, cexCol=1, keysize=1,labCol = sampleCondition)
pheatmap(scaled.mat, fontsize_row = 6, labels_col = d$condition,
         #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         color = col.pan, 
         scale = "row",
         main = "Heatmap for the top 40 differentially expressed genes",treeheight_col = 25, cluster_cols = T, 
         cluster_rows = F, clustering_distance_cols = "euclidean" )
dev.off()
graphics.off()

plotMA(res,ylim=c(-5,5))

plotDispEsts( ddsKEEP, ylim = c(1e-6, 1e1) )

hist( res$pvalue, breaks=20, col="green" )
####################################### ~ Volcano plots ~ ####
library(EnhancedVolcano)

voldu1 <- EnhancedVolcano(res,
                          lab = NA, #rownames(res),
                          x = 'log2FoldChange',
                          y = 'padj', 
                          selectLab = NULL, #rownames(res)[which(names(keyvals) %in% c('high', 'low'))],
                          #selectLab = c('OASL','CXCL10','ISG15', 'USP41','IFITM1','MX1'),
                          xlab = 'log2Fold Change', ylab = '-log10P',
                          title = NULL, #'Volcano plot - Difference in Kids compared to Adults ',
                          subtitle = NULL, #'~ without treatment effect(i.e. without Virus infection)',
                          pCutoff = 0.01, axisLabSize = 10,
                          FCcutoff = 2,
                          pointSize = 1.5,
                          labSize = 3.0,
                          labCol = 'black',
                          labFace = 'bold',
                          boxedLabels = F,
                          ylim = c(0,3.5),
                          xlim = c(-6.5,6.5),
                          colAlpha = 3/4,
                          gridlines.major = F,
                          gridlines.minor = F,
                          col = c("grey30", "palegreen3", "lightslateblue", "orangered2"),
                          legendPosition = 'bottom', #c(0.15, 0.15),
                          legendLabels = c("NS", 'log2fold change (±2)', 'p-adj < 0.01', 'log2fold change & p-adj'),
                          legendLabSize = 6,
                          legendIconSize = 1.5,
                          drawConnectors = F,
                          widthConnectors = 1.0, caption = NULL,
                          colConnectors = 'black') 

voldu1
voldu2 = voldu1 + ggplot2::theme(legend.position = "bottom", legend.text = element_text(size = 10), 
                                 axis.text.y = element_blank(), 
                                 axis.text.x = element_blank(), #axis.title = element_blank(),
                                 text = element_text(family = 'serif', face = 'bold', size = 10),
                                 axis.line = element_line(size = 0.1), 
                                 axis.ticks.x = element_blank()) + #axis.title.y = element_blank()) + 
  scale_y_continuous(breaks=NULL) 
voldu2
#graphics.off()
#setEPS()
#postscript(paste0(resultsdir, testgroup, "_", resname, "Volcano_plot.eps"), fonts = "serif")
#print(volch2)
#dev.off()
pdf(paste0(resultsdir, testgroup, "_", resname, "Volcano_plot.pdf"), width = 4, 
    height = 3.5, colormodel = "srgb",family = "serif")
print(voldu2)
dev.off()

####################################### ~ other analytic plots ~ #####
###Shrinking results to fit in the graphs

pdf(paste0(resultsdir, testgroup, "_", resname, "_rlog_PCA_plot.pdf"))
plotPCA(rld, intgroup = "condition")
dev.off()
plotPCA(rld, intgroup = "condition")
graphics.off()

####################################### ~ Data quality assessment by sample clustering ~ ####
####################################### ~ Another method of PCA plot ~ ####
pdf(paste0(resultsdir, testgroup,"_", resname, "_VSD_PCA_plot.pdf"))
vsdata <- vst(ddsKEEP, blind=FALSE)
plotPCA(vsdata, intgroup="condition")
dev.off()

plotPCA(vsdata, intgroup="condition")
head( order( rowVars( assay(rld) ), decreasing=TRUE ), 25)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing = TRUE ), 25)

heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(100),margin=c(10,8), cexRow=0.5, cexCol=1, keysize=1.5, labCol = d$sp)

pdf(paste0(resultsdir, testgroup, resname, "topVargenes_heatmap.pdf"))
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(100),margin=c(10,8), cexRow=0.5, cexCol=1, keysize=1.5, labCol = d$sp)
dev.off()



plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 )

pdf(paste0(resultsdir, testgroup,"_", resname, "_assayplot.pdf"))
plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 )
dev.off()

#Heatmap of the count matrix

select <- order(rowMeans(counts(ddsKEEP,normalized=TRUE)),
                decreasing=TRUE)[1:20]



select


####################################### ~ Heatmap of the sample-to-sample distances ~ ####
library("PoiClaClu")

poisd <- PoissonDistance(t(counts(ddsKEEP)))
(samplePoisDistMatrix <- as.matrix( poisd$dd ))

(rownames(samplePoisDistMatrix) <- paste(ddsKEEP$condition) )

colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

hc <- hclust(poisd$dd)

heatmap.2( samplePoisDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colours,
           margins=c(2,10), labCol=FALSE )


pdf(paste0(resultsdir, testgroup,"_", resname, "_sample_to_sample_distance.pdf"))

heatmap.2( samplePoisDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colours,
           margins=c(2,10), labCol=FALSE )


dev.off()


ramp <- 1:3/3

cols <- c(rgb(ramp, 0, 0),
          rgb(0, ramp, 0),
          rgb(0, 0, ramp),
          rgb(ramp, 0, ramp))
print( plotPCA( rld, intgroup = c( "condition") ) )

pdf(paste0(resultsdir, testgroup, "_", resname, "_sizefac_condition_PCA_plot.pdf"))
plotPCA( rld, intgroup = c( "condition") )
dev.off()

colData(ddsKEEP)


topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 20 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins = c(10,10),labCol = d$sp)

#sessionInfo()

####################################### ~ GO Enrichment analysis for ENSEMBL based annotations ~ #################################
res <- res1
assayed.genes <- rownames(res)
assayed.genes

de.genes <- rownames(res)[ which(res$padj < 0.01) ]

gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector) <- assayed.genes
assayed.genes[assayed.genes%in%de.genes]
(names(gene.vector))
####################################### ~ Biomart GO::terms matching*******************************########################
EG2KEGG<- getBM(mart = sp,
                values = gene.vector, attributes = c("ensembl_gene_id", "go_id"))
head(EG2KEGG, 10)

geneID2GO <- by(EG2KEGG$ensembl_gene_id,
                EG2KEGG$go_id,
                function(x) as.data.frame(x))
head(geneID2GO,15)
write.csv(as.data.frame(EG2KEGG), file = paste0(resultsdir, "/", testgroup, "_", resname, "_geneID.csv"))


GTF <- "../Anas_platyrhynchos_platyrhynchos.CAU_duck1.0.103.gtf"

library(GenomicFeatures)
txdb = makeTxDbFromGFF( GTF, format = "gtf")
txsBygene=transcriptsBy(txdb,"gene")
lengthdata=mean(width(txsBygene))
head(lengthdata)
names(gene.vector)
assayed.genes
where <- match(names(gene.vector), names(lengthdata))
where

lengthdata <- lengthdata[where]
lengthdata

library(goseq)
pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", testgroup, "_", resname, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", testgroup, "_", resname, "_enriched_go_", myspecies, ".csv"))
enriched.GO[[1]]
GOTERM=goResults$term
names(GOTERM) <- goResults$category
GOTERM[["GO:0051150"]]
length(GOTERM)
capture.output(for(go in enriched.GO[1:379]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
, file=paste0(resultsdir, "/", testgroup,"_SigGo.txt"))
goResults
goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes
#goResults[goResults$ontology == 'BP',]
#goBPRes <- goResults[goResults$ontology == 'BP',] 

####################################### ~ visualize the top 30 hits*****************************##########################
write.csv(as.data.frame(goBPRes), file = paste0(resultsdir,"/", testgroup, myspecies, "_", resname, "_go_bp.csv"))


goBP <- goBPRes %>% 
  top_n(20, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "black", fill = "coral1") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

print(goBP)

goBPdu <- goBP
goResults %>% 
  top_n(60, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

pdf(paste0(resultsdir, "/", testgroup,resname, "_", "_goseq_enrichment.pdf"), 
    width = 7.25, height = 8.25)

print(goBP)

dev.off()
####################################### ~ Pathview ~ ######################
library(pathview)
library(gage)
setwd('../duck/')
getwd()
ko <- read.table('../../../../../03.R/Ducks_DGE_Results/duck_up.txt', sep = ' ')
head(ko)

df <- ko %>% 
  group_by(V1) %>% 
  top_n(1, V3) %>%
  as.data.frame()

nrow(ko)
nrow(df)

ko <- df
colnames(ko) <- c('protein_stable_id', 'KO_numbers', 'score')
ko$protein_stable_id <- gsub('\\..*', '', ko$protein_stable_id)
head(ko)
ko  <- dplyr::select(ko, -c('score'))
head(ko)
tx2gene <- read.table("../../../../../03.R/Ducks_DGE_Results/tx2gene.txt", header = T)
head(tx2gene)
tx2gene <- dplyr::select(tx2gene, c("protein_stable_id", "gene_stable_id"))

head(tx2gene)
head(tx2gene)
nrow(tx2gene)
nrow(ko)
ko <- merge(tx2gene, ko, by = 'protein_stable_id')
head(ko)
res
rownames(ddsKEEP1)
res <- as.data.frame(results(ddsKEEP1))
res
res$gene_stable_id <- rownames(res)

head(ko)
head(res)
res <- merge(res, ko, by ='gene_stable_id')
head(res)

#res <- res[which(res$padj<0.01),]
genes <- res$log2FoldChange
names(genes) <- res$KO_numbers
#install.packages(c("boot", "cluster", "mgcv", "rtracklayer", "SparseM", "cachem"))
genes

kg.ko = kegg.gsets('ko')
kg.ko=kg.ko$kg.sets
#kg.gga = kg.gga$kg.sets

exp.fc = genes
fc.kegg.p <- gage(exp.fc, gsets = kg.ko, ref = NULL, samp = NULL)
fc.kegg.p

fc.kegg.p
sel <- fc.kegg.p$greater[, "q.val"] < 0.1 &
  + !is.na(fc.kegg.p$greater[, "q.val"])

table(sel)
path.ids <- rownames(fc.kegg.p$greater)[sel]
path.ids

path.ids2 <- substr(c(path.ids), 1, 8)
path.ids2[1:25]
path.ids2
length(path.ids2)
path.ids2

out.suffix='up_duck'
getwd()
pv.out.list <- sapply(path.ids2[1:length(path.ids2)], function(pid) pathview(
  gene.data = exp.fc, pathway.id = pid,
  species = "ko", out.suffix=out.suffix))

sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 &
  + !is.na(fc.kegg.p$less[,"q.val"])

table(sel.l)
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids.l), 1, 8)
path.ids2[1:25]
path.ids2
length(path.ids2)
path.ids2

out.suffix='down_duck'
getwd()
pv.out.list <- sapply(path.ids2[1:length(path.ids2)], function(pid) pathview(
  gene.data = exp.fc, pathway.id = pid,
  species = "ko", out.suffix=out.suffix))

#pv.out <- pathview(gene.data = genes, pathway.id = "05164", low = 'blue', mid = 'white', high = 'red',
#                   species = "ko", out.suffix = "Influenza_A_duck")

#pv.out$plot.data.gene
#pv.out <- pathview(gene.data = genes, pathway.id = "04010", low = 'blue', mid = 'white', high = 'red',
#                   species = "ko", out.suffix = "MAPK_PATHWAY_duck")

#pv.out <- pathview(gene.data = genes, pathway.id = "04622", low = 'blue', mid = 'white', high = 'red',
#                   species = "ko", out.suffix = "RIG-I_PATHWAY_duck")


#pv.out <- pathview(gene.data = genes, pathway.id = "04620", low = 'blue', mid = 'white', high = 'red',
#                   species = "ko", out.suffix = "Toll_like_receptor_PATHWAY_duck")


#pv.out <- pathview(gene.data = genes, pathway.id = "04621", low = 'blue', mid = 'white', high = 'red',
#                   species = "ko", out.suffix = "NOD-like_receptor_PATHWAY_duck")

####################################### ~ PCA plots ~ #####
MyreadCountMatrix <- read.table('../results/duck/duck_normalized_counts.txt', sep = '\t', header = T)
#read.table(paste0(resultsdir, '/', ageclass, "_hsapiens_normalized_counts.txt"), 
#                              sep = '\t', header = T)

head(MyreadCountMatrix)
MyreadCountMatrix

df <- (dplyr::select(MyreadCountMatrix, -X))

head(df)
#df
#dfmock <- df[,!c(FALSE, TRUE)]
#head(dfmock)

metadata <- data.frame(row.names = colnames(df))
head(metadata)
metadata$Group <- rep(NA, ncol(df))
metadata
metadata$Group[1:3] <- 'Control'
metadata$Group[4:6] <- 'Infected'
metadata$CRP <- sample.int(100, size=ncol(df), replace=TRUE)
metadata$ESR <- sample.int(100, size=ncol(df), replace=TRUE)

metadata
#metavirus <- metadata[metadata$Group == "Virus",]

#metadata <- metavirus
#metadata

#metadata$age <- "Adult" 
#metadata[6:10,]$age <- "Kid"
#metadata
library(PCAtools)
pca(df, metadata = metadata)
p <- pca(df, metadata = metadata)

p
p$rotated
df <- p$rotated
df
df$Group <- rep(NA, 6)
df
df[1:3,]$Group <- "Control" 
df[4:6,]$Group <- "Infected"
#df$Age
require(ggalt)
pca2df <- df
write.table(pca2df, file = paste0(resultsdir, '/pca2_data_frame.txt'), sep = '\t')
pcadu <- ggplot(df, aes(x = PC1, y = PC2, fill = Group))+
  geom_point(aes(colour = Group), size = 4) +
  theme(axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.title = element_text(colour = "black", size = 10, face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.text = element_text(colour = "black", size = 10, face = "bold")) +
  #legend.title = element_blank()) + #element_blank())+
  scale_y_continuous(breaks=NULL) +
  scale_x_continuous(breaks = NULL) +
  geom_encircle(alpha = 0.2, show.legend = FALSE) 


pcadu
library(ggpubr)
#?ggarrange
voldu2
volch2
figure1 <-ggarrange(voldu2, volch2, ncol = 2, align = "hv", legend = "bottom", 
                    font.label = list(size = 12, face = "bold"),
                    common.legend = T)
figure1

figure2 <-ggarrange(pcadu, pcach, ncol = 2, align = "hv", legend = "bottom", 
                    font.label = list(size = 12, face = "bold"),
                    common.legend = T)
figure2


figure3 <-ggarrange(figure2, figure1, nrow = 2, align = "hv", legend = "bottom", 
                    font.label = list(size = 12, face = "bold", font = "serif"),
                    common.legend = T)


figure <- figure3 + ggplot2::theme(axis.title = element_text(size = 12))

pdf(paste0(resultsdir, testgroup, "_", resname, "_figure.pdf"), width = 7.25, 
    height = 9.75, colormodel = "srgb",family = "serif")
print(figure)
dev.off()


setEPS()
postscript(paste0(resultsdir, testgroup, "_", resname, "Volcano_plot.eps"), width = 7.25, height = 9.75,
           fonts = "serif")
print(figure)
dev.off()

