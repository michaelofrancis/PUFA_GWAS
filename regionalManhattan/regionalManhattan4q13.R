#BiocManager::install("karyoploteR")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
suppressMessages(library(tidyverse))
#library(ggrepel)

#karyoploteR reference/tutorial:
#https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotManhattan/PlotManhattan.html
#https://twillis.xyz/misc/using-karyoploteR.html
#https://bioconductor.org/packages/devel/bioc/manuals/karyoploteR/man/karyoploteR.pdf


pheno<-c("FAw3", "FAw6", "DHA", "LA")

tab<-list()


# Load data -----------------------------------------------------
for (i in 1:length(pheno)){
tab[[i]]<-as_tibble(read.table(
    paste("/Users/mike/Documents/Research/PUFA-GWAS/regionalplote/UKBEURKETMET-", 
          pheno[i], ".tsv.chr4.sort", 
          sep=""), header=T))
tab[[i]]<-tab[[i]]%>%select("variant_id","chromosome",
                            "base_pair_location", "p_value")
colnames(tab[[i]])<-c("SNP", "seqnames", "start", "pval")
tab[[i]]$seqnames<-paste("chr", tab[[i]]$seqnames, sep="")
tab[[i]]$end<-tab[[i]]$start
}



gr<-list()
for (i in 1:length(pheno)){
gr[[i]]<-makeGRangesFromDataFrame(tab[[i]], keep.extra.columns = T,ignore.strand = T)
}

gr[[1]]


# Make plots ----------------------------------------------------
autotrack.margin <- 0.15
pp <- getDefaultPlotParams(plot.type=4)
kp <- plotKaryotype(genome="hg19",plot.type=4, 
                    zoom="chr4:68500000-71000000", plot.params = pp)
kpAddBaseNumbers(kp, add.units = TRUE, cex=0.7, tick.dist = 5e5, 
                 clipping=TRUE, minor.tick.dist = 100000)
#PLOT OMEGA-3 (5/5)
kpAddLabels(kp, labels = expression("-log"[10]*"(p)"), srt=90, pos=1, cex=0.8, 
    label.margin = 0.05,r0=autotrack(5,5, margin=autotrack.margin), side="left")
kpAddLabels(kp, labels = "Omega-3", srt=90, pos=1, cex=0.8, 
            label.margin = 0.02,r0=autotrack(5,5, margin=autotrack.margin), side=2)
kp <- kpPlotManhattan(kp, data=gr[[1]], points.col = "turquoise2", 
                      logp=TRUE, r0=autotrack(5,5, margin=autotrack.margin), 
                      points.cex = 0.5, ymax=15,
                      genomewide.col="red", 
                      genomewideline = -log10(2.44e-08),
)

kpAxis(kp, ymin=0, ymax=15, r0=autotrack(5,5, margin=autotrack.margin), 
       cex=0.6, numticks=4)

top1 <-which.min(gr[[1]]$pval)
top.snps<-gr[[1]][top1]
top.snps
top.snps$y<- -1 *log10(top.snps$pval)

kpPoints(kp, data = top.snps, pch=17, cex=0.6, col="red", lwd=2, 
         ymax=15, 
         r0=autotrack(5,5, margin=autotrack.margin) )
kpText(kp, data = top.snps, labels = top.snps$SNP, ymax=15, cex=0.6, 
       pos=4,col="red",
       r0=autotrack(5,5, margin=autotrack.margin) )
    

#PLOT DHA (4/5)
kpAddLabels(kp, labels = expression("-log"[10]*"(p)"), srt=90, pos=1, cex=0.8, 
            label.margin = 0.05,r0=autotrack(4,5, margin=autotrack.margin), 
            side="left")
kpAddLabels(kp, labels = "DHA", srt=90, pos=1, cex=0.8, 
            label.margin = 0.02,r0=autotrack(4,5, margin=autotrack.margin), 
            side=2)
kp <- kpPlotManhattan(kp, data=gr[[3]], points.col = "olivedrab3", 
                      logp=TRUE,
                      ymax=12, r0=autotrack(4,5, margin=autotrack.margin), 
                      points.cex = 0.5,
                      genomewide.col="red",
                      genomewideline = -log10(2.44e-08)
)

kpAxis(kp, ymin=0, ymax=12, r0=autotrack(4,5, margin=autotrack.margin), 
       cex=0.6, numticks=4)

top1 <-which.min(gr[[3]]$pval)
top.snps<-gr[[3]][top1]
top.snps
top.snps$y<- -1 *log10(top.snps$pval)

kpPoints(kp, data = top.snps, pch=17, cex=0.6, col="red", lwd=2, 
         ymax=12, 
         r0=autotrack(4,5, margin=autotrack.margin) )
kpText(kp, data = top.snps, labels = top.snps$SNP, ymax=12, cex=0.6, 
       pos=4,col="red",
       r0=autotrack(4,5, margin=autotrack.margin) )



#PLOT OMEGA-6 (3/5)
kpAddLabels(kp, labels = expression("-log"[10]*"(p)"), srt=90, pos=1, 
            cex=0.8, 
            label.margin = 0.05,r0=autotrack(3,5, margin=autotrack.margin), 
            side="left")
kpAddLabels(kp, labels = "Omega-6", srt=90, pos=1, cex=0.8, 
            label.margin = 0.02,r0=autotrack(3,5, margin=autotrack.margin), 
            side=2)
kp <- kpPlotManhattan(kp, data=gr[[2]], points.col = "tan2", 
                      logp=TRUE,
                      ymax=12, r0=autotrack(3,5, margin=autotrack.margin), 
                      points.cex = 0.5,
                      genomewide.col="red",
                      genomewideline = -log10(2.44e-08)
)

kpAxis(kp, ymin=0, ymax=12, r0=autotrack(3,5, margin=autotrack.margin), 
       cex=0.6, numticks=4)

top1 <-which.min(gr[[2]]$pval)
top.snps<-gr[[2]][top1]
top.snps
top.snps$y<- -1 *log10(top.snps$pval)

kpPoints(kp, data = top.snps, pch=17, cex=0.6, col="red", lwd=2, 
         ymax=12, 
         r0=autotrack(3,5, margin=autotrack.margin) )
kpText(kp, data = top.snps, labels = top.snps$SNP, ymax=12, cex=0.6, 
       pos=4,col="red",
       r0=autotrack(3,5, margin=autotrack.margin) )



#PLOT LA (2/5)
kpAddLabels(kp, labels = expression("-log"[10]*"(p)"), srt=90, pos=1, 
            cex=0.8, label.margin = 0.05,r0=autotrack(2,5, margin=autotrack.margin), 
            side="left")
kpAddLabels(kp, labels = "LA", srt=90, pos=1, cex=0.8, 
            label.margin = 0.02,r0=autotrack(2,5, margin=autotrack.margin), 
            side=2)
kp <- kpPlotManhattan(kp, data=gr[[4]], points.col = "gold2", 
                      logp=TRUE,
                      ymax=12, r0=autotrack(2,5,margin=autotrack.margin), 
                      points.cex = 0.5,
                      genomewide.col="red",
                      genomewideline = -log10(2.44e-08)
)

kpAxis(kp, ymin=0, ymax=12, r0=autotrack(2,5, margin=autotrack.margin), 
       cex=0.6, numticks=4)

top1 <-which.min(gr[[4]]$pval)
top.snps<-gr[[4]][top1]
top.snps
top.snps$y<- -1 *log10(top.snps$pval)
top.snps2<-top.snps
top.snps2$y<- -1 *log10(top.snps$pval)+0.7

kpPoints(kp, data = top.snps, pch=17, cex=0.6, col="red", lwd=2,
         ymax=12,
         r0=autotrack(2,5, margin=autotrack.margin) )
kpText(kp, data = top.snps2, labels = top.snps$SNP, ymax=12, cex=0.6,
       pos=4,col="red",
       r0=autotrack(2,5, margin=autotrack.margin) )

#Gene track
genes.data <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                    karyoplot = kp,
                                    plot.transcripts=TRUE,
                                    plot.transcripts.structure=TRUE 
                                    )
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)

kpPlotGenes(kp, data=genes.data, 
            transcript.margin=2,
            add.transcript.names = FALSE, 
            r1=0.18, 
            gene.name.position = "left",
            gene.name.cex=0.5,
            plot.transcripts=TRUE,
            plot.transcripts.structure=TRUE,
            clipping=FALSE
            )
