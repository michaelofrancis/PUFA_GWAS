# library(biomaRt)
#suppressMessages(library(qqman))
#BiocManager::install("karyoploteR")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
suppressMessages(library(tidyverse))
#ref: https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotManhattan/PlotManhattan.html

pheno<-c("FAw3", "DHA", "FAw6","LA")

#MAKE LD SNP ANNOTATION TABLE FROM FUMA OUTPUT

ld<-list()
for (i in 1:length(pheno)){
    ld[[i]]<-as_tibble(read.table(paste(
    "/Users/mike/Documents/Research/PUFA-GWAS/regionalplote/combineforFUMA/SNPoutput/",
    pheno[i], "/snps.txt", sep=""), header=T))
    ld[[i]]$phenotype<-pheno[i]
    ld[[i]]$SNPtype<-"LD"
    ld[[i]]<-ld[[i]]%>%select(chr,phenotype, rsID, SNPtype, r2)
}

IndSig<-as_tibble(read.csv(
    "/Users/mike/Documents/R_files/PUFA_GWAS/regionalKaryoploteLDSNPs.csv"))%>%filter(SNPtype=="IndSig")

colnames(IndSig)<-colnames(ld)

#Combine into one table
ld<-do.call(rbind, ld)

annot<-rbind(IndSig, ld)



# Set colors ---------------------------------------------

legendcol<-c("#df6022","#e78e37","#edb75c","#f4dc8b","#ffffbf",
             "#bcecb5", "#79d6b7", "#31b9c0", "#029dc2")
#0.9, 0.8, 0.7....0.2, 0.1
#reverse so you can call legendcol[9] for 0.9
legendcol<-rev(legendcol)

#pointscolor=c("turquoise2","olivedrab3", "tan2", "gold2")
#pointscolor=c("#6c8384","#838a75", "#aa9b8d", "#8a8466") #https://wtools.io/change-color-saturation 15,15,15,20
#sat: 10,8,x,15

pointscolor=c("#727d7e","#828679", "#aa9b8d", "#84806c") 



# Load data -----------------------------------------------------

tab<-list()
gr<-list()
for (i in 1:length(pheno)){
    tab[[i]]<-as_tibble(read.table(
        paste("/Users/mike/Documents/Research/PUFA-GWAS/regionalplote/combineforFUMA/UKBEURKETMET-", 
              pheno[i], ".chrs4_7", 
              sep=""), header=T))
    tab[[i]]<-tab[[i]]%>%select("variant_id","chromosome",
                                "base_pair_location", "p_value")
    colnames(tab[[i]])<-c("SNP", "seqnames", "start", "pval")
    tab[[i]]$seqnames<-paste("chr", tab[[i]]$seqnames, sep="")
    tab[[i]]$end<-tab[[i]]$start
    gr[[i]]<-makeGRangesFromDataFrame(tab[[i]], keep.extra.columns = T,ignore.strand = T)
    
}


#  --------------------------------------------------------------
# PLOT 1 chr4 ---------------------------------------------------
# ---------------------------------------------------------------

# Make plots chr4 ----------------------------------------------------

#Set params=-=-=-=-=-=-=-=-=-=-
autotrack.margin <- 0.25
n<-5
ch<-4
sigcutoff=-log10(2.44e-08)
yminval=0
ymaxval=c(15,12,12,12)
pointsize=0.3
labels=c("Omega-3", "DHA", "Omega-6", "LA")
ycorrectlabel=c(1,1,1,1)
pp <- getDefaultPlotParams(plot.type=4)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

tiff("kpchr4.tiff", 
     width = 6, height = 7, 
     units = 'in', res = 600)


kp <- plotKaryotype(genome="hg19",plot.type=4, 
                    zoom="chr4:68500000-70600000", plot.params = pp,
                    cex=1)
kpAddBaseNumbers(kp, add.units = TRUE, cex=0.7, tick.dist = 5e5, 
                 clipping=TRUE, minor.tick.dist = 100000)

kpAddLabels(kp, labels = expression("-log"[10]*"(p)"),
            srt=90, pos=1, cex=1,
            label.margin = 0.1,
            r0=autotrack(1,1, margin=0),
            side="left")


#PLOT DOTS
for (i in 1:length(pheno)){
    x=(length(pheno)+2)-i


kpAddLabels(kp, labels = labels[i],  pos=1, cex=0.8, 
                label.margin = 0.3,
                r0=autotrack(x,n, margin=autotrack.margin),
                side=2)
kp <- kpPlotManhattan(kp, data=gr[[i]], points.col = pointscolor[i],
                      logp=TRUE, 
                      r0=autotrack(x,n, margin=autotrack.margin), 
                      points.cex = 0.5, 
                      ymax=ymaxval[i],
                      genomewide.col="red", 
                      genomewideline = sigcutoff,
                    )

kpAxis(kp, ymin=0, ymax=ymaxval[i], 
       r0=autotrack(x,n, margin=autotrack.margin), 
       cex=0.6, numticks=4)


LDsnpslist<-list()
LDsnps<-list()
for (l in 1:9){
    #Get LD SNPs
    LDsnpslist[[l]]<-unlist(annot%>%filter(chr==ch, phenotype==pheno[i],
                            SNPtype=="LD", r2>=(l/10), r2<((l+1)/10))%>%select(rsID))
    LDsnps[[l]]<-gr[[i]][gr[[i]]$SNP %in% LDsnpslist[[l]]]
    LDsnps[[l]]$y<- -1 *log10(LDsnps[[l]]$pval)

    #Color LD SNPs
    if (!identical(LDsnpslist[[l]], character(0))){
        kpPoints(kp, data = LDsnps[[l]], pch=16, cex=0.5, col=legendcol[l],
                 ymax=ymaxval[i],
                 r0=autotrack(x,n, margin=autotrack.margin)
        )
    }

}

#Get Indsigsnps
IndSigsnpslist<-unlist(annot%>%filter(chr==ch, phenotype==pheno[i], SNPtype=="IndSig")%>%select(rsID))
IndSigsnpslist
IndSigsnps<-gr[[i]][gr[[i]]$SNP %in% IndSigsnpslist]
IndSigsnps$y<- -1 *log10(IndSigsnps$pval)

#Color IndSigsnps

if (!identical(IndSigsnpslist, character(0))){
    kpPoints(kp, data = IndSigsnps, pch=18, cex=0.85, col="#990000",  
             ymax=ymaxval[i], 
             r0=autotrack(x,n, margin=autotrack.margin) 
    )
}

#Get top snps
grchr<-gr[[i]][seqnames(gr[[i]]) == paste("chr", ch, sep="")]
top1 <-which.min(grchr$pval)
top.snps<-grchr[top1]
top.snps
top.snps$y<- -1 *log10(top.snps$pval)
top.snps.label<-top.snps
top.snps.label$y<-top.snps.label$y+ycorrectlabel[i]

#Color top SNP
kpPoints(kp, data = top.snps, pch=17, cex=0.85, col="red3",  
         ymax=ymaxval[i], 
         r0=autotrack(x,n, margin=autotrack.margin) 
)
#Label Top SNP
# kpText(kp, data = top.snps.label, labels = top.snps.label$SNP, 
#        ymax=ymaxval[i], cex=0.5, 
#        col="red3",pos=4,
#        r0=autotrack(x,n, margin=autotrack.margin)
# )



}

# End pheno loop ------------------------------------------------

#Gene track
genes.data <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                    karyoplot = kp,
                                    plot.transcripts=TRUE,
                                    plot.transcripts.structure=TRUE
                                    )
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)


kpPlotGenes(kp, data=genes.data,
            gene.margin=0.7,
            transcript.margin=0.7,
            add.transcript.names = FALSE, 
            r1=0.18, 
            gene.name.position = "right",
            gene.name.cex=0.55,
            plot.transcripts=TRUE,
            plot.transcripts.structure=TRUE,
            avoid.overlapping = TRUE,
            clipping=FALSE
)

dev.off()


stop()

#  --------------------------------------------------------------
# PLOT 2 chr7 ---------------------------------------------------
# ---------------------------------------------------------------

# Make plots ----------------------------------------------------

#Set params=-=-=-=-=-=-=-=-=-=-
autotrack.margin <- 0.25
n<-5
ch=7
sigcutoff=-log10(2.44e-08)
yminval=0
ymaxval=c(12,12,12,12)
pointsize=0.3
labels=c("Omega-3", "DHA", "Omega-6", "LA")
ycorrectlabel=c(1.5,1.5,1,1.5)
pp <- getDefaultPlotParams(plot.type=4)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

tiff("kpchr7.tiff", 
     width = 6, height = 7, 
     units = 'in', res = 600)

kp <- plotKaryotype(genome="hg19",plot.type=4, 
                    zoom="chr7:900000-1300000", plot.params = pp,
                    cex=1)
kpAddBaseNumbers(kp, add.units = TRUE, cex=0.7, tick.dist = 2e5, 
                 clipping=TRUE, minor.tick.dist = 100000)

# kpAddLabels(kp, labels = expression("-log"[10]*"(p)"), 
#             srt=90, pos=1, cex=1, 
#             label.margin = 0.1,
#             r0=autotrack(1,1, margin=0), 
#             side="left")


#PLOT DOTS
for (i in 1:length(pheno)){
    x=(length(pheno)+2)-i
    
    kpAddLabels(kp, labels = labels[i],  pos=1, cex=0.8, 
                label.margin = 0.3,
                r0=autotrack(x,n, margin=autotrack.margin),
                side=2)
    kp <- kpPlotManhattan(kp, data=gr[[i]], points.col = pointscolor[i],
                          logp=TRUE, 
                          r0=autotrack(x,n, margin=autotrack.margin), 
                          points.cex = 0.5, 
                          ymax=ymaxval[i],
                          genomewide.col="red", 
                          genomewideline = sigcutoff,
    )
    
    kpAxis(kp, ymin=0, ymax=ymaxval[i], 
           r0=autotrack(x,n, margin=autotrack.margin), 
           cex=0.6, numticks=4)
    
    #Get top snps
    grchr<-gr[[i]][seqnames(gr[[i]]) == paste("chr", ch, sep="")]
    top1 <-which.min(grchr$pval)
    top.snps<-grchr[top1]
    top.snps
    top.snps$y<- -1 *log10(top.snps$pval)
    top.snps.label<-top.snps
    top.snps.label$y<-top.snps.label$y+ycorrectlabel[i]
    
    #Get Indsigsnps
    IndSigsnpslist<-unlist(annot%>%filter(chr==ch, phenotype==pheno[i], SNPtype=="IndSig")%>%select(rsID))
    IndSigsnpslist
    IndSigsnps<-gr[[i]][gr[[i]]$SNP %in% IndSigsnpslist]
    IndSigsnps$y<- -1 *log10(IndSigsnps$pval)
    
    
    LDsnpslist<-list()
    LDsnps<-list()
    for (l in 1:9){
        #Get LD SNPs
        LDsnpslist[[l]]<-unlist(annot%>%filter(chr==ch, phenotype==pheno[i], 
                                               SNPtype=="LD", r2>=(l/10), r2<((l+1)/10))%>%select(rsID))
        LDsnps[[l]]<-gr[[i]][gr[[i]]$SNP %in% LDsnpslist[[l]]]
        LDsnps[[l]]$y<- -1 *log10(LDsnps[[l]]$pval)
        
        #Color LD SNPs
        
        if (!identical(LDsnpslist[[l]], character(0))){
            kpPoints(kp, data = LDsnps[[l]], pch=16, cex=0.5, col=legendcol[l], 
                     ymax=ymaxval[i], 
                     r0=autotrack(x,n, margin=autotrack.margin) 
            )
        }
        
    }
    
    #Color IndSigsnps
    
    if (!identical(IndSigsnpslist, character(0))){
        kpPoints(kp, data = IndSigsnps, pch=18, cex=0.85, col="#990000",  
                 ymax=ymaxval[i], 
                 r0=autotrack(x,n, margin=autotrack.margin) 
        )
    }
    
    #Color top SNP
    kpPoints(kp, data = top.snps, pch=17, cex=0.85, col="red3",  
             ymax=ymaxval[i], 
             r0=autotrack(x,n, margin=autotrack.margin) 
    )
    #Label Top SNP
    # kpText(kp, data = top.snps.label, labels = top.snps.label$SNP, 
    #        ymax=ymaxval[i], cex=0.5, 
    #        col="red3",pos=4,
    #        r0=autotrack(x,n, margin=autotrack.margin)
    # )
    
    
    
}

# End pheno loop ------------------------------------------------

#Gene track
genes.data <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                    karyoplot = kp,
                                    plot.transcripts=TRUE,
                                    plot.transcripts.structure=TRUE
                                    )
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)


kpPlotGenes(kp, data=genes.data,
            gene.margin=0.7,
            transcript.margin=0.7,
            add.transcript.names = FALSE, 
            r1=0.18, 
            gene.name.position = "right",
            gene.name.cex=0.55,
            plot.transcripts=TRUE,
            plot.transcripts.structure=TRUE,
            avoid.overlapping = TRUE,
            clipping=FALSE
)

dev.off()

