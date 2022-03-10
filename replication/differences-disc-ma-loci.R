library(tidyverse)
library(GenomicRanges)

disc<-as_tibble(read.csv("/Users/mike/Documents/Research/PUFA-GWAS/discovery/discovery-loci-ranges.csv"))
ma<-as_tibble(read.csv("/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/MA-loci-ranges.csv"))
disc
ma
r=500000

disc2<-disc%>%group_by(Locus.ID)%>%
    summarise(CHR=CHR, start=min(POS)-r,end=max(POS)+r)%>%distinct()


ma2<-ma%>%group_by(Locus.ID)%>%
    summarise(CHR=CHR,start=min(POS)-r, end=max(POS)+r)%>%distinct()
ma2

disc2

magr<-makeGRangesFromDataFrame(ma2,
                                 keep.extra.columns=T,
                                 ignore.strand=T)

discgr<-makeGRangesFromDataFrame(disc2,
                               keep.extra.columns=T,
                               ignore.strand=T)

discgr

#subsetByOverlaps subsets the first GRanges object 
#to include only those that overlap the second.
overlap<-as_tibble(as.data.frame(
    subsetByOverlaps(magr, discgr, ignore.strand=TRUE),
    stringsAsFactors=F))
nrow(overlap) #113

as.data.frame(overlap)
overlap$seqnames<-as.integer(overlap$seqnames)
overlap$strand<-as.character(overlap$strand)


#Find ones in MA not in discovery---------------------
joinMA<-full_join(ma2, overlap, 
                by=c("CHR"="seqnames", 
                     intersect(colnames(ma2), colnames(overlap))))

joinMA$indiscovery<-1
#The width column comes from the overlap table, so if this value is missing then they didn't overlap.
joinMA$indiscovery[is.na(joinMA$width)]<-0
# 0 = not found in this table by this method.

joinMA<-joinMA%>%select(-width, -strand)
joinMA
write.csv(joinMA, "/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/MAnotinDisc.csv",
          row.names=F, quote=F)


#Find ones in discovery not in MA---------------------

overlap<-as_tibble(as.data.frame(
    subsetByOverlaps(discgr, magr,ignore.strand=TRUE),
    stringsAsFactors=F))
nrow(overlap) #115

as.data.frame(overlap)
overlap$seqnames<-as.integer(overlap$seqnames)
overlap$strand<-as.character(overlap$strand)


joinDISC<-full_join(disc2, overlap, 
                  by=c("CHR"="seqnames", 
                       intersect(colnames(disc2), colnames(overlap))))

joinDISC$inMA<-1
#The width column comes from the overlap table, so if this value is missing then they didn't overlap.
joinDISC$inMA[is.na(joinDISC$width)]<-0
# 0 = not found in this table by this method.

joinDISC<-joinDISC%>%select(-width, -strand)

write.csv(joinDISC, "/Users/mike/Documents/Research/PUFA-GWAS/discovery/DISCnotinMA.csv",
          row.names=F, quote=F)
