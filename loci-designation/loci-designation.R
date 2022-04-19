library(readxl)
library(GenomicRanges)
library(tidyverse)

dir="/work/kylab/mike/PUFA-GWAS/loci-designation"

#Sheets that need loci merging within the sheet (all 250Kb)
#Discovery significant loci
#SMultiXcan
#Meta-analysis significant loci

#Merging / checking across sheets
#SMultiXcan to discovery
#SMultiXcan to Meta-analysis
#Discovery to meta analysis
#Meta analysis to discovery


FUMAd<-read_excel(paste(dir, "/supplementary_tables_v2.xlsx", sep=""),
                     sheet = 8, skip = 1)

FUMAm<-read_excel(paste(dir, "/supplementary_tables_v2.xlsx", sep=""),
                     sheet = 17, skip = 1)[-1]

SMul<-read_excel(paste(dir, "/supplementary_tables_v2.xlsx", sep=""),
                     sheet = 13, skip = 0)[-c(14:16)]
SMulgenepos<-as_tibble(read.csv("/work/kylab/mike/PUFA-GWAS/SMultiXcan/Huifang2/genepos/geneswithposition.csv",
		header=T))
SMul<-left_join(SMul, SMulgenepos)

colnames(SMul)[14]<-"CHR"

#---------------------------------------
###LOCI FUNCTION------------------------
#---------------------------------------
loci<-function(tab, window=250000){
#This function counts loci within one sheet and sets the same LocusID for 
#locuses within "window" bp of each other.
#Required columns in "tab" are CHR, start, and end.
#Use like this FUMA<-loci(FUMA)

tab<-tab[with(tab, order(CHR, start)),]

tab$LocusID<-1
count<-1

for (i in 2:(nrow(tab))){
    dif<- min(
	abs(tab$start[i] - tab$end[i-1]),
	abs(tab$start[i] - tab$start[i-1])
    	)
#	if (dif<500000 & dif >250000){print(paste(i, dif))}
    if (tab$CHR[i] == tab$CHR[i-1] & dif >window){ #chromosomes are same, difference more than window
        count=count+1
        tab$LocusID[i]<-count
    }
    else if (tab$CHR[i] == tab$CHR[i-1] & dif <=window){ #chromosomes are same, difference is less than window
        tab$LocusID[i]<-count
    }
    else {
        count=count+1 #chromosomes are different
        tab$LocusID[i]<-count
    }
}
tab<-tab%>%select(LocusID, CHR, start, end, everything())
return(tab)
}


FUMAd<-loci(FUMAd)
#write.csv(FUMAd, "/scratch/mf91122/PUFA-GWAS/loci-designation/FUMAd.csv", row.names=F, quote=F)

FUMAm<-loci(FUMAm)

SMul<-loci(SMul)


#---------------------------------------
###LOCI FUNCTION (BTWN 2 SHEETS)--------
#---------------------------------------

#Loci between sheets
#Make two GRanges with adding window to start and end of each....
#or half of window actually.

overlapGR<-function(tab1, tab2, w=(250000/2)){
inputparams<-as.list(sys.call())
#tab1<-FUMAd
#tab2<-FUMAm
tab1<-tab1%>%group_by(LocusID)%>%
    summarise(CHR=CHR, start=min(start)-w,end=max(end)+w)%>%distinct()

tab2<-tab2%>%group_by(LocusID)%>%
    summarise(CHR=CHR, start=min(start)-w,end=max(end)+w)%>%distinct()

tab1gr<-makeGRangesFromDataFrame(tab1,
                                 keep.extra.columns=T,
                                 ignore.strand=T)

tab2gr<-makeGRangesFromDataFrame(tab2,
                               keep.extra.columns=T,
                               ignore.strand=T)

overlap1<-as_tibble(as.data.frame(
    subsetByOverlaps(tab1gr, tab2gr, ignore.strand=TRUE),
    stringsAsFactors=F))

overlap1$seqnames<-as.numeric(levels(overlap1$seqnames))[overlap1$seqnames]
overlap1$strand<-as.character(overlap1$strand)

#tab1 vs tab2
join1<-full_join(tab1, overlap1, 
          by=c("CHR"="seqnames", 
          intersect(colnames(tab1), colnames(overlap1))))

join1$intab2<-1
#The width column comes from the overlap table, so if this value is missing then they didn't overlap.
join1$intab2[is.na(join1$width)]<-0

colnames(join1)[colnames(join1)=="intab2"]<-paste(
	inputparams[[2]], "_in_", inputparams[[3]],sep="")
join1<-join1%>%select(-width, -strand)

}

outdir="/scratch/mf91122/PUFA-GWAS/loci-designation"
write.csv(overlapGR(FUMAd, FUMAm), paste(outdir, "/FUMAd.csv", sep=""),quote=F, row.names=F)
write.csv(overlapGR(FUMAm, FUMAd), paste(outdir, "/FUMAm.csv", sep=""),quote=F, row.names=F)
write.csv(overlapGR(SMul, FUMAd),paste(outdir, "/SMul1.csv", sep=""),quote=F, row.names=F)
write.csv(overlapGR(SMul, FUMAm), paste(outdir, "/SMul2.csv", sep=""),quote=F, row.names=F)
