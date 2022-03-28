#Reference for displaying results
#https://doi.org/10.1038/s41593-020-0643-5

suppressMessages(library(plyr))
library(biomaRt)
suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyverse))

##04.SPrediXcan


phenotypes<-c("w3FA_NMR_resinv", "w3FA_NMR_TFAP_resinv", "w6FA_NMR_resinv",
              "w6FA_NMR_TFAP_resinv", "w6_w3_ratio_NMR_resinv",
              "DHA_NMR_resinv", "DHA_NMR_TFAP_resinv", "LA_NMR_resinv",
              "LA_NMR_TFAP_resinv", "PUFA_NMR_resinv", "PUFA_NMR_TFAP_resinv",
              "MUFA_NMR_resinv", "MUFA_NMR_TFAP_resinv", "PUFA_MUFA_ratio_NMR_resinv")

phenotypenames<-c("Omega-3", "Omega-3 (TFAP)", "Omega-6" , "Omega-6 (TFAP)",
                  "Omega-6:Omega-3 ratio", "DHA", "DHA (TFAP)", "LA", "LA (TFAP)",
                  "PUFAs", "PUFAs (TFAP)", "MUFAs", "MUFAs (TFAP)", "PUFAs:MUFAs ratio")



dir4="/scratch/mf91122/PUFA-GWAS/SMultiXcan/huifang2/04.SPrediXcan_tissues"
phenoresults<-as_tibble(read.csv(paste(dir4, "/significant.04.SPrediXcan.csv", sep="")))
dir5="/scratch/mf91122/PUFA-GWAS/SMultiXcan/huifang2/05.SMultiXcan"
part5out<-as_tibble(read.csv(paste(dir5, "/significant.05.SMultiXcan.csv", sep="")))

###*****************************************************************************************
###Part 5 novelty***************************************************************************
###*****************************************************************************************

r=500000

Yknown<-as_tibble(
        read.csv("/work/kylab/mike/PUFA-GWAS/novelty/Ytable/All_known_loci_02212022.csv",
        header=T))
Yknown2<-Yknown[c("rsID", "chr", "start", "end", "phenotype", "Year","PubmedID","Author")]
Yknown2$start<-Yknown2$start-r
Yknown2$end<-Yknown2$end+r
Yknown2<-Yknown2[complete.cases(Yknown2),]


#Get gene positions---------------------
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
        path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

ensembl = useDataset("hsapiens_gene_ensembl",mart=grch37)

t2g<-getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",
                'chromosome_name','start_position','end_position'), mart = ensembl)

#Side quest to get part4 gene names and positions

#file<-phenoresults
#file$ensembl_gene_id <- gsub("\\..*","", file$gene)

#my_ids.version <- merge(file, t2g, by= 'ensembl_gene_id')
#my_ids.version<-as_tibble(my_ids.version)

#dat<-my_ids.version[c("gene_name", "chromosome_name", "start_position", "end_position", "pvalue")]

##Clean up dataset
#dat<-dat[grep("HG*", dat$chromosome_name, invert=T),] #43 weird entries removed
#dat$chromosome_name<-as.numeric(dat$chromosome_name)
#dat<-dat%>%distinct()
#write.csv(dat, "/scratch/mf91122/PUFA-GWAS/SMultiXcan/huifang2/SPredictXcan_genenamepos.csv")

#End side quest-------------

#Manually curated gene positions (HG19)

geneposb<-read.csv("/work/kylab/mike/PUFA-GWAS/SMultiXcan/Huifang2/genepos/SPredictXcan_genenamepos2.csv")
colnames(geneposb)<-c("gene_name", "chromosome_name", "start_position", "end_position")


file<-part5out

file$ensembl_gene_id <- gsub("\\..*","", file$gene)

my_ids.version <- merge(file, t2g, by= 'ensembl_gene_id')
my_ids.version<-as_tibble(my_ids.version)

dat<-my_ids.version[c("gene_name", "chromosome_name", "start_position", "end_position")]

#Clean up dataset
dat<-dat[grep("HG*", dat$chromosome_name, invert=T),] #43 weird entries removed
dat$chromosome_name<-as.numeric(dat$chromosome_name)



#--------------------------------------------------

dat<-rbind(dat, geneposb)

colnames(dat)<-c("gene_name","chromosome_name","start", "end")
dat$start<-dat$start-r
dat$end<-dat$end+r

#Make GR
DATgr<-makeGRangesFromDataFrame(dat,
                         keep.extra.columns=T,
                         ignore.strand=T)

Yknowngr<-makeGRangesFromDataFrame(Yknown2,
                         keep.extra.columns=T,
                         ignore.strand=T)

#subsetByOverlaps subsets the first GRanges object 
#to include only those that overlap the second.
overlap<-as_tibble(as.data.frame(
        subsetByOverlaps(DATgr, Yknowngr, ignore.strand=TRUE),
                stringsAsFactors=F))
nrow(overlap) #1310

overlap$seqnames<-as.numeric(levels(overlap$seqnames))[overlap$seqnames]
overlap$strand<-as.character(overlap$strand)

join<-full_join(dat, overlap,
        by=c("chromosome_name"="seqnames",
                intersect(colnames(dat), colnames(overlap))))

join$foundYT<-1

#The width column comes from the overlap table, so if this value is missing then they didn't overlap.
join$foundYT[is.na(join$width)]<-0

unique(join$gene_name[join$foundYT==0])

join2<-join[c("gene_name", "chromosome_name", "start", "end", "foundYT")]
join2<-join2%>%distinct()

write.csv(join2, "/scratch/mf91122/PUFA-GWAS/SMultiXcan/huifang2/SMultiXcan_gene_novelty_03182022.csv")


#Stopped here 03-18-2022
#Find overlaps with discovery and meta-analysis loci results

###*****************************************************************************************
###Part 5 FUMA overlaps*********************************************************************
###*****************************************************************************************

FUMAdisc<-as_tibble(read.csv("/work/kylab/mike/PUFA-GWAS/FUMA-final/GenomicLoci-DISC-02102022.csv"))
FUMAma<-as_tibble(read.csv("/work/kylab/mike/PUFA-GWAS/FUMA-final/GenomicLoci-MA-UKBEUR-MET-KET-02102022.csv"))
join2<-read.csv("/scratch/mf91122/PUFA-GWAS/SMultiXcan/huifang2/SMultiXcan_gene_novelty_03182022.csv")
join2<-as_tibble(join2%>%select(-X))

#Harmonize phenotype names
phenotypes<-c("w3FA_NMR_resinv", "w3FA_NMR_TFAP_resinv", "w6FA_NMR_resinv",
                "w6FA_NMR_TFAP_resinv", "w6_w3_ratio_NMR_resinv",
                "DHA_NMR_resinv", "DHA_NMR_TFAP_resinv", "LA_NMR_resinv",
                "LA_NMR_TFAP_resinv", "PUFA_NMR_resinv", "PUFA_NMR_TFAP_resinv",
                "MUFA_NMR_resinv", "MUFA_NMR_TFAP_resinv", "PUFA_MUFA_ratio_NMR_resinv")

phenotypes2<-str_split(phenotypes, "_resinv", simplify=T)[,1]


phenotypesMA<-c("W3", "w6_METAL", "DHA", "LA_METAL", "MUFA_METAL")

phenotypenames<-c("Omega-3", "Omega-3 (TFAP)", "Omega-6" , "Omega-6 (TFAP)",
                  "Omega-6:Omega-3 ratio", "DHA", "DHA (TFAP)", "LA", "LA (TFAP)",
                  "PUFAs", "PUFAs (TFAP)", "MUFAs", "MUFAs (TFAP)", "PUFAs:MUFAs ratio")

phenotypenamesMA<-c("Omega-3","Omega-6" ,"DHA","LA","MUFAs")

FUMAdisc$Phenotype<-mapvalues(FUMAdisc$Phenotype, from=phenotypes2, to=phenotypenames)
FUMAma$Phenotype<-mapvalues(FUMAma$Phenotype, from=phenotypesMA, to=phenotypenamesMA)



r=500000
FUMAdisc$start<-FUMAdisc$start-r
FUMAdisc$end<-FUMAdisc$end+r

FUMAma$start<-FUMAma$start-r
FUMAma$end<-FUMAma$end+r

#Overlap FUMA discovery--=-==-=-=-=-=-=-=-=-=-=--=-=--=--==--=-=-=-=-=-=-=-=

dat2<-join2

#Make GR
DAT2gr<-makeGRangesFromDataFrame(dat2,
                         keep.extra.columns=T,
                         ignore.strand=T)

FUMAdiscgr<-makeGRangesFromDataFrame(FUMAdisc,
                         keep.extra.columns=T,
                         ignore.strand=T)

#subsetByOverlaps subsets the first GRanges object 
#to include only those that overlap the second.
overlap<-as_tibble(as.data.frame(
        subsetByOverlaps(DAT2gr, FUMAdiscgr, ignore.strand=TRUE),
                stringsAsFactors=F))
nrow(overlap) #497

overlap$seqnames<-as.numeric(levels(overlap$seqnames))[overlap$seqnames]
overlap$strand<-as.character(overlap$strand)

join2b<-full_join(dat2, overlap,
        by=c("chromosome_name"="seqnames",
                intersect(colnames(dat2), colnames(overlap))))

join2b$foundFUMAdisc<-1

#The width column comes from the overlap table, so if this value is missing then they didn't overlap.
join2b$foundFUMAdisc[is.na(join2b$width)]<-0

unique(join2b$gene_name[join2b$foundFUMAdisc==0])

join2b<-join2b%>%select(-width, -strand)

join2b<-join2b%>%distinct()



#Overlap FUMA Meta-analysis--=-==-=-=-=-=-=-=-=-=-=--=-=--=--==--=-=-=-=-=-=-=-=
dat3<-join2

#Make GR
DAT3gr<-makeGRangesFromDataFrame(dat3,
                         keep.extra.columns=T,
                         ignore.strand=T)

FUMAmagr<-makeGRangesFromDataFrame(FUMAma,
                         keep.extra.columns=T,
                         ignore.strand=T)

#subsetByOverlaps subsets the first GRanges object 
#to include only those that overlap the second.
overlap<-as_tibble(as.data.frame(
        subsetByOverlaps(DAT3gr, FUMAmagr, ignore.strand=TRUE),
                stringsAsFactors=F))
nrow(overlap) #459 

overlap$seqnames<-as.numeric(levels(overlap$seqnames))[overlap$seqnames]
overlap$strand<-as.character(overlap$strand)

join3b<-full_join(dat3, overlap,
        by=c("chromosome_name"="seqnames",
                intersect(colnames(dat2), colnames(overlap))))

join3b$foundFUMAma<-1

#The width column comes from the overlap table, so if this value is missing then they didn't overlap.
join3b$foundFUMAma[is.na(join3b$width)]<-0

unique(join3b$gene_name[join3b$foundFUMAma==0])

join3b<-join3b%>%select(-width, -strand)

join3b<-join3b%>%distinct()

final<-inner_join(join2b, join3b)

write.csv(final, "/scratch/mf91122/PUFA-GWAS/SMultiXcan/huifang2/SMultiXcan_overlap.03192022.csv",
	row.names=F, quote=F)
