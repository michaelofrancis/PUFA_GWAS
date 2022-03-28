#Reference for displaying results
#https://doi.org/10.1038/s41593-020-0643-5

suppressMessages(library(plyr))
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
files<-list.files(dir4, full.names=T)
files<-files[grep("significant", x=files, invert=T)]


phenoresults<-list()
#one pheno at a time
for (p in 1:length(phenotypes)){

set<-files[grepl(phenotypes[p],files)]
tis<-str_split(set, "__", simplify=T)[,3]
tis<-str_split(tis, ".csv", simplify=T)[,1]

#Read in data from many tissues for one phenotype
part4<-list()
for (i in 1:length(set)){

	part4[[i]]<-as_tibble(read.csv(set[i], header=T))
	part4[[i]]$tissue<-tis[i]
	part4[[i]]$Phenotype<-phenotypes[p]
}

part4<-do.call(rbind, part4)
part4<-part4[!is.na(part4$pvalue),]

#count number of gene-tissue combinations tested to get bf correction number---
sum<-nrow(unique(part4[c(2,15)]))
#print(paste(phenotypes[p], sum))
Pcutoff=0.05/(sum*14)
#------

part4$Phenotype<-phenotypenames[p]

part4sig<-part4[part4$pvalue<Pcutoff,]

phenoresults[[p]]<-part4sig

}#end phenotypes loop

phenoresults<-do.call(rbind, phenoresults)

phenoresults<-phenoresults%>%select(Phenotype, tissue, gene, gene_name, zscore, pvalue, var_g, 
		n_snps_used, n_snps_in_model)
#write.csv(phenoresults, paste(dir4, "/significant.04.SPrediXcan.csv", sep=""), row.names=F, quote=F)

##-----------------------------------------------------------------------------
##05.SMultiXcan----------------------------------------------------------------
##-----------------------------------------------------------------------------

dir5="/scratch/mf91122/PUFA-GWAS/SMultiXcan/huifang2/05.SMultiXcan"
files<-list.files(dir5, full.names=T)
files<-files[grep("significant", x=files, invert=T)]
filenames<-list.files(dir5)
filenames<-filenames[grep("significant", x=filenames, invert=T)]

#Read in tables
part5<-list()
for (i in 1:length(files)){
	part5[[i]]<-as_tibble(read.table(files[i], header=T))
	part5[[i]]<-part5[[i]][! is.na(part5[[i]]$pvalue) ,]
}

#Get names for tables
phen<-str_split(filenames , pattern="_resinv", simplify=T)[,1]
phen<-str_split(phen, "05.", simplify=T)[,2]
names(part5)<-filenames 

phen<-str_split(filenames , pattern="05.", simplify=T)[,2]
phen<-str_split(phen, ".M2.txt.tab_smultixcan.txt", simplify=T)[,1]

length(unique(part5[[1]]$gene_name))
#[1] 21850

#21850 unique genes to Bonferroni correction

Pcutoff<-0.05/(21850*14)
#[1] 2.28833e-06

part5sig<-list()
for (i in 1:length(files)){
	part5sig[[i]]<-part5[[i]][part5[[i]]$pvalue <Pcutoff,]

	part5sig[[i]]$Phenotype<-phen[i] 

	part5sig[[i]]<-part5sig[[i]]%>%select("Phenotype", "gene","gene_name","pvalue","n","n_indep",
		"p_i_best","t_i_best","p_i_worst","t_i_worst")
	
}

part5out<-do.call(rbind, part5sig)

part5out$Phenotype<-mapvalues(x=part5out$Phenotype, from=phenotypes, to=phenotypenames)

#write.csv(part5out, paste(dir5, "/significant.05.SMultiXcan.csv", sep=""), row.names=F, quote=F)


#Summarise pt 4 results
summary<-phenoresults%>%group_by(Phenotype)%>%summarise(countpt4 = n())
length(unique(phenoresults$gene_name))
#[1] 567

sum5<-part5out%>%group_by(Phenotype)%>%summarise(countpt5 = n())
total_summary<-full_join(summary, sum5, by="Phenotype")

part5out%>%group_by(gene_name)%>%summarise(genecount = n())%>%arrange(genecount)

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
library(biomaRt)
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
suppressMessages(library(GenomicRanges))
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

overlap$seqnames<-as.integer(overlap$seqnames)
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
FUMAdisc<-as_tibble(read.csv("/work/kylab/mike/PUFA-GWAS/FUMA-final/GenomicLoci-DISC-02102022.csv"))
FUMAma<-as_tibble(read.csv("/work/kylab/mike/PUFA-GWAS/FUMA-final/GenomicLoci-MA-UKBEUR-MET-KET-02102022.csv"))




#Harmonize phenotype names
phenotypes<-c("w3FA_NMR_resinv", "w3FA_NMR_TFAP_resinv", "w6FA_NMR_resinv",
                "w6FA_NMR_TFAP_resinv", "w6_w3_ratio_NMR_resinv",
                "DHA_NMR_resinv", "DHA_NMR_TFAP_resinv", "LA_NMR_resinv",
                "LA_NMR_TFAP_resinv", "PUFA_NMR_resinv", "PUFA_NMR_TFAP_resinv",
                "MUFA_NMR_resinv", "MUFA_NMR_TFAP_resinv", "PUFA_MUFA_ratio_NMR_resinv")


