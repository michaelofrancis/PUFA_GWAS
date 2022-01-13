library(biomaRt)
#suppressMessages(library(qqman))
library(manhattan) #https://github.com/boxiangliu/manhattan/blob/master/vignettes/manhattan.pdf
suppressMessages(library(tidyverse))
library(ggrepel)

source("/work/kylab/mike/PUFA-GWAS/SMultiXcan/plot/lociLeadSNP.R")

pheno<-c("w3FA_NMR","w3FA_NMR_TFAP",
                "w6FA_NMR", "w6FA_NMR_TFAP",
                "w6_w3_ratio_NMR",
                "DHA_NMR","DHA_NMR_TFAP",
                "LA_NMR","LA_NMR_TFAP",
                "PUFA_NMR","PUFA_NMR_TFAP",
                "MUFA_NMR", "MUFA_NMR_TFAP",
                "PUFA_MUFA_ratio_NMR")

pheno<-paste(pheno, "_resinv", sep="")


TWASdir<-"/scratch/mf91122/PUFA-GWAS/SMultiXcan/huifang2"

grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", 
	path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

ensembl = useDataset("hsapiens_gene_ensembl",mart=grch37)

t2g<-getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",
		'chromosome_name','start_position','end_position'), mart = ensembl)

for (p in 1:length(pheno)){
#p=1

#Load SMultiXcan output and get positions from gene names
file<-as_tibble(read.table(paste(TWASdir, "/05.", pheno[p],".M2.txt.tab_smultixcan.txt",sep=""),
		header=T, stringsAsFactors=F))

file$ensembl_gene_id <- gsub("\\..*","", file$gene)

t2g<-getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",
	'chromosome_name','start_position','end_position'), mart = ensembl)

my_ids.version <- merge(file, t2g, by= 'ensembl_gene_id')
my_ids.version<-as_tibble(my_ids.version)

dat<-my_ids.version%>%select(gene_name, chromosome_name, start_position, pvalue)

#Clean up dataset
dat<-dat[grep("HG*", dat$chromosome_name, invert=T),] #43 weird entries removed
dat$chromosome_name<-as.numeric(dat$chromosome_name)
dat<-dat[!is.na(dat$pvalue),]


#change 0 pvalue to minimum calculated as I don't know what to do witht hem
if (min(dat$pvalue)==0){
dat$pvalue[dat$pvalue==0]<-min(dat$pvalue[dat$pvalue != min(dat$pvalue)])
}

dat$y=-log10(dat$pvalue)


#Designate genes to label
data_frame<-dat[dat$pvalue<5e-08,]
colnames(data_frame)<-c("gene_name", "CHR", "POS", "P.value", "y")

#Run lociLeadSNP
lociLeadSNP(data_frame, setdistance=5000000)

#If multiple leadsnps in one loci, combine the labels

lls2<-locileadSNPoutput[locileadSNPoutput$LEADSNP==1,]

for (i in 1:max(lls2$LOCI)){
	lls2$gene_name[lls2$LOCI==i]<-paste(lls2$gene_name[lls2$LOCI==i], collapse=",")	
}


lls2<-lls2[!duplicated(lls2$gene_name),]

label<-lls2%>%select(gene_name, CHR, POS)
colnames(label)<-c("Label", "CHR", "POS")

dat<-left_join(dat, label, by=c("chromosome_name"="CHR", "start_position"="POS"))


#Get table into manhattan package format
colnames(dat)<-c("gene_name", "chrom", "pos", "pval", "y", "Label") 
dat$chrom<-paste("chr", dat$chrom, sep="")


#Make plot
plotoutputpath=paste("/scratch/mf91122/PUFA-GWAS/SMultiXcan/huifang2/plot/05." ,
			pheno[p], "_SMultiXcan.png", sep="")

png(filename=plotoutputpath, type="cairo",
	res=400, width=3000,height=1500)

print(

manhattan(
dat, 
build='hg19',
color1='skyblue',color2='navyblue')+

geom_hline(yintercept=-log10(5e-8),color='red')
	+geom_label_repel(
			aes(label=Label),
			colour='black',
			fill = "white",
			size=1.5,
			box.padding = 0.5,
			max.overlaps = Inf,
			arrow = arrow(length = unit(0.008, "npc")),
			ylim = c(10,NA),
			segment.size=0.3,
			segment.alpha=0.3,
			min.segment.length = 0,
			nudge_x           = 3,
			    direction         = "both"
			)



)#close print

dev.off()


#https://github.com/bernatgel/karyoploteR_questions/blob/master/biostars_470889/biostars_470889.R


}
