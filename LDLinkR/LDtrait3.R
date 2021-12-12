#Cycle through lead SNPs and find SNPs in LD in the GWAS catalog.


tok="6a30732e002f"

dir.create("/scratch/mf91122/PUFA-GWAS/LDlinkR", showWarnings=F)
setwd("/scratch/mf91122/PUFA-GWAS/LDlinkR")
library(LDlinkR)
suppressMessages(library(dplyr))

phenotypes<-c("w3FA_NMR","w3FA_NMR_TFAP",
                "w6FA_NMR", "w6FA_NMR_TFAP",
                "w6_w3_ratio_NMR",
                "DHA_NMR","DHA_NMR_TFAP",
                "LA_NMR","LA_NMR_TFAP",
                "PUFA_NMR","PUFA_NMR_TFAP",
                "MUFA_NMR", "MUFA_NMR_TFAP",
                "PUFA_MUFA_ratio_NMR")


FUMA<-as_tibble(read.csv("/work/kylab/mike/PUFA-GWAS/FUMA/GenomicLoci12112021.csv",
		header=T, stringsAsFactors=F))

#loop
for (p in 1:length(phenotypes)){
#p=1
leadvars<-FUMA[FUMA$Phenotype==phenotypes[p],]
#leadvars$rsID[!grepl("rs", leadvars$rsID)]<-paste("chr", 
#		leadvars$chr[!grepl("rs", leadvars$rsID)], 
#		":", leadvars$pos[!grepl("rs", leadvars$rsID)], sep="")

leadvars<-leadvars[grepl("rs", leadvars$rsID),]

sets<-ceiling(nrow(leadvars)/50)
outdir<-"/scratch/mf91122/PUFA-GWAS/LDlinkR/12112021/"

if (sets>1){ #multiple sets
	x=1
	while(x<sets){
        print(paste("x is ", x))

        y=x*50
        temp<-leadvars[(y-49):y,]

        leadsnps<-temp$rsID
        fileoutput<-paste(outdir, phenotypes[p],
                "-LDlinkR-", x, ".txt", sep="")

        LDtrait(leadsnps,pop = "GBR", r2d = "r2", r2d_threshold = 0.1,
        win_size = 5e+05, token = tok, file = fileoutput)

        x=x+1
	}#end while

	print("after end while")
        temp<-leadvars[ (((sets-1)*50)+1): nrow(leadvars),]

        leadsnps<-temp$rsID
        fileoutput<-paste(outdir, phenotypes[p],
               "-LDlinkR-", x, ".txt", sep="")

        LDtrait(leadsnps,pop = "GBR",r2d = "r2",r2d_threshold = 0.1,
        win_size = 5e+05, token = tok, file = fileoutput)

} else { #not multiple sets

print("else not multiple sets")
leadsnps<-leadvars$rsID

fileoutput<-paste(outdir, phenotypes[p],
                "-LDlinkR.txt", sep="")
LDtrait(leadsnps, pop = "GBR", r2d = "r2", r2d_threshold = 0.1,
        win_size = 5e+05, token = tok, file = fileoutput)
}#end else


} #end phenotype loop


#x<-read.delim("/scratch/mf91122/PUFA-GWAS/LDlinkR/w3FA_NMR-LDlinkR.txt")
#x<-as_tibble(x)
