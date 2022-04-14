#Cycle through lead SNPs and find previously reported SNPs in LD in the GWAS catalog.

tok="6a30732e002f"
outdir<-"/scratch/mf91122/PUFA-GWAS/novelty/LDtrait-again"

dir.create(outdir, showWarnings=F)
setwd(outdir)
library(LDlinkR)
suppressMessages(library(dplyr))


#Set params--
window=1000000
#------------

FUMA<-as_tibble(read.csv("/work/kylab/mike/PUFA-GWAS/FUMA-final/GenomicLoci-DISC-04092022.csv",
		header=T, stringsAsFactors=F))
FUMAm<-as_tibble(read.csv("/work/kylab/mike/PUFA-GWAS/FUMA-final/GenomicLoci-MA-UKBEUR-MET-KET-04132022.csv",
                header=T, stringsAsFactors=F))

rsID<-unique(FUMA$rsID)
rsID2<-unique(FUMAm$rsID)
rsID<-unique(c(rsID, rsID2))


for (x in rsID){

        fileoutput<-paste(outdir, "/LDtrait-", x, ".", window,".txt", sep="")

	tryCatch({
	print(x)
        LDtrait(x,pop = "EUR", r2d = "r2", r2d_threshold = 0.1,
        win_size = window, token = tok, file = fileoutput)
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}
