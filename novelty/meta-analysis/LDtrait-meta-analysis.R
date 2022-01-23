#Cycle through lead SNPs and find previously reported SNPs in LD in the GWAS catalog.

tok="6a30732e002f"
outdir<-"/scratch/mf91122/PUFA-GWAS/novelty/meta-analysis"

dir.create(outdir, showWarnings=F)
setwd(outdir)
library(LDlinkR)
suppressMessages(library(dplyr))


#Set params--
cap=25
window=500000
#------------

rsID<-FUMA$rsID
sets<-ceiling(length(rsID)/cap)

#Loop through the list of rsIDs increments of sets since the software can only process 25 at a time.
x=1
	while(x<sets){
        print(paste("x is ", x))

        y=x*cap
        snps<-rsID[(y-(cap-1)):y]

        fileoutput<-paste(outdir, "/LDlinkR-", x, ".", window,".txt", sep="")

        LDtrait(snps,pop = "EUR", r2d = "r2", r2d_threshold = 0.1,
        win_size = window, token = tok, file = fileoutput)

        x=x+1
	}#end while

	#Get the last bit
	print("Last bit...")
        snps<-rsID[ (((sets-1)*cap)+1):length(rsID)]
	fileoutput<-paste(outdir, "/LDlinkR-", x, ".", window,".txt", sep="")

        LDtrait(snps,pop = "EUR",r2d = "r2",r2d_threshold = 0.1,
        win_size = window, token = tok, file = fileoutput)
