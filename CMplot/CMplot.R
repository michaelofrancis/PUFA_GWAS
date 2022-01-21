library(manhattan) #https://github.com/boxiangliu/manhattan/blob/master/vignettes/manhattan.pdf
suppressMessages(library(CMplot)) #https://github.com/YinLiLin/CMplot/blob/master/User%20Manual%20for%20CMplot.pdf
suppressMessages(library(tidyverse))
library(ggrepel)

indir<-"/scratch/mf91122/PUFA-GWAS/METAL/M2/STDERR-postprocess"

pheno<-c("FAw3", "FAw6", "DHA", "LA", "MUFA")

file<-list()
file2<-list()
for (p in 1:length(pheno)){
#p=1

file[[p]]<-as_tibble(read.table(paste(indir, "/META_IVW_", pheno[p], ".STDERR1.POSTPROCESS.Nb.txt", 
			sep=""),
			header=T, stringsAsFactors=F))


file2[[p]]<-file[[p]]%>%select(SNP, CHR, BP, P)

if (min(file2[[p]]$P)==0){
    file2[[p]]$P[file2[[p]]$P==0]<-min(file2[[p]]$P[file2[[p]]$P != min(file2[[p]]$P)])
}

if (max(file2[[p]]$P)==1){
file2[[p]]$P[file2[[p]]$P==1]<-max(file2[[p]]$P[file2[[p]]$P != max(file2[[p]]$P)])
}

#file2[[p]]$logP=-log10(file2[[p]]$P)
colnames(file2[[p]])<-c("SNP", "CHR", "BP", paste("P_", pheno[p], sep=""))
#, paste("logP_", pheno[p], sep=""))

}

dat<-inner_join(file2[[1]], file2[[2]], by=c("SNP", "CHR", "BP"))
dat<-inner_join(dat, file2[[3]], by=c("SNP", "CHR", "BP")) 
dat<-inner_join(dat, file2[[4]], by=c("SNP", "CHR", "BP"))
dat<-inner_join(dat, file2[[5]], by=c("SNP", "CHR", "BP"))

write.csv(dat, "/scratch/mf91122/PUFA-GWAS/METAL/M2/STDERR-postprocess/plot3/CMplot.csv",
	row.names=F, quote=F)

#***************************

dat<-read.csv("/Users/mike/Documents/Research/PUFA-GWAS/CMplot.csv")
dat<-as_tibble(dat)

cutoff=100
dat$P_FAw3[-log10(dat$P_FAw3)>cutoff]=1*10^(-cutoff)
#nrow(dat[-log10(dat$P_FAw3)>cutoff,]) #[1] 140
dat$P_FAw6[-log10(dat$P_FAw6)>cutoff]=1*10^(-cutoff)
#nrow(dat[-log10(dat$P_FAw6)>cutoff,]) #[1] 2
dat$P_DHA[-log10(dat$P_DHA)>cutoff]=1*10^(-cutoff)
#nrow(dat[-log10(dat$P_DHA)>cutoff]) #[1] 132
dat$P_LA[-log10(dat$P_LA)>cutoff]=1*10^(-cutoff)
#nrow(dat[-log10(dat$P_LA)>cutoff,]) #[1] 2
dat$P_MUFA[-log10(dat$P_MUFA)>cutoff]=1*10^(-cutoff)
#nrow(dat[-log10(dat$P_MUFA)>cutoff,]) #[1] 3

start_time <- Sys.time()
CMplot(dat, 
       type="p",
       plot.type="c", #circular
       LOG10=T, #change P vals into log 
       band=0.8, #space between chromosomes
       r=2, #radius of circle
       cir.legend=T, #whether to include legend of each circle (P value labels and lines)
       outward=TRUE, #dots go out or in
       cex=0.3, #dot size
       cir.legend.col="black",
       cir.chr.h=1, #width of chromosome boundary
       chr.den.col="grey", #color of chromosome boundary
       chr.labels=paste("Chr", 1:22, sep=""), #labels for chromosomes
       threshold = c(5e-08, 5e-300), #set the second number too high to get rid of grey lines
       threshold.col = c("red", "green"),
       threshold.lty = 2, #threshold line type.
       amplify=F, #make significant points (re threshold) bigger 
       file="jpg",
       memo="",
       dpi=300,
       file.output=TRUE,
       verbose=TRUE,
       width=10,
       height=10
)

end_time <- Sys.time()
end_time - start_time
