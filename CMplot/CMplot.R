#install.packages("CMplot")

library(CMplot)
library(tidyverse)
#https://github.com/YinLiLin/CMplot/blob/master/User%20Manual%20for%20CMplot.pdf

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

#get order of SNPs
dat<-dat%>%select(SNP,CHR,BP, P_MUFA, P_LA, P_FAw6, P_DHA, P_FAw3)


novel<-as_tibble(read.csv("/Users/mike/Documents/Research/PUFA-GWAS/novelty/meta-analysis/Novel-CMplot.csv"))
novel<-novel[novel$foundYT==0,]

novel%>%group_by(Phenotype)%>%summarize(n=n())



cutoff=100
novel$P_FAw3[-log10(novel$P_FAw3)>cutoff]=1*10^(-cutoff)
novel$P_FAw6[-log10(novel$P_FAw6)>cutoff]=1*10^(-cutoff)
novel$P_DHA[-log10(novel$P_DHA)>cutoff]=1*10^(-cutoff)
novel$P_LA[-log10(novel$P_LA)>cutoff]=1*10^(-cutoff)
novel$P_MUFA[-log10(novel$P_MUFA)>cutoff]=1*10^(-cutoff)



novel<-novel%>%select(rsID,chr,pos, P_MUFA, P_LA, P_FAw6, P_DHA, P_FAw3)
colnames(novel)[1:3]<-c("SNP", "CHR", "BP")
novel


toHighlight<-list()
toHighlight[[1]]<-novel$SNP[novel$P_MUFA<5e-08]
toHighlight[[2]]<-novel$SNP[novel$P_LA<5e-08]
toHighlight[[3]]<-novel$SNP[novel$P_FAw6<5e-08]
toHighlight[[4]]<-novel$SNP[novel$P_DHA<5e-08]
toHighlight[[5]]<-novel$SNP[novel$P_FAw3<5e-08]
toHighlight 



        start_time <- Sys.time()
CMplot(dat, 
       type="p",
       plot.type="c", #circular
       LOG10=T, #change P vals into log 
       band=0.8, #space between chromosomes
       r=3, #radius of circle
       cir.legend=T, #whether to include legend of each circle (P value labels and lines)
       outward=TRUE, #dots go out or in
       cex=0.3, #dot size
       cir.legend.col="black",
       cir.chr.h=1, #width of chromosome boundary
       chr.den.col="grey85", #color of chromosome boundary
       chr.labels=paste("Chr", 1:22, sep=""), #labels for chromosomes
       threshold = c(5e-08, 5e-300), #set the second number too high to get rid of grey lines
       threshold.col = c("red4", "green"),
       threshold.lty = 2, #threshold line type.
       amplify=F, #make significant points (re threshold) bigger 
       
       highlight=toHighlight,
       highlight.col = "red3",
       highlight.pch=17,#1=circle, 2=triangle, 3 = x, 4 = x, 5= diamond, 
       #6 = down triangle, 7 = cross out circle, 8=asterisk, 9=diamondcrosshairs,
       #10=circle crosshairs, 11=starofdavid, 12=4square, 13 crosshairs,14=trianglesquare,
       #15=filled square, 16=filled circle; 17=filled triangle; 18=filledsquare
       #19=filled circle, 20=filled circle
       highlight.cex = 0.8,
       #ylim=c(0,100),
       
       file="jpg",
       memo="",
       dpi=600, #set higher when ur ready to make a nice one (600 or more?)
       col=matrix(c("purple3","mediumpurple3",  #MUFA /innermost phenotype alternating colors (MUFA)
                    "gold2","gold",  #LA
                    "tan3","tan2",  #w6
                    "chartreuse3", "olivedrab3",   #DHA 
                    "turquoise3","turquoise2"), #w3 / outermost 
                  nrow=5, byrow=T),
       file.output=TRUE,
       verbose=TRUE,
       width=10,
       height=10
)

end_time <- Sys.time()
end_time - start_time
