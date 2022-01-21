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

start_time <- Sys.time()
CMplot(dat, 
       type="p",
       plot.type="c", #circular
       band=0.8, #space between chromosomes
       r=2, #radius of circle
       cir.legend=T, #whether to include legend of each circle (P value labels and lines)
       outward=TRUE, #dots go out or in
       cex=0.3, #dot size
       cir.legend.col="black",
       cir.chr.h=1, #width of chromosome boundary
       chr.den.col="grey", #color of chromosome boundary
       threshold = 5e-08,
       threshold.col = "red",
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
