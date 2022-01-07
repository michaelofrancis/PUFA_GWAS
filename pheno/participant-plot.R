library(tidyverse)

# fileorig<-as_tibble(read.csv("Part-char-table-forplot.csv"))
# 
# file
# N<-file[1,]
# file<-fileorig[-1,]
# file
# 
# pop<-c("EUR", "AFR", "CSA", "EAS")
# 
# #Extract SEs
# 
# 
# for (i in 1:length(pop)){
#     SE<- str_split(unlist(file[i+1]), pattern = "\\(", simplify = T)[,2]
#     SE<-str_split(SE, pattern = "\\)", simplify = T)[,1]
#     file$SE<-as.double(SE)
#     colnames(file)[colnames(file)=="SE"]<-paste("UKB.", pop[i], ".SE", sep="") 
#     file[i+1]<-as.double(str_split(unlist(file[i+1]), pattern = "\\(", simplify = T)[,1])
# }
#     
# 
# UKB.EUR.SE<-str_split(file$UKB.EUR, pattern = "\\(", simplify = T)[,2]
# UKB.EUR.SE<-str_split(UKB.EUR.SE, pattern = "\\)", simplify = T)[,1]
# 
# ratio<-file[grep(":", file$X),]
# ratio
# file[!file$X %in% ratio,]
# write.csv(file, "progress.csv")
#****************************************************************
tab<-as_tibble(read.csv("progress.csv"))

lev<-c("n3FA", "n3FA/Total FA", "n6FA","n6FA/Total FA", 
                     "DHA","DHA/Total FA","LA","LA/Total FA",
                     "PUFA","PUFA/Total FA","MUFA","MUFA/Total FA",
                     "n6:n3","PUFA:MUFA")
tab$Pheno<-factor(tab$Pheno, levels=lev)
tab$Pop<-factor(tab$Pop, levels = c("EUR", "AFR", "CSA", "EAS"))


#for (i in 1:length(pop)){colnames(file)[grep(pop[1], colnames(file))][1]

#mmol/L bar plot
ggplot(tab[grep(":|/", tab$Pheno, invert=T),],
       aes(fill=Pop, 
           y=Mean, 
           x=Pheno)) + 
    geom_bar(position="dodge", stat="identity")+
        xlab("Plasma fatty acid")+
        ylab("Mean (mmol/L)")+
        geom_errorbar(
            aes(ymin=Mean-SE, ymax=Mean+SE), 
            width=.2,
                      position=position_dodge(.9)
            )+theme_bw()+
    theme(axis.text.x = element_text(angle = 30, vjust = 0.6, hjust=0.5))+
    scale_fill_manual("legend", 
                      values = c("EUR" = "lightslateblue", "AFR" = "olivedrab", 
                                 "CSA" = "tan3", "EAS"="cadetblue"))


#Ratio bar plot
ggplot(tab[grep(":|/", tab$Pheno),],
       aes(fill=Pop, 
           y=Mean, 
           x=Pheno)) + 
    geom_bar(position="dodge", stat="identity")+
    xlab("Plasma fatty acid quantity")+
    geom_errorbar(
    aes(ymin=Mean-SE, ymax=Mean+SE), 
    width=.2,
    position=position_dodge(.9)
    ) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 30, vjust = 0.6, hjust=0.5))+
    scale_fill_manual("legend", 
        values = c("EUR" = "lightslateblue", "AFR" = "olivedrab", 
                   "CSA" = "tan3", "EAS"="cadetblue"))+ 
    geom_vline(xintercept=6.5) +
    scale_y_continuous(
        # Features of the first axis
        name = "Percentage (%)",
        
        # Add a second axis and specify its features
        sec.axis = sec_axis( trans=~.*1, name="Ratio (a:b)")
    )


#*What I called SE in this script was SD the whole time.
