library(tidyverse)
library(cowplot)

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
# write.csv(file, "part.char.csv")
#****************************************************************
tab<-as_tibble(read.csv("part.char.csv"))

lev<-c("n3FA", "n3FA/Total FA",  
                     "DHA","DHA/Total FA",
       "n6FA","n6FA/Total FA",
       "LA","LA/Total FA",
                     "PUFA","PUFA/Total FA","MUFA","MUFA/Total FA",
                     "n6:n3","PUFA:MUFA")
tab$Pheno<-factor(tab$Pheno, levels=lev)
tab$Pop<-factor(tab$Pop, levels = c("EUR", "AFR", "CSA", "EAS"))


#*What I called SE in this script was SD the whole time!!!~~~~

# colors<-c("#00BFFF", "#FF4B19", "#FFC219", "#0ACC92")
# #colors<-c("#00A3D9", "#D94016", "#D9A516", "#0BD99C")
# colors<-c("#0CEBA8", "#EB4517", "#EBB217", "#00B0EB") #these look
colors<-c("#0CEBA8", "#EB2358", "#EB8823", "#17B6EB") #https://color.adobe.com/create/color-accessibility


tiff("Sup1A.tiff", 
     width = 7, height = 5, 
     units = 'in', res = 600)

#mmol/L bar plot
p1<-ggplot(tab[grep(":|/", tab$Pheno, invert=T),],
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
            )+
    theme(axis.text.x = element_text(angle = 30, vjust = 0.6, hjust=0.5))+
    theme_bw()+
    scale_fill_manual("UKB ancestry", 
                      values = c("EUR" = colors[1], "AFR" = colors[2], 
                                 "CSA" = colors[3], "EAS"=colors[4]),
                      labels = c("EUR (N=101,729)", 
                                 "AFR (N=1,564)", 
                                 "CSA (N=2,203)", 
                                 "EAS (N=633)"))

p1

dev.off()


tiff("Sup1B.tiff", 
     width = 7, height = 5, 
     units = 'in', res = 600)

#Ratio bar plot
p2<-ggplot(tab[grep(":|/", tab$Pheno),],
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
        values = c("EUR" = colors[1], "AFR" = colors[2], 
                   "CSA" = colors[3], "EAS"=colors[4]),
        labels = c("EUR (N=101,729)", 
                   "AFR (N=1,564)", 
                   "CSA (N=2,203)", 
                   "EAS (N=633)"))+ 
    geom_vline(xintercept=6.5) +
    scale_y_continuous(
        # Features of the first axis
        name = "Percentage (%)",
        # Add a second axis and specify its features
        sec.axis = sec_axis( trans=~.*1, name="Ratio (a:b)")
    )

p2

dev.off()



legend <- get_legend(
    # create some space to the left of the legend
    p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)


pgrid<-cowplot::plot_grid(p1+ theme(legend.position="none"),
                          p2+theme(legend.position="none"),
                          labels = c("S1A","S1B"), hjust=-0.01,vjust=1.2,
                          ncol=1)

tiff("FigS1-hq.tiff", width = 7, height = 10, units = 'in', res = 300)
plot_grid(pgrid, NULL, legend,NULL,ncol = 4, rel_widths = c(1,0.05,0.2,0.05))
dev.off()
