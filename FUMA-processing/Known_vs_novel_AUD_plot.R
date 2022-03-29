library(tidyverse)
library(gridExtra)
library(ggpattern)

#GGpattern reference:
#https://stackoverflow.com/questions/62393159/how-can-i-add-hatches-stripes-or-another-pattern-or-texture-to-a-barplot-in-ggp
#https://www.rdocumentation.org/packages/ggpattern/versions/0.4.2

###--------------------------------------------------------------
#MAKE PLOT- -----------------------------------------------------
###--------------------------------------------------------------
joindplot<-read.csv("alcohol_GeneSets_novelvsnon_joind_plot_v2.csv",header=T)
#joindplot<-joindplot[1:8,] #ALC USE DISORDER ONLY
joindplot$logP<-(-log10(joindplot$P))
joindplot$logadjP<-(-log10(joindplot$adjP))
joindplot$novelty<-rep(c("Novel", "Known"), 4)
joindplot$Phenotype<-str_split(joindplot$Phenotype, "\\.", simplify=T)[,1]


lev<-c("w3","DHA",
       "w6","LA")
joindplot$Phenotype<-factor(joindplot$Phenotype, levels=lev)


#joindplot$GeneSet<-wrap.labels(joindplot$GeneSet, 10)


#PLOT-----------------------==========================
w3n="turquoise2"
w3k="turquoise3"
dhan="chartreuse3"
dhak="olivedrab3"
w6n="tan2"
w6k="tan3"
lan="gold"
lak="gold2"


w3="turquoise2"
dha="chartreuse3"
w6="tan2"
la="gold"



#
p1<-ggplot(joindplot, mapping =aes(x=Phenotype,
                            y=logadjP)) +
    #geom_bar(position="dodge", stat="identity", color="black")+
    geom_col_pattern(
        aes(pattern = novelty, fill=Phenotype, pattern_fill=Phenotype), 
        position = "dodge",
        stat = "identity",
        colour          = 'black', 
        fill=c(w3,w3,dha,dha,w6,w6,la,la),
        pattern_density = 0.35 ,
        pattern_fill    = 'black',
       pattern_colour  = 'black',
       pattern_key_scale_factor = 0.3
    )+
    theme_bw()+
    theme(
          plot.title = element_text(hjust = 0.5)
          #legend.position="none"
          )+
    scale_pattern_fill_manual(values = c(w3, w3,
                                   dha, dha,
                                  w6, w6,
                                  la, la ))+
    scale_pattern_manual(values = c(Known = "stripe", Novel = "none")) +
    labs(\
         y = "-log10(adjP) for gene set enrichment", 
         pattern = "Novelty of PUFA-associated \nvariants mapped to genes") +
    ggtitle("Enrichment of gene sets mapped to 
significant variants associated with PUFA traits
for alcohol use disorder (total score)")+
   guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))

p1







