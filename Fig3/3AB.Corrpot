suppressMessages(silent <- lapply(
    c("plyr", "dplyr", "tidyverse", "ggpubr", "ggpmisc",
      "ggforce", "ggcorrplot2", 
      "data.table", "Hmisc", "corrplot", "reshape2", "viridisLite"), 
    library, character.only=T))


# PHENOTYPIC CORRELATION ----------------------------------------
df<-as_tibble( 
        read.table(
            "/Users/mike/Documents/Research/PUFA-GWAS/pheno/PUFA_GWAS_pheno_M2.txt",header=T))

df<-df[c(16:29)]
df

#For Hclust------
colnames(df)<-c("n3FA","n3FA.Total.FA","n6FA","n6FA.Total.FA","n6.n3",
                "DHA","DHA.Total.FA","LA","LA.Total.FA","PUFA","PUFA.Total.FA","MUFA",
                "MUFA.Total.FA","PUFA.MUFA")

ordercols<-c("n6.n3", "PUFA.Total.FA","PUFA.MUFA", "n6FA.Total.FA",
             "LA.Total.FA","n3FA", "DHA", "n3FA.Total.FA",
             "DHA.Total.FA", "PUFA",  "n6FA", "LA",
             "MUFA", "MUFA.Total.FA")

df<-df[ordercols]
df

#From:https://stackoverflow.com/questions/16097453/how-to-compute-p-value-and-standard-error-from-correlation-analysis-of-rs-cor
cor.test.plus <- function(x) {
    list(x, 
         Standard.Error = unname(sqrt((1 - x$estimate^2)/x$parameter)))
}

corr.results <- tibble(p1="phenotype1",p2="phenotype2",rp=0, rpse=0)
#do correlation
for (i in 1:ncol(df)){
    for (j in 1:ncol(df)){
        corr <- cor.test(df[[i]],df[[j]],method="pearson")
        corr <- corr$estimate
        corr.se <-cor.test.plus(cor.test(df[[i]],df[[j]],method="pearson"))
        corr.se <- corr.se$Standard.Error
        corr.results <- add_row(corr.results,p1=colnames(df[,i]),p2=colnames(df[,j]),rp=corr, rpse=corr.se)
    }
}

corr.results
#transfer format to matrix
pc_matrix1 <- corr.results[-1,]
pc_matrix1$p1 <- factor(pc_matrix1$p1, levels=unique(pc_matrix1$p1))
pc_matrix1$p2 <- factor(pc_matrix1$p2, levels=unique(pc_matrix1$p2))

pc_matrix <- acast(pc_matrix1,p1~p2)
pc_matrix <- as.matrix(pc_matrix)
pc_matrix #Symmetrical phenotypic correlation matrix. Should go below diag in final.
pcse_matrix<-acast(pc_matrix1,p1~p2, value.var = "rpse")
pcse_matrix<- as.matrix(pcse_matrix)
pcse_matrix
#write.csv(pcse_matrix, "phenotypic-correlation-SE-matrix.csv")


# Get p values --------------------------------------------------
#p.mat <- cor_pmat(df)


#write.csv(p.mat, "phenotypic.correlation.p.mat.csv")



# GENETIC CORRELATION -------------------------------------------


gencorr<-read.csv("/Users/mike/Documents/Research/PUFA-GWAS/discovery/ldsc_genetic_corr.csv")
gencorr
#pc_matrix

map1<-c("DHA_NMR_","DHA_NMR_TFAP_","LA_NMR_", "LA_NMR_TFAP_", 
        "MUFA_NMR_","MUFA_NMR_TFAP_","PUFA_NMR_","PUFA_NMR_TFAP_",  
        "w3FA_NMR_","w3FA_NMR_TFAP_","w6FA_NMR_","w6FA_NMR_TFAP_",  
        "w6_w3_ratio_NMR_","PUFA_MUFA_ratio_NMR_")

map2<-c("DHA","DHA.Total.FA","LA","LA.Total.FA",
        "MUFA","MUFA.Total.FA","PUFA","PUFA.Total.FA",
        "n3FA","n3FA.Total.FA","n6FA","n6FA.Total.FA","n6.n3",
        "PUFA.MUFA")
gencorr$p1<-mapvalues(gencorr$p1, from = map1, to=map2) #warning=ok
gencorr$p2<-mapvalues(gencorr$p2, from = map1, to=map2)#warning=ok

ordercols<-c("n6.n3", "PUFA.Total.FA","PUFA.MUFA", "n6FA.Total.FA",
             "LA.Total.FA","n3FA", "DHA", "n3FA.Total.FA",
             "DHA.Total.FA", "PUFA",  "n6FA", "LA",
             "MUFA", "MUFA.Total.FA")


gencorr$p1 <- factor(gencorr$p1, levels=ordercols, ordered = T)
gencorr$p2 <- factor(gencorr$p2, levels=ordercols, ordered = T)
as_tibble(gencorr)
gencorr <- acast(gencorr,p1~p2)
gencorr <- as.matrix(gencorr)

##write.csv(gencorr, "genetic_correlation.csv", quote=F)


# CORRPLOT ------------------------------------------------------

p.mat<-read.csv("GC_analysis/correlation_matrix_pvalue-fix.csv")
row.names(p.mat)<-p.mat$X
p.mat<-p.mat[-1]
p.mat<-as.matrix(p.mat)
p.mat

se.mat<-read.csv("full-correlation-SE-matrix.csv")
row.names(se.mat)<-se.mat$X
se.mat<-se.mat[-1]
se.mat<-as.matrix(se.mat)
se.mat

h<-read.csv("/Users/mike/Documents/Research/PUFA-GWAS/discovery/heatmap_data_hclustorder_02252022.csv")
h<-h[-1]
h<-as.matrix(h)
row.names(h)<-colnames(h)

labels<-c("n-6:n-3", "PUFAs pct","PUFAs:MUFAs", "n6 pct",
             "LA pct","n-3", "DHA", "n-3 pct",
             "DHA pct", "PUFAs",  "n-6", "LA",
             "MUFAs", "MUFAs pct")
colnames(h)<-labels
rownames(h)<-labels
colnames(p.mat)<-labels
rownames(p.mat)<-labels


tiff("4A.corrplot.tiff",
     width = 4.5, height = 4.5,
     units = 'in', res = 600)

ggcorrplot(h, p.mat = p.mat, show.diag = FALSE, sig.lvl = (0.05/2.98), )

dev.off()

tiff("4A.corrplot.SQUARE.tiff",
     width = 4.5, height = 4.5,
     units = 'in', res = 600)

ggcorrplot(h, p.mat = p.mat, show.diag = FALSE, 
           sig.lvl = (0.05/2.98),method="square")+
    theme(axis.text.x = element_text(angle = 45, hjust=1))


dev.off()

          
          
          
          
          

# Kaixiong's scatterplot ----------------------------------------
          
          ### convert both
          cordata <- h
          corpvalue<- p.mat
          corse<-se.mat
          
          num_factors <- nrow(cordata)
          all.names   <- rownames(cordata)
          all.names   <- as.character(all.names)
          results <- data.frame()
          
          for (i in 1:num_factors){
              for (j in 1:i){
                  if (i == j){next}
                  phe <- as.character(cordata[i,j])
                  gen <- as.character(cordata[j,i])
                  phe.p <- as.character(corpvalue[i,j])
                  gen.p <- as.character(corpvalue[j,i])
                  phe.se <- as.character(corse[i,j])
                  gen.se <- as.character(corse[j,i])
                  indicator <- 0
                  
                  if((as.numeric(phe.p) < 0.05) & (as.numeric(gen.p) < 0.05)){indicator <- 1}
                  
                  newrow <- c(all.names[i], all.names[j], phe, gen, phe.p, gen.p, phe.se, gen.se, as.character(indicator))
                  results <- rbind(results, newrow, stringsAsFactors=FALSE)
              }
          }
          
          colnames(results) <- c("trait1", "trait2", "Phe_cor", "Gen_cor", "Phe_P", "Gen_P","Phe_SE", "Gen_SE", "sig_both")
          
          #write.table(results, file="correlation.dataframe.final.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
          
          
          # ## plot the figure
          # data<-results
          # data <- read.table("GC_analysis/correlation.dataframe.final.txt", header=TRUE)
          # 
          # idx <- data$sig_both == 1
          # #pdf("comparing_phenotypic_genetic_correlation.pdf", h=4, w=4)
          # par(mar=c(2.5,2.5,1,1)+0.1, mgp=c(1.5,0.5,0))
          # 
          # plot(data$Phe_cor, data$Gen_cor, pch=20, xlab="Phenotypic Correlation", ylab="Genetic Correlation")
          # abline(a=0, b=1, lty=2)
          # abline(v=0, lty=2)
          # abline(h=0, lty=2)
          # points(data$Phe_cor[idx], data$Gen_cor[idx], pch=20, col=c("red"))
          # lr <- lm(data$Gen_cor[idx]~data$Phe_cor[idx]-1)
          # summary(lr)
          # abline(lr, col=c("red"), lty=1, lwd=2)
          # legend(x=-1.05, y=1.05, legend=c("significant in both"), pch=20, col=c("red"), bty="n")
          # legend(x=-1.08, y=1, legend=c("beta = 1.089;\np < 2.2e-16;\nR^2 = 0.93"), lty=1, lwd=2, col=c("red"), bty="n")
          # dev.off()
          # 
data<-results
data[3:8]<-sapply(data[3:8], as.numeric)

data$sig_both<-mapvalues(data$sig_both, from=c(0,1), to=c("No", "Yes"))
data$sig_both<-factor(data$sig_both, levels = c("Yes", "No"))


#https://stackoverflow.com/questions/28522831/r-scatter-plot-with-error-bars-in-both-directions
#Dummy error


tiff("4B.scatter.corr.tiff", 
     width = 5.5, height = 4, 
     units = 'in', res = 600)

ggplot(data, aes(x=Phe_cor, y=Gen_cor))+
    geom_hline(yintercept = 0,size=0.7, color="grey") +
    geom_vline(xintercept = 0, size=0.7, color="grey")+
    geom_abline(linetype='dotted', alpha=0.5,size=0.5)+
    stat_poly_line(se=FALSE, size=0.6, color="black") +
    stat_poly_eq(aes(label = paste(after_stat(eq.label), "*\", \"*", 
                                   after_stat(rr.label), "*\",  \"*",
                                after_stat(p.value.label), "*\".\"",
                                sep = "")),
                        size = 2.9)+
    geom_errorbar(aes(ymin = Gen_cor-Gen_SE, ymax = Gen_cor+Gen_SE), 
                  width = 0.02, size=0.2, alpha=0.4) + 
    # geom_errorbarh(aes(xmin = Phe_cor-Phe_SE, xmax = Phe_cor+Phe_SE), 
    #                height=0.02, size=0.2)+
    geom_point(aes(shape=sig_both), size=2, alpha=0.95)+
    theme_bw()+
    xlab("Phenotypic Correlation")+
    ylab("Genetic Correlation")+
    scale_x_continuous(limits = c(-1,1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1,1),expand = c(0, 0))+
    labs(shape="Both significant")+
        scale_shape_manual(values = c(16,7))


dev.off()    
