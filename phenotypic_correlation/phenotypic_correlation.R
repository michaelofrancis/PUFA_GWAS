suppressMessages(library(tidyverse))
library("Hmisc")
library(reshape2)


phenofile<-"/scratch/mf91122/PUFA-GWAS/pheno/PUFA_GWAS_pheno_M2.txt.pc.txt"

#read phenotypic data
df <- as_tibble(read.table(phenofile, header=T))

df<-df[c(2,17:30)]

#build a empty tibble to store result
corr_results <- tibble(p1="phenotype1",p2="phenotype2",rp=0)

#For Hclust------
colnames(df)<-c("IID", "n3FA","n3FA.Total.FA","n6FA","n6FA.Total.FA","n6.n3",
	"DHA","DHA.Total.FA","LA","LA.Total.FA","PUFA","PUFA.Total.FA","MUFA", 
	"MUFA.Total.FA","PUFA.MUFA")

ordercols<-c("IID", "n6.n3", "PUFA.Total.FA","PUFA.MUFA", "n6FA.Total.FA",
		"LA.Total.FA","n3FA", "DHA", "n3FA.Total.FA",
		"DHA.Total.FA", "PUFA",  "n6FA", "LA",   
		"MUFA", "MUFA.Total.FA")

df<-df[ordercols]
#----------------

#do correlation
for (i in 2:ncol(df)){
  for (j in 2:ncol(df)){
    corr <- cor(df[,i],df[,j],method="pearson")
    corr_results <- add_row(corr_results,p1=colnames(df[,i]),p2=colnames(df[,j]),rp=corr)
    }
}


#transfer format to matrix
pc_matrix <- corr_results[-1,]
pc_matrix$p1 <- factor(pc_matrix$p1, levels=unique(pc_matrix$p1))
pc_matrix$p2 <- factor(pc_matrix$p2, levels=unique(pc_matrix$p2))
pc_matrix <- acast(pc_matrix,p1~p2)
pc_matrix <- as.matrix(pc_matrix)

#write table
write.table(pc_matrix,"/scratch/mf91122/PUFA-GWAS/phenotypic-corr/pc_matrix.hclust.txt",
		row.names=TRUE,quote=FALSE)



#Using rcorr to get p val matrix

df2<-as.data.frame(df)[-1]

out2<-rcorr(as.matrix(df2), type ="pearson")


print.rcorr <- function(x, ..., digits = getOption("Hmisc.rcorr.digits", 17))
{
  print(round(x$r,2))
  n <- x$n
  if(all(n == n[1,1]))
    cat("\nn=", n[1,1], "\n\n")
  else {
    cat("\nn\n")
    print(n)
  }

  cat("\nP\n")
  P <- x$P
  P <- ifelse(P < .0001, 0, P)
  p <- format(round(P, digits))
  p[is.na(P)] <- ""
  print(p, quote=FALSE)
  invisible()
}
all.neg <- function(x) -1*abs(x)


lapply(out2, log10)
