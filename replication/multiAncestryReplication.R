library(plyr)
library(tidyverse)


phenotypes<-c("w3FA_NMR_resinv", "w3FA_NMR_TFAP_resinv", "w6FA_NMR_resinv",
                "w6FA_NMR_TFAP_resinv", "w6_w3_ratio_NMR_resinv",
                "DHA_NMR_resinv", "DHA_NMR_TFAP_resinv", "LA_NMR_resinv", 
                "LA_NMR_TFAP_resinv", "PUFA_NMR_resinv", "PUFA_NMR_TFAP_resinv", 
                "MUFA_NMR_resinv", "MUFA_NMR_TFAP_resinv", "PUFA_MUFA_ratio_NMR_resinv")

phenotypesb<-str_split(phenotypes, "_resinv", simplify=T)[,1]

phenotypesb<-c("Omega-3 Fatty Acids", "Omega-3 Fatty Acids to Total Fatty Acids percentage",
		"Omega-6 Fatty Acids", "Omega-6 Fatty Acids to Total Fatty Acids percentage",
		"Omega-6 Fatty Acids to Omega-3 Fatty Acids ratio", 
		"Docosahexaenoic Acid", "Docosahexaenoic Acid to Total Fatty Acids percentage",
		"Linoleic Acid", "Linoleic Acid to Total Fatty Acids percentage",
		"Polyunsaturated Fatty Acids","Polyunsaturated Fatty Acids to Total Fatty Acids percentage",
		"Monounsaturated Fatty Acids", "Monounsaturated Fatty Acids to Total Fatty Acids percentage",
		"Polyunsaturated Fatty Acids to Monounsaturated Fatty Acids ratio")

phelabels<-c("Omega-3", "Omega-3_TFAP", "Omega-6", "Omega-6_TFAP", 
	"Omega-6_Omega-3_ratio", "DHA","DHA_TFAP", "LA", "LA_TFAP", 
	"PUFAs", "PUFAs_TFAP", "MUFAs", "MUFAs_TFAP", "PUFA_MUFA_ratio")

dirM1<-"/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/M1"
dirM2<-"/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/M2"
dirmulti<-"/scratch/mf91122/PUFA-GWAS/sumstats-final/multi-ancestry"
outdir="/scratch/mf91122/PUFA-GWAS/sumstats-final/UKB"



cols2<-c(
    "NMR phenotypes",

    "Variants in UKB-EUR (unmunged)",

    "Significant variants UKB-EUR",

    "Variants in UKB-AFR",

    "Associations at P<0.05 UKB-AFR",

    "Overlapping variants UKB-EUR and UKB-AFR",

    "Significant variants UKB-EUR replicated in UKB-AFR",

   "Variants in UKB-CSA",

    "Associations at P<0.05 UKB-CSA",

    "Overlapping variants UKB-EUR and UKB-CSA",

    "Significant variants UKB-EUR replicated in UKB-CSA",

    "Variants in UKB-EAS",

    "Associations at P<0.05 UKB-EAS",

    "Overlapping variants UKB-EUR and UKB-EAS",

    "Significant variants UKB-EUR replicated in UKB-EAS"
)


tab<-as.data.frame(matrix(nrow=14,ncol=length(cols2)))
colnames(tab)<-cols2



for (p in 1:length(phenotypes)){


	UKB<-as_tibble(read_delim(
        paste(dirM2, "/", phenotypes[p], "/BOLT1-statsFile-BgenSnps-m2", sep="")
                                        ))

multi<-as_tibble(read_csv(
paste(dirmulti, "/UKB-multi-ancestry-", phelabels[p],".csv", sep="")
				))


AFR<-multi[!is.na(multi$p_AFR),][c(1:3, 4:9)]
CSA<-multi[!is.na(multi$p_CSA),][c(1:3, 10:15)]
EAS<-multi[!is.na(multi$p_EAS),][c(1:3, 16:21)]

cols<-c("CHR","SNP","BP", "A1", "A2", "Freq", "b","se","p")
colnames(AFR)<-cols
colnames(CSA)<-cols
colnames(EAS)<-cols


A<-inner_join(UKB, AFR, by=c("SNP", "CHR", "BP"))
B<-inner_join(UKB, CSA, by=c("SNP", "CHR", "BP"))
C<-inner_join(UKB, EAS, by=c("SNP", "CHR", "BP"))

tab[p,]<-c(
    phenotypesb[p],

    nrow(UKB),

    nrow(UKB%>%filter(P_BOLT_LMM<1.678e-08)),

    nrow(AFR),

    nrow(AFR%>%filter(p<0.05)),

    nrow(A),

    nrow(A%>%filter(P_BOLT_LMM<1.678e-08 & p<0.05)),

    nrow(CSA),

    nrow(CSA%>%filter(p<0.05)),

    nrow(B),

    nrow(B%>%filter(P_BOLT_LMM<1.678e-08 & p<0.05)),

    nrow(EAS),

    nrow(EAS%>%filter(p<0.05)),

    nrow(C),

    nrow(C%>%filter(P_BOLT_LMM<1.678e-08 & p<0.05))

)

}#end pheno loop
