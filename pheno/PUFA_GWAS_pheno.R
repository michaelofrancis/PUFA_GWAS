library(plyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(RNOmni)

setwd("/scratch/mf91122/UKB-pheno")

source("manyColsToDummy.R")

withdrawn<-read.csv("w48818_20210809.csv", header = FALSE)

QCids<-read.table("/work/kylab/mike/PUFA-GWAS/pheno/bd_QC-keep.txt",header=TRUE)

###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Load data=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#Load UK Biobank datasets-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
source('/work/kylab/mike/PUFA-GWAS/pheno/load_UKBphenotables.R') #20 min

#Phenotypes  ------------------------------------------------------------------------------
#Five Measured Phenotypes: Field ID 23444, Omega-3 Fatty Acids; 
#Field ID 23445, Omega-6 Fatty Acids; 
#Field ID 23446, Polyunsaturated Fatty Acids; 
#Field ID 23449, Linoleic Acid; 
#Field ID 23450, Docosahexaenoic Acid; 
#Calculated Phenotypes: Omega-6 to Omega-3 ratio; 
#Phenotype Pre-processing:  
#regress the raw phenotype on covariates; 
#the resulting residuals are transformed with rank-based inverse normal transformation (McCaw et al., 2020);
#
#Covariates 
#Model 1: sex, age, age squared, genotyping array, and assessment center indicators (sites of recruitment); 
#Model 2: Model 1 + BMI, lipid medication, socioeconomic status measured by Townsend deprivation index;  


pheno<-bd%>%select(f.eid, f.21003.0.0, f.31.0.0, 
                   f.21000.0.0, f.189.0.0,
                   f.54.0.0, f.22000.0.0
                    )

colnames(pheno)<-c("IID", "Age", "Sex",  
                   "Race", "Townsend",
                   "Assessment_center", "Geno_batch"
                    )

pheno2<-bd_join4%>%select(f.eid,
                          f.21001.0.0, f.20116.0.0, f.20160.0.0,
                          f.6177.0.0,f.6153.0.0,
                          f.23651.0.0, f.23652.0.0, 
                          f.23653.0.0, f.23654.0.0, f.23655.0.0,  
                          f.3166.0.0, f.74.0.0
                          )
#TFAP = PUFA to total fatty acid percentage

colnames(pheno2)<-c("IID",
                    "BMI", "SmokeStatus","Ever_smoked",                    
                    "lipid_med", "lipid_med_plushormones",
                    "MeasurementQualityFlag", "HighLactate",
                    "HighPyruvate", "LowGlucose","LowProtein",
                    "blood_draw_time","Fasting_time"
                    )

new<-left_join(pheno, pheno2, by="IID")
new<-as_tibble(new)

#NMR Cols
#TFAP = PUFA to total fatty acid percentage

NMRcols<-bd8%>%select(f.eid,
		f.23444.0.0,f.23451.0.0,
		f.23445.0.0,f.23452.0.0, 
		f.23459.0.0,
		f.23450.0.0,f.23457.0.0,
		f.23449.0.0,f.23456.0.0,
		f.23446.0.0,f.23453.0.0,
		f.23447.0.0,f.23454.0.0,
		f.23458.0.0)

colnames(NMRcols)<-c("IID",
		"w3FA_NMR","w3FA_NMR_TFAP",
		"w6FA_NMR", "w6FA_NMR_TFAP",
		"w6_w3_ratio_NMR",
		"DHA_NMR","DHA_NMR_TFAP",
		"LA_NMR","LA_NMR_TFAP",
		"PUFA_NMR","PUFA_NMR_TFAP",
		"MUFA_NMR", "MUFA_NMR_TFAP",
		"PUFA_MUFA_ratio_NMR")

NMRcols<-as_tibble(NMRcols)

NMRQCcols<-bd8%>%select(f.eid,
		f.23744.0.0,f.23751.0.0,
		f.23745.0.0,f.23752.0.0,
		f.23759.0.0,
		f.23750.0.0,f.23757.0.0,
		f.23749.0.0,f.23756.0.0,
		f.23746.0.0,f.23753.0.0,
		f.23747.0.0, f.23754.0.0,
		f.23758.0.0)

colnames(NMRQCcols)<-c("IID",
		"w3FA_NMR_QCflag","w3FA_NMR_TFAP_QCflag",
		"w6FA_NMR_QCflag","w6FA_NMR_TFAP_QCflag",
		"w6_w3_ratio_NMR_QCflag",
		"DHA_NMR_QCflag","DHA_NMR_TFAP_QCflag",
		"LA_NMR_QCflag", "LA_NMR_TFAP_QCflag",
		"PUFA_NMR_QCflag","PUFA_NMR_TFAP_QCflag",
		"MUFA_NMR_QCflag","MUFA_NMR_TFAP_QCflag",
		"PUFA_MUFA_ratio_NMR_QCflag")

NMRQCcols<-as_tibble(NMRQCcols)

#check if NMRcols and NMRQCcols column names are aligned:
#for(i in 1:length(colnames(NMRQCcols))){print(identical(colnames(NMRQCcols)[i],paste(colnames(NMRcols)[i],"_QCflag", sep="")))}

removeid<-NMRQCcols$IID[!is.na(NMRQCcols$w3FA_NMR_QCflag)] #[1] 4178146 Just one guy across all cols

#Merge NMR cols with main table

new<-left_join(new, NMRcols, by="IID")

#Remove this one guy because of these NMR QC cols
new<-new[!(new$IID %in% removeid), ]


#Remove withdrawn participants------------------------------------
new<-new[!(new$IID %in% withdrawn$V1), ]

#QC participants via output of UKB_participantQC.R----------------

new<-new[!(new$IID %in% QCids$V1),]

#Age squared----------------------------
new$Age2<-new$Age^2

#Make dummy 0/1 cols for each assessment center----------------------
#table(pheno$Assessment_center)
centers<-unique(new$Assessment_center)
centercols<-paste("center", 1:22, sep="")
new[centercols]<-0

for (i in 1:length(centers)){
    new[new$Assessment_center==centers[i],][centercols[i]]<-1
}

new<-new%>%select(-Assessment_center)
new

#Genotype batch
new$Geno_batch1<-0
new$Geno_batch1[new$Geno_batch>0]<-1
#sum(pheno$Geno_batch1) #[1] 438313
new$Geno_batch<-new$Geno_batch1
new<-new%>%select(-Geno_batch1)
#table(new$Geno_batch) #it worked


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#Switch sex values to numeric
new$Sex<-mapvalues(as.character(new$Sex), 
                     c("Male", "Female"), c(0,1))


#NMR QC columns-----------------------------------------------------

#First: are NA's consistent across the 5 NMR columns?
#NMRcols<-new%>%select(w3FA_NMR,w6FA_NMR,PUFA_NMR,LA_NMR, DHA_NMR)
#for (i in 1:5){print(sum(is.na(NMRcols[i])))} #yes they are

new<-new[!is.na(new$w3FA_NMR),] #only NMR participants

#Overall QC flags
new<-new[is.na(new$MeasurementQualityFlag),]



#lipid medicine columns
#f.6177 == 1 Cholesterol lowering medication
new$lipid_med[new$lipid_med!=1]<-0
#This column has NA's but I think it's okay to keep.


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Regress phenotypes, extract residuals-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#ref: https://github.com/lindgrengroup/fatdistnGWAS/commit/35ca6d65dbdb1fc3cf417432dc0575805ce80878
#Copied from CCC, not edited yet.VVVV

#CENTER 19 HAS NO 1's NO PEOPLE THERE REMOVE COLUMN
new<-new%>%select(-center19)
#CENTER22 do not include because remove one dummy variable.


phenotypes<-c("w3FA_NMR","w3FA_NMR_TFAP",
                "w6FA_NMR", "w6FA_NMR_TFAP",
                "w6_w3_ratio_NMR",
                "DHA_NMR","DHA_NMR_TFAP",
                "LA_NMR","LA_NMR_TFAP",
                "PUFA_NMR","PUFA_NMR_TFAP",
                "MUFA_NMR", "MUFA_NMR_TFAP",
                "PUFA_MUFA_ratio_NMR")
###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###MODEL 1-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
#Calculate residuals+++++++++++++++++++++++++++++
###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

resdat<-matrix(NA,nrow(new),length(phenotypes))
colnames(resdat)<-phenotypes

for (p in 1:length(phenotypes)){
    assign(paste("lm", phenotypes[p], sep="_"), 
           lm(new[[phenotypes[p]]] ~ Age + Age2 + 
            center1 + center2 + center3 + center4 + 
                center5 + center6 + center7 + center8 + 
                center9 + center10 + center11 + center12 + 
                center13 + center14 + center15 + center16 + 
                center17 + center18 + center20 + 
                center21 + Sex + Geno_batch, data=new, 
        na.action=na.exclude))
    lmname<-paste("lm", phenotypes[p], sep="_")
    lmobj<-get(lmname)
    resdat[,p]<-resid(lmobj)
}

resdat<-as_tibble(resdat)
colnames(resdat)<-paste(phenotypes, "res", sep="_")


#the resulting residuals are transformed with rank-based inverse normal transformation
###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
#Rank INT on residuals+++++++++++++++++++++++++++
###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

resdat_inv<-matrix(NA,nrow(new),length(phenotypes))
colnames(resdat_inv)<-paste(phenotypes, "resinv", sep="_")

# R code for inverse normalized transformation 
inversenormal <- function(x) { 
    # inverse normal if you have missing data 
    return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))) 
} 

for (i in 1:length(phenotypes)){
    resdat_inv[,i] <- inversenormal(resdat[,i])
}

resdat_inv<-as_tibble(as.data.frame(resdat_inv))
resdat_inv$IID<-new$IID

new<-left_join(new, resdat_inv, by="IID")

participants<-new%>%select(IID)
participants$FID<-participants$IID
participants<-participants%>%select(FID, IID)

outdir="/scratch/mf91122/PUFA-GWAS/pheno"

write.table(participants, 
	paste(outdir, "/PUFA_GWAS_phenoQC_IDS.txt",sep=""), 
	row.names=FALSE, quote=FALSE)

write.table(new, 
	paste(outdir, "PUFA_GWAS_pheno.txt", sep=""),
	row.names=FALSE, quote=FALSE)
