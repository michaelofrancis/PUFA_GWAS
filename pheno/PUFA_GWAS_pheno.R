suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(RNOmni))

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
#Phenotype Pre-processing:  
#regress the raw phenotype on covariates; 
#the resulting residuals are transformed with rank-based inverse normal transformation (McCaw et al., 2020);
#
#Covariates 
#Model 1: sex, age, age squared, genotyping array, and assessment center indicators (sites of recruitment); 
#Model 2: Model 1 + BMI, lipid medication, socioeconomic status measured by Townsend deprivation index;  


pheno<-bd%>%select(f.eid, f.21003.0.0, f.31.0.0, 
                   f.189.0.0,
                   f.54.0.0, f.22000.0.0
                    )

colnames(pheno)<-c("IID", "Age", "Sex",  
                   "Townsend",
                   "Assessment_center", "Geno_batch"
                    )

pheno2<-bd_join4%>%select(f.eid,
                          f.21001.0.0, f.20116.0.0, f.20160.0.0,
                          f.6177.0.0,f.6153.0.0,
                          f.23651.0.0, f.23652.0.0, 
                          f.23653.0.0, f.23654.0.0, f.23655.0.0  
                          )
#TFAP = PUFA to total fatty acid percentage

colnames(pheno2)<-c("IID",
                    "BMI", "SmokeStatus","Ever_smoked",                    
                    "lipid_med", "lipid_med_plushormones",
                    "MeasurementQualityFlag", "HighLactate",
                    "HighPyruvate", "LowGlucose","LowProtein"
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
new$HighPyruvate[is.na(new$HighPyruvate)]<-0
new$HighPyruvate[is.na(new$HighLactate)]<-0
new$HighPyruvate[is.na(new$LowGlucose)]<-0
new$HighPyruvate[is.na(new$LowProtein)]<-0



#lipid medicine columns
#f.6177 == 1 Cholesterol lowering medication
new$lipid_med[new$lipid_med!=1]<-0
#This column has NA's but I think it's okay to keep.

#Statins
statincols<-c(sprintf("f.20003.0.%s", 0:47))
statincodes<-c(1141146234,1141192414,1140910632,1140888594,1140864592,
	1141146138,1140861970,1140888648,1141192410,
	1141188146,1140861958,1140881748,1141200040)

manyColsToDummy(statincodes, bd_join4[,statincols], "statinoutput")
statinoutput$statins<-rowSums(statinoutput) 
statinoutput$statins[statinoutput$statins>1]<-1

statinoutput$IID<-bd_join4$f.eid

statinoutput<-statinoutput%>%select(IID, statins)

new<-left_join(new, statinoutput, by="IID")


#Clean up columns



#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Regress phenotypes, extract residuals-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#ref: https://github.com/lindgrengroup/fatdistnGWAS/commit/35ca6d65dbdb1fc3cf417432dc0575805ce80878
#Copied from CCC, not edited yet.VVVV
#regress the raw phenotype on covariates; 
#the resulting residuals are transformed with 
#rank-based inverse normal transformation (McCaw et al., 2020);


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

new1<-left_join(new, resdat_inv, by="IID")

participants1<-new1%>%select(IID)
participants1$FID<-participants1$IID
participants1<-participants1%>%select(FID, IID)


###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###MODEL 2-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

#Model 2: Model 1 + BMI, lipid medication, SES via Townsend score; 

###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
#Calculate residuals+++++++++++++++++++++++++++++
###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

resdat2<-matrix(NA,nrow(new),length(phenotypes))
colnames(resdat2)<-phenotypes

#LIPID MED IS SCREWING THIS UP

for (p in 1:length(phenotypes)){
    assign(paste("lm", phenotypes[p], sep="_"),
           lm(new[[phenotypes[p]]] ~ Age + Age2 +
            center1 + center2 + center3 + center4 +
                center5 + center6 + center7 + center8 +
                center9 + center10 + center11 + center12 +
                center13 + center14 + center15 + center16 +
                center17 + center18 + center20 +
                center21 + Sex + Geno_batch + 
		BMI + Townsend + statins, data=new,
        na.action=na.exclude))
    lmname<-paste("lm", phenotypes[p], sep="_")
    lmobj<-get(lmname)
    resdat2[,p]<-resid(lmobj)
}

resdat2<-as_tibble(resdat2)
colnames(resdat2)<-paste(phenotypes, "res", sep="_")


#the resulting residuals are transformed with rank-based inverse normal transformation
###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
#Rank INT on residuals+++++++++++++++++++++++++++
###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

resdat_inv2<-matrix(NA,nrow(new),length(phenotypes))
colnames(resdat_inv2)<-paste(phenotypes, "resinv", sep="_")

for (i in 1:length(phenotypes)){
    resdat_inv2[,i] <- inversenormal(resdat2[,i])
}

resdat_inv2<-as_tibble(as.data.frame(resdat_inv2))
resdat_inv2$IID<-new$IID

new2<-left_join(new, resdat_inv2, by="IID")

participants2<-new2%>%select(IID)
participants2$FID<-participants2$IID
participants2<-participants2%>%select(FID, IID)



###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###WRITE OUTPUT=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

outdir="/scratch/mf91122/PUFA-GWAS/pheno"

#Model 1
write.table(participants1, 
	paste(outdir, "/PUFA_GWAS_phenoQC_IDS_M1.txt",sep=""), 
	row.names=FALSE, quote=FALSE)

write.table(new1, 
	paste(outdir, "/PUFA_GWAS_pheno_M1.txt", sep=""),
	row.names=FALSE, quote=FALSE)

#Model 2
write.table(participants2,
        paste(outdir, "/PUFA_GWAS_phenoQC_IDS_M2.txt",sep=""),
        row.names=FALSE, quote=FALSE)

write.table(new2,
        paste(outdir, "/PUFA_GWAS_pheno_M2.txt", sep=""),
        row.names=FALSE, quote=FALSE)
