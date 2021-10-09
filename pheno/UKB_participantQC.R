### MF: This is a script for generating lists of UK Biobank participants which meet QC/filtering criteria.

pan<-read_tsv("../all_pops_non_eur_pruned_within_pop_pc_covs.tsv")
pan<-as_tibble(pan)
pan$s<-as.integer(pan$s)
table(pan$pop, useNA = "always")

bridge<-read.table("../T-WES-GWAS/ukb48818bridge31063.txt")
bridge<-as_tibble(bridge)
colnames(bridge)<-c("IID", "panID")

pan2<-pan%>%select(s, pop)%>%
    left_join(bridge, by=c("s"="panID"))

pan2

##===============================PART A=================================================
#Generate a list of participants who pass the following QC criteria:
#1.Genetic ethnicity = Caucasian VIA PAN UKBB
#2. Not an outlier for heterogeneity and missing genotype rate (poor quality genotype)
#3. No Sex chromosome aneuploidy
#4. Self-reported sex matches genetic sex

bd_QC<- bd %>% select(f.eid, f.31.0.0, f.22001.0.0, f.21000.0.0,
                      f.22027.0.0, f.22019.0.0, 
                      f.22021.0.0)

colnames(bd_QC)<-c("IID", "Sex", "Genetic_Sex", "Race", 
                   "Outliers_for_het_or_missing", "SexchrAneuploidy",
                   "Genetic_kinship")

#Join UKB cols with with Pan UKBB
bd_QC<-as_tibble(bd_QC) #502,527
#nrow(bd_QC) #[1] 502527
bd_QC<-bd_QC%>%inner_join(pan2, by="IID")

#Filter by Genetic ethnicity = Caucasian VIA PAN UKBB
bd_QC<-bd_QC[bd_QC$pop=="EUR",]
#nrow(bd_QC) #[1] 426881


#Test for FFQ and NMR participant overlap----------------------------
#UKB24<-read.table("UKB-took24HRparticipants-09262021.txt", header=TRUE)
#UKB24<-UKB24[UKB24$took_24HR==1,]

#test[test$IID %in% w3NMR$f.eid,] #43,754 participants with both FFQ and w3-NMR data


bd_QC<-bd_QC%>%
    filter(is.na(Outliers_for_het_or_missing) | Outliers_for_het_or_missing !="Yes") 
#nrow(bd_QC) #[1] 426433

bd_QC<-bd_QC%>%
    filter(is.na(SexchrAneuploidy) | SexchrAneuploidy != "Yes")
#nrow(bd_QC) #[1] 425854


#bd_QC<- bd_QC%>%
#    filter(is.na(Genetic_kinship) | 
#               Genetic_kinship != "Ten or more third-degree relatives identified")

#If Sex does not equal genetic sex, exclude participant
bd_QC<-bd_QC[bd_QC$Sex == bd_QC$Genetic_Sex,] 
#nrow(bd_QC) #[1] 425683


#From maximum_set_of_unrelated_individuals.MF.pl output: 
#DIDN'T DO THIS IN PUFA GWAS
#max_unrelated<-read.table("ukb48818_rel_s488282_output.dat")
#max_unrelated<-as.integer(unlist(max_unrelated))
#bd_QC<-bd_QC%>%filter(!FID %in% max_unrelated)

#Only those with NMR data
w3NMR<-bd_join4%>%select(f.eid, f.23444.0.0)
w3NMR<-w3NMR[!is.na(w3NMR$f.23444.0.0),]
bd_QC<-bd_QC[bd_QC$IID %in% w3NMR$f.eid,]



table(bd_QC$Race)

QCkeepparticipants<-bd_QC%>%select(IID)

write.table(QCkeepparticipants, file= "bd_QC-keep.txt", 
            row.names = FALSE, quote = FALSE)
