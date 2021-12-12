#Filter LD Link results

library(LDlinkR)
suppressMessages(library(dplyr))
setwd("/scratch/mf91122/T-WES-GWAS/LDlinkR")

phenotypes<-c("w3FA_NMR","w3FA_NMR_TFAP",
                "w6FA_NMR", "w6FA_NMR_TFAP",
                "w6_w3_ratio_NMR",
                "DHA_NMR","DHA_NMR_TFAP",
                "LA_NMR","LA_NMR_TFAP",
                "PUFA_NMR","PUFA_NMR_TFAP",
                "MUFA_NMR", "MUFA_NMR_TFAP",
                "PUFA_MUFA_ratio_NMR")

#read file from replicatedvars.R output


#combine = combine the LDlinkR output files
#x = count the ones who grepped fatty acid
#y = only the ones who didn't have fatty acid
list.combine<-list()
list.x<-list()
list.y<-list()

for (p in 1:length(phenotypes)){
#p=1

filelist<-list.files("/scratch/mf91122/PUFA-GWAS/LDlinkR/12112021", full.names=T)
filephen<-filelist[grepl(paste(phenotypes[p], "-", sep=""), filelist)]

list.data<-list()
for (i in 1:length(filephen)){
	list.data[[i]]<-read.delim(filephen[i], header=T, stringsAsFactors=F)
}

list.combine[[p]]<-do.call("rbind", list.data)
list.combine[[p]]<-as_tibble(list.combine[[p]])

list.combine[[p]]$fatty_acid<-0
list.combine[[p]]$fatty_acid[grepl("fatty acid", 
	list.combine[[p]]$GWAS_Trait, ignore.case = TRUE)]<-1


list.x[[p]]<-list.combine[[p]]%>%group_by(Query)%>%summarise(fatty_acid = sum(fatty_acid))
list.y[[p]]<-list.x[[p]][list.x[[p]][,2]==0,]
} #end loop p

total_number_unique_snps<-0
for (p in 1:length(phenotypes)){
total_number_unique_snps[p]<-length(unique(list.combine[[p]]$Query))
}

count_nogrep_fatty_acid<-0
for (p in 1:length(phenotypes)){
count_nogrep_fatty_acid[p]<-nrow(list.y[[p]])
}

nogrep_fatty_acid<-data.frame(phenotypes, count_nogrep_fatty_acid,total_number_unique_snps)
nogrep_fatty_acid$percentnovel<-(nogrep_fatty_acid$count_nogrep_fatty_acid/nogrep_fatty_acid$total_number_unique_snps)*100

###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#Grep fatty acid, lipid, oily fish-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

phenoterms<-c("w3FA_NMR","w3FA_NMR_TFAP",
                "w6FA_NMR", "w6FA_NMR_TFAP",
                "w6_w3_ratio_NMR",
                "DHA_NMR","DHA_NMR_TFAP",
                "LA_NMR","LA_NMR_TFAP",
                "PUFA_NMR","PUFA_NMR_TFAP",
                "MUFA_NMR", "MUFA_NMR_TFAP",
                "PUFA_MUFA_ratio_NMR")


terms<-c("fatty acid", "lipid", "oily fish")

list.combine<-list()
list.x<-list()
list.y<-list()

for (p in 1:length(phenotypes)){
#p=1

filelist<-list.files("/scratch/mf91122/PUFA-GWAS/LDlinkR/12112021", full.names=T)
filephen<-filelist[grepl(paste(phenotypes[p], "-", sep=""), filelist)]

list.data<-list()
for (i in 1:length(filephen)){
        list.data[[i]]<-read.delim(filephen[i], header=T, stringsAsFactors=F)
}

list.combine[[p]]<-do.call("rbind", list.data)
list.combine[[p]]<-as_tibble(list.combine[[p]])

list.combine[[p]]$fatty_acid<-0
list.combine[[p]]$fatty_acid[grepl(paste(terms, collapse="|"),
        list.combine[[p]]$GWAS_Trait, ignore.case = TRUE)]<-1


list.x[[p]]<-list.combine[[p]]%>%group_by(Query)%>%summarise(fatty_acid = sum(fatty_acid))
list.y[[p]]<-list.x[[p]][list.x[[p]][,2]==0,]
} #end loop p


total_number_unique_snps<-0
for (p in 1:length(phenotypes)){
total_number_unique_snps[p]<-length(unique(list.combine[[p]]$Query))
}

count_nogrep<-0
for (p in 1:length(phenotypes)){
count_nogrep[p]<-nrow(list.y[[p]])
}

nogrep<-data.frame(phenotypes, count_nogrep,total_number_unique_snps)
nogrep$percentnovel<-(nogrep$count_nogrep/nogrep$total_number_unique_snps)*100
