#Reference for displaying results
#https://doi.org/10.1038/s41593-020-0643-5


suppressMessages(library(tidyverse))

##04.SPrediXcan


phenotypes<-c("w3FA_NMR_resinv", "w3FA_NMR_TFAP_resinv", "w6FA_NMR_resinv",
              "w6FA_NMR_TFAP_resinv", "w6_w3_ratio_NMR_resinv",
              "DHA_NMR_resinv", "DHA_NMR_TFAP_resinv", "LA_NMR_resinv",
              "LA_NMR_TFAP_resinv", "PUFA_NMR_resinv", "PUFA_NMR_TFAP_resinv",
              "MUFA_NMR_resinv", "MUFA_NMR_TFAP_resinv", "PUFA_MUFA_ratio_NMR_resinv")

phenotypenames<-c("Omega-3", "Omega-3 (TFAP)", "Omega-6" , "Omega-6 (TFAP)",
                  "Omega-6:Omega-3 ratio", "DHA", "DHA (TFAP)", "LA", "LA (TFAP)",
                  "PUFAs", "PUFAs (TFAP)", "MUFAs", "MUFAs (TFAP)", "PUFAs:MUFAs ratio")



dir4="/scratch/mf91122/PUFA-GWAS/SMultiXcan/huifang2/04.SPrediXcan_tissues"
files<-list.files(dir4, full.names=T)

phenoresults<-list()
#one pheno at a time
for (p in 1:length(phenotypes)){

set<-files[grepl(phenotypes[p],files)]
tis<-str_split(set, "__", simplify=T)[,3]
tis<-str_split(tis, ".csv", simplify=T)[,1]

part4<-list()
for (i in 1:length(set)){

	part4[[i]]<-as_tibble(read.csv(set[i], header=T))

}
names(part4)<-tis

#count number of genes tested to get bf correction number---
sum<-0
for (i in 1:length(set)){
sum=sum+length(unique(part4[[i]]$gene_name))
}
print(paste(phenotypes[p], sum))
#[1] 602072
Pcutoff=0.05/sum
#------

part4sig<-list()
for (i in 1:length(set)){
	part4sig[[i]]<-part4[[i]][part4[[i]]$pvalue<Pcutoff,]
	part4sig[[i]]$tissue<-tis[i]
}
part4sig<-do.call(rbind, part4sig)
part4sig$Phenotype<-phenotypenames[p]

phenoresults[[p]]<-part4sig

}#end phenotypes loop

phenoresults<-do.call(rbind, phenoresults)

phenoresults<-phenoresults%>%select(Phenotype, tissue, gene, gene_name, zscore, pvalue, var_g, n_snps_used, n_snps_in_model)
write.csv(phenoresults, paste(dir4, "/significant.csv", sep=""), row.names=F, quote=F)

##05.SMultiXcan----------------------------------------------------------------

dir5="/scratch/mf91122/PUFA-GWAS/SMultiXcan/huifang2/05.SMultiXcan"
files<-list.files(dir5, full.names=T)
filenames<-list.files(dir5)

#Read in tables
part5<-list()
for (i in 1:length(files)){
part5[[i]]<-as_tibble(read.table(files[i], header=T))
}

#Get names for tables
phen<-str_split(filenames , pattern="_resinv", simplify=T)[,1]
phen<-str_split(phen, "05.", simplify=T)[,2]
names(part5)<-filenames 

length(unique(part5[[1]]$gene_name))
#[1] 21850

#21850 unique genes to Bonferroni correction

Pcutoff<-0.05/21850
#[1] 2.28833e-06

part5sig<-list()
for (i in 1:length(files)){
	part5sig[[i]]<-part5[[i]][part5[[i]]$pvalue <Pcutoff,]

	part5sig[[i]]$Phenotype<-phen[i] 

	part5sig[[i]]<-part5sig[[i]]%>%select("Phenotype", "gene","gene_name","pvalue","n","n_indep",
		"p_i_best","t_i_best","p_i_worst","t_i_worst")
	
}

part5out<-do.call(rbind, part5sig)

write.csv(part5out, paste(dir5, "/significant.SMultiXcan.csv", sep=""), row.names=F, quote=F)
