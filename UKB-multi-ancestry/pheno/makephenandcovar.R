suppressMessages(library(tidyverse))

phenodir="/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/pheno"
dir.create(paste(phenodir, "/phen", sep=""), showWarnings=F)

pop=c("AFR", "CSA", "EAS")

phenotypes<-c("w3FA_NMR","w3FA_NMR_TFAP",
                "w6FA_NMR", "w6FA_NMR_TFAP",
                "w6_w3_ratio_NMR",
                "DHA_NMR","DHA_NMR_TFAP",
                "LA_NMR","LA_NMR_TFAP",
                "PUFA_NMR","PUFA_NMR_TFAP",
                "MUFA_NMR", "MUFA_NMR_TFAP",
                "PUFA_MUFA_ratio_NMR")

phenotypes<-paste(phenotypes, "_resinv", sep="")
models=c(1,2)

#for testing script: p=1; j=1;m=1

for (p in 1:length(pop)){
		for (j in 1:length(phenotypes)){
			for (m in 1:length(models)){
		filename<-paste("PUFA_GWAS_pheno_M", m,".",pop[p], ".txt.pc.txt", sep="")

		#Read phenotype file
		dat<-read.table(paste(phenodir, "/", filename,sep=""), header=T)
		dat<-as_tibble(dat)
		
		#Select columns for phen file
		colstoselect<-c("IID", "IID", phenotypes[j])

		phentable<-dat[colstoselect]

		#Write phen file
		write.table(phentable, paste(phenodir, "/phen/", filename,".",
			phenotypes[j],".phen", sep=""), 
			col.names=F, row.names=F, quote=F)

		#Select columns for quantitative covariates file 
		
		#Model 1		
		if (m==1){
		quantcovarcolsM1=c("IID", "IID", "Age", "Age2", "Sex", "Geno_batch", 
				colnames(dat)[grep("center", colnames(dat))],
				sprintf("PC%s", 1:20))
		qcovartableM1<-dat[quantcovarcolsM1]
		qcovartableM1<-qcovartableM1[vapply(qcovartableM1,
                        function(x) length(unique(x)) > 1, logical(1L))]	
		#Write qcovar file
                write.table(qcovartableM1,
                        paste(phenodir, "/phen/", filename,".qcovar", sep=""),
                        col.names=F, row.names=F, quote=F)
		}

		if (m==2){

		#Model 2
		
		quantcovarcolsM2=c("IID", "IID", "Age", "Age2", "Sex", "Geno_batch",
                                colnames(dat)[grep("center", colnames(dat))],
                                sprintf("PC%s", 1:20),
				"statins","Townsend","BMI")		
		qcovartableM2<-dat[quantcovarcolsM2]
		#Remove columns that have less than one unique value.
		qcovartableM2<-qcovartableM2[vapply(qcovartableM2, 
                        function(x) length(unique(x)) > 1, logical(1L))]
                #Write qcovar file
		write.table(qcovartableM2, 
                        paste(phenodir, "/phen/", filename,".qcovar", sep=""),
                        col.names=F, row.names=F, quote=F)
		}

		}
	}#end phenotype loop
}#end pop loop
