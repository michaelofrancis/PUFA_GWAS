#install.packages("gridExtra")
#install.packages("png")
#install.packages("grobblR")
install.packages("qpdf")

library(gridExtra)
library(png)
library(grobblR)
library(qpdf)

pheno<-c("w3FA_NMR_resinv", "w3FA_NMR_TFAP_resinv", "w6FA_NMR_resinv",
         "w6FA_NMR_TFAP_resinv", "w6_w3_ratio_NMR_resinv",
         "DHA_NMR_resinv", "DHA_NMR_TFAP_resinv", "LA_NMR_resinv",
         "LA_NMR_TFAP_resinv", "PUFA_NMR_resinv", "PUFA_NMR_TFAP_resinv",
         "MUFA_NMR_resinv", "MUFA_NMR_TFAP_resinv", "PUFA_MUFA_ratio_NMR_resinv")

pop<-c("AFR", "CSA", "EAS")


#Order files---------

order_files<-function(pngfiles){
    #Use:
    #fileorder<-order_files(pngfiles)
    fileorder<-NULL
    for (j in 1:length(pop)){
        for (i in 1:length(pheno)){
            found<-pngfiles[grep(
                paste(pop[j], ".*",  pheno[i], sep=""), 
                pngfiles)]
            fileorder<-c(fileorder, found)
        }
    }
   return(fileorder)
}

order_files_nopop<-function(pngfiles){
    #Use:
    #fileorder<-order_files(pngfiles)
    fileorder<-NULL
        for (i in 1:length(pheno)){
            found<-pngfiles[grep(pheno[i], pngfiles)]
            fileorder<-c(fileorder, found)
        }
    return(fileorder)
}
prder_files_nopop


# part1:UKB-multi-ancestry --------------------------------------
dir<-"/Users/mike/Documents/Research/PUFA-GWAS/UKB-multi-ancestry/plots"
pngfiles<-list.files(dir,full.names =  T)[grep(".png", list.files(dir))]

#order files
fileorder<-order_files(pngfiles)


manpng<-fileorder[grep("Manhattan",pngfiles)]
manpng
qqpng<-fileorder[grep("QQ",pngfiles)]
qqpng
page<-list()
rows=4
bord=F
numpages<-ceiling(length(manpng)/rows)
#%% = remainder
remainder<-length(manpng) %% rows

#at 
floor(length(manpng)/rows)

for (i in 1:(floor(length(manpng)/rows))){
page[[i]]<-grob_layout(
    grob_row(
        grob_col(
            border=bord, 
            manpng[i*rows - (rows-1)]), 
        grob_col(
            qqpng[i*rows - (rows-1)], 
            border=bord)),
    grob_row(
        grob_col(
            border=bord, 
            manpng[i*rows - (rows-2)]), 
        grob_col(
            qqpng[i*rows - (rows-2)], 
            border=bord)),
    grob_row(
        grob_col(
            border=bord, 
            manpng[i*rows - (rows-3)]), 
        grob_col(
            qqpng[i*rows - (rows-3)], 
            border=bord)),
    grob_row(
        grob_col(
            border=bord, 
            manpng[i*rows - (rows-4)]), 
        grob_col(
            qqpng[i*rows - (rows-4)], 
            border=bord))
    #height = 100,
    #width = 100,
    #padding = 0
) 
}
#Do the last page
i=ceiling(length(manpng)/rows)
page[[i]]<-grob_layout(
    grob_row( #1 row
        grob_col(
            border=bord, 
            manpng[i*rows - (rows-1)]), 
        grob_col(
            qqpng[i*rows - (rows-1)], 
            border=bord)),
    grob_row( #2 rows
        grob_col(
            border=bord, 
            manpng[i*rows - (rows-2)]), 
        grob_col(
            qqpng[i*rows - (rows-2)], 
            border=bord)),
    grob_row( #3 rows
        grob_col(
            border=bord, 
            NA)),
    grob_row( #4 rows
        grob_col(
            border=bord, 
            NA))
    
)


#dir.create(paste(dir, "/pdf", sep=""))
outdir<-paste(dir, "/pdf/", sep="")

grob_to_pdf(
    page,
    file_name = paste(outdir, "UKB-MultiAncestry.pdf", sep=""),
    add_page_numbers = F,
    meta_data_title = "PUFA-GWAS supplementary fig PDF"
)


# 
# # test the for loop within the function -------------------------
# pagetest<-list()
# for (i in 1:(ceiling(length(manpng)/rows))){
#     pagetest[[i]]<-grob_layout(
#         for (j in 1:4){
#             grob_row(
#                 grob_col(
#                     border=bord, 
#                     manpng[i*rows - (rows-j)]),
#                 grob_col(
#                     qqpng[i*rows - (rows-j)], 
#                     border=bord))
#         }
#     )}
# grob_to_pdf(
#     page,
#     file_name = paste(outdir, "/list_noborder_FORLOOPTEST.pdf", sep=""),
#     add_page_numbers = TRUE,
#     meta_data_title = "Test PDF"
# )  
       


# UKB-EUR discovery M1 -----------------------------------

dir<-"/Users/mike/Documents/Research/PUFA-GWAS/discovery/plot/M1"
pngfiles<-list.files(dir,full.names =  T)[grep(".png", list.files(dir))]
fileorder<-order_files_nopop(pngfiles)

manpng<-fileorder[grep("Manhattan",pngfiles)]
manpng
qqpng<-fileorder[grep("QQ",pngfiles)]
qqpng
qqpng

page<-list()
rows=4
bord=F
numpages<-ceiling(length(manpng)/rows)
#%% = remainder
remainder<-length(manpng) %% rows

#at 
floor(length(manpng)/rows)

for (i in 1:(floor(length(manpng)/rows))){
    page[[i]]<-grob_layout(
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-1)]), 
            grob_col(
                qqpng[i*rows - (rows-1)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-2)]), 
            grob_col(
                qqpng[i*rows - (rows-2)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-3)]), 
            grob_col(
                qqpng[i*rows - (rows-3)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-4)]), 
            grob_col(
                qqpng[i*rows - (rows-4)], 
                border=bord))
        #height = 100,
        #width = 100,
        #padding = 0
    ) 
}
#Do the last page
i=ceiling(length(manpng)/rows)
page[[i]]<-grob_layout(
    grob_row( #1 row
        grob_col(
            border=bord, 
            manpng[i*rows - (rows-1)]), 
        grob_col(
            qqpng[i*rows - (rows-1)], 
            border=bord)),
    grob_row( #2 rows
        grob_col(
            border=bord, 
            manpng[i*rows - (rows-2)]), 
        grob_col(
            qqpng[i*rows - (rows-2)], 
            border=bord)),
    grob_row( #3 rows
        grob_col(
            border=bord, 
            NA)),
    grob_row( #4 rows
        grob_col(
            border=bord, 
            NA))
    
)

view_grob(page[[4]])
#dir.create(paste(dir, "/pdf", sep=""))
outdir<-paste(dir, "/pdf", sep="")

grob_to_pdf(
    page,
    file_name = paste(outdir, "/DiscoveryM1_1.pdf", sep=""),
    add_page_numbers = F,
    meta_data_title = "PUFA-GWAS supplementary fig PDF"
)

# UKB-EUR discovery M2 -----------------------------------

dir<-"/Users/mike/Documents/Research/PUFA-GWAS/discovery/plot/M2"
pngfiles<-list.files(dir,full.names =  T)[grep(".png", list.files(dir))]
fileorder<-order_files_nopop(pngfiles)

manpng<-fileorder[grep("Manhattan",pngfiles)]
manpng
qqpng<-fileorder[grep("QQ",pngfiles)]
qqpng

page<-list()
rows=4
bord=F
numpages<-ceiling(length(manpng)/rows)
#%% = remainder
remainder<-length(manpng) %% rows

#at 
floor(length(manpng)/rows)

for (i in 1:(floor(length(manpng)/rows))){
    page[[i]]<-grob_layout(
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-1)]), 
            grob_col(
                qqpng[i*rows - (rows-1)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-2)]), 
            grob_col(
                qqpng[i*rows - (rows-2)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-3)]), 
            grob_col(
                qqpng[i*rows - (rows-3)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-4)]), 
            grob_col(
                qqpng[i*rows - (rows-4)], 
                border=bord))
        #height = 100,
        #width = 100,
        #padding = 0
    ) 
}
#Do the last page
i=ceiling(length(manpng)/rows)
page[[i]]<-grob_layout(
    grob_row( #1 row
        grob_col(
            border=bord, 
            manpng[i*rows - (rows-1)]), 
        grob_col(
            qqpng[i*rows - (rows-1)], 
            border=bord)),
    grob_row( #2 rows
        grob_col(
            border=bord, 
            manpng[i*rows - (rows-2)]), 
        grob_col(
            qqpng[i*rows - (rows-2)], 
            border=bord)),
    grob_row( #3 rows
        grob_col(
            border=bord, 
            NA)),
    grob_row( #4 rows
        grob_col(
            border=bord, 
            NA))
    
)

view_grob(page[[4]])
#dir.create(paste(dir, "/pdf", sep=""))
outdir<-paste(dir, "/pdf", sep="")

grob_to_pdf(
    page,
    file_name = paste(outdir, "/DiscoveryM2.pdf", sep=""),
    add_page_numbers = F,
    meta_data_title = "PUFA-GWAS supplementary fig PDF"
)




# AFRCSAEAS -----------------------------------------------------

dir<-"/Users/mike/Documents/Research/PUFA-GWAS/UKB-multi-ancestry/meta-analyses/AFRCSAEAS/plots"
pngfiles<-list.files(dir,full.names =  T)[grep(".png", list.files(dir))]
fileorder<-order_files_nopop(pngfiles)

manpng<-fileorder[grep("Manhattan",pngfiles)]
manpng
qqpng<-fileorder[grep("QQ",pngfiles)]
qqpng

page<-list()
rows=4
bord=F
numpages<-ceiling(length(manpng)/rows)

for (i in 1:(floor(length(manpng)/rows))){
    page[[i]]<-grob_layout(
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-1)]), 
            grob_col(
                qqpng[i*rows - (rows-1)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-2)]), 
            grob_col(
                qqpng[i*rows - (rows-2)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-3)]), 
            grob_col(
                qqpng[i*rows - (rows-3)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-4)]), 
            grob_col(
                qqpng[i*rows - (rows-4)], 
                border=bord))
        #height = 100,
        #width = 100,
        #padding = 0
    ) 
}
#Do the last page
i=ceiling(length(manpng)/rows)
page[[i]]<-grob_layout(
    grob_row( #1 row
        grob_col(
            border=bord, 
            manpng[i*rows - (rows-1)]), 
        grob_col(
            qqpng[i*rows - (rows-1)], 
            border=bord)),
    grob_row( #2 rows
        grob_col(
            border=bord, 
            manpng[i*rows - (rows-2)]), 
        grob_col(
            qqpng[i*rows - (rows-2)], 
            border=bord)),
    grob_row( #3 rows
        grob_col(
            border=bord, 
            NA)),
    grob_row( #4 rows
        grob_col(
            border=bord, 
            NA))
    
)

view_grob(page[[4]])
#dir.create(paste(dir, "/pdf", sep=""))
outdir<-paste(dir, "/pdf", sep="")

grob_to_pdf(
    page,
    file_name = paste(outdir, "/AFRCSAEAS.pdf", sep=""),
    add_page_numbers = F,
    meta_data_title = "PUFA-GWAS supplementary fig PDF"
)





# SmultiXcan ----------------------------------------------------

dir<-"/Users/mike/Documents/Research/PUFA-GWAS/discovery/SMultiXCan/plot2"
pngfiles<-list.files(dir,full.names =  T)[grep(".png", list.files(dir))]
fileorder<-order_files_nopop(pngfiles)

page<-list()
rows=2
bord=F

for (i in 1:(floor(length(fileorder)/rows))){
    page[[i]]<-grob_layout(
        grob_row(
            grob_col(
                border=bord, 
                fileorder[i*rows - (rows-1)])), 
        grob_row(
            grob_col(
                border=bord, 
                fileorder[i*rows - (rows-2)]))
        # ,
        # grob_row(
        #     grob_col(
        #         border=bord, 
        #         fileorder[i*rows - (rows-3)]))
        # ,
        # grob_row(
        #     grob_col(
        #         border=bord, 
        #         fileorder[i*rows - (rows-4)]))
        #height = 100,
        #width = 90
       # padding = 5
    ) 
}

view_grob(page[[3]])
#dir.create(paste(dir, "/pdf", sep=""))
outdir<-paste(dir, "/pdf", sep="")

grob_to_pdf(
    page,
    file_name = paste(outdir, "/Smulti.pdf", sep=""),
    add_page_numbers = F,
    meta_data_title = "PUFA-GWAS supplementary fig PDF"
)


# Participant plots ----------------------------------------------------

dir<-"/Users/mike/Documents/Research/PUFA-GWAS/supplementary figures/Participant-characteristics"
pngfiles<-list.files(dir,full.names =  T)[grep(".png", list.files(dir))]
fileorder<-pngfiles

page<-list()
rows=2
bord=F

for (i in 1:(floor(length(fileorder)/rows))){
    page[[i]]<-grob_layout(
        grob_row(
            grob_col(
                border=bord, 
                fileorder[i*rows - (rows-1)])), 
        grob_row(
            grob_col(
                border=bord, 
                fileorder[i*rows - (rows-2)]))
        # ,
        # grob_row(
        #     grob_col(
        #         border=bord, 
        #         fileorder[i*rows - (rows-3)]))
        # ,
        # grob_row(
        #     grob_col(
        #         border=bord, 
        #         fileorder[i*rows - (rows-4)]))
        #height = 100,
        #width = 90
        # padding = 5
    ) 
}

view_grob(page[[3]])
#dir.create(paste(dir, "/pdf", sep=""))
outdir<-paste(dir, "/pdf", sep="")

grob_to_pdf(
    page,
    file_name = paste(outdir, "/Participant.pdf", sep=""),
    add_page_numbers = F,
    meta_data_title = "PUFA-GWAS supplementary fig PDF"
)




# Compare to opengwas -------------------------------------------

dirM1<-"/Users/mike/Documents/Research/PUFA-GWAS/discovery/compareToOpenGWAS/update03142022/M1/update03142022"
dirM2<-"/Users/mike/Documents/Research/PUFA-GWAS/discovery/compareToOpenGWAS/update03142022/M2"
pngfilesM1<-list.files(dirM1,full.names =  T)[grep(".png", list.files(dirM1))]
pngfilesM2<-list.files(dirM2,full.names =  T)[grep(".png", list.files(dirM2))]

order<-c(8,7,11,10,9,2,1,4,3,14,13,6,5,12)
pngfilesM1<-pngfilesM1[order]
pngfilesM2<-pngfilesM2[order]

page<-list()
rows=3
bord=F
numpages<-ceiling(length(pngfilesM1)/rows)

for (i in 1:(floor(length(pngfilesM1)/rows))){
    page[[i]]<-grob_layout(
        grob_row(
            grob_col(
                border=bord, 
                pngfilesM1[i*rows - (rows-1)]), 
            grob_col(
                pngfilesM2[i*rows - (rows-1)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                pngfilesM1[i*rows - (rows-2)]), 
            grob_col(
                pngfilesM2[i*rows - (rows-2)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                pngfilesM1[i*rows - (rows-3)]), 
            grob_col(
                pngfilesM2[i*rows - (rows-3)], 
                border=bord))
        # ,
        # grob_row(
        #     grob_col(
        #         border=bord, 
        #         pngfilesM1[i*rows - (rows-4)]), 
        #     grob_col(
        #         pngfilesM2[i*rows - (rows-4)], 
        #         border=bord))
        #height = 100,
        #width = 100,
        #padding = 0
    ) 
}
#Do the last page
i=ceiling(length(pngfilesM1)/rows)
page[[i]]<-grob_layout(
    grob_row( #1 row
        grob_col(
            border=bord, 
            pngfilesM1[i*rows - (rows-1)]), 
        grob_col(
            pngfilesM2[i*rows - (rows-1)], 
            border=bord)),
    grob_row( #2 rows
        grob_col(
            border=bord, 
            pngfilesM1[i*rows - (rows-2)]), 
        grob_col(
            pngfilesM2[i*rows - (rows-2)], 
            border=bord)),
    grob_row( #3 rows
        grob_col(
            border=bord, 
            NA))
    
)

#view_grob(page[[4]])
#dir.create(paste(dirM1, "/pdf", sep=""))
outdir<-paste(dirM1, "/pdf", sep="")

grob_to_pdf(
    page,
    file_name = paste(outdir, "/comparetoopenGWAS.pdf", sep=""),
    add_page_numbers = F,
    meta_data_title = "PUFA-GWAS supplementary fig PDF"
)

# EURAFRCSAEAS -----------------------------------------------------

dir<-"/Users/mike/Documents/Research/PUFA-GWAS/3MAMA/EURAFRCSAEAS/plots"
pngfiles<-list.files(dir,full.names =  T)[grep(".png", list.files(dir))]
fileorder<-order_files_nopop(pngfiles)
fileorder
manpng<-fileorder[grep("Manhattan",pngfiles)]
manpng
qqpng<-fileorder[grep("QQ",pngfiles)]
qqpng

page<-list()
rows=4
bord=F
numpages<-ceiling(length(manpng)/rows)

for (i in 1:(floor(length(manpng)/rows))){
    page[[i]]<-grob_layout(
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-1)]), 
            grob_col(
                qqpng[i*rows - (rows-1)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-2)]), 
            grob_col(
                qqpng[i*rows - (rows-2)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-3)]), 
            grob_col(
                qqpng[i*rows - (rows-3)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-4)]), 
            grob_col(
                qqpng[i*rows - (rows-4)], 
                border=bord))
        #height = 100,
        #width = 100,
        #padding = 0
    ) 
}
#Do the last page
i=ceiling(length(manpng)/rows)
page[[i]]<-grob_layout(
    grob_row( #1 row
        grob_col(
            border=bord, 
            manpng[i*rows - (rows-1)]), 
        grob_col(
            qqpng[i*rows - (rows-1)], 
            border=bord)),
    grob_row( #2 rows
        grob_col(
            border=bord, 
            manpng[i*rows - (rows-2)]), 
        grob_col(
            qqpng[i*rows - (rows-2)], 
            border=bord)),
    grob_row( #3 rows
        grob_col(
            border=bord, 
            NA)),
    grob_row( #4 rows
        grob_col(
            border=bord, 
            NA))
    
)

#view_grob(page[[4]])
#dir.create(paste(dir, "/pdf", sep=""))
outdir<-paste(dir, "/pdf", sep="")

grob_to_pdf(
    page,
    file_name = paste(outdir, "/EURAFRCSAEAS.pdf", sep=""),
    add_page_numbers = F,
    meta_data_title = "PUFA-GWAS supplementary fig PDF"
)

# MA-everything -----------------------------------------------------

dir<-"/Users/mike/Documents/Research/PUFA-GWAS/3MAMA/everything/plots"
pngfiles<-list.files(dir,full.names =  T)[grep(".png", list.files(dir))]
fileorder<-order_files_nopop(pngfiles)

manpng<-fileorder[grep("Manhattan",pngfiles)]
manpng
qqpng<-fileorder[grep("QQ",pngfiles)]
qqpng

page<-list()
rows=4
bord=F
numpages<-ceiling(length(manpng)/rows)

for (i in 1:(floor(length(manpng)/rows))){
    page[[i]]<-grob_layout(
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-1)]), 
            grob_col(
                qqpng[i*rows - (rows-1)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-2)]), 
            grob_col(
                qqpng[i*rows - (rows-2)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-3)]), 
            grob_col(
                qqpng[i*rows - (rows-3)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-4)]), 
            grob_col(
                qqpng[i*rows - (rows-4)], 
                border=bord))
        #height = 100,
        #width = 100,
        #padding = 0
    ) 
}
#Do the last page
i=ceiling(length(manpng)/rows)
page[[i]]<-grob_layout(
    grob_row( #1 row
        grob_col(
            border=bord, 
            manpng[i*rows - (rows-1)]), 
        grob_col(
            qqpng[i*rows - (rows-1)], 
            border=bord)),
   
    grob_row( #2 rows
            grob_col(
                border=bord, 
                NA)),
    grob_row( #3 rows
        grob_col(
            border=bord, 
            NA)),
    grob_row( #4 rows
        grob_col(
            border=bord, 
            NA))
    
)

view_grob(page[[4]])
#dir.create(paste(dir, "/pdf", sep=""))
outdir<-paste(dir, "/pdf", sep="")

grob_to_pdf(
    page,
    file_name = paste(outdir, "/everything.pdf", sep=""),
    add_page_numbers = F,
    meta_data_title = "PUFA-GWAS supplementary fig PDF"
)



# Combine PDF ---------------------------------------------------

files<-list.files("/Users/mike/Documents/Research/PUFA-GWAS/supplementary figures", recursive = T, full.names = T)
files<-files[grep(".pdf", files)]
files

pdf_combine(input = files, 
            output = "/Users/mike/Documents/Research/PUFA-GWAS/supplementary figures/combine.pdf")


