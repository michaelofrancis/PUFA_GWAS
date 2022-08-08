#install.packages(c("gridExtra","png","grobblR","qpdf","pdftools"))
library(gridExtra)
library(png)
library(grobblR)
library(qpdf)
library(pdftools)

##-------------------------------------------------------
##SCRIPT TO ARRANGE MANY SUPPLEMENTARY FIGURES INTO PDF##
##-------------------------------------------------------


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


# SUPPLEMENTARY FIG 2: Compare M1 and M2 -------------------------------------------

dir<-"/Users/mike/Documents/Research/PUFA-GWAS/discovery/compareM1M2"
pngfile<-list.files(dir,full.names =  T)[grep(".png", list.files(dir))]


pngfiles<-order_files_nopop(pngfile)
pngfiles

page<-list()
rows=3
bord=F
numpages<-3

j=0
while (j<2){
    page[[j+1]]<-grob_layout(
        grob_row(
            grob_col(
                border=bord, 
                pngfiles[1+(j*6)]), 
            grob_col(
                pngfiles[2+(j*6)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                pngfiles[3+(j*6)]), 
            grob_col(
                pngfiles[4+(j*6)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                pngfiles[5+(j*6)]), 
            grob_col(
                pngfiles[6+(j*6)], 
                border=bord))
    )

     j<-j+1
}

#Do the last page
page[[3]]<-grob_layout(
    grob_row( #1 row
        grob_col(
            border=bord, 
            pngfiles[13]), 
        grob_col(
            pngfiles[14], 
            border=bord)),
    grob_row( #2 rows
        grob_col(
            border=bord, 
            NA)),
     grob_row( #3 rows
        grob_col(
            border=bord, 
            NA))
)

#view_grob(page[[4]])
#dir.create(paste(dirM1, "/pdf", sep=""))
outdir<-paste(dir, "/pdf", sep="")

grob_to_pdf(
    page,
    file_name = paste(outdir, "/M2vM1.pdf", sep=""),
    add_page_numbers = F,
    meta_data_title = "PUFA-GWAS supplementary fig PDF"
)

# SUPPLEMENTARY FIG 3: UKB-EUR discovery Manhattan/QQ ---------------

dir<-"/Users/mike/Documents/Research/PUFA-GWAS/discovery/plot/newest"
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

#view_grob(page[[4]])
dir.create(paste(dir, "/pdf", sep=""))
outdir<-paste(dir, "/pdf", sep="")

grob_to_pdf(
    page,
    file_name = paste(outdir, "/Discovery-ManQQ-hq.pdf", sep=""),
    add_page_numbers = F,
    meta_data_title = "PUFA-GWAS supplementary fig PDF"
)



# SUPPLEMENTARY FIG 4: SmultiXcan ----------------------------------------------------

dir<-"/Users/mike/Documents/Research/PUFA-GWAS/discovery/SMultiXCan/plot-newest/"
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


# SUPPLEMENTARY FIG 6: UKB-multi-ancestry --------------------------------------
dir<-"/Users/mike/Documents/Research/PUFA-GWAS/UKB-multi-ancestry/plots-new"
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

# SUPPLEMENTARY FIGURE 7: LocusZoom Plots---------------------------------------------------
dir="/Users/mike/Documents/Research/PUFA-GWAS/MANUSCRIPT_FIGURES/locuszoomACSL6"
fileorder<-"/Users/mike/Documents/Research/PUFA-GWAS/MANUSCRIPT_FIGURES/locuszoomACSL6/w3-ACSL6-fromFUMA2 copy.png"

rows=2
bord=F

page<-grob_layout(
    grob_row( #1 row
        grob_col(
            border=bord, 
            fileorder[1])),
    grob_row( #2 rows
        grob_col(
            border=bord, 
            NA
        )))

outdir<-paste(dir, "/pdf", sep="")

grob_to_pdf(
    page,
    file_name = paste(outdir, "/ACSL6.LocusZoom.pdf", sep=""),
    add_page_numbers = F,
    meta_data_title = "PUFA-GWAS supplementary fig PDF"
)


# Combine PDF ---------------------------------------------------

files<-list.files("/Users/mike/Documents/Research/PUFA-GWAS/MANUSCRIPT_FIGURES", recursive = T, full.names = T)
files<-files[grep(".pdf", files)]
files<-files[c(1,2,5,6,7,8,9)]
files

pdf_combine(input = files, 
            output = "/Users/mike/Documents/Research/PUFA-GWAS/MANUSCRIPT_FIGURES/Supplementary-figs-08022022.pdf")
