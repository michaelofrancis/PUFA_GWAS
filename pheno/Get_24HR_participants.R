#Generate file to indicate whether a particpant had taken the 
#24 hour recall survey or not (n=211018)

library(plyr)
library(dplyr)
library(tidyverse)

#Load UK Biobank datasets----------------------------
source('../ukb34137_loaddata.r') #15 min


#Get participants who took 24HR coded as 0/1----------
daycols<-c("f.20080.0.0", "f.20080.1.0", "f.20080.2.0",
           "f.20080.3.0", "f.20080.4.0")
#Subset and change all the values to characters
bd1<-apply(bd[,daycols], 2, as.character)
#Change NA's to zeros and days to 1's
bd1[is.na(bd1[,daycols])]<-0
bd1[bd1[,daycols]!="0"] <-1
#Change these back to numeric
bd1<-apply(bd1[,daycols], 2, as.integer)
#Now make a new column, everyone with rowSums zero write FALSE
#and those with >0 write TRUE
sum<-apply(bd1[,daycols], 1, sum)
bd1<-as_tibble(as.data.frame(bd1))
bd1$took_24HR<-sum
bd1$took_24HR[bd1$took_24HR>0]<-1
sum(bd1$took_24HR==1) #[1] 211018 **SUCCESS**
bd1$IID<-bd$f.eid
bd1<-bd1%>%select(IID, took_24HR)
write.table(bd1, "UKB-took24HRparticipants-09262021.txt", row.names=FALSE,
                    quote=FALSE)
