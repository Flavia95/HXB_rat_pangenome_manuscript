Code adapted from:
#https://github.com/DannyArends/BXDtools
#https://github.com/mae47/Convergent_evolution/blob/main/BXD_mice/Scripts/phewas.create.geno.r

library(dplyr)
library(stringr)
library(tibble) 
##Load vcf (12 SNPs validated, called by vg-only), with splitted multiallalic sites
vcf <- read.table("12.onlyvg.nomulti.gt.txt")
#bcftools query -l file.vcf.gz > samples.txt 
samples <- read.table("samples.txt", header=FALSE, stringsAsFactors=FALSE)
#Convert the samples to a vector
samples <- as.vector(samples$V1)
#Combine the samples with column names to create the header
vcf_header <- c("Chr", "POS", samples)
colnames(vcf) <- vcf_header
#Rearrange to numerical order
#Current order
x<-colnames(vcf)
info_cols<-c("Chr","POS")
#List only HXB cols
HXB_samples<-setdiff(x, info_cols)
#rearrange
vcf<-vcf[,c(info_cols,str_sort(c(HXB_samples),numeric=TRUE))]
#Samples now in numerical order:
colnames(vcf)
##Make column "Locus" with merged chr and pos for row names
vcf<-transform(vcf, Locus=paste0(Chr,"_",POS))
#move from end of dataframe to start
vcf<-relocate(vcf, Locus, .before="Chr")

##Check all row names are unique
#number of rows in dataset, 21886 row names
nrow(vcf)
vcf[,"Locus"][duplicated(vcf[,"Locus"])] 
#keep only unique values for row name column
vcf<-vcf[!(duplicated(vcf["Locus"]) | duplicated(vcf["Locus"], fromLast = TRUE)),]
#recount number of rows in edited dataset. 21312 row names
nrow(vcf)

#fix:remove multiallelics and skip 2,3
##Replace vcf calls with numeric or alpha codes
#Do this BEFORE adding missing cM column
vcf[vcf=="0"] <- "B" #later code to -1
vcf[vcf=="1"] <- "A" #later code back to 1
vcf[vcf=="."] <- "U" #later code to Na 
vcf[vcf=="2"] <- "U" #later code to Na 
vcf[vcf=="3"] <- "U" #later code to Na 

#Add missing column cM
vcf<-add_column(vcf, cM='.', .after="Chr")
#Add column Mb as factor of POS
vcf<-transform(vcf, Mb=(POS/1000000))
#Move column Mb from end of dataframe to near start
vcf<-relocate(vcf, Mb, .before="POS")
#Remove POS column
vcf<-subset(vcf, select=-c(POS))
##remove 'chr' from in front of chr number
vcf$Chr<-gsub("chr","",vcf$Chr)

#Summary of all calls that are in the dataframe, check if these are 2 or 3
table(unlist(vcf[,-c(1:4)]), useNA="always")
#Save out file
write.table(vcf,"HXB.geno",row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
