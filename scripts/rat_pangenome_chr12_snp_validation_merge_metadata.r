# ref vs alt
df_refalt<-read.table(file="./pangenome_variant_validation_2022-12_20_500targets.csv", head=T, sep=",")

# fasta names 
df_faname<-read.table(file="./2022_D_with_complexity_pangenome_primers_meta.csv", head=T, sep=",")
df_faname$ID<-tolower(df_faname$ID) # two IDs started with "Chr" instead of "chr"
head(df_faname)
head(df_refalt)
dim(df_faname)

df_refalt$ID<-paste0(df_refalt$CHROM, ":", df_refalt$start, "-", df_refalt$end)
df1<-merge(df_faname, df_refalt, by="ID", all.x=T)
head(df1)

df_fa<-read.table(file="./pangenome_sanger_fasta_2023-01-27.tab", head=F, sep="\t")
df_fa[1,]
names(df_fa)<-c("name", "ab1name", "seq")

dfall<-merge(df1, df_fa, by="name", all=T)
subset(dfall, is.na(REF))
# "None", i.e. all fasta files have ref and alt assigned
nopcr<-subset(dfall, is.na(seq))
dim(nopcr)
# [1] 10 19
# we have five pcr failed, these samples were not sent to genewiz. We sequence from both ends for each pcr. So it is 10 rows.
write.table(file="rat_pangenome_chr12_snp_validation_all_data_2023-02-03.tab", dfall, row.names=F)


