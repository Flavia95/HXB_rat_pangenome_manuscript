library(VennDiagram)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)

######Venn diagram
ListGeno <- 258502
ListPang <- 217917 
ListInt <- 278670
#F8766D
#A3A500
#00BF7D
#00B0F6
#E76BF4
#cat.pos         = c(285, 105)
#overlap #009171
venn.plot <- draw.pairwise.venn(ListPang+ListInt, ListGeno+ListInt, ListInt,
                                category        = c("JC-only", "vg-only"),
                                fill            = c("#00B0F6", "#E76BF4"),
                                lty             = "blank",
                                cex             = 1,
                                cat.cex         = 1,
                                cat.dist        = 0.0009,
                                ext.pos         = 30,
                                ext.dist        = -0.05,
                                ext.length      = 0.85,
                                scaling.factor = 0.8,
                                family = "Arial",
                                ext.line.lwd    = 2,)
###########################################
# set the file names and corresponding region names
file_names <- c('easy.onlyjt.allsamples.gt.norm.info', 'easy.onlypg.allsamples.gt.norm.info','snps.easy.onlyjt.allsamples.gt.norm.info', 'snps.easy.onlypg.allsamples.gt.norm.info','indels.easy.onlyjt.allsamples.gt.norm.info',
                'indels.easy.onlypg.allsamples.gt.norm.info','hard.onlyjt.allsamples.gt.norm.info','hard.onlypg.allsamples.gt.norm.info','snps.hard.onlyjt.allsamples.gt.norm.info',
                'snps.hard.onlypg.allsamples.gt.norm.info','indels.hard.onlyjt.allsamples.gt.norm.info','indels.hard.onlypg.allsamples.gt.norm.info', 'easy.overlap.allsamples.gt.norm.info','indels.easy.overlap.allsamples.gt.norm.info', 'snps.easy.overlap.allsamples.gt.norm.info','hard.overlap.allsamples.gt.norm.info',  'indels.hard.overlap.allsamples.gt.norm.info',  'snps.hard.overlap.allsamples.gt.norm.info')

region_names <- c("easy", "easy", "easysnps", "easysnps", "easyindels", "easyindels", "hard","hard", "hardsnps","hardsnps","hardindels", "hardindels", "easy","easyindels", "easysnps",'hard','hardindels','hardsnps')

# create an empty list to store the processed data frames
df_list <- list()
# loop through the file names and apply the data processing steps
for (i in seq_along(file_names)) {
  file_name <- file_names[i]
  region_name <- region_names[i]
  # read in the data
  myd <- read.table(file_name, header=FALSE, sep='\t', comment.char='')
  colnames(myd) <- c('CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'AC', 'AF') 
  myd_separate <- separate_rows(myd, AF, convert=TRUE)
    # add metadata columns
  if (grepl("overlap", file_name)) {
    myd_separate$Methods <- "overlap"
  } else if (grepl("pg", file_name)) {
    myd_separate$Methods <- "Pangenomic"
  } else {
    myd_separate$Methods <- "Genomic"
  }
  myd_separate$Region <- region_name
    # add processed data frame to the list
  df_list[[i]] <- myd_separate
}

# combine the processed data frames into a single data frame
h <- do.call(rbind, df_list)
# convert Region to a factor with ordered levels
h$Region <- factor(h$Region, levels = c("easy", "easysnps", "easyindels","hard", "hardsnps",'hardindels'))
# create the plot
q <- h %>% select(CHROM, POS, AF, Methods, Region) %>% group_by(CHROM, POS, Methods, Region)

#supplementary figure
AFS_all <- ggplot(q, aes(AF, fill = Methods)) +
  geom_histogram(bins = 30, position = "dodge") +
  theme_minimal() +
  theme_light() +
  facet_grid(Region ~ ., scales = "free") +
  theme_classic() +
  labs(x = 'Alternate Frequency (SHR/OlaIpcv)', y = 'Count', fill = ' ') +
  scale_fill_manual(values = c(
    "Genomic" = "#00B0F6",
    "Pangenomic" = "#E76BF4",
    "overlap" = "#006cc2"
  ), labels = c("JC-only", "overlap", "vg-only")) +
  guides(fill = guide_legend(title = " ")) +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) 
ggsave("AFS_all.png", AFS_all, width = 7, height = 7, dpi=900)

####Only easy and complex
q_filtered <- q %>% filter(Region %in% c("easy", "hard"))
AFS_easy_hard <- ggplot(q_filtered, aes(AF, fill = Methods)) +
  geom_histogram(bins = 50, position = "dodge") + theme_classic() +
  labs(x = 'Alternate Frequency (SHR/OlaIpcv)', y = 'count', fill = ' ') +
  scale_fill_manual(values = c(
    "Genomic" = "#00B0F6",
    "Pangenomic" = "#E76BF4",
    "overlap" = "#006cc2"
  ), labels = c("jc-only", "overlap", "vg-only")) +
  guides(fill = guide_legend(title = "none")) +
  theme(legend.position = "t") +
  facet_grid(Region ~ ., scales = "free")  

#####INDELS
file_names <- c("easypg.indel.hist", "easyjt.indel.hist","easyoverlap.indel.hist","hardpg.indel.hist","hardjt.indel.hist","hardoverlap.indel.hist")
region_names <- c("easy", "easy","easy", "hard","hard","hard")

# create an empty list to store the processed data frames
df_list <- list()

# loop through the file names and apply the data processing steps
for (i in seq_along(file_names)) {
  file_name <- file_names[i]
  region_name <- region_names[i]
  
  myd <- read.table(file_name, header=T, sep='\t', comment.char='')
  myd <- subset(myd, LENGTH!=0)
  
  # add metadata columns
  #myd$Methods  <- ifelse(grepl("pg", file_name), "Pangenomic", "Genomic")
  #myd$Region  <- region_name
  # add metadata columns
  if (grepl("overlap", file_name)) {
    myd$Methods <- "overlap"
  } else if (grepl("pg", file_name)) {
    myd$Methods <- "Pangenomic"
  } else {
    myd$Methods <- "Genomic"
  }
  myd$Region <- region_name
  
  # add processed data frame to the list
  df_list[[i]] <- myd
}

# combine the processed data frames into a single data frame
h <- do.call(rbind, df_list)

# convert Region to a factor with ordered levels
h$Region <- factor(h$Region, levels = c("easy", "hard"))
# plot the data
length_indels <- ggplot(h) +
  geom_col(aes(x = LENGTH, y = COUNT, fill = Methods)) +
  xlim(-50, 50) +
  ylim(0, 6000) +
  facet_grid(Region ~ ., scales = "free") +theme_classic()  + labs(x = 'length', y = 'count', fill = ' ')  + theme_classic() +
  labs(x = 'length', y = 'count', fill = ' ') +
  scale_fill_manual(values = c(
    "Genomic" = "#00B0F6",
    "Pangenomic" = "#E76BF4",
    "overlap" = "#006cc2"
  ), labels = c("JC-only", "overlap", "vg-only")) +
  guides(fill = guide_legend(title = " "))  +
  theme(legend.position = "top")

#Panel fig2
firstrow<- plot_grid(venn.plot, length_indels, labels = c('A', 'B'), label_size = 12, ncol=2 )
c=plot_grid(firstrow, AFS_easy_hard,  labels= c('', 'C', 'D') ,nrow =2 )
ggsave('fig2.png', dpi=300, width = 30, heigh = 12, units='cm')
