library(ggplot2)
library(dplyr)
library(ComplexUpset)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

############1. upset plot
input_df <- read.delim("overalp_4methods.txt", sep="\t", header=TRUE)
indvs <- input_df %>%uncount(count) 
symptoms <- c("hall", "vg", "svim", "pav")

upset_query1 <- upset_query(intersect = c('pav','svim',"hall", "vg"), color = '#706092')
upset_query2 <- upset_query(intersect = c("svim", "vg",'hall'), color = '#706092')
upset_query3 <- upset_query(intersect = c("hall", "pav"), color = '#706092')
upset_query4 <- upset_query(intersect = c("svim", "hall"), color = '#706092')
upset_query5 <- upset_query(intersect = c("pav", "vg"), color = '#706092')
upset_query6 <- upset_query(intersect = c("svim", "hall",'vg'), color = '#706092')
upset_query7 <- upset_query(intersect = c("pav", "svim",'vg'), color = '#706092')
upset_query8 <- upset_query(intersect = c("pav", "svim"), color = '#706092')
upset_query9 <- upset_query(intersect = c("pav", "hall",'vg'), color = '#706092')
upset_query10 <- upset_query(intersect = c("svim",'vg'), color = '#706092')
upset_query11 <- upset_query(intersect = c("hall",'vg'), color = '#706092')

s4a <- upset(data = indvs, intersect = symptoms,
                    name = "chr12, SHR sample",
                    min_size = 0,
                    width_ratio = 0.125,
                    queries = list(upset_query1, upset_query2,upset_query3,upset_query4,upset_query5,upset_query6,upset_query7,upset_query8,upset_query9,upset_query10, upset_query11)
)
ggsave("upset_SHR.png", width = 11, height = 6, dpi=400)

# UPSET ALL SAMPLES
input_df <- read.delim("SV_atleast2.txt", sep="\t", header=TRUE)
indvs <- input_df %>%uncount(count) 
symptoms <- c("hall", "vg", "svim", "pav")

# Create the UpSet plot
upset_query1 <- upset_query(intersect = c('pav','svim',"hall", "vg"), color = '#706092')
upset_query2 <- upset_query(intersect = c("svim", "pav",'hall'), color = '#706092')
upset_query3 <- upset_query(intersect = c("hall", "pav"), color = '#706092')
upset_query4 <- upset_query(intersect = c("svim", "hall"), color = '#706092')
upset_query5 <- upset_query(intersect = c("pav", "vg"), color = '#706095')
upset_query6 <- upset_query(intersect = c("svim", "hall",'vg'), color = '#706092')
upset_query7 <- upset_query(intersect = c("pav", "svim",'vg'), color = '#706092')
upset_query8 <- upset_query(intersect = c("pav", "svim"), color = '#706092')
upset_query9 <- upset_query(intersect = c("pav", "hall",'vg'), color = '#706092')
upset_query10 <- upset_query(intersect = c("svim",'vg'), color = '#706092')
upset_query11 <- upset_query(intersect = c("hall",'vg'), color = '#706092')

s4b <- upset(data = indvs, intersect = symptoms,
                    name = "chr12, 31 samples",
                    min_size = 0,
                    width_ratio = 0.125,
                    queries = list(upset_query1, upset_query2,upset_query3,upset_query4,upset_query5,upset_query6,upset_query7,upset_query8,upset_query9,upset_query10, upset_query11)
)

ggsave("upset_allsamples.png", width = 5, height = 5, dpi=600)

#final plot
fa <- plot_grid(s4a,labels = c('A'), label_size = 12,nrow=1)
fc <- rasterGrob(png::readPNG("S4C.png"))
fb_c <- plot_grid(s4b,figS4c,labels = c('B','C'), label_size = 12,nrow=1)
f4 <- plot_grid(fa, fb_c,nrow = 2)

ggsave('S4.png',f4, dpi=300, width = 15, heigh = 15, units='cm')

