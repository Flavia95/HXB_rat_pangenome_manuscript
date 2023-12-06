library("RColorBrewer")
library("ggplot2")
library("tidyr") 
library("RColorBrewer") 
library(grid)
library(gridExtra)
library(cowplot)


#1. Classification type of variants
mym<-read.table('micro12_allsamples.txt', header=T , sep='\t')
mym$variant <- factor(mym$variant, levels = c("SNPs", "MNPs", "INDELs", "SNP/INDEL", "MNP/INDEL", "INDEL/CLUMPED"))
mym$type<- factor(mym$type, levels = c("simple", "complex"))
level_order <- c("SNPs", "MNPs", "INDELs", "SNP/INDEL", "MNP/INDEL", "INDEL/CLUMPED")
color_palette <- c("burlywood1", "burlywood2", "burlywood3", "darkslategray2", "darkslategray3", "darkslategray4")
mym$variant <- factor(mym$variant, levels = level_order)
ci=ggplot(mym,aes(caller, n, fill = variant)) +geom_bar(stat = 'identity',width = 0.5) +facet_grid(type ~ ., scales = 'free_y') +theme_classic() +scale_fill_manual(values = color_palette) +
  labs(x = '', y = 'count', fill = 'Variant type') +theme(legend.position = "left",rect = element_rect(fill = "transparent"))+theme(legend.position = "left",rect = element_rect(fill = "transparent")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),  # Set the x-axis label color to black
        axis.text.y = element_text(color = "black"))
ci2=ci + theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), axis.text=element_text(size=12), text = element_text(family = "Arial", size = 12), legend.position = "left", plot.tag = element_text(face = "bold"))+labs(title = " ", tag = "B",size=12)
ggsave("micro12_allsamples.png" , dpi=1000, width = 4, heigh = 3)

#2.Type classification as impact
myd=read.table("TypeCountImpact_sim_com.txt", header=T)
#mydsmall <- subset(myd, Impact %in% c("HIGH", "MODIFIER"))
mydsmall <- subset(myd, Impact %in% c("HIGH"))
myd1=ggplot(mydsmall, aes(Type, Count, fill = Variant_type)) +
  geom_col(position = "stack", width = 0.5) +
  coord_flip() +
  scale_y_log10() +
  scale_fill_manual(values = c("simple" = "burlywood3", "complex" = "darkslategray3"),
                    breaks = c("simple", "complex"), labels = c("Simple", "Complex")) +
  theme_classic() +
  labs(x = '', y = 'count (Log scale)', fill = 'Variant type') +
  theme(legend.position = "none", rect = element_rect(fill = "transparent"), plot.tag = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"),text = element_text(family = "Arial", size = 12),
        legend.position = "none") +
  ggtitle("vg-High Impact") +labs(title = " ", tag = "C",size=12)

mydsmall1 <- subset(myd, Impact %in% c("MODIFIER"))
ci=ggplot(mydsmall1, aes(Type, Count, fill = Variant_type)) +
  geom_col(position = "stack", width = 0.5) +
  coord_flip() +
  scale_y_log10() +
  scale_fill_manual(values = c("simple" = "burlywood3", "complex" = "darkslategray3"),
                    breaks = c("simple", "complex"), labels = c("Simple", "Complex")) +
  theme_classic() +
  labs(x = '', y = 'count (Log scale)', fill = 'Variant type') +
  theme(legend.position = "none", rect = element_rect(fill = "transparent"), plot.tag = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"),text = element_text(family = "Arial", size = 12),
        legend.position = "none") +
  ggtitle("vg-Modifier Impact")+labs(title = " ", tag = "D",size=12)

ggsave("typeimpact_allsamples.png" , dpi=1000, width = 5, heigh = 5)

# Read the data from the CSV file
data <- read.table("allsamples_atleast2.txt", header = TRUE, sep = " ")
color_palette <- c("#706092", "#b3aac7")

# Create a bar plot using ggplot
plot <- ggplot(data, aes(x = tool, fill = Type)) +
  geom_bar(width = 0.5) + 
  scale_fill_manual(values = color_palette, labels = c("Deletion", "Insertion"))  +theme(legend.position = "left",rect = element_rect(fill = "transparent")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), text = element_text(family = "Arial", size = 12) # Set the x-axis label color to black
  ) + labs(x = '', y = 'count', fill = 'Variant type') + theme_classic()#+ facet_grid(Type ~ ., scales = 'free_y') +theme_classic()
plot2=plot + theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), axis.text=element_text(size=12), text = element_text(size=12), legend.position = "left", plot.tag = element_text(face = "bold"))+labs(title = " ", tag = "E",size=12)
ggsave("SV12_allsamples.png" , plot2, dpi=1000, width = 4, heigh = 3)

#2.Type classification as impact
myd2=read.table("TypeCountImpact_SV_ins_del.txt", header=T)
#myd3 <- subset(myd2, Impact %in% c("HIGH", "MODIFIER"))
myd3 <- subset(myd2, Impact %in% c("HIGH"))
p=ggplot(myd3, aes(Type, Count, fill = Variant_type)) +
  geom_col(position = "stack", width = 0.5) +
  coord_flip() +
  scale_y_log10() +
  scale_fill_manual(values = c("deletion" = "#706092", "insertion" = "#b3aac7"),
                    breaks = c("deletion", "insertion"), labels = c("Deletion", "Insertion")) +
  theme_classic() +
  labs(x = '', y = 'count(Log scale)', fill = 'Variant type') +
  theme(legend.position = "none", rect = element_rect(fill = "transparent"), plot.tag = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"),text = element_text(family = "Arial", size = 12),
        legend.position = "none") +
  ggtitle("High Impact")+labs(title = " ", tag = "F",size=12)

myd4 <- subset(myd2, Impact %in% c("MODIFIER"))
p2=ggplot(myd4, aes(Type, Count, fill = Variant_type)) +
  geom_col(position = "stack", width = 0.5) +
  coord_flip() +
  scale_y_log10() +
  scale_fill_manual(values = c("deletion" = "#706092", "insertion" = "#b3aac7"),
                    breaks = c("deletion", "insertion"), labels = c("Deletion", "Insertion")) +
  theme_classic() +
  labs(x = '', y = 'count (Log scale)', fill = 'Variant type') +
  theme(legend.position = "none", rect = element_rect(fill = "transparent"), plot.tag = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"),text = element_text(family = "Arial", size = 12),
        legend.position = "none") +
  ggtitle("Modifier Impact") +labs(title = " ", tag = "G",size=12)

ggsave("typeimpactSV_allsamples.png" ,myd2, dpi=1000, width = 5, heigh = 5)

fig1a <- rasterGrob(png::readPNG("chr12.pan+ref.fa.gz.f6a1f45.c2fac19.eec7a56.smooth.final.sort.png"))
# fig1a in ggplot
fig1a_plot <- ggplot() + 
  annotation_custom(fig1a, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+theme_bw() +  # Or any other theme you prefer
  theme(
    plot.background = element_rect(fill = "white"),
    panel.border = element_blank()  # Remove panel borders
  )

# Paired plots 
a1 = (ci2 / plot2) 
a2 = (myd1 / p) 
a3 = (ci / p2)

figure <- plot_grid(
  fig1a_plot,  
  plot_grid(a1,a2, a3,ncol = 3), 
  ncol = 1  
)

figure <- figure +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')

ggsave("Figure1.tiff" ,figure, dpi=700, width = 13, heigh = 10)



#firstrow <- plot_grid(fig1a, labels = c('A'), label_size = 12,ncol=1)
#secondrow <- plot_grid(ci2,myd1,ci, labels = c('B', 'C','D'), label_size = 12, ncol = 3)
#thirdrow <- plot_grid(plot2, p,p2, labels = c('E', 'F','G'), label_size = 12, ncol = 3)
#final_plot <- plot_grid(firstrow, secondrow,thirdrow, nrow = 3)
#ggsave('fig2.eps', dpi=300, width = 32, heigh = 18, units='cm')
#. Number of effect by region
#myd1=read.table("effectbyregion_sim_com.txt", header=T)
#ggplot(myd1, aes(Type, Count, fill = Variant_type)) + geom_col(position = "stack", width = 0.5) +scale_y_log10() +scale_fill_brewer(palette = "Set2")+ 
#scale_fill_manual(values = c("simple" = "burlywood3", "complex" = "darkslategray3"),breaks = c("simple", "complex"),labels = c("Simple", "Complex")) +
#theme_classic()  +   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +coord_flip()+ labs(x = '', y = 'count', fill = 'Variant type') +theme(legend.position = "top",rect = element_rect(fill = "transparent"))
#ggsave("effectregion_allsamples.png" , dpi=1000, width = 5, heigh = 4)