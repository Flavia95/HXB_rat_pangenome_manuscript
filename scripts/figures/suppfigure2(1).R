library(ggplot2)
library(gridExtra)
library(VennDiagram)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(grid)
library(cowplot)
### S Fig. B
myd <- read.table('qc', header = TRUE, sep = '\t')
region_colors <- c("#452049", "#009E73")

# Create the Precision plot
pp <- ggplot(myd, aes(x = Variant.type, y = Precision...., color = Genomic.Region)) +
  geom_violin() +
  geom_jitter(width = 0, height = 0, alpha = 0) +  # Remove points
  geom_text(aes(label = Chromosomes, color = Genomic.Region), vjust = -0.7, size = 3) +
  theme_classic() +
  labs(title = '', x = '', y = 'Precision (%)') +
  scale_color_manual(values = region_colors) +
  theme(legend.position = "none") + ggtitle("Genome wide (SHR/OlaIpcv)- Microvariants")+theme(legend.position = "left",rect = element_rect(fill = "transparent"))+theme(legend.position = "none",rect = element_rect(fill = "transparent")) +
  theme(axis.text.x = element_text( hjust = 1, vjust = 0.5, color = "black"),  # Set the x-axis label color to black
        axis.text.y = element_text(color = "black"))

# Create the Sensitivity plot
ps <- ggplot(myd, aes(x = Variant.type, y = Sensitivity...., color = Genomic.Region)) +
  geom_violin() +
  geom_jitter(width = 0, height = 0, alpha = 0) +  # Remove points
  geom_text(aes(label = Chromosomes, color = Genomic.Region), vjust = -0.7, size = 3) +
  theme_classic() +
  labs(title = '', x = '', y = 'Sensitivity (%)') +
  theme(legend.position = "none") +
  scale_color_manual(values = region_colors)  + theme(axis.text.x = element_text( hjust = 1, vjust = 0.5, color = "black"),  # Set the x-axis label color to black
                                                      axis.text.y = element_text(color = "black"))

# Create the F-Measure plot
pf <- ggplot(myd, aes(x = Variant.type, y = F.Measure...., color = Genomic.Region)) +
  geom_violin() +
  geom_jitter(width = 0, height = 0, alpha = 0) +  # Remove points
  geom_text(aes(label = Chromosomes, color = Genomic.Region), vjust = -0.7, size = 3) +
  theme_classic() +
  labs(title = '', x = '', y = 'F-Measure (%)') +
  theme(legend.position = "none") +
  scale_color_manual(values = region_colors) + theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, color = "black"),  # Set the x-axis label color to black
                                                     axis.text.y = element_text(color = "black"))
# Custom function to extract legend
get_legend <- function(plot) {
  g <- ggplotGrob(plot)
  leg <- g$grobs[[which(g$layout$name == "guide-box")]]
  legend <- leg$grobs[[1]]
  legend
}

# Create a fake plot to extract the legend
fake <- ggplot(myd, aes(x = Variant.type, y = Precision...., color = Genomic.Region)) +
  geom_violin() +
  labs(color = 'Genomic\nRegion') +
  scale_color_manual(values = region_colors)

combined_legend <- get_legend(fake)

# Add a title to the combined legend
combined_legend <- combined_legend 

# Combine the plots with the legend
fb <- grid.arrange(pp, ps, pf, combined_legend, ncol = 4, widths = c(2, 2, 2, 0.5))
ggsave("qc.png", fb, dpi = 500, width = 10, height = 4)


## S Fig C
df <- read.table("allvariants_chr.txt", header = TRUE)
df$counts_K <- df$count / 1000

chromosome_order <- c(paste0("chr", 1:20), "chrX", "chrY")
df$chromosomes <- factor(df$chromosomes, levels = chromosome_order)
microvariant_colors <- c("complex" = "darkslategray4", "simple" = "burlywood2")

fc <- ggplot(df, aes(x = counts_K, y = reorder(methods, -as.numeric(Microvariant)), fill = Microvariant)) +
  geom_col(position = "stack",width = 0.4) +
  facet_wrap(~ chromosomes, ncol = 8, scales = "free") +  # Arrange in 3 rows and 8 columns
  labs(x = "Count (K)", y = " ", fill = "Variant type") +
  scale_fill_manual(values = microvariant_colors, name = "Variant type") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "top"
  ) +
  ggtitle("Chromosome 12 - SHR/Olalpcv")

ggsave("varianttype_allchr.png", height = 6, width = 8, dpi=500)

ListGeno <- 12089
ListPang <- 6682 
ListInt <- 61630

venn.plot <- draw.pairwise.venn(
  ListGeno + ListInt,
  ListPang + ListInt, 
  ListInt,
  category = c("vg-only","JC-only"),
  fill = c("#00B0F6","#E76BF4"),
  lty = "solid",
  cex = 1.5, 
  cat.cex = 3.5,
  cat.dist = 0.0009,
  ext.pos = 30,
  ext.dist = -0.05, 
  ext.length = 0.85,
  inverted = TRUE,
  family = "Times",
  ext.line.lwd = 2,
  col = "black",
  scaling.factor = 0.8
)

ggsave("vennplot_shr.png", venn.plot, width = 8, height = 6, bg='transparent') #figd

##S Fig E
png(file = "pie.png", res = 100, bg = "transparent") #fige
par(family = "Times")
#labels = c("SNPs:2,627", "MNPs", "Others", "Indels")
labels = c(expression(bold("SNPs: 2,627")), "MNPs", "Others", "Indels")
Prop = c(39.3, 5.5, 5.6, 49.6)
colors <- c("#E76BF4", "#f0a6f8", "#f5c3fa", "#ee97f7")

order <- c(4, 1, 2, 3)  # Imposta "SNPs" all'ultimo posto
labels <- labels[order]
Prop <- Prop[order]

pie(Prop, labels = labels, radius = 1, col = colors, cex = 1.5)
label_radii <- 0.7  # Regola questa distanza dal centro
label_positions <- cumsum(Prop) - 0.5 * Prop + 0.1  # Posizioni delle etichette
text(label_radii * cos(label_positions * 2 * pi / sum(Prop)),
     label_radii * sin(label_positions * 2 * pi / sum(Prop)),
     labels = paste(Prop, "%"),
     cex = 1.5)

dev.off()

fa <- rasterGrob(png::readPNG("S2A.png"))
f1 <- plot_grid(fa,fb, labels = c('A','B'), label_size = 12, nrow=2)
fc <- plot_grid(fc,labels = c('C'), label_size = 12, nrow=1)
fd <- rasterGrob(png::readPNG("vennplot_shr.png"))
fe <- rasterGrob(png::readPNG("pie.png"))
ff <- rasterGrob(png::readPNG("SF.png"))
secondrow <- plot_grid(fd,fe, ff,labels = c('D', 'E', 'F'), label_size = 12, nrow=1, ncol=3)
final_plot <- plot_grid(f1, fc, secondrow, nrow = 3)

ggsave('figS2.png',final_plot ,dpi=800, width = 48, heigh = 45, units='cm')