# Load the required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(cowplot) 
library(grid)
library(png) 


# Define a function to create the scatter plot with correlation line and annotations for Spearman's rank correlation
create_scatter_plot <- function(data, trait_name) {
  # Filter out rows with missing values (".") in Value or Variant columns
  data <- data %>%
    filter(Value != "." & Variant != ".")
  
  # Convert Value column to numeric
  data$Value <- as.numeric(data$Value)
  
  # Convert Variant column to a factor with custom labels
  data$Variant <- factor(data$Variant, levels = c("0", "1"))
  
  # Check if there are enough data points for correlation analysis
  if (nrow(data) < 2) {
    cat("Not enough data points for correlation analysis in", trait_name, "\n")
    return(NULL)
  }
  
  result <- cor.test(data$Value, as.numeric(data$Variant), method = "spearman")
  cor_coef <- result$estimate
  p_value <- result$p.value
  
  # Fit a linear model to calculate the slope and intercept for the correlation line
  lm_model <- lm(Value ~ as.numeric(Variant), data = data)
  slope <- coef(lm_model)[2]
  intercept <- coef(lm_model)[1]
  
  # Replace "HRP_" with "GN_" in the title
  plot_title <- sub("HRP_", "GN_", trait_name)
  
  p <- ggplot(data, aes(x = Variant, y = Value)) +
    geom_point() + # Create a scatter plot with points
    geom_abline(intercept = intercept, slope = slope, color = "gray", linetype = "dashed") + # Add the correlation line
    #annotate("text", x = "0", y = Inf, 
             #label = paste("p-value:", format(p_value, scientific = TRUE, digits = 2)), 
            # hjust = 0, vjust = 3, size = 5) + # Add the p-value as an annotation
    labs(title = plot_title, # Use the modified title
         x = "Genotype",
         y = "Value") +
    theme_classic() +
    scale_x_discrete(labels = c("0" = "B", "1" = "H")) +
    theme(
      axis.text.x = element_text(color = "black"),   # Set x-axis tick text color to black
      axis.text.y = element_text(color = "black"),   # Set y-axis tick text color to black
      axis.ticks.x = element_line(color = "black"),  # Set x-axis tick line color to black
      axis.ticks.y = element_line(color = "black")   # Set y-axis tick line color to black
    )
  
  return(p)
}


# Define a function to read and preprocess data from a file
read_and_preprocess <- function(file_path, trait_name) {
  data <- read.table(file_path, header = TRUE, na.strings = ".")
  return(list(data = data, trait_name = trait_name))
}

# Read and preprocess data for each trait
data_files <- c(
  'HRP_10153.txt',
  'HRP_10191.txt',
  'HRP_10006.txt'
)

plots <- list()

for (file_path in data_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  data_info <- read_and_preprocess(file_path, file_name)
  
  if (!is.null(data_info$data)) {
    plots[[file_name]] <- create_scatter_plot(data_info$data, data_info$trait_name)
  }
}

plot_10153 <- NULL
plot_10191 <- NULL
plot_10006 <- NULL

for (file_path in data_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  data_info <- read_and_preprocess(file_path, file_name)
  
  if (!is.null(data_info$data)) {
    if (file_name == 'HRP_10153') {
      plot_10153 <- create_scatter_plot(data_info$data, data_info$trait_name)
    } else if (file_name == 'HRP_10191') {
      plot_10191 <- create_scatter_plot(data_info$data, data_info$trait_name)
    } else if (file_name == 'HRP_10006') {
      plot_10006 <- create_scatter_plot(data_info$data, data_info$trait_name)
    }
  }
}

# Create the plot for 'HRP_10191' with modified title aligned to the left
fe <- plot_10191 +
  labs(
    x = "Genotype at chr12:18797475",
    y = "Chromogranin concentration\n(ng/mL)"
  ) +
  theme(
    plot.title = element_text(color = "#cc0c00", size = 10, hjust = 0), # Align title to the left
    plot.subtitle = element_text(hjust = 0.5),plot.tag = element_text(face = "bold")
  ) +labs(tag = "E",size=12)+
  ggtitle("CLASS: CNS/Gene expression\nPHENOTYPE: Chromogranin\n in hippocampus\nr: 0.12 p-value < 0.05")

fd <- plot_10006 +
  labs(
    x = "Genotype at chr12:18797475",
    y = "Insulin Concentrations\n(mmol/L)"
  ) +
  theme(
    plot.title = element_text(color = "#674ea7", size = 10, hjust = 0), # Align title to the left
    plot.subtitle = element_text(hjust = 0.5),plot.tag = element_text(face = "bold")
  ) +labs(tag = "D",size=12)+
  ggtitle("CLASS: Insulin/Glucose Ratio\nPHENOTYPE: Insulin Concentrations\nr: -0.83 p-value < 0.05")

fc <- plot_10153 +
  labs(
    x = "Genotype at chr12:4347739",
    y = "Blood glucose \n(mmol/L)"
  ) +
  theme(
    plot.title = element_text(color = "#517524", size = 10, hjust = 0), # Align title to the left
    plot.subtitle = element_text(hjust = 0.5),plot.tag = element_text(face = "bold")
  ) +  labs(tag = "C",size=12)+
  ggtitle("CLASS: Blood Chemistry-Glucose\nPHENOTYPE: Blood Glucose\nr: 0.65 p-value < 0.05")

ta=expression(paste('Chr12:4347739 (', italic('Zfp958l1'), ') ')) 
pa <- rasterGrob(png::readPNG("/home/flavia/Desktop/4347739_edit.png"))
fa <- ggplot() + annotation_custom(pa, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +theme_bw() +  # Or any other theme you prefer
  theme(
    plot.background = element_rect(fill = "white"),
    panel.border = element_blank()  # Remove panel borders
  ) +
   labs(title= ta, x='Phenosome', y='LRS', tag = "A",size=12) +
   theme(plot.margin = unit(c(0,0,0,0), "null"),plot.tag = element_text(face = "bold"))

tb=expression(paste('Chr12:18797475 (', italic('LOC685157'), ')'))
pb <- rasterGrob(png::readPNG("/home/flavia/Desktop/18797475_edit.png"))
fb <- ggplot() + annotation_custom(pb, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +theme_bw() +  # Or any other theme you prefer
  theme(
    plot.background = element_rect(fill = "white"),
    panel.border = element_blank()  # Remove panel borders
  )+labs(title= tb, x='Phenosome', y='LRS', tag = "B",size=12) +
   theme(plot.margin = unit(c(0,0,0,0), "null"),plot.tag = element_text(face = "bold"))

ff <- rasterGrob(png::readPNG("ff.png"))

figff_plot <- ggplot() + 
  annotation_custom(ff, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+theme_bw() +  # Or any other theme you prefer
  theme(
    plot.background = element_rect(fill = "white"),
    panel.border = element_blank(),plot.tag = element_text(face = "bold")  # Remove panel borders
  ) + labs(tag='F')


fig3<-grid.arrange(arrangeGrob(fa, fb , ncol = 2, bottom=''), 
             arrangeGrob(fc, fd, fe, ncol = 3), 
             arrangeGrob(figff_plot, ncol = 1), 
             heights = c(1.5, 1, 1))



ggsave('fig3.tiff', fig3, height=10, width=10, dpi=600)
