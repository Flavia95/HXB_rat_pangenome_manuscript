#Trait metadata - https://genenetwork.org/api/v_pre1/traits/HXBBXHPublish.csv (or .json)
#Sample data - https://genenetwork.org/api/v_pre1/sample_data/HXBBXHPublish.csv (or .json)
require(BXDtools)
#devtools::install_github("dannyarends/bxdtools")

my_marker="chr12_4347739"
annotation_file_loc="./HXBBXHPublish_traits.csv"
phenotypes_file_loc="./HXBBXHPublish_sample_data.csv"

############GENERATING GENOTYPE FILE
# reload this edited function with appropriate filename
download.BXD.genotypes <- function() {
  bxd.geno <- read.table("HXB.geno", sep = "\t",
                         header = TRUE,row.names=1, colClasses = "character", stringsAsFactors = FALSE)
  return(bxd.geno)
}

#reload this edited function. mantain only samples columns
only.BXD.genotypes <- function(bxd.geno){
  genotypes <- bxd.geno[,4:ncol(bxd.geno)]
  return(genotypes)
}

bxd.geno <- download.BXD.genotypes()
bxd.genotypes <- only.BXD.genotypes(bxd.geno)

# Return a bxd.genotypes recoded as B = -1, A = 1, and U = NA
#This function makes the whole thing numeric, and seems to be absolutely necessary
recode.BXD.genotypes <- function(bxd.genotypes, coding = rbind(c("A", "1"), c("B", "-1"), c("U", "NA"))) {
  for(x in 1:nrow(coding)){
    bxd.genotypes[bxd.genotypes == coding[x, 1]] <- coding[x, 2]
  }
  bxd.numeric <- apply(bxd.genotypes, 2, function(x){
    return(as.numeric(x))
  })
  rownames(bxd.numeric) <- rownames(bxd.genotypes)
  attr(bxd.numeric, "map") <- attr(bxd.genotypes, "map")
  return(bxd.numeric)
}

#Recode from letters to numbers
bxd.genotypes <- recode.BXD.genotypes(bxd.genotypes)
table(unlist(bxd.genotypes), useNA="always")

########################GENERATING PHENOTYPE FILE
###reload edited function
download.BXD.phenotypes <- function(verbose = FALSE){
#tmp.annot <- tempfile()
#tmp.data <- tempfile()
#download.file("https://genenetwork.org/api/v_pre1/traits/HXBBXHPublish.csv", tmp.annot)
#download.file("https://genenetwork.org/api/v_pre1/sample_data/HXBBXHPublish.csv", tmp.data)

# read in the annotation and phenotype data
annotation <- read.csv(annotation_file_loc, sep = ",", header=TRUE)
phenotypes <- read.csv(phenotypes_file_loc, sep = ",", header=TRUE, na.strings=c("x"), row.names=1)

# format the annotation and phenotype data
annotation <- annotation[!duplicated(annotation),]
row.names(annotation) <- paste0("GN_", annotation[,"Id"])

#this renames the columns of the phenotypes dataframe, replacing the prefix "HRP_" with "GN_"
colnames(phenotypes) <- gsub("^HRP_", "GN_", colnames(phenotypes))

#removes any rows or columns that are present in one dataframe but not the other.
phenotypes <- phenotypes[, which(colnames(phenotypes) %in% rownames(annotation))]
annotation <- annotation[which(rownames(annotation) %in% colnames(phenotypes)),]

#this transposes the phenotypes dataframe so that the phenotypes are in rows and samples are in columns
phenotypes <- t(phenotypes[, rownames(annotation)])

#create a new dataframe called phenotype.descriptions that contains the phenotype names, descriptions, and years
#The column and row names are updated to match the annotation dataframe
phenotype.descriptions <- annotation[,c("Id", "Description", "Year")]
colnames(phenotype.descriptions) <- c("name", "description", "year")
rownames(phenotype.descriptions) <- rownames(annotation)

#this updates the description column in the phenotype.descriptions dataframe to remove a specific string "Original post publication description: " that appears at the beginning of each description.
phenotype.descriptions[,"description"] <- gsub("Original post publication description: ", "", phenotype.descriptions[,"description"])
#weird <- rownames(phenotype.descriptions)[grep("un-named trait", tolower(phenotype.descriptions[,"description"]))]
#phenotype.descriptions <- phenotype.descriptions[-which(rownames(phenotype.descriptions) %in% weird),]
#phenotypes <- phenotypes[-which(rownames(phenotypes) %in% weird),]

# add annotation to phenotype data
attr(phenotypes, "annotation") <- phenotype.descriptions
return(phenotypes)
}

#load phenotype data
bxd.pheno<-download.BXD.phenotypes()

#Modify bxd.genotypes column names so match bxd.pheno file
colnames(bxd.genotypes) <- gsub("BNLxCub", "BN-Lx/Cub", colnames(bxd.genotypes))
colnames(bxd.genotypes) <- gsub("SHROlaIpcv", "SHR/OlaIpcv", colnames(bxd.genotypes))
##remove any columns that are present in genotypes but not in phenotypes
ncol(bxd.genotypes)
bxd.genotypes <- bxd.genotypes[, which(colnames(bxd.genotypes) %in% colnames(bxd.pheno))]
ncol(bxd.genotypes)

#run matrix function
bxd.phenotypes<-as.phenotype.matrix(bxd.genotypes,bxd.pheno)
# Don't need to add a phenotype class (extracted from the description column) to each phenotype, as this is contained within the only.phenosomes function below
#save out annotation attribute
aa<-attr(bxd.phenotypes, "annotation")

#subset
pheno_samples<-colnames(bxd.genotypes)
bxd.phenotypes<-bxd.phenotypes[,pheno_samples]

#####Following code deals with marker errors and AMOVA perfect fit warning (no variation)
#print names of any phenotype (to remove) that has 1 or fewer unique phenotype values, or 2 phenotype values represented with only a single genotype each
marker = my_marker
pheno_row_remove<-c()
nrow(bxd.phenotypes)
#for every row/phenotype
for (r in 1:(nrow(bxd.phenotypes))) {
  y<-row.names(bxd.phenotypes)[r] #save row/phenotype name
  vals<-unique(na.omit(bxd.phenotypes[r,])) #list unique phenotype values
  #1.
  if (length(vals)<=1) {
    print(paste0("only 0 or 1 phenotype value so removing row: ",r))
    pheno_row_remove<-append(pheno_row_remove,y) } 
  #2.
  else if (length(vals)==2) {
    #only 2 phenotype values so checking there is some genotype variation
    for (i in vals) {
      samples<-colnames(bxd.phenotypes)[which(bxd.phenotypes[r,]==i)] #select all the samples that have this same value for a single row/phenotype
      subset_genotype<-bxd.genotypes[marker,samples] #select the same samples from the genotype file (all rows)
      #if number of different genotypes in these samples > 1, break out for loop
      if (length(unique(na.omit(subset_genotype))) > 1) {
        break }
    }
    #if have cycled through both unique phenotype values having only 1 genotype, add row/phenotype name to 'remove' list 
    #else at least one value has more than 1 genotype, so keep row/phenotype and move on to test next one
    if (length(unique(na.omit(subset_genotype))) ==1) {
      print(paste0("only 2 phenotype values and no genotype variation so removing row: ",r))
      pheno_row_remove<-append(pheno_row_remove,y) } else { }
  } #closing else if loop providing values was = 2
  #3.
  else {} #closing if, else if, else loop
} #closing row/phenosome for loop
length(pheno_row_remove)
#remove these phenotypes
bxd.phenotypes<-bxd.phenotypes[!(row.names(bxd.phenotypes) %in% pheno_row_remove),]

#Add annotation attribute back in
attr(bxd.phenotypes, "annotation")<-aa


# RELOAD only.phenosomes function to remove same phenotypes from phenotype.description object
# (This function gets the phenotypes from phenosomes with more than a minimum number of phenotypes)
only.phenosomes <- function(bxd.phenotypes, minimum = 1, splitCNS = TRUE, verbose = FALSE) {
  bxd.phenotypes <- phenotype.add.class(bxd.phenotypes, splitCNS = splitCNS)
  phenotype.descriptions <- attr(bxd.phenotypes, "annotation")
  #######ADDED SECTION
  #if did cut some phenotypes above (due to too many NA values) then do the following loop
  if (length(pheno_row_remove)>0) { 
  print("cutting same rows from phenotype.descriptions object as cut from phenotypes file. Row count before and after:")
  print(nrow(phenotype.descriptions))
  phenotype.descriptions<-phenotype.descriptions[setdiff(rownames(phenotype.descriptions),pheno_row_remove),]
  print(nrow(phenotype.descriptions)) } else { 
  }
  ########
  classes <- phenotype.descriptions[, "class"]
  phenosomes <- names(which(table(classes) > minimum))
  useable.phenotypes <- which(classes %in% phenosomes)
  if(verbose) cat("Useable phenotypes:", length(useable.phenotypes), "before:", nrow(bxd.phenotypes), "\n")
  bxd.phenosomes <- bxd.phenotypes[useable.phenotypes, ]
  phenotype.descriptions <- phenotype.descriptions[useable.phenotypes, ]
  
  ordering <- sort(phenotype.descriptions[,"class"], index.return = TRUE)
  bxd.phenosomes <- bxd.phenosomes[ordering$ix, ]
  attr(bxd.phenosomes, "annotation") <- phenotype.descriptions[ordering$ix,]
  return(bxd.phenosomes)
}

#Group phenosomes together
bxd.phenosomes <- only.phenosomes(bxd.phenotypes, verbose=TRUE)
#Run PHEWAS with my marker
scores <- do.BXD.phewas(bxd.genotypes, bxd.phenosomes, marker = my_marker, LRS=TRUE)
#Plot function
plot.phewas <- function(scores, bxd.phenosomes, do.sort = FALSE, decreasing = FALSE,
                        main = "BxD PheWAS results", pch = 19, cex = 0.4, significance = 16, type = "h",
                        colorSeed = 1, colorRange = c("darkslateblue", "hotpink1", "forestgreen", "orange", "black", "firebrick1"), ...) {
  classes <- attr(bxd.phenosomes, "annotation")[,"class"]
  names(classes) <- rownames(bxd.phenosomes)
  marker <- attr(scores, "marker")
  ylab <- "-log10(P)"
  if(attr(scores, "LRS")) ylab = "LRS"
  if(do.sort) {
    ordering <- c()
    for(pheclass in unique(classes)){
      indexes <- which(classes == pheclass)
      # Use Radix sorting and put the NA scores last
      sorting <- sort(scores[indexes], decreasing = decreasing, index.return = TRUE, method = "quick", na.last=NA)
      ordering <- c(ordering, names(sorting$x))
    }
    scores <- scores[ordering]
    classes.inOrder <- classes[ordering]
  }else{
    classes.inOrder <- classes
  }
  #return(classes.inOrder)
  # Set the colorSeed to NA, to cycle colors to groups each plot
  if(!is.na(colorSeed)) set.seed(colorSeed)
  # Generate colors from colorRange and mix them up for beter visability
  mcolors <- sample(colorRampPalette(colorRange)(length(unique(classes))))
  names(mcolors) <- unique(classes)
  ######ADDED LINE:
  par(oma=c(0,0,0,10))
  
  plot(scores, col = mcolors[classes.inOrder], pch = pch, 
       cex = cex, type = type, xaxt='n', xlab = "Phenosome", ylab = ylab, main=main)
  #legend("bottomright", xpd=NA, inset=c(-0.6,0), paste0(names(mcolors), " (N = ", table(classes)[names(mcolors)], ")"), 
  legend("topright", xpd=NA, inset=c(-0.5,0), paste0(names(mcolors), " (N = ", table(classes)[names(mcolors)], ")"), 
         col = mcolors, pch = pch, cex = cex, ncol = 1)
  abline(h = significance, lty=2)
  results <- cbind(classes.inOrder, scores)
  significant <- results[which(as.numeric(results[,2]) > significance),,drop=FALSE]
  legend("topleft", paste0(rownames(significant), " ", significant[,1], " ", ylab, "=", round(as.numeric(significant[,2]),1)), cex = cex, fill = mcolors[significant[,1]])
  #legend("topright", xpd=NA, inset=c(-0.69,0), paste0(rownames(significant), " ", significant[,1], " ", ylab, "=", round(as.numeric(significant[,2]),1)), cex = cex, fill = mcolors[significant[,1]])
  invisible(results)
}

pdf("myplot.pdf", width=7, height=7) # Open PDF file for writing
plot.phewas(scores, bxd.phenosomes, TRUE, main = paste0("BXD PheWAS at ", my_marker))
dev.off() 