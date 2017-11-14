# Online script to generate cohesion metrics for a set of samples 
# CMH 29Nov16; cherren@wisc.edu
# Editted for Danczak et al.
# RED 2017; danczak.6@osu.edu

# User instructions: read in a sample table (in absolute or relative abundance) as object "b". Run the entire script, and the 4 vectors (2 of connectedness and 2 of cohesion) are generated for each sample at the end.
# Parameters that can be adjusted include pers.cutoff (persistence cutoff for retaining taxa in analysis) and iter (number of iterations for the null model)


well = c("Greene 0.2") # Select which dataset to examine
iter <- 200 # Decide the number of iterations to run for each taxon. (>= 200 is recommended)
pers.cutoff <- 0.10 # Choose a persistence cutoff (min. fraction of taxon presence across samples) for retaining taxa in the analysis

#--------------------- Do not edit below this line ---------------------#

library(picante) # Loaded for rarefy (also in vegan)
library(doSNOW) # Used for multithreading

cl <- makeCluster(10)
registerDoSNOW(cl)

zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)
} #find the number of zeroes in a vector

neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  return(n.mean)
} #create function that averages only negative values in a vector


pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  return(p.mean)
} #create function that averages only positive values in a vector

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol # Need this to fix the fact that you can't exactly have half of a read created from my averaging step

###################################################################

# Read in dataset
## Data should be in a matrix where each row is a sample.
setwd("/home/rdanczak/EPA_16S_Analysis/Cohesion")
b <- read.table("Raw_EPA_no_pynast.txt", header = T, row.names = 1)
factors = read.csv("OTU_table_factors.csv")

row.names(factors) = factors[,1] # This and the next three lines are just processing the data into more useful formats
factors = factors[,-1]
factors = factors[-which(row.names(factors) %in% c("YB04","BL01","BL02")),]

########################################
#### Averaging my duplicate results ####
########################################

x = as.character(unique(factors[,1])) # Setting the unique factors to a variable
data.avg = matrix(nrow = length(b[,1]), ncol = length(x)) # Creating an empty matrix to spit the averaged data into

for(i in 1:length(x)){
  w = which(factors[,1] %in% x[i])
  if(length(w)==1){
    data.avg[,i] = b[,w]
  } else {
    data.avg[,i] = apply(b[,w], 1, mean)
  }
} # This loop goes through by column and averages the two columns with matching factors

colnames(data.avg) = x
row.names(data.avg) = row.names(b)

b = data.avg
b = as.data.frame(t(b))

# Selecting the well subset
b = b[grep(well, row.names(b)),]

# Need to round all 0.5 values up as 0.5 means that sequences were present
b = as.matrix(b)
b[which(is.wholenumber(b) == FALSE)] = b[which(is.wholenumber(b) == FALSE)]+0.5

# Rarefying the community to the lowest sequence count
b = as.data.frame(rrarefy(b, min(rowSums(b))))

# Removing OTUs which have zero abundance after rarefaction
zero.members = apply(b, 2, max) # Finding the maximum number of sequence per OTU
b = b[,-which(zero.members==0)] # Removing those OTUs which have 0 sequences after rarefaction

# Suggested steps to re-format data. At the end of these steps, the data should be in a matrix "c" where there are no empty samples or blank taxon columns. 
c <- as.matrix(b)
c <- c[rowSums(c) > 0, colSums(c) > 0]

# Optionally re-order dataset to be in chronological order. Change date format for your data. 
#c <- c[order(as.Date(rownames(c), format = "%m/%d/%Y")), ]

# Save total number of individuals in each sample in the original matrix. This will be 1 if data are in relative abundance, but not if matrix c is count data
rowsums.orig <- rowSums(c)

# Based on persistence cutoff, define a cutoff for the number of zeroes allowed in a taxon's distribution
zero.cutoff <- ceiling(pers.cutoff * dim(c)[1])

# Remove taxa that are below the persistence cutoff
d <- c[ , apply(c, 2, zero) < (dim(c)[1]-zero.cutoff) ]
# Remove any samples that no longer have any individuals, due to removing taxa
d <- d[rowSums(d) > 0, ]

# Create relative abundance matrix.  
rel.d <- d / rowsums.orig

# Create observed correlation matrix
cor.mat.true <- cor(rel.d)

# Create vector to hold median otu-otu correlations for initial otu
med.tax.cors <- vector()

# Run this loop with each taxon as the focal taxon
for(which.taxon in 1:dim(rel.d)[2]){
  
  #create vector to hold correlations from every permutation for each single otu
  ## perm.cor.vec.mat stands for permuted correlations vector matrix
  perm.cor.vec.mat <- vector()
  
  perm.cor.vec.mat = foreach(i = 1:iter, .verbose = T, .combine = cbind) %dopar% {
    # Permute the abundances of each taxon
    perm.rel.d <- apply(rel.d, 2, sample)
    row.names(perm.rel.d) = row.names(rel.d)
    
    # Replace focal taxon abundance with original distribution
    perm.rel.d[,which.taxon] <- rel.d[ ,which.taxon]
    
    # Calculate correlation matrix of permuted matrix
    cor.mat.null <- cor(perm.rel.d)
    
    # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
    cor.mat.null[,which.taxon]
  }
  # Save the median correlations between the focal taxon and all other taxa  
  med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
  
  # For large datasets, this can be helpful to know how long this loop will run
  print(sprintf("Finished taxon number %s at %s", which.taxon, date()))
}

# Save observed minus expected correlations
obs.exp.cors.mat <- cor.mat.true - med.tax.cors
diag(obs.exp.cors.mat) <- 0

#### 
#### Produce desired vectors of connectedness and cohesion 

# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)

# Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

####
#### Combine vectors into one list and print 
output1 <- cbind(connectedness.neg, connectedness.pos)
output2 <- cbind(cohesion.neg, cohesion.pos)
colnames(output1) <- c("Negative Connectedness", "Positive Connectedness")
row.names(output1) <- colnames(rel.d)
colnames(output2) <- c("Negative Cohesion", "Positive Cohesion")

write.csv(output1, sprintf("%s_Connectedness.csv", well), quote = F)
write.csv(output2, sprintf("%s_Cohesion.csv", well), quote = F)

stopCluster(cl)
