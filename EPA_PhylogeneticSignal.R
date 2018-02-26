# Separated this from the Ecological Modeling script because it was a bit too much for the computer to handle
setwd("/home/rdanczak/EPA_16S_Analysis/No_Pynast_bNTI")
# setwd("/Users/danczak.6/Documents/Wilkins Lab/Ohio EPA Project/Sequencing Data/16S Data/")

library(vegan)
library(picante)

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol # Need this to fix the fact that you can't exactly have half of a read created from my averaging step

### Phylogenetic Signal ###

tree = read.tree("Raw_EPA_OTU_Tree.tre")
data = read.delim("Raw_EPA_no_pynast.txt", row.names = 1)
factors = read.csv("OTU_table_factors.csv", row.names = 1) # Importing factors to make the averaging of replicates possible
table = read.csv("ORP_Table.csv", row.names = 1)

# Doing some data cleaning
data = data[,-which(colnames(data) %in% c("YB04","BL01","BL02"))] 
factors = factors[-which(row.names(factors) %in% c("YB04","BL01","BL02")),]

### Averaging my duplicate results
x = as.character(unique(factors[,1])) # Setting the unique factors to a variable
data.avg = matrix(nrow = length(data[,1]), ncol = length(x)) # Creating an empty matrix to spit the averaged data into

for(i in 1:length(x)){
  w = which(factors[,1] %in% x[i])
  if(length(w)==1){
    data.avg[,i] = data[,w]
  } else {
    data.avg[,i] = apply(data[,w], 1, mean)
  }
} # This loop goes through by column and averages the two columns with matching factors

colnames(data.avg) = x
row.names(data.avg) = row.names(data)

data = as.data.frame(t(data.avg))

rm('data.avg')

### Removing the 0.1um filter data from the dataset
w = grep(" 0.1", row.names(data))
data = data[-w,]

table = as.matrix(table)
table = as.data.frame(table[-w,])

### Need to rarefy or else cophenetic isn't happy
data = as.matrix(data)
data[which(is.wholenumber(data) == FALSE)] = data[which(is.wholenumber(data) == FALSE)]+0.5 # Rounding the 0.5s up as half a sequence is a sequence - it was present at least once

data = as.data.frame(rrarefy(data, 8000))

# Determining OTUs which have zero abundance after rarefaction
zero.members = apply(data, 2, max) # Finding the maximum number of sequence per OTU
data = data[,-which(zero.members==0)] # Removing those OTUs

# Matching the tree to the newly rarefied OTU dataset
phylo = match.phylo.data(tree,t(data))

data = t(phylo$data)
tree = phylo$phy

print("Finished rarefying...hopefully")
### Figuring out the abundance-weighted water table mean, this part was picked up from Stegen 2012
r = apply(data, 2, function(x) x*table) # Weighting ORP by abudances
x = matrix(unlist(r), ncol = length(r), byrow = F)
y = apply(x, 2, mean, na.rm = T) # Determing the mean of the column (i.e. the mean of abundance-weighted depths)
n = abs(outer(y, y, '-')) # Pairwise comparison of the abundance-weighted depths

### Reading the tree for phylogenetic distance information
print("Starting cophenetic")
x = cophenetic.phylo(tree)

### Trying to plot the Mantel Correlogram
print("We've gotten this far! Starting Mantel!")
q = mantel.correlog(n,x)

pdf("PhyloSignal.pdf")
plot(q)
dev.off()

print("Done")