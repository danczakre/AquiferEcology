# Writing a script to both rarefy and calculate bNTI for OTU tables
# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/bNTI_Local_Machine.r; Stegen et al, 2013
# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r; Stegen et al, 2013
# RED 2017; danczak.6@osu.edu


local = 1 # Whether on the server or not
micron = 1 # Whether or not to include the 0.1um fraction (1 is include, 0 is exclude)

library(vegan)
library(picante)
library(doMC)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol # Need this to fix the fact that you can't exactly have half of a read created from my averaging step

rc.reps = 9999
registerDoMC(cores = 1)

###################################
#### Data Loading and cleaning ####
###################################
if(local==1){
  setwd("/Users/danczak.6/Documents/Wilkins Lab/Ohio EPA Project/Sequencing Data/16S Data/") # Working directory for 16S comm on my local computer
  data = read.delim("Raw_EPA_no_pynast.txt") # Importing the raw OTU data
  tree = read.tree("rep_set.tre")
  factors = read.csv("OTU_table_factors.csv") # Importing factors to make the averaging of replicates possible
  
} else {
  setwd("/home/rdanczak/EPA_16S_Analysis/No_Pynast_bNTI") # Working directory for 16S comm on the server
  data = read.delim("Raw_EPA_no_pynast.txt") # Importing the raw OTU data
  tree = read.tree("Raw_EPA_OTU_Tree.tre") # Importing the raw OTU tree
  factors = read.csv("OTU_table_factors.csv") # Importing factors to make the averaging of replicates possible
  
}

row.names(data) = data[,1] # Setting the the first column to the row.names and then removing them
data = data[,-1]

data = data[,-which(colnames(data) %in% c("YB04","BL01","BL02"))] # Removing YB04 becuase it was a failed sequencing run

# Averaging the replicates for our sequencing endeavors
row.names(factors) = factors[,1] # This and the next three lines are just processing the data into more useful formats
factors = factors[,-1]
factors = factors[-which(row.names(factors) %in% c("YB04","BL01","BL02")),]

########################################
#### Averaging my duplicate results ####
########################################

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

############################
#### Running Raup-Crick ####
############################

# Removing the 0.1um filters from the dataset
if(micron==0){
  data.avg = data.avg[,-grep("0.1", colnames(data.avg))]
  
  zero.members = apply(data.avg, 1, max) # Finding the maximum number of sequence per OTU
  data.avg = data.avg[-which(zero.members==0),]
}

# Setting up the data to mesh with the Stegen et al. code
spXsite = t(data.avg)
spXsite[which(is.wholenumber(spXsite) == FALSE)] = spXsite[which(is.wholenumber(spXsite) == FALSE)]+0.5 # Rounding the 0.5s up as half a sequence is a sequence - it was present at least once

# Count number of sites and total species richness across all plots (gamma)
n_sites = nrow(spXsite)
gamma = ncol(spXsite)

# Build a site by site matrix for the results, with the names of the sites in the row and col names:
results = matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))

# Make the spXsite matrix into a new, pres/abs. matrix:
spXsite.inc = ceiling(spXsite/max(spXsite))

# Create an occurrence vector- used to give more weight to widely distributed species in the null model
occur = apply(spXsite.inc, MARGIN=2, FUN=sum)

# Create an abundance vector- used to give more weight to abundant species in the second step of the null model
abundance = apply(spXsite, MARGIN=2, FUN=sum)

# Loops through every pairwise comparison, generating null results

for(null.one in 1:(nrow(spXsite)-1)){
  for(null.two in (null.one+1):nrow(spXsite)){
    
    null_bray_curtis<-NULL
    null_bray_curtis = foreach(i=1:rc.reps, .packages = c("vegan","picante")) %dopar% {
      
      # Generates two empty communities of size gamma
      com1<-rep(0,gamma)
      com2<-rep(0,gamma)
      
      # Add observed number of species to com1, weighting by species occurrence frequencies
      com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
      com1.samp.sp = sample(which(com1>0), (sum(spXsite[null.one,])-sum(com1)), replace=TRUE, prob=abundance[which(com1>0)]);
      com1.samp.sp = cbind(com1.samp.sp,1);
      com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2], com1.samp.sp[,1], FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
      com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts));
      com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts;
      #sum(com1) - sum(spXsite[null.one,]); # This should be zero if everything worked properly
      rm('com1.samp.sp','com1.sp.counts');			
      
      # Again for com2
      com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
      com2.samp.sp = sample(which(com2>0), (sum(spXsite[null.two,])-sum(com2)), replace=TRUE, prob=abundance[which(com2>0)]);
      com2.samp.sp = cbind(com2.samp.sp,1);
      com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2], com2.samp.sp[,1], FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
      com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts));
      com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts;
      # sum(com2) - sum(spXsite[null.two,]); # This should be zero if everything worked properly
      rm('com2.samp.sp','com2.sp.counts');
      
      null.spXsite = rbind(com1,com2); # Null.spXsite
      
      # Calculates the null Bray-Curtis
      null_bray_curtis[i] = vegdist(null.spXsite, method='bray');
      
    }; # End of the null loop
    
    # Unlisting the parallel list
    null_bray_curtis = unlist(null_bray_curtis)
    
    # Calculates the observed Bray-Curtis
    obs.bray = vegdist(spXsite[c(null.one,null.two),], method='bray');
    
    # How many null observations is the observed value tied with?
    num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
    
    # How many null values are smaller than the observed *dissimilarity*?
    num_less_than_in_null = sum(null_bray_curtis<obs.bray);
    
    rc = ((num_less_than_in_null + (num_exact_matching_in_null)/2)/rc.reps) # This variation of rc splits ties
    
    rc = (rc-.5)*2 # Adjusts the range of the  Raup-Crick caclulation to -1 to 1
    
    results[null.two,null.one] = round(rc,digits=2); # Stores rc into the results matrix
    
    print(c(null.one,null.two,date())); # Keeps track of position
    
  }; # End of inner loop
  
}; # End of outer loop

rc.results = as.dist(results) # Converts results into a distance matrix

write.csv(as.matrix(rc.results), "RC_BC_MC_Results.csv", quote = F)

rm('spXsite.inc')

####################################
#### Beginning the bNTI Process ####
####################################

# Transposing data because the rrarefy command doesn't like samples in columns
# alt.data = as.data.frame(t(data.avg))
# alt.data[which(is.wholenumber(alt.data) == FALSE)] = alt.data[which(is.wholenumber(alt.data) == FALSE)]+0.5

# Rarefying the data because bNTI calculations don't seem to tolerate too many OTUs over 20000; specifically, the error is produced in cophenetic
data.rare = as.data.frame(rrarefy(spXsite, 8000)) # Rarefying to the sample with the lowest sequence number

# Determining OTUs which have zero abundance
zero.members = apply(data.rare, 2, max) # Finding the maximum number of sequence per OTU
data.rare = data.rare[,-which(zero.members==0)] # Removing those OTUs which have 0 sequences after rarefaction

# Matching the tree to the newly rarefied OTU dataset
phylo = match.phylo.data(tree,t(data.rare))

# Calculating bMNTD for my samples
bMNTD = as.matrix(comdistnt(t(phylo$data), cophenetic(phylo$phy), abundance.weighted = T))

# Calculating the bMNTD for 999 random distributions
bMNTD.rand = array(c(-999),dim=c(ncol(phylo$data),ncol(phylo$data),999)) # Creating 999 'dummy' matrices to put random results into

for(i in 1:999){
  bMNTD.rand[,,i] = as.matrix(comdistnt(t(phylo$data),taxaShuffle(cophenetic(phylo$phy)),abundance.weighted = T,exclude.conspecifics = F))
  print(c(date(),i)) # Measuring how far the loop has gotten
} # Performing the calculations on using the OTU table but with randomized taxonomic affiliations

# Calculating the bNTI between these randomized communities and the empirical dataset
bNTI = matrix(c(NA),nrow=ncol(phylo$data),ncol=ncol(phylo$data))

for(i in 1:(ncol(phylo$data)-1)){ 
  for(j in (i+1):ncol(phylo$data)){
    m = bMNTD.rand[j,i,] # Just setting all the randomizations for a given comparison to a matrix
    bNTI[j,i] = ((bMNTD[j,i]-mean(m))/sd(m)) # The bNTI calculation
  }
}

rownames(bNTI) = colnames(phylo$data)
colnames(bNTI) = colnames(phylo$data)

write.csv(bNTI, "MC_EPA_bNTI.csv", quote = F)