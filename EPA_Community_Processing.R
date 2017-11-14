# Processing the Merged EPA OTU Table
# The reason for the reupdate was to generate static color vectors so my figures always come out the same - I was getting too easily confused
# RED 2017; danczak.6@osu.edu


#---------------- Switches ----------------#
Var.Num = 5 # Change if the number of variables/factors in the first few columns changes
Var = 4 # Switch controlling which column the variable of interest exists in
PW = 1 # Controls whether I want to remove the Licking-PW from the analysis as it is an outlier (1 = remove, 0 = keep)
WGCNA = 0 # Controls whether I want to generate a WGCNA analysis for network analysis (1 = yes, 0 = no)
switch.otu = 0 # Controls whether or not I want to work with the raw OTU table (1 = yes, 0 = no)
rfy = 0 # Controls rarefaction for Faith's PD (1 = rarefy, 0 = don't rarefy)
SIMP = FALSE

#-------------- Loading packages --------------#
library(vegan) # For general ecology functions
library(GUniFrac) # For UniFrac
library(picante) # For bNTI and Faith's PD
library(ade4) # For advanced tree building
library(Hmisc) # I don't remember...
library(reshape2) # For editing data to better fit the ggplot2 inputs
library(WGCNA) # Looking at networking potential
library(yarrr) # For pirate plots - sounds silly, but is a slightly more advanced boxplot

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol # Need this to fix the fact that you can't exactly have half of a read created from my averaging step

############################################  Do not edit below this line ##############################################################

#######################################
#### Importing and setting up data ####
#######################################

setwd("~/Documents/Wilkins Lab/Ohio EPA Project/Sequencing Data/16S Data") # Working directory for 16S comm
data = read.csv("Family_Level_OTU_edit.csv")
otu = read.delim("Raw_EPA_no_pynast.txt") # Loading in the raw OTU table for UniFrac analyses
geo = read.csv("EPA_Geochem_16S.csv")
# otu.meta=read.csv("Merged EPA L5 OTU Table (Feb 22 2016).csv")
tree = read.tree("rep_set.tre") # Loading in the tree for UniFrac analyses

factors = t(data[1:Var.Num,])
colnames(factors) = factors[1,]
factors = as.data.frame(factors[-1,])

data = data[-1:-Var.Num,]
data[,-1] = apply(data[,-1], 2, function(x) as.numeric(as.character(x)))
data = aggregate(. ~ OTU.ID, data = data, FUN = sum) # Need to aggregate the newly assigned taxonomies (I replaced OD1 things with Parcubacteria)

row.names(data) = data[,1] # Setting row names
data = data[,-1] # Removing the vector in which the row names were stored

# data = data[,-which(colnames(data) %in% c("YB04","BL01","BL02"))] # Pulling out data that I am not interested in (YB04 because it failed and the blanks until I figure out how to deal with them)

data = t(data) # Transposing data so it is actually useful
# factors = as.data.frame(data[,1:Var.Num]) # Hardcoded to the number of factors
# data = data[,-1:-Var.Num] # Removes the factor lengths

data = apply(data, 2, function(x) as.numeric(as.character(x))) # Reconverts the data to numbers from characters
row.names(data) <- row.names(factors)

######################################
#### Averaging duplicated samples ####
######################################

x = unique(factors[,1])
data.avg = matrix(nrow = length(x), ncol = length(data[1,]))

# Generating the averaged Family level table
for(i in 1:length(x)){
  w = which(factors[,1] %in% x[i])
  if(length(w)==1){
    data.avg[i,] = data[w,]
  } else{
  data.avg[i,] = apply(data[w,], 2, mean)
  }
}

row.names(data.avg) = x
colnames(data.avg) = colnames(data)

# Generating the corresponding factors sheet
w = which(factors[,1] %in% x[1])[1]
factors.avg = factors[w,]

for(i in 2:length(x)){
  w = which(factors[,1] %in% x[i])[1]
  factors.avg = rbind(factors.avg, factors[w,])
}

row.names(factors.avg) = factors.avg[,1]
factors.avg = factors.avg[,-1]

#### Removing species with a max abundance of 0 ####
if(length(which(apply(data.avg, 2, max) == 0)) > 0){
  data.avg = data.avg[,-which(apply(data.avg, 2, max) == 0)]
}

#### Removing the Licking-PW data ####
if(PW ==1){
  data.avg = data.avg[-grep("Licking-PW", row.names(data.avg)),]
  factors.avg = factors.avg[-grep("Licking-PW", row.names(factors.avg)),]
}

rm('data')

######################################
#### Performing a SIMPER analysis ####
######################################

if(SIMP){
  sim = simper(data.avg, factors.avg$`Well/Filter`, permutations = 999, trace = T)
  x = summary(sim, ordered = T) # Needed to convert the simper results from a list of lists to a list of data frames
  
  x = mapply(cbind, x, Well = as.list(names(x)), SIMPLIFY = F) # Appends the comparison as a column (might not be needed with the do.call option)
  x = lapply(x, function(y) head(y, n = 20L)) # Finding the top 20 contributors
  names(x) = NULL # This will prevent do call from chaning the names making my for-loop super redudant
  xx = do.call(rbind, x) # Making the final matrix for all comparisons to export
  
  # The for-loop is an alternate method which I like because it doesn't mess with the row names (now redundant)
  # y = NULL
  # 
  # for(i in 1:length(x)){
  #   y = rbind(y, x[[i]])
  # }
  
  sim = xx
}

#########################################################
#### Calculating alpha diversity for the all filters ####
#########################################################

unique.factor = unique(factors.avg[,Var]) # Sets the desired variable as a trackable factor

alpha = as.data.frame(matrix(data = NA, nrow = length(data.avg[,1]), ncol = 6))
row.names(alpha) = row.names(data.avg)

alpha[,1] = as.character(factors.avg$`Well/Filter`)
alpha[,2] = factors.avg$Date
alpha[,3] = diversity(data.avg, index = "shannon")
alpha[,4] = diversity(data.avg, index = "simpson")
alpha[,5] = diversity(data.avg, index = "shannon")/log(specnumber(data.avg))
alpha[,6] = specnumber(data.avg)

colnames(alpha) = c("Well", "Date", "Shannon", "Simpson", "Pielou", "Species_Richness")

alpha$Date = gsub("_","/1/",alpha$Date)
alpha$Date = as.Date(alpha$Date, "%m/%d/%Y")

ggplot(data = alpha, aes(x = Date, y = Shannon))+
  geom_point(aes(colour = Well))+
  geom_line(aes(colour = Well))+
  theme_bw()+
  ggtitle("Shannon's H for 0.2um filters")

ggplot(data = alpha, aes(x = Date, y = Simpson))+
  geom_point(aes(colour = Well))+
  geom_line(aes(colour = Well))+
  theme_bw()+
  ggtitle("Simpson's 1-D for 0.2um filters")

ggplot(data = alpha, aes(x = Date, y = Pielou))+
  geom_point(aes(colour = Well))+
  geom_line(aes(colour = Well))+
  theme_bw()+
  ggtitle("Pielou's J for 0.2um filters")

ggplot(data = alpha, aes(x = Date, y = Species_Richness))+
  geom_point(aes(colour = Well))+
  geom_line(aes(colour = Well))+
  theme_bw()+
  ggtitle("SR for 0.2um filters")

boxplot(Simpson~Well, data = alpha, main = "Simpson Diversity for 0.2um filters")
boxplot(Shannon~Well, data = alpha, main = "Shannon Diversity for 0.2um filters")
boxplot(Pielou~Well, data = alpha, main = "Pielou's J for 0.2um filters")
boxplot(Species_Richness~Well, data = alpha, main = "SR for 0.2um filters")


############################################################
#### Calculating alpha diversity for the 0.2 um filters ####
############################################################

unique.factor = unique(factors.avg[,Var]) # Sets the desired variable as a trackable factor

w = which(factors.avg$Filter %in% 0.2)

alpha = as.data.frame(matrix(data = NA, nrow = length(data.avg[w,1]), ncol = 4))
row.names(alpha) = row.names(data.avg[w,])

alpha[,1] = as.character(factors.avg$`Well/Filter`[w])
alpha[,2] = factors.avg$Date[w]
alpha[,3] = diversity(data.avg[w,], index = "shannon")
alpha[,4] = diversity(data.avg[w,], index = "simpson")

colnames(alpha) = c("Well", "Date", "Shannon", "Simpson")

alpha$Date = gsub("_","/1/",alpha$Date)
alpha$Date = as.Date(alpha$Date, "%m/%d/%Y")

ggplot(data = alpha, aes(x = Date, y = Shannon))+
  geom_point(aes(colour = Well))+
  geom_line(aes(colour = Well))+
  theme_bw()+
  ggtitle("Shannon's H for 0.2um filters")

ggplot(data = alpha, aes(x = Date, y = Simpson))+
  geom_point(aes(colour = Well))+
  geom_line(aes(colour = Well))+
  theme_bw()+
  ggtitle("Simpson's 1-D for 0.2um filters")

boxplot(Simpson~Well, data = alpha, main = "Simpson Diversity for 0.2um filters")
boxplot(Shannon~Well, data = alpha, main = "Shannon Diversity for 0.2um filters")

#######################################################
#### Determing beta diversity for the 0.2um filter ####
#######################################################

w = which(factors.avg$Filter %in% 0.2) # Selecting only the 0.2 filter data

# Generating temporary NMDS data
nms.data = data.avg[w,] 
nms.factors = factors.avg[w,]
nms.factors$`Well/Filter` = factor(nms.factors$`Well/Filter`, levels = unique(nms.factors$`Well/Filter`))
nms.uniq = unique(nms.factors[, Var])

# Setting col and shape vectors for figure generation
col = c(rgb(190,38,37, maxColorValue = 255), rgb(0,97,28, maxColorValue = 255), rgb(13,79,139, maxColorValue = 255)) # Generates vector of consistent colors for the NMDS
shape = c(15,16,17) # Generating consistent shapes for use in the NMDS plot

# Calculating the NMDS space
nmds = metaMDS(nms.data, distance = "bray", k=2, autotransform=F)

# Plotting the NMDS
nmds.plot = ordiplot(nmds,type = "n",xlim = c(min(nmds$points[,1] -.1),max(nmds$points[,1]+.1)), main = "0.2 Filter 16S Data")

for (i in 1:length(nms.uniq)){
  points(nmds.plot$sites[which(nms.factors[, Var] %in% nms.uniq[i]), 1], nmds.plot$sites[which(nms.factors[, Var] %in% nms.uniq[i]), 2],
         col = col[i], pch = shape[i])
} # For loop allows for dynamic plotting of variables; also avoids hard coding

legend("topright",legend = nms.uniq, col = col, pch = shape) # Creates legend from the random points

text(x = nmds.plot$sites[,1]-.3, y = nmds.plot$sites[,2], labels = nms.factors$Date) # Plots dates on top of points

text(x = max(nmds.plot$sites[,1])+.2, y = min(nmds.plot$sites[,2])-.2, labels = sprintf("Stress:%s",nmds$stress)) # Plots the stress value

# Determining the beta dispersion for the 0.2 um filters
bd = betadisper(vegdist(nms.data, method = "bray"), nms.factors$`Well/Filter`)
boxplot(bd, main = "Bray Curtis Beta Dispersion (0.2um)")

x = as.data.frame(cbind(as.character(nms.factors$`Well/Filter`), bd$distances))
x[,2] = as.numeric(as.character(x[,2]))
pirateplot(V2~V1, data = as.data.frame(x), , main = "Bray Curtis Beta Dispersion Distances (0.2um)")

##############################################################
#### Calculating the PERMANOVA groupings for the EPA data ####
##############################################################

x = adonis(formula = data.avg~factors.avg$`Well/Filter`, method = "bray") # Distance needs to be bray because this is ecological/community data

data.perm = cbind("Overall", x$aov.tab$F.Model[1], x$aov.tab$`Pr(>F)`[1]) # Creating seed data

for(i in 1:length(unique.factor)){
  for(n in i:length(unique.factor)){
    if(i==n){
      print(sprintf("Skipped comparison between %s and %s",unique.factor[i],unique.factor[n])) # Skipping comparisons bewtween identical groupings
    } else {
      w = c(which(factors.avg$`Well/Filter` %in% unique.factor[i]), which(factors.avg$`Well/Filter` %in% unique.factor[n])) # Finding the locations for unique factor 'i' and then appending the locations for 'n' to it
      x = adonis(data.avg[w,]~factors.avg$`Well/Filter`[w], method="bray") # Calculating PERMANOVA for those statistics
      w = cbind(sprintf("%s -> %s", unique.factor[i], unique.factor[n]), x$aov.tab$F.Model[1], x$aov.tab$`Pr(>F)`[1]) # Creating a vector of my data
      
      data.perm = rbind(data.perm,w) # Adding my vector of data to the existing seed data
    }
  }
}


#### Finding the maxixums for each sample ####
maxi = matrix(ncol = 5, nrow = length(data.avg[,1])) # Creating the empty matrix to house the top values for each sample
maxab = matrix(ncol = 5, nrow = length(data.avg[,1]))
row.names(maxi) = row.names(data.avg) # Adding an identifier

for (i in 1:length(data.avg[,1])){
  maxi[i,] = names(tail(sort(data.avg[i,]),5))
  maxab[i,] = tail(sort(data.avg[i,]),5)
  
} # The for loop goes through and finds the names corresponding to the top 5 values in each row and writes it to the empty matrix

maxi = cbind(as.character(factors$`Sample Name`),maxi) # Adding further identifiers to the matrix to make it easier to interpret


############################################
#### Organism co-occurence calculations ####
############################################

x = c("Athens 0.2", "Greene 0.2", "Licking 0.2")

corr.data = data.avg # Tried various kinds of different transformations - no change noticable

for(i in 1:length(x)){
  if(i == 1){
    w = rcorr(corr.data[which(factors.avg[,4] %in% x[i]),], type = c("pearson")) # Selected rcorr from package 'Hmisc' because it generates a significance value
    
    m = melt(w$r) # Need to melt the data so that the data is organized in a bit more "user-friendly" manner; doing that with r-values and p-values
    n = melt(w$P)
    
    cor.mat = cbind(x[i], m, n[,3]) # Adding the sample name as a column so I can sort my data
    
    colnames(cor.mat)[1] = c("Sample")
    colnames(cor.mat)[4:5] = c("R-value","p-value") # Naming the columns so it's easier to read
  } else {
    w = rcorr(corr.data[which(factors.avg[,4] %in% x[i]),], type = c("pearson"))
    
    m = melt(w$r)
    n = melt(w$P)
    
    w = cbind(x[i], m, n[,3])
    
    colnames(w)[1] = c("Sample")
    colnames(w)[4:5] = c("R-value","p-value")
    
    cor.mat = rbind(cor.mat, w) # Merging newly generated data to the previously generated data
  }
}

cor.mat = cor.mat[which(cor.mat$`p-value` < 0.05),] # Filtering out those 'insignificant' values, i.e. those with p-values lower than 0.05


########################################
#### Organism-Geochem. Correlations ####
########################################

# Removing data with NAs
geo = geo[,colSums(is.na(geo))==0]

# Data organization for geochem. data
row.names(geo) = geo[,1]

# Removing the Licking-PW data
if(PW == 1){
  geo = geo[-grep("PW", row.names(geo)),]
}

# Setting up the factor list
geo.factors = geo[,2:3]
geo = geo[,-1:-3]

data.geo = data.avg[grep("0.2",factors.avg$`Well/Filter`),] # Picking only the necessary data for these comparisons

# Ensuring that both of the datasets are organized the same way
data.geo = data.geo[order(row.names(data.geo)),] 
geo = geo[order(row.names(geo)),]
geo.factors = geo.factors[order(row.names(geo.factors)),]

# Setting up the for-loop
x = c("AT-6", "GR-13", "LI-4")

for(i in 1:length(x)){
  if(i == 1){
    q = data.geo[which(geo.factors$Well %in% x[i]),]
    r = as.matrix(geo[which(geo.factors$Well %in% x[i]),])
    
    w = rcorr(q, r, type = c("pearson")) # Selected rcorr from package 'Hmisc' because it generates a significance value
    
    m = melt(w$r) # Need to melt the data so that the data is organized in a bit more "user-friendly" manner; doing that with r-values and p-values
    n = melt(w$P)
    
    geo.cor = cbind(x[i], m, n[,3]) # Adding the sample name as a column so I can sort my data
    
    colnames(geo.cor)[1] = c("Sample")
    colnames(geo.cor)[4:5] = c("R-value","p-value") # Naming the columns so it's easier to read
  } else {
    q = data.geo[which(geo.factors$Well %in% x[i]),]
    r = as.matrix(geo[which(geo.factors$Well %in% x[i]),])
    
    w = rcorr(q, r, type = c("pearson"))
    
    m = melt(w$r)
    n = melt(w$P)
    
    w = cbind(x[i], m, n[,3])
    
    colnames(w)[1] = c("Sample")
    colnames(w)[4:5] = c("R-value","p-value")
    
    geo.cor = rbind(geo.cor, w) # Merging newly generated data to the previously generated data
  }
}

geo.cor = geo.cor[which(geo.cor$`p-value` < 0.05),] # Remove non-signficant correlations
geo.cor = geo.cor[which(geo.cor$`R-value` > 0.85 | geo.cor$`R-value` < -0.85),] # Keep only strong relationships


#############################################################
#### Generating a PCA only for the dates we sampled only ####
#############################################################

# Removing the uninteresting/problematic datasets
geo = geo[,-which(names(geo) %in% c("Water_Level","Copper","Zinc","Calcium","Sodium","pH","Temp"))] # I know where the NAs are due to errors in running things or things simply not being measured
geo = scale(geo)

### Generating PCA
pca = prcomp(geo)
pca.plot = ordiplot(pca, type = "n")

for (i in 1:length(unique(geo.factors$Well))){
  points(pca.plot$sites[which(geo.factors$Well %in% unique(geo.factors$Well)[i]), 1],
         pca.plot$sites[which(geo.factors$Well %in% unique(geo.factors$Well)[i]), 2],
         col = col[i], pch = shape[i])
} # For loop will automatically process my data into a much prettier and specific graphs

arrows(x0 = 0, y0 = 0, x1 = pca.plot$species[,1]*3, y1 = pca.plot$species[,2]*3)
text(pca.plot$species[,1]*3, pca.plot$species[,2]*3, row.names(pca.plot$species))
text(pca.plot$sites[,1], pca.plot$sites[,2], row.names(pca.plot$sites))
legend("bottomright", legend = unique(geo.factors$Well), col = col, pch = shape)


########################################
#### Testing out the WGCNA commands ####
########################################
x = unique(factors.avg$`Well/Filter`)

if(WGCNA == 1){ # Switch controlling whether WCGNA will be run or not
  
  for(i in 1:length(x)){
    w = which(factors.avg$`Well/Filter` %in% x[i])
    alt.data = data.avg[w,] # Generating an alternate dataset for only a given location
    
    if(length(which(apply(alt.data, 2, max) == 0)) > 0){
      alt.data = alt.data[,-which(apply(alt.data, 2, max) == 0)]
    } # This is to remove taxa that don't appear in a given sample; given that this is a network, this won't be an issue
    
    powers = c(1:20) # Generates a series of powers to test for future soft thresholding; soft thresholding allows for continuous data to be described
    sft = pickSoftThreshold(alt.data, powerVector = powers, verbose = 5) # Tests the powers against the data to see which power fits the best
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2]) # Plots the test to see if they worked (where it flattens is where it 'worked')
    
    p = -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
    p = sft$fitIndices[,1][which(p %in% max(-sign(sft$fitIndices[,3])*sft$fitIndices[,2]))] # Automatically generating the most 'optimal' power
    
    
    block = blockwiseModules(alt.data, power = p, TOMType = "unsigned", minModuleSize = 10,
                             numericLabels = TRUE) # Assigns the data to modules by transforming the data with the calculated power
    
    moduleLabels = block$colors # Saving labels from the modules
    moduleColors = labels2colors(block$colors) # Saves colors of the modules
    MEs = block$MEs # Saving module eigenvalues which contain the information used to designate modules
    taxTree = block$dendrograms[[1]] # Saving the dendrograms which were generated from TOM data created during module assignment
    
    plotDendroAndColors(taxTree, moduleColors[block$blockGenes[[1]]], "Module colors",
                        dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05) # Plots modules and dendrograms
    
    plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", plotDendrograms = FALSE, xLabelsAngle = 90) # Plots relationships between eigenvectors to see what modules are most closely related
    
    plot(hclust(as.dist(1-cor(MEs)), method = "average")) # Just a plot which can help to see if any module merging can occur
    
    TOM = TOMsimilarity(adjacency(alt.data, power = p)) # Calculates the "Topographcial Overlap Matrix"; a measure of how 'related' two nodes are, what their connections are, how adjacent, etc.
    dissTOM = 1 - TOM # Changes similarity matrix to dissimilarity matrix
    plotTOM = dissTOM^7 # Raises the dis. matrix to a power in order to bring out the effects of 'moderate' effects
    diag(plotTOM) = NA # Sets the diagonal of the data to NA's so that there isn't any weird scaling (they are normally one due to self-comparisons)
    TOMplot(plotTOM, taxTree, moduleColors) # Generates the plot which looks at relationships between each constituent
    
    w = apply(alt.data, 2, mean)
    w = cbind(w, moduleColors)
    
    dimnames(TOM) = list(colnames(alt.data), colnames(alt.data)) # Sets the dimension names to over this matrix to the taxa it was derived from
    cyt = exportNetworkToCytoscape(TOM, edgeFile = sprintf("Class-CytoscapeEdges_%s.txt", x[i]), 
                                   nodeFile = sprintf("Class-CytoscapeNodes_%s.txt", x[i]), threshold = 0.1,
                                   nodeAttr = w) # Threshold is the minimum allowable value for in the TOM matrix
    
    write.csv(cbind(moduleColors, colnames(alt.data)), sprintf("Class-ModuleMembership_%s.csv", x[i]), row.names = F, quote = F)
    
  } # End of the loop running through the samples
  
} # End of the WGCNA switch

###########################################################
#### Processing the raw OTU table and running analyses ####
###########################################################

if(switch.otu == 1){ # Switch to proceed with raw OTU table analyses

  ### Cleaning up the data
  row.names(otu) = otu[,1] # Setting the the first column to the row.names and then removing them
  otu = otu[,-1]
  
  otu = otu[,-which(colnames(otu) %in% c("YB04","BL01","BL02"))] # Removing YB04 becuase it was a failed sequencing run
  
  #otu = apply(otu,2,function(x) {x/sum(x)}) # Turning abundances to relative abundances
  
  otu = t(otu) # Transposing data so it is actually useful
  otu = apply(otu, 2, function(x) as.numeric(as.character(x))) # Reconverts the data to numbers from characters
  row.names(otu) <- row.names(factors)
  
  #### Averaging duplicated OTU sequences ####
  x = unique(factors[,1])
  otu.avg = matrix(nrow = length(x), ncol = length(otu[1,]))
  
  # Generating the averaged raw OTU table
  for(i in 1:length(x)){
    w = which(factors[,1] %in% x[i])
    if(length(w)==1){
      otu.avg[i,] = otu[w,]
    } else {
      otu.avg[i,] = apply(otu[w,], 2, mean)
    }
  }
  
  # Generating the corresponding factors sheet
  w = which(factors[,1] %in% x[1])[1]
  factors.avg = factors[w,]
  
  for(i in 2:length(x)){
    w = which(factors[,1] %in% x[i])[1]
    factors.avg = rbind(factors.avg, factors[w,])
  }
  
  row.names(factors.avg) = factors.avg[,1]
  factors.avg = factors.avg[,-1]
  
  row.names(otu.avg) = row.names(factors.avg)
  colnames(otu.avg) = colnames(otu)
  
  otu = otu.avg # Setting otu equal to otu.avg to make my downstream coding a bit simpler
  
  if(PW == 1){
    otu = otu[-grep("Licking-PW", row.names(otu)),]
    factors.avg = factors.avg[-grep("Licking-PW", row.names(factors.avg)),]
  }
  
  otu.meta = factors.avg
  
  rm('otu.avg')
  
  # Performing some ecological transformations to ensure the tree data is matching the OTU data
  tree = root(tree, 1, r = TRUE) # Rooting the tree because UniFrac requires it
  
  phylo = match.phylo.data(tree, t(otu)) # Need to ensure that the OTUs from the OTU table match those in the tree
  otu = phylo$data # Writing the pruned results for both the OTU table and the tree
  tree = phylo$phy
  
  ##################################
  ### OTU-based WGCNA generation ###
  ##################################
  # Using the point roughly where the graph plateaus in order to generate 
  # the WGCNA plots, rather than the highest point - maybe not
  
  temp.otu = otu # Storing OTU as a temp to restore after the WGCNA because I don't want to recode it
  
  # Converting the alt.data into relative abundance to filter low-abundance OTUs (Cytoscape can't handle all OTUs)
  temp = apply(otu, 2, function(x) x/sum(x))
  p = apply(temp, 1, max)
  otu = otu[-which(p < 0.00005),]
  
  otu = decostand(otu, method = "hellinger")
  otu = log(otu+1)
  
  x = unique(otu.meta$`Well/Filter`)
  y = unique(geo.factors$Well)
  
  TOM.stats = NULL
  
  for(i in 1:length(x)){
    w = which(otu.meta$`Well/Filter` %in% x[i]) # Selecting the well at hand
    
    # Data selection and sorting
    alt.data = t(otu[,w]) # Generating an alternate dataset for only a given location
    alt.meta = otu.meta[w,]
    alt.geo = geo[which(geo.factors$Well %in% y[i]),]
    
    alt.data = alt.data[order(row.names(alt.data)),]
    alt.geo = alt.geo[order(row.names(alt.geo)),]
    
    # Removing absent taxa
    if(length(which(apply(alt.data, 2, max) == 0)) > 0){
      alt.data = alt.data[,-which(apply(alt.data, 2, max) == 0)]
    } # This is to remove taxa that don't appear in a given sample; given that this is a network, this won't be an issue
    
    # Starting WGCNA commands
    powers = c(seq(1,20,by=0.5)) # Generates a series of powers to test for future soft thresholding; soft thresholding allows for continuous data to be described
    sft = pickSoftThreshold(alt.data, powerVector = powers, verbose = 5) # Tests the powers against the data to see which power fits the best
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2]) # Plots the test to see if they worked (where it flattens is where it 'worked')
    
    p = -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
    p = sft$fitIndices[,1][which(p %in% max(-sign(sft$fitIndices[,3])*sft$fitIndices[,2]))] # Automatically generating the most 'optimal' power
    
    
    block = blockwiseModules(alt.data, power = p, TOMType = "signed", minModuleSize = 70,
                             numericLabels = TRUE) # Assigns the data to modules by transforming the data with the calculated power
    
    moduleLabels = block$colors # Saving labels from the modules
    moduleColors = labels2colors(block$colors) # Saves colors of the modules
    MEs = block$MEs # Saving module eigenvalues which contain the information used to designate modules
    taxTree = block$dendrograms[[1]] # Saving the dendrograms which were generated from TOM data created during module assignment
    
    plotDendroAndColors(taxTree, moduleColors[block$blockGenes[[1]]], "Module colors",
                        dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05) # Plots modules and dendrograms
    
    plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", plotDendrograms = FALSE, xLabelsAngle = 90) # Plots relationships between eigenvectors to see what modules are most closely related
    
    plot(hclust(as.dist(1-cor(MEs)), method = "average")) # Just a plot which can help to see if any module merging can occur
    
    TOM = TOMsimilarity(adjacency(alt.data, power = p)) # Calculates the "Topographcial Overlap Matrix"; a measure of how 'related' two nodes are, what their connections are, how adjacent, etc.
    
    dissTOM = 1 - TOM # Changes similarity matrix to dissimilarity matrix
    plotTOM = dissTOM^7 # Raises the dis. matrix to a power in order to bring out the effects of 'moderate' effects
    diag(plotTOM) = NA # Sets the diagonal of the data to NA's so that there isn't any weird scaling (they are normally one due to self-comparisons)
    TOMplot(plotTOM, taxTree, moduleColors) # Generates the plot which looks at relationships between each constituent
    
    w = apply(alt.data, 2, mean)
    w = cbind(w, moduleColors)
    
    dimnames(TOM) = list(colnames(alt.data), colnames(alt.data)) # Sets the dimension names to over this matrix to the taxa it was derived from
    cyt = exportNetworkToCytoscape(TOM, edgeFile = sprintf("OTU(hellinger,03)-CytoscapeEdges_%s.txt", x[i]), 
                                   nodeFile = sprintf("OTU(hellinger,03)-CytoscapeNodes_%s.txt", x[i]), threshold = 0.3,
                                   nodeAttr = w) # Threshold is the minimum allowable value for in the TOM matrix
    
    write.csv(cbind(moduleColors, colnames(alt.data)), sprintf("OTU(hellinger)-ModuleMembership_%s.csv", x[i]), row.names = F, quote = F)
    
  } # End of the loop running through the samples
  

  otu = temp.otu
  
  ##########################################
  ### Comparing OTU data to geochem data ###
  ##########################################
  
  # Looking on at the 0.2 portion of the data
  w = which(otu.meta$Filter %in% 0.2) # Selecting only 0.2
  otu.2 = t(otu[,w]) # Continuning selecting
  
  # Rearranging data to match the organization in the geochem.
  otu.2.geo = otu.2[order(row.names(otu.2)),] # Rearranging to match geo matrix
  w = which(colSums(otu.2.geo) == 0) # Finding OTUs with a max of zero
  otu.2.geo = otu.2.geo[,-w] # Removing said OTUs
  
  # Conducting the procrustes permuted and mantel tests
  prot = protest(otu.2.geo, geo, scale = T, permutations = 999)
  mant = mantel(as.matrix(vegdist(otu.2.geo, method = 'bray')), as.matrix(vegdist(geo, method = 'euclidean')),  permutations = 999)
  
  ##########################################
  ### Analysing average OTU distribution ###
  ##########################################
  for(i in 1:length(x)){
    w = which(otu.meta$`Well/Filter` %in% x[i])
    alt.data = t(otu[,w]) # Generating an alternate dataset for only a given location
    
    # Converting the alt.data into relative abundance and removing OTUs absent from samples
    alt.data = t(apply(alt.data, 1, function(x) x/sum(x)))
    
    if(length(which(apply(alt.data, 2, max) == 0)) > 0){
      alt.data = alt.data[,-which(apply(alt.data, 2, max) == 0)]
    } # This is to remove taxa that don't appear in a given sample; given that this is a network, this won't be an issue
    
    barplot(sort(apply(alt.data,2,mean) ,decreasing = T)[1:250], main = sprintf("%s", x[i]))
    
    p = apply(alt.data, 2, max)
    legend(x = "topright", legend = c(length(which(p < 0.00005))/length(p)))
  }
  
  
  ###########################################
  ### Calculating and plotting Faith's PD ###
  ###########################################
  faith = read.csv("Faith_PD_Total.csv", row.names = 1) # For reanalysis
  # faith = pd(samp = t(otu), tree = tree)
  
  if(rfy ==1){ # Rarefying the OTU table to determine if differences in sequence depth change PD
  
  rare = t(otu)
  rare = rare[which(otu.meta$Filter %in% "0.2"),]
  rare[which(is.wholenumber(rare) == FALSE)] = rare[which(is.wholenumber(rare) == FALSE)]+0.5
  rare = as.data.frame(rrarefy(rare, 8000)) # Rarefy to 8000 which is the value that I used in bNTI

  zero.members = apply(rare, 2, max) # Finding the maximum number of sequence per OTU
  rare = rare[,-which(zero.members==0)]

  phylo.rare = match.phylo.data(tree, t(rare))
  rare = phylo.rare$data
  rare.tree = phylo.rare$phy

  faith = pd(samp = t(rare), tree = rare.tree)
  }
  
  # Data processing
  faith$Well = as.character(factors.avg$`Well/Filter`)
  faith$Date = as.character(factors.avg$Date)
  
  faith$Date = gsub("_","/1/",faith$Date)
  faith$Date = as.Date(faith$Date, "%m/%d/%Y")
  
  # Removing the 0.1 filter data
  w = grep(pattern = "0.1", x = faith$Well)
  faith = faith[-w,]
  
  ggplot(data = faith, aes(x = Date, y = PD))+
    geom_point(aes(colour = Well))+
    geom_line(aes(colour = Well))+
    theme_bw()+
    scale_color_manual(values = col)+
    ggtitle("Faith's PD")
  
  ggplot(data = faith, aes(x = Date, y = SR))+
    geom_point(aes(colour = Well))+
    geom_line(aes(colour = Well))+
    theme_bw()+
    scale_color_manual(values = col)+
    ggtitle("Species Richness")
  
  boxplot(PD~Well, data = faith, main = "Faith's PD")
  pirateplot(PD~Well, data = faith, pal = c(rgb(190,38,37, maxColorValue = 255),
                                                       rgb(0,97,28, maxColorValue = 255),
                                                       rgb(13,79,139, maxColorValue = 255)),
             point.o = 1, inf.f.o = 0.2, bean.b.o = 0, bean.f.o = 0.75, main = "Faith's PD")
  boxplot(SR~Well, data = faith, main = "Species Richness")
  pirateplot(SR~Well, data = faith, pal = c(rgb(190,38,37, maxColorValue = 255),
                                            rgb(0,97,28, maxColorValue = 255),
                                            rgb(13,79,139, maxColorValue = 255)),
             point.o = 1, inf.f.o = 0.2, bean.b.o = 0, bean.f.o = 0.75, main = "Species Richness")
  
  ####################################################################
  ### Other OTU-based alpha diversity metrics to make things equal ###
  ####################################################################
  
  unique.factor = unique(factors.avg[,Var]) # Sets the desired variable as a trackable factor
  
  w = grep("0.2", factors.avg$`Well/Filter`)
  
  alpha = as.data.frame(matrix(data = NA, nrow = length(t(otu)[w,1]), ncol = 6))
  row.names(alpha) = row.names(t(otu)[w,])
  
  alpha[,1] = as.character(factors.avg$`Well/Filter`[w])
  alpha[,2] = factors.avg$Date[w]
  alpha[,3] = diversity(t(otu)[w,], index = "shannon")
  alpha[,4] = diversity(t(otu)[w,], index = "simpson")
  alpha[,5] = diversity(t(otu)[w,], index = "shannon")/log(specnumber(t(otu)[w,])) # This is divided by the log because vegan defined Shannon's H with log, rather than ln
  alpha[,6] = specnumber(t(otu)[w,])
  
  colnames(alpha) = c("Well", "Date", "Shannon", "Simpson", "Pielou", "Species_Richness")
  
  alpha$Date = gsub("_","/1/",alpha$Date)
  alpha$Date = as.Date(alpha$Date, "%m/%d/%Y")
  
  # Time series
  ggplot(data = alpha, aes(x = Date, y = Shannon))+
    geom_point(aes(colour = Well))+
    geom_line(aes(colour = Well))+
    theme_bw()+
    scale_color_manual(values = col)+
    ggtitle("Shannon's H for 0.2um filters")
  
  ggplot(data = alpha, aes(x = Date, y = Simpson))+
    geom_point(aes(colour = Well))+
    geom_line(aes(colour = Well))+
    theme_bw()+
    scale_color_manual(values = col)+
    ggtitle("Simpson's 1-D for 0.2um filters")
  
  ggplot(data = alpha, aes(x = Date, y = Pielou))+
    geom_point(aes(colour = Well))+
    geom_line(aes(colour = Well))+
    theme_bw()+
    scale_color_manual(values = col)+
    ggtitle("Pielou's J for 0.2um filters")
  
  ggplot(data = alpha, aes(x = Date, y = Species_Richness))+
    geom_point(aes(colour = Well))+
    geom_line(aes(colour = Well))+
    theme_bw()+
    scale_color_manual(values = col)+
    ggtitle("SR for 0.2um filters")
  
  # Pirate plots
  pirateplot(Shannon~Well, data = alpha, pal = c(rgb(190,38,37, maxColorValue = 255),
                                            rgb(0,97,28, maxColorValue = 255),
                                            rgb(13,79,139, maxColorValue = 255)),
             point.o = 1, inf.f.o = 0.2, bean.b.o = 0, bean.f.o = 0.75, main = "Shannon's H")
  
  pirateplot(Simpson~Well, data = alpha, pal = c(rgb(190,38,37, maxColorValue = 255),
                                            rgb(0,97,28, maxColorValue = 255),
                                            rgb(13,79,139, maxColorValue = 255)),
             point.o = 1, inf.f.o = 0.2, bean.b.o = 0, bean.f.o = 0.75, main = "Simpsons 1-D")
  
  pirateplot(Pielou~Well, data = alpha, pal = c(rgb(190,38,37, maxColorValue = 255),
                                            rgb(0,97,28, maxColorValue = 255),
                                            rgb(13,79,139, maxColorValue = 255)),
             point.o = 1, inf.f.o = 0.2, bean.b.o = 0, bean.f.o = 0.75, main = "Pielou's J")
  
  pirateplot(Species_Richness~Well, data = alpha, pal = c(rgb(190,38,37, maxColorValue = 255),
                                            rgb(0,97,28, maxColorValue = 255),
                                            rgb(13,79,139, maxColorValue = 255)),
             point.o = 1, inf.f.o = 0.2, bean.b.o = 0, bean.f.o = 0.75, main = "Species Richness")
  
  
  ######################################################################
  ### Generating, processing and plotting Unifrac Data for all data  ###
  ######################################################################
  
  unifracs = GUniFrac(t(otu),tree) # Calculating the UniFrac distances; transposing is important because all of Picante's scripts look for the OTUs differently than I have them (it seems to be phylomatching causes it)
  
  uni.c = c(rgb(190,38,37, maxColorValue = 255), rgb(0,97,28, maxColorValue = 255), rgb(13,79,139, maxColorValue = 255),
            rgb(240,166,166, maxColorValue = 255), rgb(105,144,182, maxColorValue = 255), rgb(109,164,125, maxColorValue = 255)) # Generates vector of consistent colors for the NMDS
  uni.s = c(15, 16, 17, 15, 17, 16) # Generating consistent shapes for use in the NMDS plot
  
  # Plotting the total weighted UniFrac data
  uni.plot = ordiplot(cmdscale(unifracs$unifracs[, , "d_1"],k=2), type="n",main="Weighted UniFrac (All)") # d_1 is the Weighted UniFrac
  for(i in 1:length(unique(otu.meta$`Well/Filter`))){
    w = which(otu.meta$`Well/Filter` %in% unique(otu.meta$`Well/Filter`)[i])
    points(uni.plot$sites[w,1], uni.plot$sites[w,2], col = uni.c[i], pch = uni.s[i])
  }
  
  legend("bottomright", legend = unique(otu.meta$`Well/Filter`), col = uni.c, pch = uni.s) # Creates legend from the random points
  text(x = uni.plot$sites[,1], y = uni.plot$sites[,2], labels = factors.avg$Date)
  
  # Environmental fit for species in the weighted, total UniFrac data
  sp.fit = envfit(uni.plot, data.avg) # Fitting the family data to the graph to see which families are the most correlated with the UniFrac data
  
  # Figuring out the top 20 highest scoring species
  x = scores(sp.fit, display = "vectors")
  y = ((x[,1]-0)^2) + ((x[,2]-0)^2)
  z = names(head(sort(y, decreasing = T), n = 20))
  y = which(names(y) %in% z)
  
  # Attempting to plot the species as size based circles
  points(x = x[y,1]*0.25, y = x[y,2]*0.25, cex = apply(data.avg[,y], 2, sum)*5) # Scaling the sizes of the circles up
  text(x = x[y,1]*0.25, y = x[y,2]*0.25, labels = z)
  
  # Hierarchical clustering
  hc = hclust(as.dist(unifracs$unifracs[,,"d_1"]), method = "average")
  plot(hc)
  
  # Plotting the total unweighted UniFrac data
  uni.plot = ordiplot(cmdscale(unifracs$unifracs[, , "d_UW"],k=2), type="n", main="Unweighted UniFrac (All)") # d_UW is the Unweighted UniFrac
  for(i in 1:length(unique(otu.meta$`Well/Filter`))){
    w = which(otu.meta$`Well/Filter` %in% unique(otu.meta$`Well/Filter`)[i])
    points(uni.plot$sites[w,1], uni.plot$sites[w,2], col = uni.c[i], pch = uni.s[i])
  }
  
  legend("bottomright", legend = unique(otu.meta$`Well/Filter`), col = uni.c, pch = uni.s) # Creates legend from the random points
  text(x = uni.plot$sites[,1], y = uni.plot$sites[,2], labels = factors.avg$Date)
  
  # Environmental fit for the unweighted, total UniFrac data
  sp.fit = envfit(uni.plot, data.avg) # Fitting the family data to the graph to see which families are the most correlated with the UniFrac data
  
  # Figuring out the top 20 highest scoring species
  x = scores(sp.fit, display = "vectors")
  y = ((x[,1]-0)^2) + ((x[,2]-0)^2)
  z = names(head(sort(y, decreasing = T), n = 20))
  y = which(names(y) %in% z)
  
  # Attempting to plot the species as size based circles
  points(x = x[y,1]*0.25, y = x[y,2]*0.25, cex = apply(data.avg[,y], 2, sum)*10) # Scaling the sizes of the circles up
  text(x = x[y,1]*0.25, y = x[y,2]*0.25, labels = z)
  
  # Hierarchical clust
  hc = hclust(as.dist(unifracs$unifracs[,,"d_UW"]), method = "average")
  plot(hc)
  
  #### Determining the permanova relationships between UniFrac clustering
  
  phy.fac = otu.meta$`Well/Filter`
  unique.factor = unique(as.character(otu.meta$`Well/Filter`))
  uni.names = dimnames(unifracs$unifracs)[[3]]
  
  uni.perm = NULL
  
  for(h in 1:length(uni.names)){
    
    uni.state = unifracs$unifracs[,, uni.names[h]]
    
    if(length(which(is.nan(uni.state) %in% TRUE)) < 1){
      x = adonis(formula = uni.state~phy.fac, method = "bray") # Distance needs to be bray because this is ecological/community data
      y = cbind("Overall", x$aov.tab$F.Model[1], x$aov.tab$`Pr(>F)`[1]) # Creating seed data
      uni.perm = rbind(uni.perm, y)
      
      for(i in 1:length(unique.factor)){
        for(n in i:length(unique.factor)){
          if(i==n){
            print(sprintf("Skipped comparison between %s and %s",unique.factor[i],unique.factor[n])) # Skipping comparisons bewtween identical groupings
          } else {
            w = c(which(phy.fac %in% unique.factor[i]), which(phy.fac %in% unique.factor[n])) # Finding the locations for unique factor 'i' and then appending the locations for 'n' to it
            x = adonis(as.dist(uni.state[w,w])~phy.fac[w], method="bray") # Calculating PERMANOVA for those statistics
            w = cbind(sprintf("UniFrac %s - %s -> %s", uni.names[h], unique.factor[i], unique.factor[n]), x$aov.tab$F.Model[1], x$aov.tab$`Pr(>F)`[1]) # Creating a vector of my data
            
            uni.perm = rbind(uni.perm, w) # Adding my vector of data to the existing seed data
          }
        }
      }
    }
  }
  
  uni.stat = NULL
  uni.stat[[1]] = uni.perm
  
  ### Looking at the UniFrac beta dispersion for the complete dataset
  
  # Calculating wUniFrac (Total) BD
  bd = betadisper(as.dist(unifracs$unifracs[,,"d_1"]), otu.meta$`Well/Filter`)
  boxplot(bd, main = "wUniFrac (Total) Beta Dispersion")
  
  x = as.data.frame(cbind(as.character(otu.meta$`Well/Filter`), bd$distances))
  x[,2] = as.numeric(as.character(x[,2]))
  pirateplot(V2~V1, data = as.data.frame(x), pal = c(rgb(240,166,166, maxColorValue = 255), rgb(190,38,37, maxColorValue = 255),
                                                     rgb(109,164,125, maxColorValue = 255), rgb(0,97,28, maxColorValue = 255),
                                                     rgb(105,144,182, maxColorValue = 255), rgb(13,79,139, maxColorValue = 255)),
             main = "wUniFrac (Total) Beta Dispersion Distances")
  
  # Calculating uwUniFrac (Total) BD
  bd = betadisper(as.dist(unifracs$unifracs[,,"d_UW"]), otu.meta$`Well/Filter`)
  boxplot(bd, main = "uwUniFrac (Total) Beta Dispersion")
  
  x = as.data.frame(cbind(as.character(otu.meta$`Well/Filter`), bd$distances))
  x[,2] = as.numeric(as.character(x[,2]))
  pirateplot(V2~V1, data = as.data.frame(x), pal = c(rgb(240,166,166, maxColorValue = 255), rgb(190,38,37, maxColorValue = 255),
                                                     rgb(109,164,125, maxColorValue = 255), rgb(0,97,28, maxColorValue = 255),
                                                     rgb(105,144,182, maxColorValue = 255), rgb(13,79,139, maxColorValue = 255)),
             main = "uwUniFrac (Total) Beta Dispersion Distances")
  
  
  ###############################################################################
  ### Generating, processing and plotting Unifrac Data for 0.2um filters only ###
  ###############################################################################
  # So, the order for my geochemistry doesn't match the order of the UniFrac data - can't guarantee it, 
  # but I expect that will cause issues so I'm fixing it for 0.2 for now.
  env.geo = geo
  row.names(env.geo) = row.names(data.geo)
  env.geo = env.geo[match(row.names(data.avg[which(factors.avg$Filter %in% 0.2),]), row.names(env.geo)),]
  
  
  # Recalculating the UniFrac for only the 0.2um filters
  phy.temp = match.phylo.data(tree,otu[,which(otu.meta$Filter %in% 0.2)]) # Need to ensure that the OTUs from the OTU table match those in the tree
  otu.temp = phy.temp$data # Writing the pruned results for both the OTU table and the tree
  tree.temp = phy.temp$phy
  
  phy.fac = otu.meta$`Well/Filter`[which(otu.meta$Filter %in% 0.2)]
  
  col = c(rgb(190,38,37, maxColorValue = 255), rgb(0,97,28, maxColorValue = 255), rgb(13,79,139, maxColorValue = 255)) # Generates vector of consistent colors for the NMDS
  shape = c(15,16,17) # Generating consistent shapes for use in the NMDS plot
  
  unifracs = GUniFrac(t(otu.temp),tree.temp)
  
  ### Weighted UniFrac for 0.2
  #uni.plot = ordiplot(cmdscale(unifracs$unifracs[, , "d_1"],k=2), type="n", main="Weighted UniFrac (0.2)") # d_1 is the Weighted UniFrac
  uni.pca = princomp(unifracs$unifracs[,,"d_1"])
  uni.plot = ordiplot(uni.pca, type = "n", main = "Weighted UniFrac (0.2um)")
  
  for(i in 1:length(unique(otu.meta$Well))){
    n = which(otu.meta$Filter %in% 0.2)
    w = which(otu.meta$Well[n] %in% unique(otu.meta$Well)[i])
    points(uni.plot$sites[w,1], uni.plot$sites[w,2], col = col[i], pch = shape[i])
  }
  
  legend("bottomleft", legend = unique(otu.meta$Well), col = col, pch = shape) # Creates legend from the random points
  text(x = uni.plot$sites[,1], y = uni.plot$sites[,2], labels = factors.avg$Date[which(otu.meta$Filter %in% 0.2)])
  
  # Environmental fit for the weighted UniFrac data from the 0.2 filters
  sp.fit = envfit(uni.plot, data.avg[which(factors.avg$Filter %in% 0.2),]) # Fitting the family data to the graph to see which families are the most correlated with the UniFrac data
  geo.fit = envfit(uni.plot, env.geo) # Fitting the geochemistry to the UniFrac for the 0.2 data
  
  # Figuring out the top 20 highest scoring species
  x = scores(sp.fit, display = "vectors") # Pulls out the scores enabling me to customize the plot
  y = ((x[,1]-0)^2) + ((x[,2]-0)^2) # Finding the distances for the scores on the ordination
  z = names(head(sort(y, decreasing = T), n = 30)) # Finding the longest 20 arrows - might need to recode this to exclude insignificnat p-values
  y = which(names(y) %in% z) # Finding the location in scores vector to plot only these values
  
  # Attempting to plot the species as size based circles
  points(x = x[y,1]*0.25, y = x[y,2]*0.25, cex = apply(data.avg[,y], 2, sum)*10) # Scaling the sizes of the circles up
  text(x = x[y,1]*0.25, y = x[y,2]*0.25, labels = z)
  
  # Attempting to plot geochemistry as vectors
  plot(geo.fit)
  
  ### Unweighted Unifrac for 0.2
  #uni.plot = ordiplot(cmdscale(unifracs$unifracs[, , "d_UW"],k=2), type="n", main="Unweighted UniFrac (0.2)") # d_UW is the Unweighted UniFrac
  uni.pca = princomp(unifracs$unifracs[,,"d_UW"])
  uni.plot = ordiplot(uni.pca, type = "n", main = "Uweighted UniFrac (0.2um)")
  
  for(i in 1:length(unique(otu.meta$Well))){
    n = which(otu.meta$Filter %in% 0.2)
    w = which(otu.meta$Well[n] %in% unique(otu.meta$Well)[i])
    points(uni.plot$sites[w,1], uni.plot$sites[w,2], col = col[i], pch = shape[i])
  }
  
  legend("bottomright", legend = unique(otu.meta$Well), col = col, pch = shape) # Creates legend from the random points
  text(x = uni.plot$sites[,1], y = uni.plot$sites[,2], labels = factors.avg$Date[which(otu.meta$Filter %in% 0.2)])
  
  # Environmental fit for the unweighted UniFrac data from the 0.2 filters
  sp.fit = envfit(uni.plot, data.avg[which(factors.avg$Filter %in% 0.2),]) # Fitting the family data to the graph to see which families are the most correlated with the UniFrac data
  geo.fit = envfit(uni.plot, env.geo) # Fitting the geochemistry to the UniFrac for the 0.2 data
  
  # Figuring out the top 20 highest scoring species in 0.2 filters
  x = scores(sp.fit, display = "vectors") # Pulls out the scores enabling me to customize the plot
  y = ((x[,1]-0)^2) + ((x[,2]-0)^2) # Finding the distances for the scores on the ordination
  z = names(head(sort(y, decreasing = T), n = 30)) # Finding the longest 20 arrows - might need to recode this to exclude insignificnat p-values
  y = which(names(y) %in% z) # Finding the location in scores vector to plot only these values
  
  # Attempting to plot the species as size based circles
  points(x = x[y,1]*0.25, y = x[y,2]*0.25, cex = apply(data.avg[,y], 2, sum)*10) # Scaling the sizes of the circles up
  text(x = x[y,1]*0.25, y = x[y,2]*0.25, labels = z)
  
  # Attempting to plot geochemistry as vectors
  plot(geo.fit)
  
  
  #### Determining the permanova relationships between UniFrac clustering in 0.2 filters
  unique.factor = unique(as.character(phy.fac))
  uni.names = dimnames(unifracs$unifracs)[[3]]
  
  uni.perm = NULL
  
  for(h in 1:length(uni.names)){
    
    uni.state = unifracs$unifracs[,, uni.names[h]]
    
    if(length(which(is.nan(uni.state) %in% TRUE)) < 1){
      x = adonis(formula = uni.state~phy.fac) # Distance needs to be bray because this is ecological/community data
      y = cbind(sprintf(" %s Overall", uni.names[h]), x$aov.tab$F.Model[1], x$aov.tab$`Pr(>F)`[1]) # Creating seed data
      uni.perm = rbind(uni.perm, y)
      
      for(i in 1:length(unique.factor)){
        for(n in i:length(unique.factor)){
          if(i==n){
            print(sprintf("Skipped comparison between %s and %s",unique.factor[i],unique.factor[n])) # Skipping comparisons bewtween identical groupings
          } else {
            w = c(which(phy.fac %in% unique.factor[i]), which(phy.fac %in% unique.factor[n])) # Finding the locations for unique factor 'i' and then appending the locations for 'n' to it
            x = adonis(as.dist(uni.state[w,w])~phy.fac[w]) # Calculating PERMANOVA for those statistics
            w = cbind(sprintf("UniFrac %s - %s -> %s", uni.names[h], unique.factor[i], unique.factor[n]), x$aov.tab$F.Model[1], x$aov.tab$`Pr(>F)`[1]) # Creating a vector of my data
            
            uni.perm = rbind(uni.perm, w) # Adding my vector of data to the existing seed data
          }
        }
      }
    }
  }
  
  uni.stat[[2]] = uni.perm
  
  ### Looking at the UniFrac beta dispersion for the complete dataset
  
  # Calculating wUniFrac (0.2um) BD
  bd = betadisper(as.dist(unifracs$unifracs[,,"d_1"]), phy.fac)
  boxplot(bd, main = "wUniFrac (0.2um) Beta Dispersion")
  
  x = as.data.frame(cbind(as.character(phy.fac), bd$distances))
  x[,2] = as.numeric(as.character(x[,2]))
  pirateplot(V2~V1, data = as.data.frame(x), pal = c(rgb(190,38,37, maxColorValue = 255),
                                                     rgb(0,97,28, maxColorValue = 255),
                                                     rgb(13,79,139, maxColorValue = 255)),
             main = "wUniFrac (0.2um) Beta Dispersion Distances")
  
  # Calculating uwUniFrac (0.2um) BD
  bd = betadisper(as.dist(unifracs$unifracs[,,"d_UW"]), phy.fac)
  boxplot(bd, main = "uwUniFrac (0.2um) Beta Dispersion")
  
  x = as.data.frame(cbind(as.character(phy.fac), bd$distances))
  x[,2] = as.numeric(as.character(x[,2]))
  pirateplot(V2~V1, data = as.data.frame(x), pal = c(rgb(190,38,37, maxColorValue = 255),
                                                     rgb(0,97,28, maxColorValue = 255),
                                                     rgb(13,79,139, maxColorValue = 255)),
             main = "uwUniFrac (0.2um) Beta Dispersion Distances")
  
  #### Cleaning up unecessary files
  rm("phy.temp", "otu.temp", "tree.temp", "uni.perm")
} # End of the raw OTU table switch