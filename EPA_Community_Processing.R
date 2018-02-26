### This script processing the OTU table to generate some ecological information ###
# RED 2017; danczak.6@osu.edu

#---------------- Switches ----------------#
Var.Num = 5 # Change if the number of variables/factors in the first few columns changes
Var = 4 # Switch controlling which column the variable of interest exists in
PW = 1 # Controls whether I want to remove the Licking-PW from the analysis as it is an outlier (1 = remove, 0 = keep)
WGCNA = 0 # Controls whether I want to generate a WGCNA analysis for network analysis (1 = yes, 0 = no)
rfy = 0 # Controls rarefaction for Faith's PD (1 = rarefy, 0 = don't rarefy)

#-------------- Loading packages --------------#
library(vegan) # For general ecology functions
library(GUniFrac) # For UniFrac
library(picante) # For bNTI and Faith's PD
library(ade4) # For advanced tree building
library(Hmisc) # For rcorr
library(reshape2) # For editing data to better fit the ggplot2 inputs
library(WGCNA) # Looking at networking
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
tree = read.tree("rep_set.tre") # Loading in the tree for UniFrac analyses

### Importing factors information from the Family level OTU table
factors = t(data[1:Var.Num,])
colnames(factors) = factors[1,]
factors = as.data.frame(factors[-1,])

rm("data") # Removing the Family level OTU table as it is no longer needed

### Cleaning up the data
row.names(otu) = otu[,1] # Setting the the first column to the row.names and then removing them
otu = otu[,-1]

otu = otu[,-which(colnames(otu) %in% c("YB04","BL01","BL02"))] # Removing YB04 becuase it was a failed sequencing run

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

# Creating the color and shapes vectors to ensure everything is consistent
col = c(rgb(190,38,37, maxColorValue = 255), rgb(13,79,139, maxColorValue = 255), rgb(0,97,28, maxColorValue = 255)) # Generates vector of consistent colors for the NMDS
shape = c(15,17,16) # Generating consistent shapes for use in the NMDS plot


##################################
### OTU-based WGCNA generation ###
##################################
# Using the point roughly where the graph plateaus in order to generate 
# the WGCNA plots, rather than the highest point - maybe not
if(WGCNA == 1){
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
    
    # q = fundamentalNetworkConcepts(TOM)
    # 
    # TOM.stats[i,] =
    
    # dissTOM = 1 - TOM # Changes similarity matrix to dissimilarity matrix
    # plotTOM = dissTOM^7 # Raises the dis. matrix to a power in order to bring out the effects of 'moderate' effects
    # diag(plotTOM) = NA # Sets the diagonal of the data to NA's so that there isn't any weird scaling (they are normally one due to self-comparisons)
    # TOMplot(plotTOM, taxTree, moduleColors) # Generates the plot which looks at relationships between each constituent
    
    # Determining geochemical and module relationships
    MEs = moduleEigengenes(alt.data, moduleColors)$eigengenes
    MEs = orderMEs(MEs)
    moduleTraitCor = cor(MEs, alt.geo, use = "p")
    moduleTraitPval = corPvalueStudent(moduleTraitCor, nrow(alt.data))
    
    textMatrix = signif(moduleTraitPval, 1)
    dim(textMatrix) = dim(moduleTraitCor)
    
    labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(alt.geo), yLabels = names(MEs), ySymbols = names(MEs), 
                   colorLabels = F, colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = F,
                   cex.text = 0.5, zlim = c(-1,1))
    
    # Sorting the module membership into a list of dataframes
    q = unique(moduleColors)
    
    modules = NULL
    
    for(f in 1:length(q)){
      modules[[f]] = alt.data[,which(moduleColors %in% q[f])]
    } # Loops through each module and adds it to the list; allows maniuplation/examination of each modules
    
    names(modules) = q
    
    w = apply(alt.data, 2, mean)
    w = cbind(w, moduleColors)
    
    dimnames(TOM) = list(colnames(alt.data), colnames(alt.data)) # Sets the dimension names to over this matrix to the taxa it was derived from
    cyt = exportNetworkToCytoscape(TOM, edgeFile = sprintf("OTU(hellinger,03)-CytoscapeEdges_%s.txt", x[i]), 
                                   nodeFile = sprintf("OTU(hellinger,03)-CytoscapeNodes_%s.txt", x[i]), threshold = 0.3,
                                   nodeAttr = w) # Threshold is the minimum allowable value for in the TOM matrix
    
    write.csv(cbind(moduleColors, colnames(alt.data)), sprintf("OTU(hellinger)-ModuleMembership_%s.csv", x[i]), row.names = F, quote = F)
    
  } # End of the loop running through the samples
  
  otu = temp.otu
  
} # End of if-loop for WGCNA
##########################################
### Comparing OTU data to geochem data ###
##########################################

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

rm("otu.2", "otu.2.geo")

##########################################
### Analysing average OTU distribution ###
##########################################
x = unique(factors$`Well/Filter`)

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

##############################################################
#### Calculating the PERMANOVA groupings for the OTU data ####
##############################################################
otu.temp = t(otu[,which(otu.meta$Filter %in% 0.2)]) # Looking only the 0.2 fraction
meta.temp = otu.meta[which(otu.meta$Filter %in% 0.2),] # Matching the metadata
temp.uniq = unique(meta.temp$Well) # Just generating a temporary unique vector

# Removing zeroes from the temp. otu dataset
w = which(colSums(otu.temp) == 0) # Finding OTUs with a max of zero
otu.temp = otu.temp[,-w] # Removing said OTUs

# Moving onto actually running PERMANOVA
x = adonis(formula = otu.temp~meta.temp$Well, method = "bray") # Distance needs to be bray because this is ecological/community data

otu.perm = cbind("Overall", x$aov.tab$F.Model[1], x$aov.tab$`Pr(>F)`[1]) # Creating seed data

for(i in 1:length(temp.uniq)){
  for(n in i:length(temp.uniq)){
    if(i==n){
      print(sprintf("Skipped comparison between %s and %s",temp.uniq[i],temp.uniq[n])) # Skipping comparisons bewtween identical groupings
    } else {
      w = c(which(meta.temp$Well %in% temp.uniq[i]), which(meta.temp$Well %in% temp.uniq[n])) # Finding the locations for unique factor 'i' and then appending the locations for 'n' to it
      x = adonis(otu.temp[w,]~meta.temp$Well[w], method="bray") # Calculating PERMANOVA for those statistics
      w = cbind(sprintf("%s -> %s", temp.uniq[i], temp.uniq[n]), x$aov.tab$F.Model[1], x$aov.tab$`Pr(>F)`[1]) # Creating a vector of my data
      
      otu.perm = rbind(otu.perm,w) # Adding my vector of data to the existing seed data
    }
  }
}

rm("otu.temp", "meta.temp", "temp.uniq")

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
uni.pca = prcomp(unifracs$unifracs[,,"d_1"])
uni.plot = ordiplot(uni.pca, type = "n", main = "Weighted UniFrac (0.2um)")

for(i in 1:length(unique(otu.meta$Well))){
  n = which(otu.meta$Filter %in% 0.2)
  w = which(otu.meta$Well[n] %in% unique(otu.meta$Well)[i])
  points(uni.plot$sites[w,1], uni.plot$sites[w,2], col = col[i], pch = shape[i])
}

legend("bottomleft", legend = unique(otu.meta$Well), col = col, pch = shape) # Creates legend from the random points
text(x = uni.plot$sites[,1], y = uni.plot$sites[,2], labels = factors.avg$Date[which(otu.meta$Filter %in% 0.2)])

# Environmental fit for the weighted UniFrac data from the 0.2 filters
geo.fit = envfit(uni.plot, env.geo) # Fitting the geochemistry to the UniFrac for the 0.2 data

# Attempting to plot geochemistry as vectors
plot(geo.fit)

### Unweighted Unifrac for 0.2
#uni.plot = ordiplot(cmdscale(unifracs$unifracs[, , "d_UW"],k=2), type="n", main="Unweighted UniFrac (0.2)") # d_UW is the Unweighted UniFrac
uni.pca = prcomp(unifracs$unifracs[,,"d_UW"])
uni.plot = ordiplot(uni.pca, type = "n", main = "Uweighted UniFrac (0.2um)")

for(i in 1:length(unique(otu.meta$Well))){
  n = which(otu.meta$Filter %in% 0.2)
  w = which(otu.meta$Well[n] %in% unique(otu.meta$Well)[i])
  points(uni.plot$sites[w,1], uni.plot$sites[w,2], col = col[i], pch = shape[i])
}

legend("bottomright", legend = unique(otu.meta$Well), col = col, pch = shape) # Creates legend from the random points
text(x = uni.plot$sites[,1], y = uni.plot$sites[,2], labels = factors.avg$Date[which(otu.meta$Filter %in% 0.2)])

# Environmental fit for the unweighted UniFrac data from the 0.2 filters
geo.fit = envfit(uni.plot, env.geo) # Fitting the geochemistry to the UniFrac for the 0.2 data

# Attempting to plot geochemistry as vectors
plot(geo.fit)

# Hierarchical clustering of the 0.2 filters
hc = hclust(as.dist(unifracs$unifracs[,,"d_UW"]), method = "average")
plot(hc, main = "Unweighted UniFrac")

hc = hclust(as.dist(unifracs$unifracs[,,"d_1"]), method = "average")
plot(hc, main = "Weighted UniFrac")

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
