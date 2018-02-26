#### Code to analyze ecological modeling results ####
# RED 2017; danczak.6@osu.edu

library(reshape2) # Needed to reorganize the data
library(ecodist) # Primarily needed for the Mantel test (so far...)
library(ggplot2) # For plotting
library(Hmisc) # For rcorr
library(yarrr) # For the wonderful pirate plots

fac.name = c("Well.Filter") # This decides what the unique factors will be generated from - will test things out

corr = T # Switch controlling correlations
PW = T # Switch controlling whether I'm including the Licking-PW data or not (usually not, to simplify my life)

#----------------------------------------------------------#
setwd("~/Documents/Wilkins Lab/Ohio EPA Project/Sequencing Data/16S Data/Ecological Modeling Files/")

bNTI = read.csv("MC_EPA_bNTI.csv")
raup.bc = read.csv("RC_BC_MC_Results.csv")
geo = read.csv("Eco_Model_Geochem.csv")

if(all(bNTI[,1:5] == raup.bc[,1:5] & bNTI[,1:5] == geo[,1:5]) != TRUE){
  stop("bNTI and Raup-Crick orders do not match each other.")
} else {
  print("bNTI and Raup-Crick match.")
} # Ensure that my data is arranged the same or else the splitting won't work correctly

row.names(bNTI) = bNTI[,1]
row.names(raup.bc) = raup.bc[,1] # Setting row names
row.names(geo) = geo[,1]

factors = bNTI[,2:5] # Creating factor lists
row.names(factors) = row.names(bNTI) # Giving the factors list row names

# Deleting the Licking-PW data (if the switch is active)
if(PW){
  bNTI = bNTI[-grep("PW",row.names(bNTI)),]
  bNTI = bNTI[,-grep("PW",colnames(bNTI))]
  
  raup.bc = raup.bc[-grep("PW",row.names(raup.bc)),]
  raup.bc = raup.bc[,-grep("PW",colnames(raup.bc))]
  
  geo = geo[-grep("PW",row.names(geo)),]
  factors = factors[-grep("PW",row.names(factors)),]
  
  # Factors get wonky - need to reassign them
  factors = as.data.frame(apply(factors, 2, function(x) factor(x, levels = unique(x))))
}

# Establishing the unique factors vector which enables me to build loops easily
fac.name = which(colnames(factors) %in% fac.name) # Finding column where the desired factors exists
unique.factors = unique(factors[,fac.name]) # Generating unique list

# Finalizing some cleaning steps
bNTI = bNTI[,-1:-5]
raup.bc = raup.bc[,-1:-5] # Removing row names and factors from main data
geo = geo[,-1:-5]

#---------------------------------------------------------#
#--- Generating the necessary geochemistry comparisons ---#
#---------------------------------------------------------#

### Generating pairwise geochemistry comparisons ###
geo.melt = lapply(geo, function(x) abs(outer(x,x,'-'))) # Calculating the pairwise comparison for each variable outputting a list

geo.dist = geo.melt # Storing the "dissimilarity" matrix for the geochemical values

geo.melt = lapply(geo.melt, function(x) {x[upper.tri(x, diag = T)] = NA; x}) # Changing the upper triangle to NA's
geo.melt = lapply(geo.melt, function(x) {row.names(x) <- row.names(geo); x}) # Assigns row names to matrices in list
geo.melt = lapply(geo.melt, function(x) {colnames(x) <- row.names(geo); x}) # Assigns column names to the matrices in the list

geo.melt = lapply(geo.melt, function(x) as.data.frame(x))
geo.melt = lapply(geo.melt, function(x) {x$names <- factors[,fac.name]; x})
geo.melt = lapply(geo.melt, function(x) melt(x, id.vars = "names"))
geo.melt = lapply(geo.melt, function(x) {x <- x[!is.na(x$value),]; x})


#------------------------------------------#
#--- Looking at within well comparisons ---#
#------------------------------------------#

### Looking at within well geochemistry comparisons ###
geo.alt = NULL

for(i in 1:length(geo.melt)){
  
  temp = geo.melt[[i]]
  
  for(j in 1:length(unique.factors)){
    if(j == 1){
      w = temp[which(temp$names %in% unique.factors[j]),]
      temp.alt = w[grep(unique.factors[j], w$variable),] # Only examining the within well comparisons
    } else {
      w = temp[which(temp$names %in% unique.factors[j]),]
      x = w[grep(unique.factors[j], w$variable),]
      temp.alt = rbind(temp.alt, x) # Adding the within well comparisons of other wells to the previous comparisons
      
      rm('w','x')
    }
  }
  geo.alt[[i]] = temp.alt
  rm('temp','temp.alt')
}

names(geo.alt) = names(geo.melt)

### Reorganizing the bNTI data ###
bNTI$names = factors[,fac.name] # Adding a column according to which I can melt the data
bNTI.melt = melt(bNTI, id.vars = "names") # Melting by the newly added column
bNTI.melt = bNTI.melt[!is.na(bNTI.melt$value),] # Removing those NA values derived from the upper triangle

bNTI.melt$variable = gsub(".0.2", " 0.2 ", bNTI.melt$variable)
bNTI.melt$variable = gsub(".0.1", " 0.1 ", bNTI.melt$variable) # These two commands are fixing some random aberrations generated during the melting step

for(i in 1:length(unique.factors)){
  if(i == 1){
    w = bNTI.melt[which(bNTI.melt$names %in% unique.factors[i]),]
    bNTI.alt = w[grep(unique.factors[i], w$variable),] # Only examining the within well comparisons
  } else {
    w = bNTI.melt[which(bNTI.melt$names %in% unique.factors[i]),]
    x = w[grep(unique.factors[i], w$variable),]
    bNTI.alt = rbind(bNTI.alt, x) # Adding the within well comparisons of other wells to the previous comparisons
    
    rm('w','x')
  }
}

# Generating a boxplot for the within well bNTI values
boxplot(value~names, data = bNTI.alt)
abline(h = c(-2,2), col = "red", lty = 2, lwd = 2)

# Generating pirate plot for within well bNTI values
pirateplot(value~names, data = bNTI.alt, pal = c(rgb(240,166,166, maxColorValue = 255), rgb(190,38,37, maxColorValue = 255),
                                                 rgb(109,164,125, maxColorValue = 255), rgb(0,97,28, maxColorValue = 255),
                                                 rgb(105,144,182, maxColorValue = 255), rgb(13,79,139, maxColorValue = 255)),
           point.o = 1, inf.f.o = 0.2, bean.b.o = 0, bean.f.o = 0.75, main = "Within-well bNTI", xlab = NULL, ylab = "bNTI") # Be aware that pirate plot detects when a factor level is missing values, unlike 'boxplot' in R
abline(h = c(-2,2), col = "red", lty = 2, lwd = 2)

### Parsing out the Raup-Crick data ###
raup.bc[upper.tri(raup.bc, diag = T)]=NA
raup.bc$names = factors[,fac.name] # Adding a column according to which I can melt the data
raup.melt = melt(raup.bc, id.vars = "names") # Melting by the newly added column
raup.melt = raup.melt[!is.na(raup.melt$value),] # Removing those NA values derived from the upper triangle

raup.melt$variable = gsub(".0.2", " 0.2 ", raup.melt$variable)
raup.melt$variable = gsub(".0.1", " 0.1 ", raup.melt$variable) # These two commands are fixing some random aberrations generated during the melting step

for(i in 1:length(unique.factors)){
  if(i == 1){
    w = raup.melt[which(raup.melt$names %in% unique.factors[i]),]
    raup.alt = w[grep(unique.factors[i], w$variable),] # Only examining the within well comparisons
  } else {
    w = raup.melt[which(raup.melt$names %in% unique.factors[i]),]
    x = w[grep(unique.factors[i], w$variable),]
    raup.alt = rbind(raup.alt, x) # Adding the within well comparisons of other wells to the previous comparisons
    
    rm('w','x')
  }
}

# Selecting those Raup-Crick values where the bNTI values are insignificant
w = (bNTI.alt$value > -2)&(bNTI.alt$value < 2)
raup.alt = raup.alt[which(w %in% TRUE),]

# Boxplot for Raup-Crick within-well values
boxplot(value~names, data = raup.alt) # Plotting the Raup-Crick data for those non-significant bNTI values
abline(h = c(-0.95,0.95), col = "red", lty = 2, lwd = 2)

# Pirate plot for Raup-Crick within-well values
pirateplot(value~names, data = raup.alt, pal = c(rgb(240,166,166, maxColorValue = 255), rgb(190,38,37, maxColorValue = 255),
                                                 rgb(109,164,125, maxColorValue = 255),
                                                 rgb(105,144,182, maxColorValue = 255), rgb(13,79,139, maxColorValue = 255)),
           point.o = 1, inf.f.o = 0.2, bean.b.o = 0, bean.f.o = 0.75, main = "Within-well Raup-Crick(BC)", xlab = NULL, ylab = "RC-BC")
abline(h = c(-0.95,0.95), col = "red", lty = 2, lwd = 2)

### Generating Mantel correlations ###
mant = NULL
sig.names = NULL

for(i in 1:length(geo.alt)){
  
  temp = geo.dist[[i]] # Storing the geochem matrix at hand into a temporary variable to save space
  
  if(length(temp[,1]) < length(bNTI[,1])){
    
    sprintf("%s correlation was skipped due to missing values", names(geo.alt[i])) # Printing the skipped geochem variable
    
  } else {
    
    for(j in 1:length(unique.factors)){
      w = grep(unique.factors[j], bNTI$names) # Finding the locations for the unique factors at hand
      
      if(length(w) > 0){
        if(sd(temp[w,w], na.rm = T) == 0){ # This if switch checks to make sure that the data is not all 0
          
          mant = rbind(mant, c("0", "0", "0", "0", "0", "0")) # Dummy vector to fill abent data
          
        } else {
          
        m = mantel(as.dist(temp[w,w])~as.dist(bNTI[w,w]), nperm = 9999, mrank = T) # Running spearman correlation, permuted
        mant = rbind(mant, m)
        
        }
      }
      
      sig.names = rbind(sig.names,
                        sprintf("%s correlation in %s", names(geo.alt)[i], unique.factors[j]))
    }
  }
}

mant = apply(mant, 2, as.numeric)

w = which(mant[,4] <= 0.05) # Filtering out insignificant comparisons
mant = mant[w,]
sig.names = sig.names[w,]

row.names(mant) = sig.names # Setting names

#-------------------------------------------#
#--- Looking at between well comparisons ---#
#-------------------------------------------#

# So, I think I'm going to ensure that it is only 0.2 vs 0.2 and 0.1 vs 0.1, at least for now.
# My reasoning is that the 0.1 are strange and distinct enough from even the same well 0.2 filters
# that they will artificially inflate variable selection.

### Looking at between well geochemistry comparisons ###
geo.alt = NULL

for(i in 1:length(geo.melt)){
  
  temp = geo.melt[[i]]
  temp.alt = NULL
  
  for(j in 1:length(unique.factors)){
    if(length(grep("0.2", unique.factors[j])) >= 1){
      # Find the factor in the "being compared" column
      v = temp[which(temp$names %in% unique.factors[j]),]
      
      # Find the factor in the "compared to" column - necessary to ensure we get every comparison
      w = temp[grep(unique.factors[j], temp$variable),]
      w = w[,c("variable","names","value")] # Flipping the columns around to ensure that I know the "being compared" part
      w$variable = unique.factors[j] # Making it consistent with the previous data
      colnames(w) = c("names","variable","value") # Again, just for consistency
      
      # Merge the two into a final variable
      w = rbind(v,w)
      
      # Removing self comparisons and any 0.1 filters
      x = c(grep(unique.factors[j], w$variable), grep("e 0.1|s 0.1|g 0.1|W 0.1", w$variable))
      y = w[-x,]
      
      # Setting final matrix
      temp.alt = rbind(temp.alt, y) # Only examining the between well comparisons (i.e. removing within well comparisons)
      
      rm('v','w','x','y')
    } else {
      
      # Find the factor in the "being compared" column
      v = temp[which(temp$names %in% unique.factors[j]),]
      
      # Find the factor in the "compared to" column - necessary to ensure we get every comparison
      w = temp[grep(unique.factors[j], temp$variable),]
      w = w[,c("variable","names","value")] # Flipping the columns around to ensure that I know the "being compared" part
      w$variable = unique.factors[j] # Making it consistent with the previous data
      colnames(w) = c("names","variable","value") # Again, just for consistency
      
      # Merge the two into a final variable
      w = rbind(v,w)
      
      # Removing self comparisons and any 0.1 filters
      x = c(grep(unique.factors[j], w$variable), grep("e 0.2|s 0.2|g 0.2|W 0.2", w$variable))
      y = w[-x,]
      
      # Setting final matrix
      temp.alt = rbind(temp.alt, y) # Only examining the between well comparisons (i.e. removing within well comparisons)
      
      rm('v','w','x','y')
    }
  }
  geo.alt[[i]] = temp.alt
  rm('temp','temp.alt')
}

names(geo.alt) = names(geo.melt)

### Reorganizing the bNTI data ###

bNTI.alt = NULL

for(i in 1:length(unique.factors)){
  if(length(grep("0.2", unique.factors[i])) >= 1){
    
    # Find the factor in the "being compared" column
    v = bNTI.melt[which(bNTI.melt$names %in% unique.factors[i]),]
    
    # Find the factor in the "compared to" column - necessary to ensure we get every comparison
    w = bNTI.melt[grep(unique.factors[i], bNTI.melt$variable),]
    w = w[,c("variable","names","value")] # Flipping the columns around to ensure that I know the "being compared" part
    w$variable = unique.factors[i] # Making it consistent with the previous data
    colnames(w) = c("names","variable","value") # Again, just for consistency
    
    # Merge the two into a final variable
    w = rbind(v,w)
    
    # Removing self comparisons and any 0.1 filters
    x = c(grep(unique.factors[i], w$variable), grep("e 0.1|s 0.1|g 0.1|W 0.1", w$variable))
    y = w[-x,]
    
    # Setting final matrix
    bNTI.alt = rbind(bNTI.alt, y) # Only examining the between well comparisons (i.e. removing within well comparisons)
    
    rm('v','w','x','y')
  } else {
    
    # Find the factor in the "being compared" column
    v = bNTI.melt[which(bNTI.melt$names %in% unique.factors[i]),]
    
    # Find the factor in the "compared to" column - necessary to ensure we get every comparison
    w = bNTI.melt[grep(unique.factors[i], bNTI.melt$variable),]
    w = w[,c("variable","names","value")] # Flipping the columns around to ensure that I know the "being compared" part
    w$variable = unique.factors[i] # Making it consistent with the previous data
    colnames(w) = c("names","variable","value") # Again, just for consistency
    
    # Merge the two into a final variable
    w = rbind(v,w)
    
    # Removing self comparisons and any 0.2 filter samples
    x = c(grep(unique.factors[i], w$variable), grep("e 0.2|s 0.2|g 0.2|W 0.2", w$variable))
    y = w[-x,]
    
    # Setting final matrix
    bNTI.alt = rbind(bNTI.alt, y) # Only examining the between well comparisons (i.e. removing within well comparisons)
    
    rm('v','w','x','y')
  }
}

boxplot(value~names, data = bNTI.alt)
abline(h = c(-2,2), col = "red", lty = 2, lwd = 2)

# Generating pirate plot for between well bNTI values
pirateplot(value~names, data = bNTI.alt, pal = c(rgb(240,166,166, maxColorValue = 255), rgb(190,38,37, maxColorValue = 255),
                                                 rgb(109,164,125, maxColorValue = 255), rgb(0,97,28, maxColorValue = 255),
                                                 rgb(105,144,182, maxColorValue = 255), rgb(13,79,139, maxColorValue = 255)),
           point.o = 1, inf.f.o = 0.2, bean.b.o = 0, bean.f.o = 0.75, main = "Between-well bNTI", xlab = NULL, ylab = "bNTI") # Be aware that pirate plot detects when a factor level is missing values, unlike 'boxplot' in R
abline(h = c(-2,2), col = "red", lty = 2, lwd = 2)

### Parsing out the Raup-Crick data ###

raup.alt = NULL

for(i in 1:length(unique.factors)){
  if(length(grep("0.2", unique.factors[i])) >= 1){
    
    # Find the factor in the "being compared" column
    v = raup.melt[which(raup.melt$names %in% unique.factors[i]),]
    
    # Find the factor in the "compared to" column - necessary to ensure we get every comparison
    w = raup.melt[grep(unique.factors[i], bNTI.melt$variable),]
    w = w[,c("variable","names","value")] # Flipping the columns around to ensure that I know the "being compared" part
    w$variable = unique.factors[i] # Making it consistent with the previous data
    colnames(w) = c("names","variable","value") # Again, just for consistency
    
    # Merge the two into a final variable
    w = rbind(v,w)
    
    # Removing self comparisons and any 0.1 filters
    x = c(grep(unique.factors[i], w$variable), grep("e 0.1|s 0.1|g 0.1|W 0.1", w$variable))
    y = w[-x,]
    
    # Setting final matrix
    raup.alt = rbind(raup.alt, y) # Only examining the between well comparisons (i.e. removing within well comparisons)
    
    rm('v','w','x','y')
  } else {
    
    # Find the factor in the "being compared" column
    v = raup.melt[which(raup.melt$names %in% unique.factors[i]),]
    
    # Find the factor in the "compared to" column - necessary to ensure we get every comparison
    w = raup.melt[grep(unique.factors[i], bNTI.melt$variable),]
    w = w[,c("variable","names","value")] # Flipping the columns around to ensure that I know the "being compared" part
    w$variable = unique.factors[i] # Making it consistent with the previous data
    colnames(w) = c("names","variable","value") # Again, just for consistency
    
    # Merge the two into a final variable
    w = rbind(v,w)
    
    # Removing self comparisons and any 0.1 filters
    x = c(grep(unique.factors[i], w$variable), grep("e 0.2|s 0.2|g 0.2|W 0.2|W 0.1", w$variable))
    y = w[-x,]
    
    # Setting final matrix
    raup.alt = rbind(raup.alt, y) # Only examining the between well comparisons (i.e. removing within well comparisons)
    
    rm('v','w','x','y')
  }
}

w = (bNTI.alt$value > -2)&(bNTI.alt$value < 2)
raup.alt = raup.alt[which(w %in% TRUE),]

# Boxplot for Raup-Crick, between-well values
boxplot(value~names, data = raup.alt) # Plotting the Raup-Crick data for those non-significant bNTI values
abline(h = c(-0.95,0.95), col = "red", lty = 2, lwd = 2)

# Pirate plot for Raup-Crick between-well values
pirateplot(value~names, data = raup.alt, pal = c(rgb(240,166,166, maxColorValue = 255), rgb(190,38,37, maxColorValue = 255),
                                                 rgb(109,164,125, maxColorValue = 255),
                                                 rgb(105,144,182, maxColorValue = 255), rgb(13,79,139, maxColorValue = 255)),
           point.o = 1, inf.f.o = 0.2, bean.b.o = 0, bean.f.o = 0.75, main = "Between-well Raup-Crick(BC)", xlab = NULL, ylab = "RC-BC")
abline(h = c(-0.95,0.95), col = "red", lty = 2, lwd = 2)


#----------------------------------#
#--- Looking at total fractions ---#
#----------------------------------#
# So, this time, instead of looking at the between well comparisons exclusively, I'm looking at them with
# the self comparisons intact. This is mostly just concerning the correlations.

### Looking at between well geochemistry comparisons ###
geo.alt = NULL

for(i in 1:length(geo.melt)){
  
  temp = geo.melt[[i]]
  temp.alt = NULL
  
  for(j in 1:length(unique.factors)){
    if(length(grep("0.2", unique.factors[j])) >= 1){
      # Find the factor in the "being compared" column
      v = temp[which(temp$names %in% unique.factors[j]),]
      
      # Find the factor in the "compared to" column - necessary to ensure we get every comparison
      w = temp[grep(unique.factors[j], temp$variable),]

      # Merge the two into a final variable
      w = rbind(v,w)
      
      # Removing any 0.1 filter samples
      x = c(grep("e 0.1|s 0.1|g 0.1|W 0.1", w$names), grep("e 0.1|s 0.1|g 0.1|W 0.1", w$variable))
      y = w[-x,]
      
      # Setting final matrix
      temp.alt = rbind(temp.alt, y) # Only examining the between well comparisons (i.e. removing within well comparisons)
      
      rm('v','w','x','y')
    } else {
      
      # Find the factor in the "being compared" column
      v = temp[which(temp$names %in% unique.factors[j]),]
      
      # Find the factor in the "compared to" column - necessary to ensure we get every comparison
      w = temp[grep(unique.factors[j], temp$variable),]
      
      # Merge the two into a final variable
      w = rbind(v,w)
      
      # Removing self comparisons and any 0.1 filters
      x = c(grep("e 0.2|s 0.2|g 0.2|W 0.2", w$names), grep("e 0.2|s 0.2|g 0.2|W 0.2", w$variable))
      y = w[-x,]
      
      # Setting final matrix
      temp.alt = rbind(temp.alt, y) # Only examining the between well comparisons (i.e. removing within well comparisons)
      
      rm('v','w','x','y')
    }
  }
  geo.alt[[i]] = temp.alt
  rm('temp','temp.alt')
}

names(geo.alt) = names(geo.melt)

### Reorganizing the bNTI data ###

bNTI.alt = NULL

for(i in 1:length(unique.factors)){
  if(length(grep("0.2", unique.factors[i])) >= 1){
    
    # Find the factor in the "being compared" column
    v = bNTI.melt[which(bNTI.melt$names %in% unique.factors[i]),]
    
    # Find the factor in the "compared to" column - necessary to ensure we get every comparison
    w = bNTI.melt[grep(unique.factors[i], bNTI.melt$variable),]
    
    # Merge the two into a final variable
    w = rbind(v,w)
    
    # Removing any 0.1 filters but keeping self comparisons
    x = c(grep("e 0.1|s 0.1|g 0.1|W 0.1", w$names), grep("e 0.1|s 0.1|g 0.1|W 0.1", w$variable))
    y = w[-x,]
    
    # Setting final matrix
    bNTI.alt = rbind(bNTI.alt, y) # Only examining the between well comparisons (i.e. removing within well comparisons)
    
    rm('v','w','x','y')
  } else {
    
    # Find the factor in the "being compared" column
    v = bNTI.melt[which(bNTI.melt$names %in% unique.factors[i]),]
    
    # Find the factor in the "compared to" column - necessary to ensure we get every comparison
    w = bNTI.melt[grep(unique.factors[i], bNTI.melt$variable),]
    
    # Merge the two into a final variable
    w = rbind(v,w)
    
    # Removing any 0.2 filter samples but keeping self comparisons
    x = c(grep("e 0.2|s 0.2|g 0.2|W 0.2", w$names), grep("e 0.2|s 0.2|g 0.2|W 0.2", w$variable))
    y = w[-x,]
    
    # Setting final matrix
    bNTI.alt = rbind(bNTI.alt, y) # Only examining the between well comparisons (i.e. removing within well comparisons)
    
    rm('v','w','x','y')
  }
}

boxplot(value~names, data = bNTI.alt)
abline(h = c(-2,2), col = "red", lty = 2, lwd = 2)

# Generating pirate plot for between well bNTI values
pirateplot(value~names, data = bNTI.alt, pal = c(rgb(240,166,166, maxColorValue = 255), rgb(190,38,37, maxColorValue = 255),
                                                 rgb(109,164,125, maxColorValue = 255), rgb(0,97,28, maxColorValue = 255),
                                                 rgb(105,144,182, maxColorValue = 255), rgb(13,79,139, maxColorValue = 255)),
           point.o = 1, inf.f.o = 0.2, bean.b.o = 0, bean.f.o = 0.75, main = "Between-well w/ self bNTI", xlab = NULL, ylab = "bNTI") # Be aware that pirate plot detects when a factor level is missing values, unlike 'boxplot' in R
abline(h = c(-2,2), col = "red", lty = 2, lwd = 2)

### Parsing out the Raup-Crick data ###

raup.alt = NULL

for(i in 1:length(unique.factors)){
  if(length(grep("0.2", unique.factors[i])) >= 1){
    
    # Find the factor in the "being compared" column
    v = raup.melt[which(raup.melt$names %in% unique.factors[i]),]
    
    # Find the factor in the "compared to" column - necessary to ensure we get every comparison
    w = raup.melt[grep(unique.factors[i], bNTI.melt$variable),]
    
    # Merge the two into a final variable
    w = rbind(v,w)
    
    # Removing any 0.1 filters but retaining self comparisons
    x = c(grep("e 0.1|s 0.1|g 0.1|W 0.1", w$names), grep("e 0.1|s 0.1|g 0.1|W 0.1", w$variable))
    y = w[-x,]
    
    # Setting final matrix
    raup.alt = rbind(raup.alt, y) # Only examining the between well comparisons (i.e. removing within well comparisons)
    
    rm('v','w','x','y')
  } else {
    
    # Find the factor in the "being compared" column
    v = raup.melt[which(raup.melt$names %in% unique.factors[i]),]
    
    # Find the factor in the "compared to" column - necessary to ensure we get every comparison
    w = raup.melt[grep(unique.factors[i], bNTI.melt$variable),]
    
    # Merge the two into a final variable
    w = rbind(v,w)
    
    # Removing any 0.1 filters but retaining self comparisons
    x = c(grep("e 0.2|s 0.2|g 0.2|W 0.2", w$names), grep("e 0.2|s 0.2|g 0.2|W 0.2", w$variable))
    y = w[-x,]
    
    # Setting final matrix
    raup.alt = rbind(raup.alt, y) # Only examining the between well comparisons (i.e. removing within well comparisons)
    
    rm('v','w','x','y')
  }
}

w = (bNTI.alt$value > -2)&(bNTI.alt$value < 2)
raup.alt = raup.alt[which(w %in% TRUE),]

# Boxplot for Raup-Crick, between-well values
boxplot(value~names, data = raup.alt) # Plotting the Raup-Crick data for those non-significant bNTI values
abline(h = c(-0.95,0.95), col = "red", lty = 2, lwd = 2)

# Pirate plot for Raup-Crick between-well values
pirateplot(value~names, data = raup.alt, pal = c(rgb(240,166,166, maxColorValue = 255), rgb(190,38,37, maxColorValue = 255),
                                                 rgb(109,164,125, maxColorValue = 255),
                                                 rgb(105,144,182, maxColorValue = 255), rgb(13,79,139, maxColorValue = 255)),
           point.o = 1, inf.f.o = 0.2, bean.b.o = 0, bean.f.o = 0.75, main = "Between-well w/ self Raup-Crick(BC)", xlab = NULL, ylab = "RC-BC")
abline(h = c(-0.95,0.95), col = "red", lty = 2, lwd = 2)

### Locking at between well Mantel comparisons ###
w = grep(0.2, bNTI$names) # Ensuring that only the 0.2 fraction is being examined
bNTI.bet = bNTI[w,w] # Selecting those subsets of data
geo.bet = lapply(geo.dist, function(x) x[w,w])

mant.bet = NULL

for(i in 1:length(geo.alt)){
  mant.bet = rbind(mant.bet, mantel(as.dist(geo.bet[[i]])~as.dist(bNTI.bet), nperm = 9999, mrank = T))
}

row.names(mant.bet) = names(geo.alt)

rm("bNTI.bet", "raup.alt", "bNTI.alt", "w")

#------------------------------#
#--- Looking at the heatmap ---#
#------------------------------#

### Coding a heatmap to look at those comparisons which were the most signficant ###
bNTI = bNTI[,-length(bNTI[1,])] # Given the current code, I need to remove the names first and then reimplement them - might debug that later
raup.bc = raup.bc[,-length(raup.bc[1,])]

bNTI[upper.tri(bNTI)] = t(bNTI)[upper.tri(bNTI)] # Mirroring the comparisons so that the data looks prettier in a heat map
raup.bc[upper.tri(raup.bc)] = t(raup.bc)[upper.tri(raup.bc)]

bNTI = as.matrix(bNTI) # Need to do some silly matrix work to make the changes I want...
raup.bc = as.matrix(raup.bc)

w = (bNTI>-2)&(bNTI<2)
raup.bc[which(w %in% FALSE)] = NA # Replacing the values that have significant bNTI values with NA to only look at the insignificant portion

bNTI = as.data.frame(bNTI)
raup.bc = as.data.frame(raup.bc) # Need to reconvert back to a data frame or else the next command causes problems

bNTI$names = row.names(bNTI) # Adding the data to melt by
raup.bc$names = row.names(raup.bc)


### bNTI values first ###

bNTI.heat = melt(bNTI, id.vars = "names") # Melting the data by the newly appended row names

names.ord = c("Athens 0.1 7/8/14","Athens 0.1 1/19/15","Athens 0.1 10/27/15","Athens 0.1 4/15/15","Athens 0.1 7/22/15",
              "Athens 0.1 2/2016","Athens 0.1 4/2016","Athens 0.1 7/2016","Athens 0.2 7/8/14","Athens 0.2 10/7/14",
              "Athens 0.2 1/19/15","Athens 0.2 4/15/15","Athens 0.2 7/22/15","Athens 0.2 10/27/15","Athens 0.2 2/2016",
              "Athens 0.2 4/2016","Athens 0.2 7/2016","Greene 0.1 7/8/14","Greene 0.1 10/7/14","Greene 0.1 1/19/15",
              "Greene 0.1 4/15/15","Greene 0.1 7/22/15","Greene 0.1 10/27/15","Greene 0.1 2/2016","Greene 0.1 4/2016",
              "Greene 0.1 7/2016","Greene 0.2 7/8/14","Greene 0.2 10/7/14","Greene 0.2 1/19/15","Greene 0.2 4/15/15",
              "Greene 0.2 7/22/15","Greene 0.2 10/27/15","Greene 0.2 2/2016","Greene 0.2 4/2016","Greene 0.2 7/2016",
              "Licking 0.1 7/8/14","Licking 0.1 1/19/15","Licking 0.1 10/27/15","Licking 0.1 4/15/15","Licking 0.1 7/22/15",
              "Licking 0.1 2/2016","Licking 0.1 4/2016","Licking 0.1 7/2016","Licking 0.2 7/8/14","Licking 0.2 10/7/14",
              "Licking 0.2 1/19/15","Licking 0.2 4/15/15","Licking 0.2 7/22/15","Licking 0.2 10/27/15","Licking 0.2 2/2016",
              "Licking 0.2 4/2016","Licking 0.2 7/2016","Licking-PW 0.1","Licking-PW 0.2")

var.ord = c("Athens.0.1.7.8.14","Athens.0.1.1.19.15","Athens.0.1.10.27.15","Athens.0.1.4.15.15","Athens.0.1.7.22.15",
            "Athens.0.1.2.2016","Athens.0.1.4.2016","Athens.0.1.7.2016","Athens.0.2.7.8.14","Athens.0.2.10.7.14",
            "Athens.0.2.1.19.15","Athens.0.2.4.15.15","Athens.0.2.7.22.15","Athens.0.2.10.27.15","Athens.0.2.2.2016",
            "Athens.0.2.4.2016","Athens.0.2.7.2016","Greene.0.1.7.8.14","Greene.0.1.10.7.14","Greene.0.1.1.19.15",
            "Greene.0.1.4.15.15","Greene.0.1.7.22.15","Greene.0.1.10.27.15","Greene.0.1.2.2016","Greene.0.1.4.2016",
            "Greene.0.1.7.2016","Greene.0.2.7.8.14","Greene.0.2.10.7.14","Greene.0.2.1.19.15","Greene.0.2.4.15.15",
            "Greene.0.2.7.22.15","Greene.0.2.10.27.15","Greene.0.2.2.2016","Greene.0.2.4.2016","Greene.0.2.7.2016",
            "Licking.0.1.7.8.14","Licking.0.1.1.19.15","Licking.0.1.10.27.15","Licking.0.1.4.15.15","Licking.0.1.7.22.15",
            "Licking.0.1.2.2016","Licking.0.1.4.2016","Licking.0.1.7.2016","Licking.0.2.7.8.14","Licking.0.2.10.7.14",
            "Licking.0.2.1.19.15","Licking.0.2.4.15.15","Licking.0.2.7.22.15","Licking.0.2.10.27.15","Licking.0.2.2.2016",
            "Licking.0.2.4.2016","Licking.0.2.7.2016","Licking.PW.0.1","Licking.PW.0.2")

bNTI.heat$names = factor(bNTI.heat$names, levels = names.ord)
bNTI.heat$variable = factor(bNTI.heat$variable, levels = var.ord)

# Looking at the complete picture

p = ggplot(bNTI.heat, aes_string(x = "names", y = "variable"))+ # Heat map magic
  geom_tile(aes(fill = value))+
  scale_fill_gradient2(low = "steelblue", high = "red", mid = "gray95")+
  scale_y_discrete(limits = rev(levels(bNTI.heat$variable)))+ # Reversing the y-axis to make the order make sense
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p

# Looking at only the 0.2 filters

w = bNTI.heat[grep("0\\.2", bNTI.heat$names),] # Selecting the 0.2 filters in the first column
w = w[grep("e.0.2|s.0.2|g.0.2|W.0.2", w$variable),] # Selecting the 0.2 filter comparisons
w$variable = factor(w$variable, levels = var.ord[grep("e.0.2|s.0.2|g.0.2|W.0.2", var.ord)]) # Adjusting the factors to the fact that a lot of data is now missing

p = ggplot(w, aes_string(x = "names", y = "variable"))+ # Heat map magic
  geom_tile(aes(fill = value))+
  scale_fill_gradient2(low = "steelblue", high = "red", mid = "gray95")+
  scale_y_discrete(limits = rev(levels(w$variable)))+ # Just matching the order of labels on the x-axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p

### Raup-Crick values second ###

raup.heat = melt(raup.bc, id.vars = "names") # Melting the data by the newly appended row names

raup.heat$names = factor(raup.heat$names, levels = names.ord)
raup.heat$variable = factor(raup.heat$variable, levels = var.ord)

# Looking at the complete picture

p = ggplot(raup.heat, aes_string(x = "names", y = "variable"))+ # Heat map magic
  geom_tile(aes(fill = value))+
  scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "gray95")+
  scale_y_discrete(limits = rev(levels(raup.heat$variable)))+ # Just matching the order of labels on the x-axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p

# Looking at only the 0.2 filters

w = raup.heat[grep("0\\.2", raup.heat$names),] # Selecting the 0.2 filters in the first column
w = w[grep("e.0.2|s.0.2|g.0.2|W.0.2", w$variable),] # Selecting the 0.2 filter comparisons
w$variable = factor(w$variable, levels = var.ord[grep("e.0.2|s.0.2|g.0.2|W.0.2", var.ord)]) # Adjusting the factors to the fact that a lot of data is now missing

p = ggplot(w, aes_string(x = "names", y = "variable"))+ # Heat map magic
  geom_tile(aes(fill = value))+
  scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "gray95")+
  scale_y_discrete(limits = rev(levels(w$variable)))+ # Just matching the order of labels on the x-axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p

############################################################
### Generating percentages from bNTI and Raup-Crick data ###
############################################################
## Assigning a percent explained by selection/neutral processes to each sample

# Data sorting and whatnot
x = bNTI[,-length(bNTI[1,])] # bNTI without the names vector tagged on
y = raup.bc[,-length(raup.bc[1,])] # Raup-Crick without the names vector tagged on

w = grep("e.0.2|s.0.2|g.0.2", row.names(x))
x = x[w,w]
y = y[w,w]

# Assigning names to begin calculating percentages
model = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model[x > 2] = "Variable Selection"
model[x < -2] = "Homogenizing Selection"
model[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model) = dimnames(x)

rm("x", "y", "w")

# Compiling model results for all comparisons
model.results = NULL
model.results = cbind(model.results, apply(model, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.results = cbind(model.results, apply(model, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.results = cbind(model.results, apply(model, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.results = cbind(model.results, apply(model, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.results = cbind(model.results, apply(model, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.results) = row.names(model)
colnames(model.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

# Compiling model results for self comparisons only
x = c("Athens", "Greene", "Licking")
model.self = NULL

for(i in 1:length(x)){
  w = grep(x[i], row.names(model))
  temp = model[w,w]
  
  model.temp = NULL
  model.temp = cbind(model.temp, apply(temp, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
  model.temp = cbind(model.temp, apply(temp, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
  model.temp = cbind(model.temp, apply(temp, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
  model.temp = cbind(model.temp, apply(temp, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
  model.temp = cbind(model.temp, apply(temp, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))
  
  model.self = rbind(model.self, model.temp)
  rm("temp", "model.temp")
}

colnames(model.self) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")
