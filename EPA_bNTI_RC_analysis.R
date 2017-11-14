#### Code to analyze ecological modeling results ####
# RED 2017; danczak.6@osu.edu

library(reshape2) # Needed to reorganize the data
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

### Setting up bNTI/geochemistry correlations ###
if(corr){
  
  sig = NULL
  sig.names = NULL # Creating null variables which will come to store my signifcance data
  
  for(i in 1:length(geo.alt)){
    
    temp = geo.alt[[i]] # Storing the geochem matrix at hand into a temporary variable to save space
    
    if(length(temp[,1]) < length(bNTI.alt[,1])){
      
      sprintf("%s correlation was skipped due to missing values", names(geo.alt[i])) # Printing the skipped geochem variable
      
    } else {
      
      for(j in 1:length(unique.factors)){
        w = which(bNTI.alt$names %in% unique.factors[j]) # Finding the locations for the unique factors at hand
        
        if(length(w) > 0){
          q = rcorr(temp$value[w], bNTI.alt$value[w]) # Determining the linear r/p-values
          q = c(q$r[1,2], q$P[1,2]) # Storing those values
          sig = rbind(sig, q) # Storing those values, permanently
          
          sig.names = rbind(sig.names,
                            sprintf("%s correlation in %s", names(geo.alt)[i], unique.factors[j])) # Needed to save the names so I could apply them as row names later
          if(is.nan(q[2])){
            
            print("Writing error codes is annoying...")  
            
          } else if(q[2] < 0.05){
            
            pdf(sprintf("W_%s_%s.pdf",names(geo.alt)[i], unique.factors[j]), width = 6.5, height = 3.5)
            
            plot(temp$value[w], bNTI.alt$value[w], ylab = "bNTI", xlab = names(geo.alt)[i], 
                 main = sprintf("Within-well: %s correlation in %s", names(geo.alt)[i], unique.factors[j])) # Generation of the actual plot
            
            if(sum(temp$value[w]) == 0 | sum(bNTI.alt$value[w]) == 0){
              legend("bottomright", 
                     sprintf("No r/p-value for %s v %s", names(geo.alt)[i], unique.factors[j]),
                     bty = "n")
            } else {
              legend("bottomright", sprintf("r: %s  p: %s", q[1], q[2]), bty = "n") # Printing the r/p-values onto the graphs
              abline(lm(bNTI.alt$value[w]~temp$value[w]), col = "blue", lwd = 2) # Plotting the linear model
            } # Needed to code a failsafe in the even that a variable only had 0's
            
            abline(h = c(-2,2), col = "red", lty = 2, lwd = 2) # Plotting significance boundaries
            
            dev.off()
          } # This if-switch will only generate figures that have signficant p-values
          
          rm('w', 'q')
        } # If loop to ensure only the right correlations go through, i.e. LiPW don't have any values as they were point measurments
      } # End of the inner for-loop controlling the plotting and correlations
    } # Ending of the out if-loop which controls which geo chem variables are going to proceed
    rm('temp')
  } # Wrapping up the whole thing; outer for-loop looping through the different geochemical parameters
  
  row.names(sig) = sig.names
  rm('sig.names')
  
  sig[,2] = p.adjust(sig[,2], method = "fdr")
  write.csv(sig, "Significance_Within-Well.csv", quote = F)
  
} # Switch to perform correlations or not - tired of the figures regenerating

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

### Setting up bNTI/geochemistry correlations ###
if(corr){
  
  bw.sig = NULL
  bw.sig.names = NULL # Creating null variables which will come to store my signifcance data
  
  for(i in 1:length(geo.alt)){
    
    temp = geo.alt[[i]] # Storing the geochem matrix at hand into a temporary variable to save space
    
    if(length(temp[,1]) < length(bNTI.alt[,1])){
      
      sprintf("%s correlation was skipped due to missing values", names(geo.alt[i])) # Printing the skipped geochem variable
      
    } else {
      
      for(j in 1:length(unique.factors)){
        w = which(bNTI.alt$names %in% unique.factors[j]) # Finding the locations for the unique factors at hand
        
        if(length(w) > 0){
          q = rcorr(temp$value[w], bNTI.alt$value[w]) # Determining the linear r/p-values
          q = c(q$r[1,2], q$P[1,2]) # Storing those values
          bw.sig = rbind(bw.sig, q) # Storing those values, permanently
          
          bw.sig.names = rbind(bw.sig.names,
                               sprintf("%s correlation in %s", names(geo.alt)[i], unique.factors[j])) # Needed to save the names so I could apply them as row names later
          if(is.nan(q[2])){
            print("I don't want to write anything else right now...")
          } else if (q[2] < 0.05){
            pdf(sprintf("BW_%s_%s.pdf", names(geo.alt)[i], unique.factors[j]), width = 6.5, height = 3.5)
            
            plot(temp$value[w], bNTI.alt$value[w], ylab = "bNTI", xlab = names(geo.alt)[i], 
                 main = sprintf("Between-well: %s correlation in %s", names(geo.alt)[i], unique.factors[j]), type = "n") # Generation of the actual plot
            
            if(length(grep("Athens", bNTI.alt$variable[w])) >= 1){
              r = grep("Athens", bNTI.alt$variable[w])
              points(temp$value[w][r], bNTI.alt$value[w][r], pch = 16, col = rgb(190,38,37, maxColorValue = 255))
            }
            
            if(length(grep("Greene", bNTI.alt$variable[w])) >= 1){
              r = grep("Greene", bNTI.alt$variable[w])
              points(temp$value[w][r], bNTI.alt$value[w][r], pch = 16, col = rgb(0,97,28, maxColorValue = 255))
            }
            
            if(length(grep("Licking", bNTI.alt$variable[w])) >= 1){
              r = grep("Licking", bNTI.alt$variable[w])
              points(temp$value[w][r], bNTI.alt$value[w][r], pch = 16, col = rgb(13,79,139, maxColorValue = 255))
            }
            
            if(length(grep("PW", bNTI.alt$variable[w])) >= 1){
              r = grep("PW", bNTI.alt$variable[w])
              points(temp$value[w][r], bNTI.alt$value[w][r], pch = 16, col = rgb(13,79,139, maxColorValue = 255))
            }
            
            if(sum(temp$value[w]) == 0 | sum(bNTI.alt$value[w]) == 0){
              legend("bottomright", 
                     sprintf("No r/p-value for %s v %s", names(geo.alt)[i], unique.factors[j]),
                     bty = "n")
            } else {
              legend("bottomright", sprintf("r: %s  p: %s", q[1], q[2]), bty = "n") # Printing the r/p-values onto the graphs
              abline(lm(bNTI.alt$value[w]~temp$value[w]), col = "blue", lwd = 2) # Plotting the linear model
            } # Needed to code a failsafe in the even that a variable only had 0's
            
            abline(h = c(-2,2), col = "red", lty = 2, lwd = 2) # Plotting significance boundaries
            
            dev.off()
          } # This if switch ensures only significant figures are plotted
          
          rm('w', 'q', 'r')
        } # If loop to ensure only the right correlations go through, i.e. LiPW don't have any values as they were point measurments
      } # End of the inner for-loop controlling the plotting and correlations
    } # Ending of the out if-loop which controls which geo chem variables are going to proceed
    rm('temp')
  } # Wrapping up the whole thing; outer for-loop looping through the different geochemical parameters
  
  row.names(bw.sig) = bw.sig.names
  rm('bw.sig.names')
  
  bw.sig[,2] = p.adjust(bw.sig[,2], method = "fdr")
  write.csv(bw.sig, "Significance_Between-Well.csv", quote = F)
  
} # Same reason as above - tired of correlations running every time I test the completeness of the script


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

### Setting up bNTI/geochemistry correlations ###
if(corr){
  
  tot.sig = NULL
  tot.sig.names = NULL # Creating null variables which will come to store my signifcance data
  
  for(i in 1:length(geo.alt)){
    
    temp = geo.alt[[i]] # Storing the geochem matrix at hand into a temporary variable to save space
    
    if(length(temp[,1]) < length(bNTI.alt[,1])){
      
      sprintf("%s correlation was skipped due to missing values", names(geo.alt[i])) # Printing the skipped geochem variable
      
    } else {
      
      w = grep("0.2", bNTI.alt$names) # Finding the locations for the unique factors at hand
      
      if(length(w) > 0){
        q = rcorr(temp$value[w], bNTI.alt$value[w]) # Determining the linear r/p-values
        q = c(q$r[1,2], q$P[1,2]) # Storing those values
        tot.sig = rbind(tot.sig, q) # Storing those values, permanently
        
        tot.sig.names = rbind(tot.sig.names,
                              sprintf("%s correlation", names(geo.alt)[i])) # Needed to save the names so I could apply them as row names later
        if(is.nan(q[2])){
          print("I don't want to write anything else right now...")
        } else if (q[2] < 0.05){
          pdf(sprintf("Total_%s.pdf", names(geo.alt)[i]), width = 6.5, height = 3.5)
          
          plot(temp$value[w], bNTI.alt$value[w], ylab = "bNTI", xlab = names(geo.alt)[i], 
               main = sprintf("Total %s correlation", names(geo.alt)[i]), type = "n") # Generation of the actual plot
          
          if(length(grep("Athens", bNTI.alt$names[w])) >= 1){
            r = grep("Athens", bNTI.alt$names[w])
            points(temp$value[w][r], bNTI.alt$value[w][r], pch = 16, col = rgb(190,38,37, maxColorValue = 255))
          }
          
          if(length(grep("Greene", bNTI.alt$names[w])) >= 1){
            r = grep("Greene", bNTI.alt$names[w])
            points(temp$value[w][r], bNTI.alt$value[w][r], pch = 16, col = rgb(0,97,28, maxColorValue = 255))
          }
          
          if(length(grep("Licking", bNTI.alt$names[w])) >= 1){
            r = grep("Licking", bNTI.alt$names[w])
            points(temp$value[w][r], bNTI.alt$value[w][r], pch = 16, col = rgb(13,79,139, maxColorValue = 255))
          }
          
          if(length(grep("PW", bNTI.alt$names[w])) >= 1){
            r = grep("PW", bNTI.alt$names[w])
            points(temp$value[w][r], bNTI.alt$value[w][r], pch = 16, col = rgb(13,79,139, maxColorValue = 255))
          }
          
          if(sum(temp$value[w]) == 0 | sum(bNTI.alt$value[w]) == 0){
            legend("bottomright", 
                   sprintf("No r/p-value for %s v %s", names(geo.alt)[i], unique.factors[j]),
                   bty = "n")
          } else {
            legend("bottomright", sprintf("r: %s  p: %s", q[1], q[2]), bty = "n") # Printing the r/p-values onto the graphs
            abline(lm(bNTI.alt$value[w]~temp$value[w]), col = "blue", lwd = 2) # Plotting the linear model
          } # Needed to code a failsafe in the even that a variable only had 0's
          
          abline(h = c(-2,2), col = "red", lty = 2, lwd = 2) # Plotting significance boundaries
          
          dev.off()
        } # This if switch ensures only significant figures are plotted
        
        rm('w', 'q', 'r')
      } # If loop to ensure only the right correlations go through, i.e. LiPW don't have any values as they were point measurments
    } # Ending of the out if-loop which controls which geo chem variables are going to proceed
    rm('temp')
  } # Wrapping up the whole thing; outer for-loop looping through the different geochemical parameters
  
  row.names(tot.sig) = tot.sig.names
  rm('tot.sig.names')
  
  tot.sig[,2] = p.adjust(tot.sig[,2], method = "fdr")
  write.csv(tot.sig, "Significance_Total.csv", quote = F)
  
} # Same reason as above - tired of correlations running every time I test the completeness of the script


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

# Looking only at the 0.1 filters

w = bNTI.heat[grep("0\\.1", bNTI.heat$names),]
w = w[grep("e.0.1|s.0.1|g.0.1|W.0.1", w$variable),]
w$variable = factor(w$variable, levels = var.ord[grep("e.0.1|s.0.1|g.0.1|W.0.1", var.ord)])

p = ggplot(w, aes_string(x = "names", y = "variable"))+
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

# Looking only at the 0.1 filters

w = raup.heat[grep("0\\.1", raup.heat$names),]
w = w[grep("e.0.1|s.0.1|g.0.1|W.0.1", w$variable),]
w$variable = factor(w$variable, levels = var.ord[grep("e.0.1|s.0.1|g.0.1|W.0.1", var.ord)])

p = ggplot(w, aes_string(x = "names", y = "variable"))+
  geom_tile(aes(fill = value))+
  scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "gray95")+
  scale_y_discrete(limits = rev(levels(w$variable)))+ # Just matching the order of labels on the x-axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p