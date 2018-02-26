### Processing the cohesion data

# Loading in packages
library(reshape2)
library(ggplot2)
library(vegan)
library(ecodist)
library(Hmisc)

####################################
### Loading and sorting the data ###
####################################

# Setting work directory and loading data
setwd("~/Documents/Wilkins Lab/Ohio EPA Project/Sequencing Data/16S Data/Ecological Modeling Files/")
data = read.csv("Total_0.2_Cohesion.csv", row.names = 1)
comm = read.table("Raw_EPA_no_pynast.txt", row.names = 1, header = T)
fac.file = read.csv("Family_Level_OTU_edit.csv")
w.uni = read.csv("../Weighted_Uni_output.csv", row.names = 1)
uw.uni = read.csv("../Unweighted_Uni_output.csv", row.names = 1)
bnti = read.csv("MC_EPA_bNTI.csv", row.names = 1)

# Storing the factors
factors = t(fac.file[1:5,])
colnames(factors) = factors[1,]
factors = factors[-1,]

# General data scrubbing
comm = comm[,-which(colnames(comm) %in% c("YB04","BL01","BL02"))] # Removing YB04 becuase it was a failed sequencing run
colnames(comm) = colnames(fac.file[,-1])

# # Removing the factors from the community matrix and bnti dataset
# comm = comm[-1:-5,]
bnti = bnti[,-1:-4]

# # Fixing the data to be numbers
# comm[,-1] = apply(comm[,-1], 2, function(x) as.numeric(as.character(x)))
# 
# # Aggregating the data because duplicates need to be removed
# comm = aggregate(.~OTU.ID, data = comm, FUN = sum)
# 
# # Setting row names on the newly aggregated data
# row.names(comm) = comm[,1]
# comm = comm[,-1]
# 
# # Transposing the community matrix
comm = t(comm)

# Fixing Dates
data$Time = gsub("_","/1/",data$Time)
data$Time = as.Date(data$Time, "%m/%d/%Y")

rm("fac.file")

######################################
#### Averaging duplicated samples ####
######################################

x = unique(factors[,1])
comm.avg = matrix(nrow = length(x), ncol = length(comm[1,]))

# Generating the averaged Family level table
for(i in 1:length(x)){
  w = which(factors[,1] %in% x[i])
  if(length(w)==1){
    comm.avg[i,] = comm[w,]
  } else{
    comm.avg[i,] = apply(comm[w,], 2, mean)
  }
}

row.names(comm.avg) = x
colnames(comm.avg) = colnames(comm)
comm.avg = as.data.frame(comm.avg)

comm = comm.avg

# Generating the corresponding factors sheet
w = which(factors[,1] %in% x[1])[1]
factors.avg = factors[w,]

for(i in 2:length(x)){
  w = which(factors[,1] %in% x[i])[1]
  factors.avg = rbind(factors.avg, factors[w,])
}

row.names(factors.avg) = factors.avg[,1]
factors.avg = factors.avg[,-1]

factors = as.data.frame(factors.avg)

rm('comm.avg', 'factors.avg')

##################################
### Plotting the cohesion data ###
##################################

# Plotting the data by time
ggplot(data = data, aes(x = Time, y = Positive_Cohesion))+
  geom_point(aes(colour = Well))+
  geom_line(aes(colour = Well))+
  theme_bw()+
  ggtitle("Positive Cohesion thru Time")

ggplot(data = data, aes(x = Time, y = Negative_Cohesion))+
  geom_point(aes(colour = Well))+
  geom_line(aes(colour = Well))+
  theme_bw()+
  ggtitle("Negative Cohesion thru Time")

# Generating summary boxplots for Cohesion data
ggplot(data = data, aes(x = Well, y = Positive_Cohesion))+
  geom_boxplot(aes(group = Well, fill = Well))+
  theme_bw()+
  ggtitle("Positive Cohesion by Well")

ggplot(data = data, aes(x = Well, y = Negative_Cohesion))+
  geom_boxplot(aes(group = Well, fill = Well))+
  theme_bw()+
  ggtitle("Negative Cohesion by Well")

#########################################################
### Relating the cohesion value to community turnover ###
#########################################################

# Removing the 0.1 fraction of the data
w = c(which(factors$Filter %in% "0.1"), which(factors$Well %in% "Licking-PW"))
comm = comm[-w,]
factors = factors[-w,]
bnti = bnti[-w,-w]

# Generating a Bray-Curtis distance matrix
comm = as.matrix(vegdist(x = comm, method = "bray"))

#--------------------------------------------------------------------#

### Attempting to plot Bray-Curtis and bNTI by Cohesion

# Plotting Positive
x = NULL
y = NULL

w = grep("Athens", row.names(data))
ww = grep("Athens", row.names(comm))

for(i in 1:length(ww)){
  x = rbind(x, cbind(data$Positive_Cohesion[w], comm[ww,ww[i]]))
  y = rbind(y, cbind(data$Positive_Cohesion[w], bnti[ww,ww[i]]))
}

w = grep("Greene", row.names(data))
ww = grep("Greene", row.names(comm))

for(i in 1:length(ww)){
  x = rbind(x, cbind(data$Positive_Cohesion[w], comm[ww,ww[i]]))
  y = rbind(y, cbind(data$Positive_Cohesion[w], bnti[ww,ww[i]]))
}

w = grep("Licking", row.names(data))
ww = grep("Licking", row.names(comm))

for(i in 1:length(ww)){
  x = rbind(x, cbind(data$Positive_Cohesion[w], comm[ww,ww[i]]))
  y = rbind(y, cbind(data$Positive_Cohesion[w], bnti[ww,ww[i]]))
}

x = x[-which(x[,2] == 0),]
y = y[-is.na(y[,2]),]

plot(x[,1], x[,2], main = "Bray-Curtis by Positive Cohesion")
abline(lm(x[,2]~x[,1]))
corr1 = rcorr(x[,1], x[,2])
legend("bottomright", sprintf("r-value: %s, p-value: %s", corr1$r[1,2], corr1$P[1,2]))

plot(y[,1], y[,2], main = "bNTI by Positive Cohesion")
abline(lm(y[,2]~y[,1]))
corr2 = rcorr(y[,1], y[,2])
legend("bottomright", sprintf("r-value: %s, p-value: %s", corr2$r[1,2], corr2$P[1,2]))

boxplot(x[,2]~x[,1])
boxplot(y[,2]~y[,1])

# Plotting Negative
x = NULL
y = NULL

w = grep("Athens", row.names(data))
ww = grep("Athens", row.names(comm))

for(i in 1:length(ww)){
  x = rbind(x, cbind(data$Negative_Cohesion[w], comm[ww,ww[i]]))
  y = rbind(y, cbind(data$Negative_Cohesion[w], bnti[ww,ww[i]]))
}

w = grep("Greene", row.names(data))
ww = grep("Greene", row.names(comm))

for(i in 1:length(ww)){
  x = rbind(x, cbind(data$Negative_Cohesion[w], comm[ww,ww[i]]))
  y = rbind(y, cbind(data$Negative_Cohesion[w], bnti[ww,ww[i]]))
}

w = grep("Licking", row.names(data))
ww = grep("Licking", row.names(comm))

for(i in 1:length(ww)){
  x = rbind(x, cbind(data$Negative_Cohesion[w], comm[ww,ww[i]]))
  y = rbind(y, cbind(data$Negative_Cohesion[w], bnti[ww,ww[i]]))
}

x = x[-which(x[,2] == 0),]
y = y[-which(is.na(y[,2]) %in% TRUE),]

plot(x[,1], x[,2], main = "Bray-Curtis by Negative Cohesion")
abline(lm(x[,2]~x[,1]))
corr1 = rcorr(x[,1], x[,2])
legend("bottomright", sprintf("r-value: %s, p-value: %s", corr1$r[1,2], corr1$P[1,2]))

plot(y[,1], y[,2], main = "bNTI by Negative Cohesion")
abline(lm(y[,2]~y[,1]))
corr2 = rcorr(y[,1], y[,2])
legend("bottomright", sprintf("r-value: %s, p-value: %s", corr2$r[1,2], corr2$P[1,2]))

boxplot(x[,2]~x[,1])
boxplot(y[,2]~y[,1])

rm("x", "y")

#######################################################
### Working on some different statistics approaches ###
#######################################################

### Attempting to use Mantel to solve my problem ###
# Using difference in cohesion between samples to relate it to bNTI and Bray-Curtis
bnti[upper.tri(bnti)]=t(bnti)[upper.tri(bnti)]
data = data[match(row.names(bnti), row.names(data)),] # Reordering data to match bNTI and Bray-Curtis

pos.coh = abs(data$Positive_Cohesion) # Absolute value for positive cohesion values (shouldn't matter but doing this anyway)
pos.coh = abs(outer(pos.coh, pos.coh, "-")) # Generating matrix of differences
row.names(pos.coh) = row.names(data)
colnames(pos.coh) = row.names(data)

neg.coh = abs(data$Negative_Cohesion) # Absolute value for positive cohesion values (shouldn't matter but doing this anyway)
neg.coh = abs(outer(neg.coh, neg.coh, "-")) # Generating matrix of differences
row.names(neg.coh) = row.names(data)
colnames(neg.coh) = row.names(data)

# Total Mantel analysis
mant = mantel(formula = as.dist(bnti)~as.dist(neg.coh), nperm = 9999, mrank = T)
mant = rbind(mant, mantel(formula = as.dist(bnti)~as.dist(pos.coh), nperm = 9999, mrank = T))
mant = rbind(mant, mantel(formula = as.dist(uw.uni)~as.dist(neg.coh), nperm = 9999, mrank = T))
mant = rbind(mant, mantel(formula = as.dist(uw.uni)~as.dist(pos.coh), nperm = 9999, mrank = T))
mant = rbind(mant, mantel(formula = as.dist(w.uni)~as.dist(neg.coh), nperm = 9999, mrank = T))
mant = rbind(mant, mantel(formula = as.dist(w.uni)~as.dist(pos.coh), nperm = 9999, mrank = T))
mant = rbind(mant, mantel(formula = as.dist(comm)~as.dist(neg.coh), nperm = 9999, mrank = T))
mant = rbind(mant, mantel(formula = as.dist(comm)~as.dist(pos.coh), nperm = 9999, mrank = T))

# By-well Mantel analysis
n = c("Athens", "Greene", "Licking")

r.names = NULL
for(i in 1:length(n)){
  w = grep(n[i], row.names(bnti))
  
  temp = mantel(formula = as.dist(bnti[w,w])~as.dist(neg.coh[w,w]), nperm = 9999, mrank = T)
  temp = rbind(temp, mantel(formula = as.dist(bnti[w,w])~as.dist(pos.coh[w,w]), nperm = 9999, mrank = T))
  temp = rbind(temp, mantel(formula = as.dist(uw.uni[w,w])~as.dist(neg.coh[w,w]), nperm = 9999, mrank = T))
  temp = rbind(temp, mantel(formula = as.dist(uw.uni[w,w])~as.dist(pos.coh[w,w]), nperm = 9999, mrank = T))
  temp = rbind(temp, mantel(formula = as.dist(w.uni[w,w])~as.dist(neg.coh[w,w]), nperm = 9999, mrank = T))
  temp = rbind(temp, mantel(formula = as.dist(w.uni[w,w])~as.dist(pos.coh[w,w]), nperm = 9999, mrank = T))
  temp = rbind(temp, mantel(formula = as.dist(comm[w,w])~as.dist(neg.coh[w,w]), nperm = 9999, mrank = T))
  temp = rbind(temp, mantel(formula = as.dist(comm[w,w])~as.dist(pos.coh[w,w]), nperm = 9999, mrank = T))
  
  mant = rbind(mant, temp)
  r.names = c(r.names, paste(n[i],"bNTI-Neg"), paste(n[i],"bNTI-Pos"), paste(n[i],"UWU-Neg"), paste(n[i],"UWU-Pos"),
              paste(n[i],"WU-Neg"), paste(n[i],"WU-Pos"), paste(n[i],"BC-Neg"), paste(n[i],"BC-Pos"))
}


row.names(mant) = c("Total bNTI-Neg", "Total bNTI-Pos", "Total UWU-Neg", "Total UWU-Pos", 
                    "Total WU-Neg", "Total WU-Pos", "Total BC-Neg", "Total BC-Pos",
                    r.names)

### Multivariate comparisons between cohesion and turnover measurements ###
col = c(rgb(190,38,37, maxColorValue = 255), rgb(0,97,28, maxColorValue = 255), rgb(13,79,139, maxColorValue = 255)) # Generates vector of consistent colors for the NMDS
shape = c(15,17,16) # Generating consistent shapes for use in the NMDS plot

n = c("Athens", "Greene", "Licking")

# Plotting BC
y = prcomp(as.dist(comm)) # Generating principle component analysis for Bray-Curtis
ordiplot(y, type = "n", xlab = paste("PC1 (", (summary(y)$importance[2,1]*100), "% )", sep = ""), 
         ylab = paste("PC2 (", (summary(y)$importance[2,2]*100), "%)", sep = "")) # Blank plot

for(i in 1:length(n)){
  w = grep(n[i], row.names(y$x))
  points(x = y$x[w,1], y =y$x[w,2], col = col[i], pch = shape[i])
  # text(x = y$x[w,1], y =y$x[w,2], labels = row.names(y$x)[w]) # Comment this out usually, just making sure my points are colored correctly
} # Plotting graph

x = envfit(y, data[,3:4], permutations = 999) # Generating the environmental fit for cohesion values
plot(x, p.max = 0.05) # Plotting only significant envfit results

multi.coh = cbind(x$vectors$r, x$vectors$pvals) # Storing all envfit results

# Plotting Unweighted UniFrac
y = prcomp(as.dist(uw.uni)) # Generating principle component analysis for Bray-Curtis
ordiplot(y, type = "n", xlab = paste("PC1 (", (summary(y)$importance[2,1]*100), "% )", sep = ""), 
         ylab = paste("PC2 (", (summary(y)$importance[2,2]*100), "%)", sep = "")) # Blank plot

for(i in 1:length(n)){
  w = grep(n[i], row.names(y$x))
  points(x = y$x[w,1], y =y$x[w,2], col = col[i], pch = shape[i])
  # text(x = y$x[w,1], y =y$x[w,2], labels = row.names(y$x)[w])
} # Plotting graph

x = envfit(y, data[,3:4], permutations = 999) # Generating the environmental fit for cohesion values
plot(x, p.max = 0.05) # Plotting only significant envfit results

temp = cbind(x$vectors$r, x$vectors$pvals) # Storing all envfit results
multi.coh = rbind(multi.coh, temp) # Merging them with the previous results

# Plotting Weighted UniFrac
y = prcomp(as.dist(w.uni)) # Generating principle component analysis for Bray-Curtis
ordiplot(y, type = "n", xlab = paste("PC1 (", (summary(y)$importance[2,1]*100), "% )", sep = ""), 
         ylab = paste("PC2 (", (summary(y)$importance[2,2]*100), "%)", sep = "")) # Blank plot

for(i in 1:length(n)){
  w = grep(n[i], row.names(y$x))
  points(x = y$x[w,1], y =y$x[w,2], col = col[i], pch = shape[i])
  # text(x = y$x[w,1], y =y$x[w,2], labels = row.names(y$x)[w])
} # Plotting graph

x = envfit(y, data[,3:4], permutations = 999) # Generating the environmental fit for cohesion values
plot(x, p.max = 0.05) # Plotting only significant envfit results

temp = cbind(x$vectors$r, x$vectors$pvals) # Storing all envfit results
multi.coh = rbind(multi.coh, temp) # Merging them with the previous results

# Plotting bNTI
y = prcomp(as.dist(bnti)) # Generating principle component analysis for Bray-Curtis
ordiplot(y, type = "n", xlab = paste("PC1 (", (summary(y)$importance[2,1]*100), "% )", sep = ""), 
         ylab = paste("PC2 (", (summary(y)$importance[2,2]*100), "%)", sep = "")) # Blank plot

for(i in 1:length(n)){
  w = grep(n[i], row.names(y$x))
  points(x = y$x[w,1], y =y$x[w,2], col = col[i], pch = shape[i])
  # text(x = y$x[w,1], y =y$x[w,2], labels = row.names(y$x)[w])
} # Plotting graph

x = envfit(y, data[,3:4], permutations = 999) # Generating the environmental fit for cohesion values
plot(x, p.max = 0.05) # Plotting only significant envfit results

temp = cbind(x$vectors$r, x$vectors$pvals) # Storing all envfit results
multi.coh = rbind(multi.coh, temp) # Merging them with the previous results

row.names(multi.coh) = c("Bray-Curtis Negative", "Bray-Curtis Postive", "Unweighted Negative", "Unweighted Postive",
                         "Weighted Negative", "Weighted Postive", "bNTI Negative", "bNTI Postive")
colnames(multi.coh) = c("r-value", "pval")

### Selecting the previous timepoint comparisons to correlate, similar to Herren and McMahon ###

data = data[order(data$Well, data$Time),]
bnti = bnti[match(row.names(data), row.names(bnti)), match(row.names(data), row.names(bnti))]
comm = comm[match(row.names(data), row.names(comm)), match(row.names(data), row.names(comm))]

n = c("Athens", "Greene", "Licking")
bnti.vec = NULL
comm.vec = NULL

for(i in 1:length(n)){
  w = grep(n[i], row.names(bnti))
  
  temp.bnti = matrix(nrow = length(w)-1, ncol = 3)
  temp.comm = matrix(nrow = length(w)-1, ncol = 3)
  
  for(j in 1:length(w)-1){
    temp.bnti[j,1] = as.numeric(bnti[j+1,j])
    temp.bnti[j,2] = as.numeric(data[j+1,3])
    temp.bnti[j,3] = as.numeric(data[j+1,4])
    
    temp.comm[j,1] = as.numeric(comm[j+1,j])
    temp.comm[j,2] = as.numeric(data[j+1,3])
    temp.comm[j,3] = as.numeric(data[j+1,4])
  }
  
  bnti.vec = rbind(bnti.vec, temp.bnti)
  comm.vec = rbind(comm.vec, temp.comm)
  
  
  
  rm("temp.comm", "temp.bnti")
}

comm.mod = lm(comm.vec[,1]~comm.vec[,2])
bnti.mod = lm(bnti.vec[,1]~bnti.vec[,2])

### Looking at differences between cohesion distributions
sum.stat = NULL
r.names = NULL

temp = kruskal.test(data$Positive_Cohesion ~ data$Well)
sum.stat = rbind(sum.stat, c(temp$statistic, temp$p.value))
temp = kruskal.test(data$Negative_Cohesion ~ data$Well)
sum.stat = rbind(sum.stat, c(temp$statistic, temp$p.value))

for(i in 1:(length(n)-1)){
  for(j in 1:length(n)){
    if(i ==j){
      print("Skipping self-comparison")
    } else {
      w = c(grep(n[i], row.names(data)), grep(n[j], row.names(data)))
      
      temp = wilcox.test(data$Positive_Cohesion[w] ~ data$Well[w])
      sum.stat = rbind(sum.stat, c(temp$statistic, temp$p.value))
      temp = wilcox.test(data$Negative_Cohesion[w] ~ data$Well[w])
      sum.stat = rbind(sum.stat, c(temp$statistic, temp$p.value))
      
      r.names = c(r.names, paste(n[i], "-", n[j], "Pos"), paste(n[i], "-", n[j], "Neg"))
    }
  }
}

row.names(sum.stat) = c("Overall Pos", "Overall Neg", r.names)
