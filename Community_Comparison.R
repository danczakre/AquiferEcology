### Thinking about different ways to analyze the core community ###
# RED 2017; danczak.6@osu.edu

Var.Num = 5
filter = TRUE # Switch to determine whether or not I'm removing low abundance OTUs

library(plyr)
library(data.table)

####################################
#### Data Import and Processing ####
####################################

# Setting up the variables
setwd("~/Documents/Wilkins Lab/Ohio EPA Project/Sequencing Data/16S Data") # Working directory for 16S comm
data = read.csv("Family_Level_OTU_edit.csv") # Loading a previoulsy condensed OTU table at the family level
otu = read.delim("Raw_EPA_no_pynast.txt") # Loading in the raw OTU table for OTU level analyses
tax = read.delim("rep_set_tax_assignments.txt", header = F) # Loading in taxonomy data to assign taxonmies to OTUs

# Creating the factor list
factors = t(data[1:Var.Num,])
colnames(factors) = factors[1,]
factors = as.data.frame(factors[-1,])
data = data[-1:-Var.Num,]

## Data Variable

# Aggregating the newly assigned taxonomies for the 'data' variable
data[,-1] = apply(data[,-1], 2, function(x) as.numeric(as.character(x))) # Converting everything except the OTU.IDs into numbers
data = aggregate(. ~ OTU.ID, data = data, FUN = sum) # Need to aggregate the newly assigned taxonomies (I replaced OD1 things with Parcubacteria)

# Setting up row names for the 'data' variable
row.names(data) = data[,1] # Setting row names
data = data[,-1] # Removing the vector in which the row names were stored

data = t(data) # Transposing data so it is actually useful
data = apply(data, 2, function(x) as.numeric(as.character(x))) # Reconverts the data to numbers from characters
row.names(data) = row.names(factors)

# Cleaning the 'otu' varaible
otu = otu[,-which(colnames(otu) %in% c("YB04","BL01","BL02"))] # Removing YB04 becuase it was a failed sequencing run; the others are blanks

## OTU Variable

# Setting up the row names for the 'otu' variable
row.names(otu) = otu[,1] # Setting the the first column to the row.names and then removing them
otu = otu[,-1]

otu = t(otu) # Transposing data so it is actually useful
otu = apply(otu, 2, function(x) as.numeric(as.character(x))) # Reconverts the data to numbers from characters
row.names(otu) = row.names(factors)

############################################################
#### Averaging the duplicates from both of the datasets ####
############################################################

## Processing the 'data' table
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

# Removing members with a max abundance of 0
if(length(which(apply(data.avg, 2, max) == 0)) > 0){
  data.avg = data.avg[,-which(apply(data.avg, 2, max) == 0)]
}


## Processing the 'otu' table
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

row.names(otu.avg) = row.names(factors.avg)
colnames(otu.avg) = colnames(otu)

rm('data','otu') # Clearing unecessary data


##################################
#### Shared community members ####
##################################

# Calculating relative abundances
otu.avg = apply(otu.avg, 1, function(x) {x/sum(x)}) # Turning abundances to relative abundances
otu.avg = t(otu.avg)

# Removing OTUs with fewer than 0.05% relative abundance
if(filter == TRUE){
  w = which(apply(otu.avg, 2, max) > 5e-4)
  otu.avg = otu.avg[,w]
}

# Working with the 'otu' table first
otu.agg = as.data.frame(otu.avg)
otu.agg$Well = factors.avg$`Well/Filter` # Setting the well/filter classifications as a column in the variable
otu.agg = aggregate(.~Well, data = otu.agg, FUN = mean) # Aggregating it based on well/filter

# Reseting row names to the new ones
row.names(otu.agg) = otu.agg$Well
otu.agg = otu.agg[,-1]


## Looking at the comparison between 0.1 and 0.2 filters
filter = NULL
filter.rn = NULL
filter.uniq = unique(factors$Well)

for(i in 1:length(filter.uniq)){
  if(row.names(otu.agg)[i] == c("Li-PW 0.1") | row.names(otu.agg)[i] == c("Li-PW 0.2")){
    w = which(factors.avg$Well %in% filter.uniq[i])
    x = as.data.frame(otu.avg[w,])
    
    x = x[,which(apply(x, 2, max) != 0)]
    x = apply(x, 2, function(x) as.numeric(x))
    
    row.names(x) = x[,1]
    x = x[,-1]
    
    x = rbind(x, apply(x, 2, diff))
    
    filter[[i]] = x
    names(filter)[[i]] = as.character(filter.uniq[i])
    
  } else {
    w = which(factors.avg$Well %in% filter.uniq[i])
    x = as.data.frame(otu.avg[w,])
    x = data.table(x)
    x$Filter = factors.avg$Filter[w]
    x = melt(x, id.vars = "Filter")
    x = dcast.data.table(x, Filter~variable, value.var = "value", fun.aggregate = mean, na.rm = T)
    x = as.data.frame(x)
    
    x = x[,which(apply(x, 2, max) != 0)]
    x = apply(x, 2, function(x) as.numeric(x))
    
    row.names(x) = x[,1]
    x = x[,-1]
    
    x = rbind(x, apply(x, 2, diff))
    
    filter[[i]] = x
    names(filter)[[i]] = as.character(filter.uniq[i])
  }
}

# Writing results

# for(i in 1:length(filter)){
#   write.csv(t(filter[[i]]), sprintf("Filter_Comparison_%s.csv", names(filter)[i]), quote = F)
# }

## Looking at common membership through time
timed = NULL
timed.rn = NULL

for(i in 1:length(row.names(otu.agg))){
  if(row.names(otu.agg)[i] == c("Li-PW 0.1") | row.names(otu.agg)[i] == c("Li-PW 0.2")){
    print(sprintf("Skipped %s due to one time point", row.names(otu.agg)[i])) # Skipping comparisons bewtween identical groupings
  } else {
    x = otu.avg[which(factors.avg$`Well/Filter` %in% row.names(otu.agg)[i]),]
    y = as.data.frame(rbind((colnames(x[,apply(x, 2, min) != 0]))))
    timed.rn = c(timed.rn, row.names(otu.agg)[i])
    if(length(y) == 0){
      y = as.data.frame(rbind(c("None")))
      timed = rbind.fill(timed, y)
    } else {
      timed = rbind.fill(timed, y)
    }
  }
}

timed = apply(timed, 2, function(x) as.character(x))
row.names(timed) = timed.rn

## Looking at community composition between wells

# All we are doing here is finding the OTUs which are greater than zero
w = grep(0.2, row.names(otu.agg)) # Only interested in the 0.2 for now; can change this at any moment though
OTUs = as.data.frame(rbind(colnames(otu.agg[w[1], which(otu.agg[w[1],] > 0)])))
OTUs.leng = length(OTUs)

for(i in 2:length(w)){
  x = as.data.frame(rbind(colnames(otu.agg[w[i], which(otu.agg[w[i],] > 0)])))
  OTUs.leng = c(OTUs.leng, length(x))
  OTUs = rbind.fill(OTUs, x)
}

OTUs = apply(OTUs, 2, function(x) as.character(x))
row.names(OTUs) = row.names(otu.agg[w,])

# Going through and finding the 'intersect' of all combinations
inter = NULL
inter.rn = NULL

for(i in 1:length(OTUs[,1])){
  for(n in i:length(OTUs[,1])){
    if(i==n){
      
      print(sprintf("Skipped comparison between %s and %s", row.names(OTUs)[i], row.names(OTUs)[n])) # Skipping comparisons bewtween identical groupings
    
    } else {
      
      x = OTUs[i,!is.na(OTUs[i,])] # One-liners got too long for me to continue bothering with them
      y = OTUs[n,!is.na(OTUs[n,])]
      z = as.data.frame(rbind(intersect(x, y))) # Finding the locations for unique factor 'i' and then appending the locations for 'n' to it

      inter.rn = c(inter.rn, sprintf("%s -> %s", row.names(OTUs)[i], row.names(OTUs)[n]))
      inter = rbind.fill(inter, z)
      
      for(m in n:length(OTUs[,1])){
        if(n==m){

          print(sprintf("Skipped comparison between %s, %s and %s", row.names(OTUs)[i], row.names(OTUs)[n], row.names(OTUs)[m])) # Skipping comparisons bewtween identical groupings

        } else {

          x = OTUs[i,!is.na(OTUs[i,])] # One-liners got too long for me to continue bothering with them
          y = OTUs[n,!is.na(OTUs[n,])]
          z = OTUs[m,!is.na(OTUs[m,])]

          a = as.data.frame(rbind(intersect(intersect(x, y),z))) # Finding the locations for unique factor 'i' and then appending the locations for 'n' to it
          b = as.data.frame(rbind(c(sprintf("%s", row.names(OTUs)[i]), sprintf("%s", row.names(OTUs)[n]), sprintf("%s", row.names(OTUs)[m]))))
          
          inter.rn = c(inter.rn, sprintf("%s -> %s -> %s", row.names(OTUs)[i], row.names(OTUs)[n], row.names(OTUs)[m]))
          inter = rbind.fill(inter, a)
        } # Inner 'if' loop
      } # End of the three comparison loop
    } # Outer 'if' loop
  } # Second comparison for-loop
} # First comparison for-loop

rm('x','y','z','a')

inter = apply(inter, 2, function(x) as.character(x))
row.names(inter) = inter.rn

######################################
#### Assigning taxonomies to OTUs ####
######################################

# Giviing taxonomies to the between well comparisons
inter.tax = NULL
inter.abn = NULL
inter.comp = strsplit(row.names(inter), split = " -> ")
inter.comp = t(sapply(inter.comp, "[", i = seq(max(sapply(inter.comp, length)))))

for(i in 1:length(inter[,1])){
  w = which(row.names(otu.agg) %in% inter.comp[i,]) # Selecting the right filter/well
  x = sapply(inter[i,!is.na(inter[i,])], function(x) which(tax[,1] %in% x)) # Finding the location linked to the OTU_ID
  y = as.data.frame(rbind(as.character(tax[x, 2]))) # Writing the taxonomy linked to the above location
  z = otu.agg[w, inter[i,!is.na(inter[i,])]] # Writing the averaged abundances for the OTUs
  yy = as.data.frame(rbind(apply(z, 2, sd)))
  inter.abn[[i]] = z # Storing the abundances
  
  colnames(y) = colnames(z)
  colnames(yy) = colnames(z)
  
  zz = rbind(z,yy,y) # Writing the taxonomies such that they are linked to abundances
  row.names(zz)[(length(z[,1])+1):length(zz[,1])] <- c("SD", "Tax")
  inter.tax[[i]] = zz
}

names(inter.tax) = apply(inter.comp, 1, function(x) paste(x[!is.na(x)], collapse = "-"))
names(inter.abn) = apply(inter.comp, 1, function(x) paste(x[!is.na(x)], collapse = "-"))

# Saving the between well comparisons
for(i in 1:length(names(inter.tax))){
  write.csv(t(inter.tax[[i]]), sprintf("Tax_between_%s.csv", names(inter.tax)[i]), quote = F)
}


# Assigning taxonomies to the through time results
timed.tax = NULL
timed.abn = NULL

for(i in 1:length(timed[,1])){
  if(length(which(timed[i,] %in% "None")) == 1){
    print(sprintf("Skipped %s due to no OTUs", row.names(timed)[i]))
    timed.abn[[i]] = matrix(data = sample(1:10), nrow = 10, ncol = 10) # Need to create a dummy matrix or else downstream analyses freak out a bit
  } else {
    x = grep(row.names(timed)[i], factors.avg$`Well/Filter`) # Variable to include only the well/filter for i
    xx = sapply(timed[i,!is.na(timed[i,])], function(x) which(tax[,1] %in% x)) # Getting the desired taxonomies
    y = as.data.frame(rbind(as.character(tax[xx, 2])))
    y = cbind(y, c("Taxonomy"))
    z = as.data.frame(otu.avg[x, timed[i, !is.na(timed[i,])]]) # Pulling the abundances for the desired OTUs
    yy = as.data.frame(rbind(apply(z, 2, sd)))
    yy = cbind(yy, c("StDev"))
    
    timed.abn[[i]] = z
    z$Date = factors.avg$Date[x]
    
    colnames(y) = colnames(z)
    colnames(yy) = colnames(z)
    
    zz = rbind(z,yy,y)
    row.names(zz)[(length(z[,1])+1):length(zz[,1])] <- c("SD", "Tax")
    timed.tax[[i]] = zz
  }
}

names(timed.tax) = row.names(timed)
names(timed.abn) = row.names(timed)

# Saving the timed results
for(i in 1:length(names(timed.tax))){
  if(is.null(timed.tax[[i]]) != T){
    write.csv(t(timed.tax[[i]]), sprintf("Tax_thru_time_%s.csv", names(timed.tax)[i]), quote = F)
  }
}

rm('x','xx','y','z','zz')