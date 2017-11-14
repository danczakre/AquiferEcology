### Processing the cohesion data
# RED 2017; danczak.6@osu.edu

# Loading in packages
library(reshape2)
library(ggplot2)
library(vegan)
library(yarrr)

####################################
### Loading and sorting the data ###
####################################

# Setting work directory and loading data
setwd("~/Documents/Wilkins Lab/Ohio EPA Project/Sequencing Data/16S Data/Ecological Modeling Files/")
data = read.csv("Total_0.2_Cohesion.csv", row.names = 1)
comm = read.csv("Family_Level_OTU_edit.csv")
bnti = read.csv("MC_EPA_bNTI.csv", row.names = 1)

# Storing the factors
factors = t(comm[1:5,])
colnames(factors) = factors[1,]
factors = factors[-1,]

# Removing the factors from the community matrix and bnti dataset
comm = comm[-1:-5,]
bnti = bnti[,-1:-4]

# Fixing the data to be numbers
comm[,-1] = apply(comm[,-1], 2, function(x) as.numeric(as.character(x)))

# Aggregating the data because duplicates need to be removed
comm = aggregate(.~OTU.ID, data = comm, FUN = sum)

# Setting row names on the newly aggregated data
row.names(comm) = comm[,1]
comm = comm[,-1]

# Transposing the community matrix
comm = t(comm)

# Fixing Dates
data$Time = gsub("_","/1/",data$Time)
data$Time = as.Date(data$Time, "%m/%d/%Y")

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

pirateplot(Positive_Cohesion~Well, data = data, pal = c(rgb(190,38,37, maxColorValue = 255),
                                                        rgb(0,97,28, maxColorValue = 255),
                                                        rgb(13,79,139, maxColorValue = 255)),
           point.o = 1, inf.f.o = 0.2, bean.b.o = 0, bean.f.o = 0.75, main = "Positive Cohesion by Well")


ggplot(data = data, aes(x = Well, y = Negative_Cohesion))+
  geom_boxplot(aes(group = Well, fill = Well))+
  theme_bw()+
  ggtitle("Negative Cohesion by Well")

pirateplot(Negative_Cohesion~Well, data = data, pal = c(rgb(190,38,37, maxColorValue = 255),
                                                        rgb(0,97,28, maxColorValue = 255),
                                                        rgb(13,79,139, maxColorValue = 255)),
           point.o = 1, inf.f.o = 0.2, bean.b.o = 0, bean.f.o = 0.75, main = "Negative Cohesion by Well")


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

