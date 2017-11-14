#### Testing out a simple dodged bar chart to demonstrate the distribution
# RED 2017; danczak.6@osu.edu

library(ggplot2)
library(reshape2)

col = c(rgb(190,38,37, maxColorValue = 255), rgb(0,97,28, maxColorValue = 255), rgb(13,79,139, maxColorValue = 255)) # Colour vector for plots

# Loading in data
setwd("~/Documents/Wilkins Lab/Ohio EPA Project/Sequencing Data/16S Data/Community Comparisons/")
data = read.csv("Tax_between_Athens 0.2-Greene 0.2-Licking 0.2.csv")

# Cleaning up the data to make it completely numbers, dumbing taxonomies into another vector, removing SD
tax = as.matrix(data[,length(data[1,])])
data = data[,-length(data[1,])] # Remove tax

data = data[,-length(data[1,])] # Removes SD

#
# Sorting the data by Athens
#
tax = as.matrix(tax[order(data[,2], decreasing = T),])
data = data[order(data[,2], decreasing = T),]

temp = data[1:10,]
temp = melt(data = temp, id.vars = "X")
temp$X = factor(temp$X, levels = temp$X)

ggplot(data = temp, aes(x = X, y = value))+
  geom_bar(stat = "identity", aes(fill = variable), position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_manual(values = col)+
  ggtitle("Top 10 Shared OTUs from Athens")

#
# Sorting the data by Greene
#
tax = as.matrix(tax[order(data[,3], decreasing = T),])
data = data[order(data[,3], decreasing = T),]

temp = data[1:10,]
temp = melt(data = temp, id.vars = "X")
temp$X = factor(temp$X, levels = temp$X)

ggplot(data = temp, aes(x = X, y = value))+
  geom_bar(stat = "identity", aes(fill = variable), position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_manual(values = col)+
  ggtitle("Top 10 Shared OTUs from Greene")

#
# Sorting the data by Licking
#
tax = as.matrix(tax[order(data[,4], decreasing = T),])
data = data[order(data[,4], decreasing = T),]

temp = data[1:10,]
temp = melt(data = temp, id.vars = "X")
temp$X = factor(temp$X, levels = temp$X)

ggplot(data = temp, aes(x = X, y = value))+
  geom_bar(stat = "identity", aes(fill = variable), position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_manual(values = col)+
  ggtitle("Top 10 Shared OTUs from Licking")
