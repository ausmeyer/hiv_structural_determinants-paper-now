rm(list = ls())

map <- read.table('aas.fasta_full_length_monomer1.map', sep=',', head=F, stringsAsFactors=F)
rsa <- read.table('integrase_full_length.rsa', head=T, stringsAsFactors=F)$RSA
rsa <- rsa[map$V3!='-']
distances <- read.table('distances_integrase_full_length_monomer1.dat', head=F, sep=',', stringsAsFactors=T)
distances <- distances[map$V3!='-', map$V3!='-']


rates <- read.table('sites.dat', head=T, sep=',', stringsAsFactors=F)
keep.alignment.sites <- map$V3[map$V3 != '-']
dN.dS <- rates$dN.dS[as.numeric(keep.alignment.sites)]

best.site <- c()
best.r <- c()
rfree <- c()

print(dim(map))
print(length(dN.dS))
print(dim(distances))

for(i in 1:100) {
  small.set <- sample(1:length(dN.dS), length(dN.dS)*0.75)
  fit.dN.dS <- dN.dS[small.set]
  fit.distances <- distances[small.set, ]
  fit.rsa <- rsa[small.set]
  
  fit.correlations <- as.vector(sapply(1:ncol(fit.distances), function(x) cor(fit.dN.dS, 
                                                                              predict(lm(fit.dN.dS ~ fit.rsa + fit.distances[, x])))))
  
  best <- which(fit.correlations^2 == max(fit.correlations^2))
  best.site <- append(best.site, best)
  
  free.dN.dS <- dN.dS[-small.set]
  free.distances <- distances[-small.set, best]
  free.rsa <- rsa[-small.set]
  
  free.correlations <- cor(free.dN.dS, 
                           predict(lm(free.dN.dS ~ free.rsa + free.distances)))
  
  best.r <- append(best.r, max(fit.correlations^2))
  rfree <- append(rfree, max(free.correlations^2))
}

print(mean(best.r))
print(mean(rfree))
print(cor(dN.dS, rsa)^2)
fit.rsa <- lm(dN.dS ~ rsa)
print(summary(fit.rsa)$r.squared)

print(table(best.site))

fit.site <- as.numeric(names(sort(-table(best.site)))[1])
fit <- lm(dN.dS ~ rsa + distances[, fit.site])
#print(summary(fit))

fit.distance <- lm(dN.dS ~ distances[, fit.site])
print(summary(fit.distance))

correlations <- as.vector(sapply(1:ncol(distances), function(x) cor(dN.dS, 
                                                                    distances[, x])))
new.correlations <- rep(median(correlations), length(map$V3))
new.correlations[as.numeric(map$V2[map$V3 != '-'])] <- correlations
write.table(data.frame(new.correlations), file= 'distance_model.correlations', row.names=F, col.names=F)

correlations <- as.vector(sapply(1:ncol(distances), function(x) cor(dN.dS, 
                                                                    predict(lm(dN.dS ~ rsa + distances[, x])))))
new.correlations <- rep(median(correlations), length(map$V3))
new.correlations[as.numeric(map$V2[map$V3 != '-'])] <- correlations
write.table(data.frame(new.correlations), file= 'combined_model.correlations', row.names=F, col.names=F)

library(ggplot2)
library(cowplot)

p <- ggplot(data.frame(x=1:length(dN.dS), y=dN.dS), aes(x=x, y=y)) + geom_point()
show(p)
