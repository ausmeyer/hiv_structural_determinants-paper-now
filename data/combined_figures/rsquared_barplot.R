rm(list = ls())

library(ggplot2)
library(cowplot)

setwd('~/Google Drive/Documents/PostDoc/HIVStructuralDeterminants/hiv_structural_determinants/data/combined_figures/')

get.rates <- function(location) {
  rates <- read.table(paste(location, 'sites.dat', sep=''), sep=',', head=T)
  map <- read.table(paste(location, 'aas.fasta.map', sep=''), sep=',', head=F, stringsAsFactors=F)
  keep.alignment.sites <- map$V3[map$V3 != '-']
  return(rates$dN.dS[as.numeric(keep.alignment.sites)])
}

fit.gamma.distribution <- function(dat, par.pass) {
  fr <- function(par.pass) {
    expected <- diff(pgamma(seq(0, max(dat), by = 0.05), shape = par.pass[1], rate = par.pass[2]))*length(dat)
    observed <- as.numeric(table(cut(dat, breaks = seq(0, max(dat), by = 0.05))))
    return(sum((observed - expected)^2/expected))
  }
  
  fit <- optim(par.pass, fr, method="Nelder-Mead", control = list(maxit = 10000))
  return(fit)
}

fit.beta.distribution <- function(dat, par.pass) {
  fr <- function(par.pass) {
    expected <- diff(pbeta(seq(0, 1, by = 0.05), shape1 = par.pass[1], shape2 = par.pass[2]))*length(dat)
    observed <- as.numeric(table(cut(dat, breaks = seq(0, 1, by = 0.05))))
    return(sum((observed - expected)^2/expected))
  }
  
  fit <- optim(par.pass, fr, method="Nelder-Mead", control = list(maxit = 10000))
  return(fit)
}

plot.gamma.histogram <- function(data) {
  #data <- data[data <= 1]
  fit <- fit.gamma.distribution(data, c(1, 1)) 
  
  #Calculate chi-squared statistic
  expected <- diff(pgamma(seq(0, max(data), by = 0.05), shape = fit$par[1], rate = fit$par[2]))*length(data)
  observed <- as.numeric(table(cut(data, breaks = seq(0, max(data), by = 0.05))))
  chisq.stat <- sum((observed - expected)^2/expected)

  #Calculate degrees of freedom
  dof <- ceiling(max(data)/0.05) - 1 - 2
  
  print("Gamma")
  print(pchisq(q = chisq.stat, df = dof, lower.tail=FALSE))
  print(fit$par)

  p.tmp <- ggplot(data = data.frame(x=data), aes(x=x)) + 
    geom_histogram(aes(y=..density..), color='gray', fill='gray', binwidth = 0.05) +
    stat_function(fun=dgamma, args=list(shape = fit$par[1], rate = fit$par[2]), size=0.4) +
    scale_x_continuous(limits = c(0, 3)) +
    scale_y_continuous(limits = c(0, 12)) + 
    xlab('dN/dS') +
    ylab('Count') +
    geom_vline(xintercept = 1, linetype = "longdash")
  
  return(list(plot = p.tmp, fit = fit))
}

plot.beta.histogram <- function(data) {
  #data <- data[data <= 1]
  fit <- fit.beta.distribution(data, c(1, 1)) 
  
  #Calculate chi-squared statistic
  expected <- diff(pbeta(seq(0, 1, by = 0.05), shape1 = fit$par[1], shape2 = fit$par[2]))*length(data)
  observed <- as.numeric(table(cut(data, breaks = seq(0, 1, by = 0.05))))
  chisq.stat <- sum((observed - expected)^2/expected)
  
  #Calculate degrees of freedom
  dof <- ceiling(1/0.05) - 1 - 2
  
  print("Beta")
  print(pchisq(q = chisq.stat, df = dof, lower.tail=FALSE))
  print(fit$par)
  
  p.tmp <- ggplot(data = data.frame(x=data), aes(x=x)) + 
    geom_histogram(aes(y=..density..), color='gray', fill='gray', binwidth = 0.05) +
    stat_function(fun=dgamma, args=list(shape = fit$par[1], rate = fit$par[2]), size=0.4) +
    scale_x_continuous(limits = c(0, 3)) +
    scale_y_continuous(limits = c(0, 12)) + 
    xlab('dN/dS') +
    ylab('Count') +
    geom_vline(xintercept = 1, linetype = "longdash")
  
  return(list(plot = p.tmp, fit = fit))
}

r.value.1 <- 0.1001862
r.free.1 <- 0.09242819
r.rsa.1 <- 0.05076719
rates.1 <- get.rates('../capsid/capsid_structure/predictors/')
mean.dN.dS.1 <- mean(rates.1)
tmp.gamma <- plot.gamma.histogram(rates.1)
p.1 <- tmp.gamma$plot
tmp.beta <- plot.beta.histogram(rates.1)
p.1.beta <- tmp.beta$plot
fit.1 <- tmp.gamma$fit
r.1 <- "Capsid"

r.value.2 <- 0.3716488
r.free.2 <- 0.3637618
r.rsa.2 <- 0.1812471
rates.2 <- get.rates('../gp120/gp120_structure/predictors/')
mean.dN.dS.2 <- mean(rates.2)
tmp.gamma <- plot.gamma.histogram(rates.2)
p.2 <- tmp.gamma$plot
tmp.beta <- plot.beta.histogram(rates.2)
p.2.beta <- tmp.beta$plot
fit.2 <- tmp.gamma$fit
r.2 <- "gp120"

r.value.3 <- 0.1325805
r.free.3 <- 0.09904008
r.rsa.3 <- 0.0609476
rates.3 <- get.rates('../matrix/matrix_structure/predictors/')
mean.dN.dS.3 <- mean(rates.3)
tmp.gamma <- plot.gamma.histogram(rates.3)
p.3 <- tmp.gamma$plot
tmp.beta <- plot.beta.histogram(rates.3)
p.3.beta <- tmp.beta$plot
fit.3 <- tmp.gamma$fit
r.3 <- "Matrix"

r.value.4 <- 0.01859244
r.free.4 <- 0.03031527
r.rsa.4 <- 0.00638955
rates.4 <- get.rates('../integrase/integrase_structures/predictors/')
mean.dN.dS.4 <- mean(rates.4)
tmp.gamma <- plot.gamma.histogram(rates.4)
p.4 <- tmp.gamma$plot
tmp.beta <- plot.beta.histogram(rates.4)
p.4.beta <- tmp.beta$plot
fit.4 <- tmp.gamma$fit
r.4 <- "Integrase"

r.value.5 <- 0.3104171
r.free.5 <- 0.3069456
r.rsa.5 <- 0.04737184
rates.5 <- get.rates('../protease/protease_structure/predictors/')
mean.dN.dS.5 <- mean(rates.5)
tmp.gamma <- plot.gamma.histogram(rates.5)
p.5 <- tmp.gamma$plot
tmp.beta <- plot.beta.histogram(rates.5)
p.5.beta <- tmp.beta$plot
fit.5 <- tmp.gamma$fit
r.5 <- "Protease"

r.value.6 <- 0.06046455
r.free.6 <- 0.05723965
r.rsa.6 <- 0.04350685
rates.6 <- get.rates('../reverse_transcriptase/rt_structures/predictors/')
mean.dN.dS.6 <- mean(rates.6)
tmp.gamma <- plot.gamma.histogram(rates.6)
p.6 <- tmp.gamma$plot
tmp.beta <- plot.beta.histogram(rates.6)
p.6.beta <- tmp.beta$plot
fit.6 <- tmp.gamma$fit
r.6 <- "Reverse Transcriptase"

#print(mean(r.value.1, r.value.2, r.value.3, r.value.4, r.value.5, r.value.6))
#print(mean(r.free.1, r.free.2, r.free.3, r.free.4, r.free.5, r.free.6))
#print(mean(r.rsa.1, r.rsa.2, r.rsa.3, r.rsa.4, r.rsa.5, r.rsa.6))

rs <- c(r.value.1, r.free.1, r.rsa.1, 
        r.value.2, r.free.2, r.rsa.2, 
        r.value.3, r.free.3, r.rsa.3, 
        r.value.4, r.free.4, r.rsa.4, 
        r.value.5, r.free.5, r.rsa.5, 
        r.value.6, r.free.6, r.rsa.6)

r.names <- c(r.1, r.1, r.1,
             r.2, r.2, r.2,
             r.3, r.3, r.3,
             r.4, r.4, r.4, 
             r.5, r.5, r.5,
             r.6, r.6, r.6)
r.names <- factor(r.names, levels=r.names)

r.type <- rep(c('Combined - Training', 'Combined - Test', 'RSA - only'), 6)
r.type <- factor(r.type, levels=r.type)

df <- data.frame(r.square = rs, names = r.names, r.type = r.type, stringsAsFactors = F)

p <- ggplot(aes(x = names, y = r.square, fill=r.type, colour=NULL), data = df) +
  geom_bar(stat = 'identity', position=position_dodge()) +
  scale_fill_hue(name="Model Type") +
  ylab(expression(paste("Variance Explained (R"^"2", ')', sep=''))) +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(limits = c(0, 0.4))

ggsave("r_squared.png", p, width=7.5, height=7.5)
ggsave("r_squared.pdf", p, width=7.5, height=7.5)

df.comp <- data.frame(x = df$r.square[df$r.type == 'RSA - only'], y = df$r.square[df$r.type == 'Combined - Test'], Protein = df$names[df$r.type == 'Combined - Test'])

p <- ggplot(aes(x = x, y = y, colour=Protein), data = df.comp) +
  geom_point(size=4) +
  ylab(expression(paste("Combined Model - Test (R"^"2", ')', sep=''))) +
  xlab(expression(paste("RSA - only (R"^"2", ')', sep=''))) +
  scale_x_continuous(limits = c(0, 0.2)) +
  scale_y_continuous(limits = c(0, 0.4))

ggsave("combined_RSA.png", p, width=7.5, height=5.5)
ggsave("combined_RSA.pdf", p, width=7.5, height=5.5)

p <- plot_grid(p.1, p.2, p.3, p.4, p.5, p.6, ncol = 2, labels = c('A', 'B', 'C', 'D', 'E', 'F'))
ggsave("rate_distribution.png", p, width=8, height=12)
ggsave("rate_distribution.pdf", p, width=8, height=12)

p.beta <- plot_grid(p.1.beta, p.2.beta, p.3.beta, p.4.beta, p.5.beta, p.6.beta, ncol = 2, labels = c('A', 'B', 'C', 'D', 'E', 'F'))
ggsave("rate_distribution_beta.png", p.beta, width=8, height=12)
ggsave("rate_distribution_beta.pdf", p.beta, width=8, height=12)