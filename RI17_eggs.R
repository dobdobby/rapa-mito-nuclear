rm(list=ls())
library(lme4)
library(MASS)
library(car)
library(ggh4x)
library(lme4)
library(MASS)
library(car)

## import data
data2 <- read.csv("RI17_eggs.csv", header=TRUE)
data2 <- na.omit(data2)

names(data2)[names(data2) == 'genotype'] <- 'line'

# create a geno variable
data2$geno <- gsub('[[:digit:]]+', '', data2$line)

#recode rapa
data2$rapa <- as.factor(data2$rapa)
levels(data2$rapa) <- c("Control", "+Rapamycin")

# create a Number of eggs per fly variable
data2$Neggsflies <- with(data2, Neggs / Nfemales)

## PLOT N EGGS

# create a Number of eggs per fly variable
data2$Neggsflies <- with(data2, Neggs / Nfemales)

p1 <- ggplot(data=data2, aes(y=Neggsflies+1, x=rapa))
set.seed(1984)
p2 <- p1 + 
	geom_boxplot(outlier.shape=NA) + 
	geom_jitter(width=0.2, aes(col=rapa), alpha=0.5) +
	scale_y_continuous(trans='log10') +
	facet_grid(. ~ geno) + 
	theme_bw() +
	scale_color_manual(values=c("grey", "darksalmon"), name="") +
	ylab("1 + N. eggs / fly") +
	xlab("")

pdf(h=2, w=8, "eggs.pdf")
print(p2)
dev.off()

negbinomial3 <- MASS::glm.nb(Neggs ~ rapa*mito*nuclear, data= data2, na.action="na.fail", contrasts=list(rapa="contr.sum", mito="contr.sum", nuclear="contr.sum"), control = glm.control(maxit = 5000))
negbinomial3
par(mfrow=c(2,2))
plot(negbinomial3)
Anova(negbinomial3, type="III")