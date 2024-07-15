rm(list=ls())

library(segmented)
library(car)
library(ggplot2)
library(emmeans)
 library(survival) # 3.5-5
  library(coxme) # 2.2-18.1
  library(car) #3.1-2
  library(dplyr) #1.1.2
  library(emmeans) #1.8.6
  library(ggplot2) #3.4.2
  library(ggh4x) #0.2.4
library(effectsize)
library(knitr)
library(ciTools)
library(here)
library(survminer)
library(rms)
library(see)

options(na.action="na.fail", contrasts=c("contr.sum", "contr.poly"))

AAc <- subset(read.csv("RI12_AA_ctrl.csv"), !is.na(ln.ux))
AAr <- subset(read.csv("RI12_AA_rapa.csv"), !is.na(ln.ux))

ABc <- subset(read.csv("RI12_AB_ctrl.csv"), !is.na(ln.ux))
ABr <- subset(read.csv("RI12_AB_rapa.csv"), !is.na(ln.ux))

BAc <- subset(read.csv("RI12_BA_ctrl.csv"), !is.na(ln.ux))
BAr <- subset(read.csv("RI12_BA_rapa.csv"), !is.na(ln.ux))

BBc <- subset(read.csv("RI12_BB_ctrl.csv"), !is.na(ln.ux))
BBr <- subset(read.csv("RI12_BB_rapa.csv"), !is.na(ln.ux))


AAc$geno <- "AA"
AAc$mito <- "A"
AAc$nuclear <- "A"
AAc$trt <- "control"
AAc$percTime <- 100 * with(AAc, time / max(time))

AAr$geno <- "AA"
AAr$mito <- "A"
AAr$nuclear <- "A"
AAr$trt <- "rapa"
AAr$percTime <- 100 * with(AAr, time / max(time))

ABc$geno <- "AB"
ABc$mito <- "A"
ABc$nuclear <- "B"
ABc$trt <- "control"
ABc$percTime <- 100 * with(ABc, time / max(time))

ABr$geno <- "AB"
ABr$mito <- "A"
ABr$nuclear <- "B"
ABr$trt <- "rapa"
ABr$percTime <- 100 * with(ABr, time / max(time))

BAc$geno <- "BA"
BAc$mito <- "B"
BAc$nuclear <- "A"
BAc$trt <- "control"
BAc$percTime <- 100 * with(BAc, time / max(time))

BAr$geno <- "BA"
BAr$mito <- "B"
BAr$nuclear <- "A"
BAr$trt <- "rapa"
BAr$percTime <- 100 * with(BAr, time / max(time))

BBc$geno <- "BB"
BBc$mito <- "B"
BBc$nuclear <- "B"
BBc$trt <- "control"
BBc$percTime <- 100 * with(BBc, time / max(time))

BBr$geno <- "BB"
BBr$mito <- "B"
BBr$nuclear <- "B"
BBr$trt <- "rapa"
BBr$percTime <- 100 * with(BBr, time / max(time))

head(AAc)
head(AAr)

head(ABc)
head(ABr)

head(BAc)
head(BAr)

head(BBc)
head(BBr)


d <- rbind(AAc, AAr, ABc, ABr, BAc, BAr, BBc, BBr)

str(d)

dim(d)

ggplot(d, aes(x=time, y=ln.ux, group=trt)) +
  geom_point() +
  facet_grid(rows=vars(mito), cols=vars(nuclear)) +
  geom_smooth()


##############################
##	segmented regression		##
##############################

	#rationale for choosing npsi ==> which of 1-3 gives lowest AIC?
pdf("segmented_regression_percTime.pdf")
par(mfrow=c(2,2))

##AAc
#AAcPlot
AAcM1 <- with(AAc, lm(ln.ux ~ percTime, weights=no.surv))
summary(AAcM1)
for(i in 1:3){
	print(
		AIC(
				segmented(AAcM1, seg.Z = ~percTime,npsi=i)
				)
			)
}
AAcM1_segmented <- segmented(AAcM1, seg.Z = ~percTime,npsi=2)
AAcM1_segmented
summary(AAcM1_segmented)
AAcM1_segmented_CI <- confint(AAcM1_segmented)
AAcM1_segmented_CI
plot(AAcM1_segmented, col=1, ylim=c(-8.5,0), main="AA", ylab="logN mortality rate", xlab="Relative age (% max. lifespan)")
lines(AAcM1_segmented, shift=T)
#abline(v=AAcM1_segmented_CI[,2:3], lty=3)
#abline(v=AAcM1_segmented_CI[,1], lty=2)
points(x=AAc$percTime, y=AAc$ln.ux, pch=16)
newx <- seq(min(AAc$percTime), max(AAc$percTime), length.out=100)
preds <- predict(AAcM1_segmented, newdata = data.frame(percTime=newx), interval="confidence")
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = alpha("black", 0.25), border = NA)
lines(x=newx, y=preds[,2], col=1, lty=2, lwd=0.5)
lines(x=newx, y=preds[,3], col=1, lty=2, lwd=0.5)
lines(x=newx, y=preds[,1], col=1, lty=1, lwd=1)

##AAr
#AArPlot
#AArM1 <- with(AAr, lm(ln.ux ~ percTime, weights=no.surv))
#plot(AArM1)	#shows that row 7 has undue weight
AArM1 <- with(AAr[c(1:6, 8:nrow(AAr)),], lm(ln.ux ~ percTime, weights=no.surv))
summary(AArM1)
for(i in 1:3){
	print(
		AIC(
				segmented(AArM1, seg.Z = ~percTime,npsi=i)
				)
			)
}

AArM1_segmented <- segmented(AArM1, seg.Z = ~percTime, npsi=1)
AArM1_segmented
summary(AArM1_segmented)
AArM1_segmented_CI <- confint(AArM1_segmented)
AArM1_segmented_CI
plot(AArM1_segmented, col="salmon", add=T)
lines(AArM1_segmented, shift=T, col="salmon")
#abline(v=AArM1_segmented_CI[,2:3], lty=3, col="salmon")
#abline(v=AArM1_segmented_CI[,1], lty=2, col="salmon")
points(x=AAr$percTime, y=AAr$ln.ux, pch=16, col="salmon")
newx <- seq(min(AAr$percTime), max(AAr$percTime), length.out=100)
preds <- predict(AArM1_segmented, newdata = data.frame(percTime=newx), interval="confidence")
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = alpha("salmon", 0.25), border = NA)
lines(x=newx, y=preds[,2], col="salmon", lty=2, lwd=0.5)
lines(x=newx, y=preds[,3], col="salmon", lty=2, lwd=0.5)
lines(x=newx, y=preds[,1], col="salmon", lty=1, lwd=1)

text(y=-7, x=70, labels=paste("Rsq =", signif(summary(AAcM1_segmented)$r.squared, 2)), col=1)
text(y=-7.5, x=70, labels=paste("Rsq =", signif(summary(AArM1_segmented)$r.squared, 2)), col="salmon")
text(x=AAr$percTime[7]+15, y=AAr$ln.ux[7], labels="Excluded", col="salmon")


##ABc
#ABcPlot
ABcM1 <- with(ABc, lm(ln.ux ~ percTime, weights=no.surv))
summary(ABcM1)

for(i in 1:3){
	print(
		AIC(
				segmented(ABcM1, seg.Z = ~percTime,npsi=i)
				)
			)
}

ABcM1_segmented <- segmented(ABcM1, seg.Z = ~percTime,npsi=2)	#very marginal Rsq increase with npsi=3. 
ABcM1_segmented
summary(ABcM1_segmented)
ABcM1_segmented_CI <- confint(ABcM1_segmented)
ABcM1_segmented_CI
plot(ABcM1_segmented, col=1, ylim=c(-8.5,0), main="AB", ylab="logN mortality rate", xlab="Relative age (% max. lifespan)")
lines(ABcM1_segmented, shift=T)
#abline(v=ABcM1_segmented_CI[,2:3], lty=3)
#abline(v=ABcM1_segmented_CI[,1], lty=2)
points(x=ABc$percTime, y=ABc$ln.ux, pch=16)
newx <- seq(min(ABc$percTime), max(ABc$percTime), length.out=100)
preds <- predict(ABcM1_segmented, newdata = data.frame(percTime=newx), interval="confidence")
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = alpha("black", 0.25), border = NA)
lines(x=newx, y=preds[,2], col=1, lty=2, lwd=0.5)
lines(x=newx, y=preds[,3], col=1, lty=2, lwd=0.5)
lines(x=newx, y=preds[,1], col=1, lty=1, lwd=1)

##ABr
#ABrPlot
ABrM1 <- with(ABr, lm(ln.ux ~ percTime, weights=no.surv))
summary(ABrM1)

for(i in 1:3){
	print(
		AIC(
				segmented(ABrM1, seg.Z = ~percTime,npsi=i)
				)
			)
}

ABrM1_segmented <- segmented(ABrM1, seg.Z = ~percTime, npsi=1)	#increase at nspi=2 so marginal it's silly.
ABrM1_segmented
summary(ABrM1_segmented)
ABrM1_segmented_CI <- confint(ABrM1_segmented)
ABrM1_segmented_CI
plot(ABrM1_segmented, col="salmon", add=T)
lines(ABrM1_segmented, shift=T, col="salmon")
#abline(v=ABrM1_segmented_CI[,2:3], lty=3, col="salmon")
#abline(v=ABrM1_segmented_CI[,1], lty=2, col="salmon")
points(x=ABr$percTime, y=ABr$ln.ux, pch=16, col="salmon")
newx <- seq(min(ABr$percTime), max(ABr$percTime), length.out=100)
preds <- predict(ABrM1_segmented, newdata = data.frame(percTime=newx), interval="confidence")
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = alpha("salmon", 0.25), border = NA)
lines(x=newx, y=preds[,2], col="salmon", lty=2, lwd=0.5)
lines(x=newx, y=preds[,3], col="salmon", lty=2, lwd=0.5)
lines(x=newx, y=preds[,1], col="salmon", lty=1, lwd=1)

text(y=-7, x=70, labels=paste("Rsq =", signif(summary(ABcM1_segmented)$r.squared, 2)), col=1)
text(y=-7.5, x=70, labels=paste("Rsq =", signif(summary(ABrM1_segmented)$r.squared, 2)), col="salmon")


##BAc
#BAcPlot
BAcM1 <- with(BAc, lm(ln.ux ~ percTime, weights=no.surv))
summary(BAcM1)

for(i in 1:3){
	print(
		AIC(
				segmented(BAcM1, seg.Z = ~percTime,npsi=i)
				)
			)
}

BAcM1_segmented <- segmented(BAcM1, seg.Z = ~percTime,npsi=2)
BAcM1_segmented
summary(BAcM1_segmented)
BAcM1_segmented_CI <- confint(BAcM1_segmented)
BAcM1_segmented_CI
plot(BAcM1_segmented, col=1, ylim=c(-8.5,0), main="BA", ylab="logN mortality rate", xlab="Relative age (% max. lifespan)")
lines(BAcM1_segmented, shift=T)
#abline(v=BAcM1_segmented_CI[,2:3], lty=3)
#abline(v=BAcM1_segmented_CI[,1], lty=2)
points(x=BAc$percTime, y=BAc$ln.ux, pch=16)
newx <- seq(min(BAc$percTime), max(BAc$percTime), length.out=100)
preds <- predict(BAcM1_segmented, newdata = data.frame(percTime=newx), interval="confidence")
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = alpha("black", 0.25), border = NA)
lines(x=newx, y=preds[,2], col=1, lty=2, lwd=0.5)
lines(x=newx, y=preds[,3], col=1, lty=2, lwd=0.5)
lines(x=newx, y=preds[,1], col=1, lty=1, lwd=1)

##BAr
#BArPlot
BArM1 <- with(BAr, lm(ln.ux ~ percTime, weights=no.surv))
summary(BArM1)

for(i in 1:3){
	print(
		AIC(
				segmented(BArM1, seg.Z = ~percTime,npsi=i)
				)
			)
}

BArM1_segmented <- segmented(BArM1, seg.Z = ~percTime, npsi=2)
BArM1_segmented
summary(BArM1_segmented)
BArM1_segmented_CI <- confint(BArM1_segmented)
BArM1_segmented_CI
plot(BArM1_segmented, col="salmon", add=T)
lines(BArM1_segmented, shift=T, col="salmon")
#abline(v=BArM1_segmented_CI[,2:3], lty=3, col="salmon")
#abline(v=BArM1_segmented_CI[,1], lty=2, col="salmon")
points(x=BAr$percTime, y=BAr$ln.ux, pch=16, col="salmon")
newx <- seq(min(BAr$percTime), max(BAr$percTime), length.out=100)
preds <- predict(BArM1_segmented, newdata = data.frame(percTime=newx), interval="confidence")
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = alpha("salmon", 0.25), border = NA)
lines(x=newx, y=preds[,2], col="salmon", lty=2, lwd=0.5)
lines(x=newx, y=preds[,3], col="salmon", lty=2, lwd=0.5)
lines(x=newx, y=preds[,1], col="salmon", lty=1, lwd=1)

text(y=-7, x=60, labels=paste("Rsq =", signif(summary(BAcM1_segmented)$r.squared, 2)), col=1)
text(y=-7.5, x=60, labels=paste("Rsq =", signif(summary(BArM1_segmented)$r.squared, 2)), col="salmon")


##BBc
#BBcPlot
BBcM1 <- with(BBc, lm(ln.ux ~ percTime, weights=no.surv))
summary(BBcM1)

for(i in 1:3){
	print(
		AIC(
				segmented(BBcM1, seg.Z = ~percTime,npsi=i)
				)
			)
}

BBcM1_segmented <- segmented(BBcM1, seg.Z = ~percTime,npsi=2)
BBcM1_segmented
summary(BBcM1_segmented)
BBcM1_segmented_CI <- confint(BBcM1_segmented)
BBcM1_segmented_CI
plot(BBcM1_segmented, col=1, ylim=c(-8.5,0), main="BB", ylab="logN mortality rate", xlab="Relative age (% max. lifespan)")
lines(BBcM1_segmented, shift=T)
#abline(v=BBcM1_segmented_CI[,2:3], lty=3)
#abline(v=BBcM1_segmented_CI[,1], lty=2)
points(x=BBc$percTime, y=BBc$ln.ux, pch=16)
newx <- seq(min(BBc$percTime), max(BBc$percTime), length.out=100)
preds <- predict(BBcM1_segmented, newdata = data.frame(percTime=newx), interval="confidence")
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = alpha("black", 0.25), border = NA)
lines(x=newx, y=preds[,2], col=1, lty=2, lwd=0.5)
lines(x=newx, y=preds[,3], col=1, lty=2, lwd=0.5)
lines(x=newx, y=preds[,1], col=1, lty=1, lwd=1)

##BBr
#BBrPlot
BBrM1 <- with(BBr, lm(ln.ux ~ percTime, weights=no.surv))
summary(BBrM1)

for(i in 1:3){
	print(
		AIC(
				segmented(BBrM1, seg.Z = ~percTime,npsi=i)
				)
			)
}

BBrM1_segmented <- segmented(BBrM1, seg.Z = ~percTime, npsi=2)
BBrM1_segmented
summary(BBrM1_segmented)
BBrM1_segmented_CI <- confint(BBrM1_segmented)
BBrM1_segmented_CI
plot(BBrM1_segmented, col="salmon", add=T)
lines(BBrM1_segmented, shift=T, col="salmon")
#abline(v=BBrM1_segmented_CI[,2:3], lty=3, col="salmon")
#abline(v=BBrM1_segmented_CI[,1], lty=2, col="salmon")
points(x=BBr$percTime, y=BBr$ln.ux, pch=16, col="salmon")
newx <- seq(min(BBr$percTime), max(BBr$percTime), length.out=100)
preds <- predict(BBrM1_segmented, newdata = data.frame(percTime=newx), interval="confidence")
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = alpha("salmon", 0.25), border = NA)
lines(x=newx, y=preds[,2], col="salmon", lty=2, lwd=0.5)
lines(x=newx, y=preds[,3], col="salmon", lty=2, lwd=0.5)
lines(x=newx, y=preds[,1], col="salmon", lty=1, lwd=1)

text(y=-7, x=70, labels=paste("Rsq =", signif(summary(BBcM1_segmented)$r.squared, 2)), col=1)
text(y=-7.5, x=70, labels=paste("Rsq =", signif(summary(BBrM1_segmented)$r.squared, 2)), col="salmon")
dev.off()

##############################
##	plot inflection points	##
##############################

AAcM1_segmented_CI <- data.frame(AAcM1_segmented_CI)
AArM1_segmented_CI <- data.frame(AArM1_segmented_CI)

ABcM1_segmented_CI <- data.frame(ABcM1_segmented_CI)
ABrM1_segmented_CI <- data.frame(ABrM1_segmented_CI)

BAcM1_segmented_CI <- data.frame(BAcM1_segmented_CI)
BArM1_segmented_CI <- data.frame(BArM1_segmented_CI)

BBcM1_segmented_CI <- data.frame(BBcM1_segmented_CI)
BBrM1_segmented_CI <- data.frame(BBrM1_segmented_CI)


AAcM1_segmented_CI$rapa <- "-"
AArM1_segmented_CI$rapa <- "+"

ABcM1_segmented_CI$rapa <- "-"
ABrM1_segmented_CI$rapa <- "+"

BAcM1_segmented_CI$rapa <- "-"
BArM1_segmented_CI$rapa <- "+"

BBcM1_segmented_CI$rapa <- "-"
BBrM1_segmented_CI$rapa <- "+"


AAcM1_segmented_CI$mito <- "A"
AArM1_segmented_CI$mito <- "A"

ABcM1_segmented_CI$mito <- "A"
ABrM1_segmented_CI$mito <- "A"

BAcM1_segmented_CI$mito <- "B"
BArM1_segmented_CI$mito <- "B"

BBcM1_segmented_CI$mito <- "B"
BBrM1_segmented_CI$mito <- "B"


AAcM1_segmented_CI$nuclear <- "A"
AArM1_segmented_CI$nuclear <- "A"

ABcM1_segmented_CI$nuclear <- "B"
ABrM1_segmented_CI$nuclear <- "B"

BAcM1_segmented_CI$nuclear <- "A"
BArM1_segmented_CI$nuclear <- "A"

BBcM1_segmented_CI$nuclear <- "B"
BBrM1_segmented_CI$nuclear <- "B"


segmented_CIs <- rbind(
	AAcM1_segmented_CI,
	AArM1_segmented_CI,
	ABcM1_segmented_CI,
	ABrM1_segmented_CI,
	BAcM1_segmented_CI,
	BArM1_segmented_CI,
	BBcM1_segmented_CI,
	BBrM1_segmented_CI)
segmented_CIs

segmented_CIs$psi1 <- substr(rownames(segmented_CIs), 1, 13)
segmented_CIs
rownames(segmented_CIs) <- 1:nrow(segmented_CIs)
colnames(segmented_CIs)[1:3] <- c("time", "CIlow", "CIhigh")
segmented_CIs$genotype <- with(segmented_CIs, factor(paste(mito, nuclear, sep="")))

p1 <- ggplot(
	data = segmented_CIs, aes(y = interaction(rapa, genotype), x = time, col=genotype, group=genotype)) +
	scale_y_discrete(NULL, guide="axis_nested") +
        geom_errorbar(aes(xmin=CIlow, xmax=CIhigh, colour=genotype), width=0.6) +
        geom_point(aes(shape=rapa)) + 
        theme_bw()+
#        ylim(c(0, 0.9)) +
        theme(axis.text = element_text(size = 12)) +
        labs(x="%age lifespan elapsed at mortality rate inflection ± CI") +
        scale_colour_manual(labels= c("AA", "AB", "BA", "BB"),
                            values = c("darkolivegreen3", "steelblue2", "coral2","gold")) +
                            facet_grid(cols=vars(psi1)) +
                            xlim(min(segmented_CIs$CIlow),100)

p1

pdf("psi1.pdf", h=4, w=4)
print(p1)
dev.off()

##############################################
##	model mortality rate at specific points	##
##############################################

ind <- data.frame(
	rapa=rep(c("-", "+"), 4),
	mito=rep(c("A", "B"), each=4),
	nuclear=rep(rep(c("A", "B"), each=2),2))

mortRate50 <- cbind(ind, rbind(
predict(AAcM1_segmented, newdata = data.frame(percTime=50), interval="confidence"),
predict(AArM1_segmented, newdata = data.frame(percTime=50), interval="confidence"),

predict(ABcM1_segmented, newdata = data.frame(percTime=50), interval="confidence"),
predict(ABrM1_segmented, newdata = data.frame(percTime=50), interval="confidence"),

predict(BAcM1_segmented, newdata = data.frame(percTime=50), interval="confidence"),
predict(BArM1_segmented, newdata = data.frame(percTime=50), interval="confidence"),

predict(BBcM1_segmented, newdata = data.frame(percTime=50), interval="confidence"),
predict(BBrM1_segmented, newdata = data.frame(percTime=50), interval="confidence")))


mortRate30 <- cbind(ind, rbind(
predict(AAcM1_segmented, newdata = data.frame(percTime=30), interval="confidence"),
predict(AArM1_segmented, newdata = data.frame(percTime=30), interval="confidence"),

predict(ABcM1_segmented, newdata = data.frame(percTime=30), interval="confidence"),
predict(ABrM1_segmented, newdata = data.frame(percTime=30), interval="confidence"),

predict(BAcM1_segmented, newdata = data.frame(percTime=30), interval="confidence"),
predict(BArM1_segmented, newdata = data.frame(percTime=30), interval="confidence"),

predict(BBcM1_segmented, newdata = data.frame(percTime=30), interval="confidence"),
predict(BBrM1_segmented, newdata = data.frame(percTime=30), interval="confidence")))


mortRate70 <- cbind(ind, rbind(
predict(AAcM1_segmented, newdata = data.frame(percTime=70), interval="confidence"),
predict(AArM1_segmented, newdata = data.frame(percTime=70), interval="confidence"),

predict(ABcM1_segmented, newdata = data.frame(percTime=70), interval="confidence"),
predict(ABrM1_segmented, newdata = data.frame(percTime=70), interval="confidence"),

predict(BAcM1_segmented, newdata = data.frame(percTime=70), interval="confidence"),
predict(BArM1_segmented, newdata = data.frame(percTime=70), interval="confidence"),

predict(BBcM1_segmented, newdata = data.frame(percTime=70), interval="confidence"),
predict(BBrM1_segmented, newdata = data.frame(percTime=70), interval="confidence")))

mortRate50$time <- 50
mortRate30$time <- 30
mortRate70$time <- 70

mortRate <- rbind(mortRate50, mortRate30, mortRate70)
mortRate$genotype <- with(mortRate, factor(paste(mito, nuclear, sep="")))
#colnames(mortRate)[c(4)] <- "Mortality rate"

mortRateP1 <- ggplot(data=mortRate,  aes(x=interaction(rapa, mito, nuclear), y = fit, col=genotype, group=genotype)) +
	scale_x_discrete(NULL, guide="axis_nested") +
        geom_errorbar(aes(ymin=lwr, ymax=upr, colour=genotype), width=0.6) +
        geom_point(aes(shape=rapa), size=5) + 
        theme_bw()+
#        ylim(c(0, 0.9)) +
        theme(axis.text = element_text(size = 12)) +
        labs(y="Mortality rate ± CI95") +
        scale_colour_manual(labels= c("AA", "AB", "BA", "BB"),
                            values = c("darkolivegreen3", "steelblue2", "coral2","gold")) +
                            facet_grid(cols=vars(time))
                           
pdf("mortalityRates.pdf", h=4, w=8)
print(mortRateP1)
dev.off()