#########################################################################################################
######## Research question: Is there a difference in response to rapamycin across mito-nuclear backgrounds?
######## R Studio version: 2023.06.0+421
######## Date: 11/2023 
######## Author: RI and AJD

rm(list=ls())

#Load packages
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

#options(na.action="na.fail", contrasts=c("contr.sum", "contr.poly"))
options(es.use_symbols = TRUE)

##########################################################################################################
##	MAIN ANALYSIS	########################################################################################
##########################################################################################################

## RUN CODE TO PROCESS DEATHS AND CENSOR DATA---------------------------------------------
    jmp2r <- function(x){
      for(i in 1:nrow(x)){
        #		print(paste("Processing row", i))
        xx <- x[i,]
        if(xx$event[1] > 0){
          xxx <- matrix(rep(as.matrix(xx), xx[1,"event"]), ncol=ncol(xx), byrow=T)
          if(exists("a")){a <- rbind(a, xxx)}else{(a <- xxx)}
        }
      }
      a <- as.data.frame(a)
      for(i in 1:ncol(a)){
        colnames(a)[i] <- colnames(x)[i]
        if(class(a[,i]) != class(x[,i])) {a[,i] <- as.character(a[,i]); class(a[,i]) <- class(x[,i])}
      }	
      #    a$censor <- ifelse(a$censor==0, 1, 0)  #run this line if you are working with truly jmp format data. but not in this case, where we have control so we can define things much more sensibly than jmp does...
      a <- a[,!colnames(a) %in% c("count", "deathIsZero_censorIsOne", "event")]
      a
    }

## INPUT DATA INTO EXCEL-------------------------------------------------------------------
    d_A <- read.csv("RI12_survData_JMP.csv")

## PREPARE DATA FOR ANALYSIS --------------------------------------------------------------

  d_A <- jmp2r(d_A)

  #Subset of data without ancestral haplotypes: australia and benin
    d_A <- droplevels(subset(d_A, Bxed=="yes"))

  #Create mitoNuclear variable
    d_A$mitoNuclear <- factor(paste(d_A$mito, d_A$nuclear))

  #Change experimental variables to factors 
    d_A$mito <- as.factor(d_A$mito)
    d_A$nuclear <- as.factor(d_A$nuclear)
    d_A$line <- as.factor(d_A$line)
    d_A$rapa <- as.factor(d_A$rapa)
      
    

## VISUALISE LIFESPAN CURVES (KAPLAIN-MEIER) FOR EACH MITONUCLEOGENOTYPE ---------------------

  d_A$geno <- gsub('[[:digit:]]+', '', d_A$line)
  d_A$cond <- with(d_A, paste(line, rapa))
  d_A$geno <- as.factor(d_A$geno)
  d_A$cond <- as.factor(d_A$cond)

  # BBs--------

    BBs_A <- droplevels(subset(d_A, geno =="BB"))
    BBs_S_A <- Surv(BBs_A$time, BBs_A$censor)
    BBs_survfit_A <- survfit(BBs_S_A ~ BBs_A$rapa)
    bbs <- plot(BBs_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))

      # Create BB as a variable that I can later merge with BB1, BB2 and BB3 to create a plot that
      # has BB1, BB2 and BB3 plotted under control and rapamycin treatments
      # PLUS BB (which is the average lifespan of BB1, BB2 and BB3)
        BBs_A$line <- BBs_A$geno

      # Now plot BB1, BB3 and BB3 plus BB (their average)
    allBBs_A <-  droplevels(subset(d_A, line %in% c("BB1","BB2", "BB3")))
    allBBs_A2 <- rbind(allBBs_A, BBs_A)
    allBBs_A2$cond <- with(allBBs_A2, paste(line, rapa))
    allBBs_A2$cond <- as.factor(allBBs_A2$cond)
    allBBs_S_A <- Surv(allBBs_A2$time, allBBs_A2$censor)
    allBBs_survfit_A <- survfit(allBBs_S_A ~ allBBs_A2$cond)
    plot(allBBs_survfit_A, lty=c(1,1,1,1,1,1,1,1), col=c("grey","salmon","grey","salmon","grey","salmon","grey","salmon"), lwd=c(1,1,0.2, 0.2,0.2,0.2,0.2,0.2),
         cex.axis=1.5, ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))

   # ABs--------
    ABs_A <- droplevels(subset(d_A, geno =="AB"))
    ABs_S_A <- Surv(ABs_A$time, ABs_A$censor)
    ABs_survfit_A <- survfit(ABs_S_A ~ ABs_A$rapa)
    abs <- plot(ABs_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))
      
    # Create AB as a variable that I can later merge with AB1, AB2 and AB3 to create a plot that
    # has AB1, AB2 and AB3 plotted under control and rapamycin treatments
    # PLUS AB (which is the average lifespan of AB1, AB2 and AB3)
      ABs_A$line <- ABs_A$geno
    
    allABs_A <-  droplevels(subset(d_A, line %in% c("AB1","AB2", "AB3")))
    allABs_A2 <- rbind(allABs_A, ABs_A)
    allABs_A2$cond <- with(allABs_A2, paste(line, rapa))
    allABs_A2$cond <- as.factor(allABs_A2$cond)
    allABs_S_A <- Surv(allABs_A2$time, allABs_A2$censor)
    allABs_survfit_A <- survfit(allABs_S_A ~ allABs_A2$cond)
    plot(allABs_survfit_A, lty=c(1,1,1,1,1,1,1,1), col=c("grey","salmon","grey","salmon","grey","salmon","grey","salmon"), lwd=c(1,1,0.2, 0.2,0.2,0.2,0.2,0.2),
         cex.axis=1.5, ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90)) 


    # BAs--------
    BAs_A <- droplevels(subset(d_A, geno =="BA"))
    BAs_S_A <- Surv(BAs_A$time, BAs_A$censor)
    BAs_survfit_A <- survfit(BAs_S_A ~ BAs_A$rapa)
    bas <- plot(BAs_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))
      
      # Create BA as a variable that I can later merge with BA1, BA2 and BA3 to create a plot that
      # has BA1, BA2 and BA3 plotted under control and rapamycin treatments
      # PLUS BA (which is the average lifespan of BA1, BA2 and BA3)
        BAs_A$line <- BAs_A$geno
      
    allBAs_A <-  droplevels(subset(d_A, line %in% c("BA1","BA2", "BA3")))
    allBAs_A2 <- rbind(allBAs_A, BAs_A)
    allBAs_A2$cond <- with(allBAs_A2, paste(line, rapa))
    allBAs_A2$cond <- as.factor(allBAs_A2$cond)
    allBAs_S_A <- Surv(allBAs_A2$time, allBAs_A2$censor)
    allBAs_survfit_A <- survfit(allBAs_S_A ~ allBAs_A2$cond)
    plot(allBAs_survfit_A, lty=c(1,1,1,1,1,1,1,1), col=c("grey","salmon","grey","salmon","grey","salmon","grey","salmon"), lwd=c(1,1,0.2, 0.2,0.2,0.2,0.2,0.2),
           cex.axis=1.5, ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90)) 

    # AAs--------
    AAs_A <- droplevels(subset(d_A, geno =="AA"))
    AAs_S_A <- Surv(AAs_A$time, AAs_A$censor)
    AAs_survfit_A <- survfit(AAs_S_A ~ AAs_A$rapa)
    aas <- plot(AAs_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))
        
      # Create AA as a variable that I can later merge with AA1, AA2 and AA3 to create a plot that
      # has AA1, AA2 and AA3 plotted under control and rapamycin treatments
      # PLUS AA (which is the average lifespan of AA1, AA2 and AA3)
        AAs_A$line <- AAs_A$geno
      
    allAAs_A <-  droplevels(subset(d_A, line %in% c("AA1","AA2", "AA3")))
    allAAs_A2 <- rbind(allAAs_A, AAs_A)
    allAAs_A2$cond <- with(allAAs_A2, paste(line, rapa))
    allAAs_A2$cond <- as.factor(allAAs_A2$cond)
    allAAs_S_A <- Surv(allAAs_A2$time, allAAs_A2$censor)
    allAAs_survfit_A <- survfit(allAAs_S_A ~ allAAs_A2$cond)
    plot(allAAs_survfit_A, lty=c(1,1,1,1,1,1,1,1), col=c("grey","salmon","grey","salmon","grey","salmon","grey","salmon"), lwd=c(1,1,0.2, 0.2,0.2,0.2,0.2,0.2),
           cex.axis=1.5, ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90)) 
  

##################################
##	parametric survival models	##
##################################

d_500 <- droplevels(subset(d_A, rapa=="no"))
dd_500 <- datadist(d_500)
options(datadist = 'dd_500')
psm500 <- psm(Surv(time, censor) ~ mito * nuclear, 
	data=d_500,
	dist="weibull")
psm5000 <- psm(Surv(time, censor) ~ line, 
	data=d_500,
	dist="weibull")
AIC(psm500)
AIC(psm5000)

cph500 <- cph(Surv(time, censor) ~ mito * nuclear, 
	data=d_500)
cph5000 <- cph(Surv(time, censor) ~ line, 
	data=d_500)
AIC(cph500)
AIC(cph5000)

dd_A <- datadist(d_A)
options(datadist = 'dd_A')
#psm0 <- psm(Surv(time, censor) ~ rapa * line, 
#	data=d_A,
#	dist="weibull")
#pairs(emmeans(psm0, ~ rapa * line), by="line")

psm1 <- psm(Surv(time, censor) ~ rapa * mito * nuclear, 
	data=d_A,
	dist="weibull")
psm1
psm1_exponential <- update(psm1, dist="exponential")
psm1_gaussian <- update(psm1, dist="gaussian")
psm1_logistic <- update(psm1, dist="logistic")
psm1_logNormal <- update(psm1, dist="lognormal")
psm1_loglogistic <- update(psm1, dist="loglogistic")

psm1_exponential
psm1_gaussian
psm1_logistic
psm1_logNormal
psm1_loglogistic

AIC(psm1)
AIC(psm1_exponential)
AIC(psm1_gaussian)
AIC(psm1_logistic)
AIC(psm1_logNormal)
AIC(psm1_loglogistic)

	#we accept the logistic model as the best fitting.
psm1 <- psm1_logistic
psm1
Anova(psm1, type="3")
psm1_jointTests <- joint_tests(psm1)
psm1_jointTests
fbw <- fastbw(psm1)
fbw
plot(effectsize(psm1))
psm1_ES <- with(psm1_jointTests, F_to_eta2(f=F.ratio, df=df1, df_error=df2, alternative="two"))
plot(psm1_ES)

	#logistic has highest r^2
emmip(psm1, ~ rapa | mito | nuclear, CIs=T)
psm1_emm <- emmeans(psm1, ~ rapa | mito | nuclear)
joint_tests(psm1_emm)
joint_tests(psm1_emm, by=c("mito", "nuclear"))
pairs(psm1_emm)
pairs(psm1_emm, by=c("mito", "nuclear"))
pairsDf <- data.frame(pairs(psm1_emm, by=c("mito", "nuclear")))
pairsDf$genotype <- with(pairsDf, factor(paste(mito, nuclear, sep="")))


	#plot emmeans
psm1Df <- data.frame(psm1_emm)
psm1Df$genotype <- with(psm1Df, factor(paste(mito, nuclear, sep="")))
levels(psm1Df$rapa) <- c("-", "+")

p1 <- ggplot(
	data = psm1Df, aes(x = interaction(rapa, genotype), y = emmean, col=genotype, group=genotype)) +
	scale_x_discrete(NULL, guide="axis_nested") +
        geom_errorbar(aes(ymin=upper.CL, ymax=lower.CL, colour=genotype), width=0.6) +
        geom_point(aes(shape=rapa), size=5) + 
        theme_bw()+
#        ylim(c(0, 0.9)) +
        theme(axis.text = element_text(size = 12)) +
        labs(y="Emmean ± CI") +
        scale_colour_manual(labels= c("AA", "AB", "BA", "BB"),
                            values = c("darkolivegreen3", "steelblue2", "coral2","gold")) 

p2 <- ggplot(
	data = pairsDf, aes(x = genotype, y = estimate, col=genotype, group=genotype)) +
        geom_errorbar(aes(ymin=estimate-SE, ymax=estimate+SE, colour=genotype), width=.6)+
        geom_point(size=5) + 
        theme_bw()+
#        ylim(c(0, 0.9)) +
        theme(axis.text = element_text(size = 12)) +
        labs(y="Rapamycin effect (coefficient ± SE)", x="genotype") +
        scale_colour_manual(labels= c("AA", "AB", "BA", "BB"),
                            values = c("darkolivegreen3", "steelblue2", "coral2","gold")) +
                            ylim(-4, -9)                       

pdf(h=4, w=4, "emmeans.pdf")
print(p1)
dev.off()

pdf(h=4, w=4, "emmeans_effects.pdf")
print(p2)
dev.off()

psm2 <- update(psm1, . ~ rapa * line)
psm2
joint_tests(psm2, by="line")