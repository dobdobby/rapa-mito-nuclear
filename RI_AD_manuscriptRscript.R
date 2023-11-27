#########################################################################################################
######## Research question: Is there a difference in starvation resistance across mito-nuclear backgrounds?
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

  d_A <- read.csv("RI12_RapaLifespan_datasetforjmp2r_130622.csv")

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
      

## MODEL THE DATA ---------------------------------------------------------------------------

  #Create  "survival" object - this is the Y variable for the model 
    S_A <- with(d_A, Surv(time, censor))
    #the number gives time, "+" means "this individual was censored at this time

## ANALYSIS OF MODEL ------------------------------------------------------------------------
    
  #ANOVA type III
    m1_A <- coxme(S_A ~ rapa * mitoNuclear + (1|line), d_A)
    Anova(m1_A, type="III")

  #Does the effect of Rapa differ in the AAs, ABs, BAs, BBs?
  #emmeans: post hoc comparisons among groups after fitting a model.
    m1emm_A <- emmeans(m1_A, ~ rapa | mitoNuclear)
    pairs(m1emm_A, simple="each")
    
    joint_tests(m1_A, by="mitoNuclear")
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
  
# ---------------------------------------------------------------------------------------------
# Plot data to see how different lifespans are under: 
# (A) Control treatment
# (B) Rapamycin treatment
# (C) In response to rapamycin

    ########################-------------------------------------------------------------------
    ##	Controls lifespan	##-------------------------------------------------------------------
    ########################-------------------------------------------------------------------

    #make m1emm into a dataframe and subset by controls
      m1emm_df_A <- data.frame(m1emm_A)
      m1emm_controls_A <- subset(m1emm_df_A, rapa=="no")

    #create an emmeans controls variable
      emmeans_controls_A <- m1emm_controls_A$emmean

    #create an SEs controls variable
      SEs_controls_A <- m1emm_controls_A$SE

    genotype <- c("AA", "AB", "BA", "BB") 
    controls_df_A <- data.frame(emmeans_controls_A, SEs_controls_A, genotype)
      
    # plot the data under control treatment
      ggplot(data = controls_df_A , aes(x = genotype, y = emmeans_controls_A, col=genotype, group=genotype)) +
        geom_errorbar(aes(ymin=emmeans_controls_A-SEs_controls_A, ymax=emmeans_controls_A+SEs_controls_A), colour="black", width=.1)+
        geom_point(size=10) + 
        theme_bw()+
        ylim(-1,1.5)+
        theme(axis.text = element_text(size = 10)) +
        geom_hline(yintercept=0) +
        labs(y="Controls Lifespan(emmeans)", x="genotype", title="Lifespan across genotypes controls") +
        scale_colour_manual(labels= c("AA", "AB", "BA", "BB"),
                              values = c("darkolivegreen3", "steelblue2", "coral2","gold"))

    ##########################-------------------------------------------------------------------
    ##	Rapamycin lifespan	##-------------------------------------------------------------------
    ##########################-------------------------------------------------------------------

    #subset m1emm by rapamycin condition
      m1emm_rapa_A <- subset(m1emm_df_A, rapa=="yes")

    #create an emmeans rapa variable
      emmeans_rapa_A <- m1emm_rapa_A$emmean

    #create an SEs rapa variable
      SEs_rapa_A <- m1emm_rapa_A$SE

    genotype <- c("AA", "AB", "BA", "BB")
    rapa_df_A <- data.frame(emmeans_rapa_A, SEs_rapa_A, genotype)

    #plot the data under control treatment
      ggplot(data = rapa_df_A , aes(x = genotype, y = emmeans_rapa_A, col=genotype, group=genotype)) +
        geom_errorbar(aes(ymin=emmeans_rapa_A-SEs_rapa_A, ymax=emmeans_rapa_A+SEs_rapa_A), colour="black", width=.1)+
        geom_point(size=10) + 
        theme_bw()+
        ylim(-1.4,1.5)+
        theme(axis.text = element_text(size = 10)) +
        geom_hline(yintercept=0) +
        labs(y="Rapamycin Lifespan(emmeans)", x="genotype", title="Lifespan across genotype rapamycin") +
        scale_colour_manual(labels= c("AA", "AB", "BA", "BB"),
                            values = c("darkolivegreen3", "steelblue2", "coral2","gold")) 
      

    ######################################------------------------------------------------------
    ##	Effect of Rapamycin on lifespan	##------------------------------------------------------
    ######################################------------------------------------------------------
      
      m1pairs_A <- pairs(m1emm_A, by="mitoNuclear")
      m1pairs_A
      m1pairsDF_A <- data.frame(m1pairs_A)

    #Plot the data showing the magnitude of the response to rapamycin treatment
      ggplot(data = m1pairsDF_A, aes(x = genotype, y = estimate, col=genotype, group=genotype)) +
        geom_errorbar(aes(ymin=estimate-SE, ymax=estimate+SE), colour="black", width=.1)+
        geom_point(size=10) + 
        theme_bw()+
        ylim(c(0, 0.9)) +
        theme(axis.text = element_text(size = 12)) +
        labs(y="Rapamycin effect (Estimate value)", x="genotype", title="Rapamycin effect on lifespan across genotypes") +
        scale_colour_manual(labels= c("AA", "AB", "BA", "BB"),
                            values = c("darkolivegreen3", "steelblue2", "coral2","gold")) 
      

## MORTALITY RATE-------------------------------------------------------------------------------
      
  d <- read.csv("RI12_Mortalitydata.csv")

  strip1 <- strip_themed(background_x = elem_list_rect(fill = c("darkolivegreen4", "steelblue3", "gold", "coral1")))
  p0 <- ggplot(d, aes(x=time, y=mort_values))
  p0  + 
  geom_line(aes(group=condition), alpha=0.2, color="black")+ geom_smooth(aes(group=rapa, color=rapa)) +
  facet_grid(. ~ MN, scales="free_x") +
  theme_bw() +
  theme(axis.text = element_text(size = 12)) +
  scale_colour_manual(labels= c("Control", "Rapamycin"),
                      values = c("black", "darksalmon")) +
  facet_wrap2(~MN, strip=strip1)  +
  labs(y="Mortality rate", x="Time (Days)") 

  
######################################################################################################
##	SUPPLEMENTARY FIGURES: POPULATIONS	########################################################################
######################################################################################################

  ######################------------------------------------------------------
  ##	AA1, AA2, AA3, ##------------------------------------------------------
  ######################------------------------------------------------------

    #AA1--------

      AA1_A <- droplevels(subset(d_A, cond %in% c("AA1 no", "AA1 yes")))
      AA1_S_A <- Surv(AA1_A$time, AA1_A$censor)
      AA1_survfit_A <- survfit(AA1_S_A ~ AA1_A$cond)
      plot(AA1_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))

    #AA2--------
      
      AA2_A <- droplevels(subset(d_A, cond %in% c("AA2 no", "AA2 yes")))
      AA2_S_A <- Surv(AA2_A$time, AA2_A$censor)
      AA2_survfit_A <- survfit(AA2_S_A ~ AA2_A$cond)
      plot(AA2_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))

    #AA3-------- 

      AA3_A <- droplevels(subset(d_A, cond %in% c("AA3 no", "AA3 yes")))
      AA3_S_A <- Surv(AA3_A$time, AA3_A$censor)
      AA3_survfit_A <- survfit(AA3_S_A ~ AA3_A$cond)
      plot(AA3_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))

  ######################------------------------------------------------------
  ##	BB1, BB2, BB3	##------------------------------------------------------
  ######################------------------------------------------------------

    #BB1-------- 

      BB1_A <- droplevels(subset(d_A, cond %in% c("BB1 no", "BB1 yes")))
      BB1_S_A <- Surv(BB1_A$time, BB1_A$censor)
      BB1_survfit_A <- survfit(BB1_S_A ~ BB1_A$cond)
      plot(BB1_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))

    #BB2-------- 

      BB2_A <- droplevels(subset(d_A, cond %in% c("BB2 no", "BB2 yes")))
      BB2_S_A <- Surv(BB2_A$time, BB2_A$censor)
      BB2_survfit_A <- survfit(BB2_S_A ~ BB2_A$cond)
      plot(BB2_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))
      legend("bottomleft", legend=levels(BB2_A$cond), lty=c(1,1), col=c(1,2), lwd=3)

    #BB3--------

      BB3_A <- droplevels(subset(d_A, cond %in% c("BB3 no", "BB3 yes")))
      BB3_S_A <- Surv(BB3_A$time, BB3_A$censor)
      BB3_survfit_A <- survfit(BB3_S_A ~ BB3_A$cond)
      plot(BB3_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))
      legend("bottomleft", legend=levels(BB3_A$cond), lty=c(1,1), col=c(1,2), lwd=3)

  ######################------------------------------------------------------
  ##	AB1, AB2, AB3	##------------------------------------------------------
  ######################------------------------------------------------------

    #AB1-------- 
      
      AB1_A <- droplevels(subset(d_A, cond %in% c("AB1 no", "AB1 yes")))
      AB1_S_A <- Surv(AB1_A$time, AB1_A$censor)
      AB1_survfit_A <- survfit(AB1_S_A ~ AB1_A$cond)
      plot(AB1_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))
            
    #AB2-------- 
      
      AB2_A <- droplevels(subset(d_A, cond %in% c("AB2 no", "AB2 yes")))
      AB2_S_A <- Surv(AB2_A$time, AB2_A$censor)
      AB2_survfit_A <- survfit(AB2_S_A ~ AB2_A$cond)
      plot(AB2_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))

    #AB3--------

      AB3_A <- droplevels(subset(d_A, cond %in% c("AB3 no", "AB3 yes")))
      AB3_S_A <- Surv(AB3_A$time, AB3_A$censor)
      AB3_survfit_A <- survfit(AB3_S_A ~ AB3_A$cond)
      plot(AB3_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))

  ######################------------------------------------------------------
  ##	BA1, BA2, BA3	##------------------------------------------------------
  ######################------------------------------------------------------

    #BA1-------- 

      BA1_A <- droplevels(subset(d_A, cond %in% c("BA1 no", "BA1 yes")))
      BA1_S_A <- Surv(BA1_A$time, BA1_A$censor)
      BA1_survfit_A <- survfit(BA1_S_A ~ BA1_A$cond)
      plot(BA1_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))

    #BA2--------

      BA2_A <- droplevels(subset(d_A, cond %in% c("BA2 no", "BA2 yes")))
      BA2_S_A <- Surv(BA2_A$time, BA2_A$censor)
      BA2_survfit_A <- survfit(BA2_S_A ~ BA2_A$cond)
      plot(BA2_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))
      legend("bottomleft", legend=levels(BA2_A$cond), lty=c(1,1), col=c(1,2), lwd=3)

    #BA3-------- 

      BA3_A <- droplevels(subset(d_A, cond %in% c("BA3 no", "BA3 yes")))
      BA3_S_A <- Surv(BA3_A$time, BA3_A$censor)
      BA3_survfit_A <- survfit(BA3_S_A ~ BA3_A$cond)
      plot(BA3_survfit_A, lty=c(1,1), col=c("grey", "darksalmon"), lwd=3, xlab= substitute(paste(bold("Time (Days)"))), ylab= substitute(paste(bold("Survival"))), cex.lab=1.6, mgp=c(2.3,1,0), xlim=c(0,90))
      legend("bottomleft", legend=levels(BA3_A$cond), lty=c(1,1), col=c(1,2), lwd=3)

#############################################------------------------------------------------------
#########MORTALITY RATE POPULATIONS##########------------------------------------------------------
#############################################------------------------------------------------------

  p1 <- ggplot(d, aes(x=time, y=mort_values)) +
    facet_grid(rapa ~ MN)
  p1  + geom_line(aes(group=condition, color=rapa), alpha=1) +
    theme_bw() +
    theme(axis.text = element_text(size = 12)) +
    scale_colour_manual(labels= c("Control", "Rapamycin"),
                        values = c("black", "darksalmon")) +
    facet_wrap2(~MN, strip=strip1)  +
    labs(y="Mortality rate", x="Time (Days)") 

