library(ggplot2)
library(vegan)
library(dplyr)
library(magrittr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
#install.packages("randomForest")
library(randomForest)
library(knitr)
library(qiime2R)
library(tidyr) #for separate function
library(naniar)# Ffor replace all function
library(ggpubr)

## Random Forest for 16S sequencing data
setwd("~/Desktop/eunice/PhD/PhDProject/FungalBRD/Qiime/QiimeFungalComplete/exported/NetworkBacFun/")
metadata <-read.csv("~/Desktop/eunice/PhD/PhDProject/FungalBRD/Qiime/QiimeFungalComplete/exported/BRD.metadadaFungal.csv", na.strings = c("","NA"), header=TRUE)
#OTU table (shared file)
#The OTU table as exported from qiime has a pound sign before the header row. You need to delete that pound sign in a text editor.
str(metadata)
metadata$BRD <- factor(metadata$BRD) 
metadata$PenCode <- factor(metadata$PenCode)
metadata$season <- factor(metadata$season)
levels(metadata$season) <- list("Summer"="Summer", "Fall"="Fall", "Winter"="Winter")
order_groups <- metadata$ID
row.names(metadata) = metadata[,1]
metadata %>% tally()
metadata %>% count(season)
metadata %>% count(BRD)

##test the fingi-fungi interaction
S16tax <- read.csv("NewTax16SCooccur.csv")
rownames(S16tax) <- S16tax$X
S16tax <- S16tax[-1]
ITStax <- read.csv("NewTaxFungiCooccur.csv")
rownames(ITStax) <- ITStax$X
ITStax <- ITStax[-1]

#subsett for the BRD pathogens
myco <- subset(S16tax, Species == "Mycoplasma bovis")
histo<- subset(S16tax, Species == "Histophilus")
mann <- subset(S16tax, Species == "Mannheimia haemolytica")
paste <- subset(S16tax, Species == "Pasteurella multocida")

##test the fingi-bacteria interaction
ITStax.G <- ITStax[,-c(7:8)] #only genus
S16tax.G <- S16tax[,-c(7:8)]
AllBRDassos <- read.csv("prob.tableBRD.csv", sep = ",", header = T) # a total of 913725 associations (fungi-fungi, fungi-bacteria)
pos <- subset(AllBRDassos, p_gt < 0.05) #correct 96004 
neg <- subset(AllBRDassos, p_lt < 0.05) #correct 20483 

##this is only the table selecting for the bacteria and fungal associations
BRD_bacfun <- read.csv("BRD_bacfun.csv", sep = ",", header = T) #total associations 205095
tail(BRD_bacfun$sp1_name)
#p_gt < 0.05 #positive associations
#p_lt < 0.05 #negative associations
posBF <- subset(BRD_bacfun, p_gt < 0.05) ##otal positive 5593 
posBF <- subset(posBF, obs_cooccur > 37) ### probability that the species will be in at least 60% of the total samples
posBF <- subset(posBF, prob_cooccur > 0.75) ### probability that the two species will be in the same site

#MERGE taxonomy
posBF <- merge(posBF, ITStax.G, by.x = "sp1_name", by.y = 0)
str(posBF)
colnames(posBF) <- c("sp1_name", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                     "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "sp2_name", "PhylumITS", "ClassITS","OrderITS", "FamilyITS","GenusITS", "SpeciesITS")
posBF <- merge(posBF,S16tax.G, by.x = "sp2_name", by.y = 0)
str(posBF)
colnames(posBF) <- c("sp2_name","sp1_name", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                     "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "PhylumITS", "ClassITS","OrderITS", "FamilyITS","GenusITS", "SpeciesITS",
                     "Phylum16S", "Class16S","Order16S", "Family16S","Genus16S", "Species16S")
write.csv(posBF, "posBFBRD.csv")

##test if there is fungi associated with BRD-pathogens
mycoF <- merge(BRD_bacfun, myco, by.x = "sp2_name", by.y = 0)
str(mycoF)
mycoF$sp1_name <- as.factor(mycoF$sp1_name) #227 ASVs
unique(mycoF$sp1_name)
mycoFsP <-subset(mycoF, p_gt < 0.05) #positive
mycoFsP$sig <- c("Positive")
mycoFsP <- mycoFsP[,-c(12:16)]
mycoFsP <- merge(mycoFsP,ITStax.G, by.x = "sp1_name", by.y = 0)

mycoFsN <-subset(mycoF, p_lt < 0.05) #negative
mycoFsN %>% count(sp1_name)
mycoFsN$sig <- c("Negative")
mycoFsN <- mycoFsN[,-c(12:16)]
mycoFsN <- merge(mycoFsN,ITStax.G, by.x = "sp1_name", by.y = 0)

mycoSig <- rbind(mycoFsP,mycoFsN)
write.csv(mycoSig, "mycoSigBRD.csv")

pasteF <- merge(BRD_bacfun, paste, by.x = "sp2_name", by.y = 0)
str(pasteF)
pasteF$sp1_name <- as.factor(pasteF$sp1_name) #227 ASVs
unique(pasteF$sp1_name)
pasteFsN <-subset(pasteF, p_lt < 0.05) #negative
pasteFsP <-subset(pasteF, p_gt < 0.05) #positive, no significant

histoF <- merge(BRD_bacfun, histo, by.x = "sp2_name", by.y = 0)
unique(histoF$sp1_name)
histoF$sp1_name <- as.factor(histoF$sp1_name) #227 ASVs
unique(histoF$sp1_name)
histoFsP <-subset(histoF, p_gt < 0.05) #positive
histoFsP$sig <- c("Positive")
histoFsP <-histoFsP[,-c(12:16)]
histoFsP <- merge(histoFsP,ITStax.G, by.x = "sp1_name", by.y = 0)

histoFsN <-subset(histoF, p_lt < 0.05) #negative
histoFsN$sig <- c("Negative")
histoFsN <- histoFsN[,-c(12:16)]
histoFsN <- merge(histoFsN,ITStax.G, by.x = "sp1_name", by.y = 0)
histoSig <- rbind(histoFsP,histoFsN)
write.csv(histoSig, "histoSigBRD.csv")

mannF <- merge(BRD_bacfun, mann, by.x = "sp2_name", by.y = 0)

##negative associations
negBF <-subset(BRD_bacfun, p_lt < 0.05) #negative associations 5414
random <- (205095 - 5593 - 5414)
random
negBF <-subset(negBF, obs_cooccur <1)
negBF <-subset(negBF, prob_cooccur <0.05)
negBF <- merge(negBF, ITStax.G, by.x = "sp1_name", by.y = 0)
str(negBF)
colnames(negBF) <- c("sp1_name", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                     "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "sp2_name", "PhylumITS", "ClassITS","OrderITS", "FamilyITS","GenusITS", "SpeciesITS")
negBF <- merge(negBF,S16tax.G, by.x = "sp2_name", by.y = 0)
str(negBF)
colnames(negBF) <- c("sp2_name","sp1_name", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                     "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "PhylumITS", "ClassITS","OrderITS", "FamilyITS","GenusITS", "SpeciesITS",
                     "Phylum16S", "Class16S","Order16S", "Family16S","Genus16S", "Species16S")
write.csv(negBF, "negBFBRD.csv")

### healthy samples bacteria-fungi
Healthy_bacfun <- read.csv("Healthy_bacfun.csv", sep = ",", header = T) #total associations 281420
posBF.H <- subset(Healthy_bacfun, p_gt < 0.05) ##positive 8734
posBF.H <- subset(posBF.H, obs_cooccur > 44) ### probability that the species will be in at least 60% of the total samples
posBF.H<- subset(posBF.H, prob_cooccur > 0.75) ### probability that the two species will be in the same site

#MERGE taxonomy
posBF.H <- merge(posBF.H, ITStax.G, by.x = "sp1_name", by.y = 0)
str(posBF.H)
colnames(posBF.H) <- c("sp1_name", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                     "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "sp2_name", "PhylumITS", "ClassITS","OrderITS", "FamilyITS","GenusITS", "SpeciesITS")
posBF.H <- merge(posBF.H,S16tax.G, by.x = "sp2_name", by.y = 0)
str(posBF.H)
colnames(posBF.H) <- c("sp2_name","sp1_name", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                     "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "PhylumITS", "ClassITS","OrderITS", "FamilyITS","GenusITS", "SpeciesITS",
                     "Phylum16S", "Class16S","Order16S", "Family16S","Genus16S", "Species16S")
write.csv(posBF.H, "posBFHealthy.csv")

##test if there is fungi associated with BRD-pathogens
mycoFH <- merge(Healthy_bacfun, myco, by.x = "sp2_name", by.y = 0)
str(mycoFH)
mycoFH$sp1_name <- as.factor(mycoFH$sp1_name) #196 ASVs
unique(mycoFH$sp1_name)
mycoFHsP <-subset(mycoFH, p_gt < 0.05) #positive
mycoFHsP %>% count(sp1_name)
mycoFHsP$sig <- c("Positive")
mycoFHsP <-mycoFHsP[,-c(12:16)]
mycoFHsP <- merge(mycoFHsP,ITStax.G, by.x = "sp1_name", by.y = 0)

mycoFHsN <-subset(mycoFH, p_lt < 0.05) #negative
mycoFHsN %>% count(sp1_name)
mycoFHsN$sig <- c("Negative")
mycoFHsN <-mycoFHsN[,-c(12:16)]
mycoFHsN <- merge(mycoFHsN,ITStax.G, by.x = "sp1_name", by.y = 0)

mycoSigH <- rbind(mycoFHsP,mycoFHsN)
write.csv(mycoSigH, "mycoSigHealthy.csv")

pasteFH <- merge(Healthy_bacfun, paste, by.x = "sp2_name", by.y = 0)
str(pasteFH)
pasteFH$sp1_name <- as.factor(pasteFH$sp1_name) #308 ASVs
unique(pasteFH$sp1_name)
pasteFHsP <-subset(pasteFH, p_gt < 0.05) #positive
pasteFHsN <-subset(pasteFH, p_lt < 0.05) #negative

histoFH <- merge(Healthy_bacfun, histo, by.x = "sp2_name", by.y = 0)
unique(histoFH$sp1_name)
histoFH$sp1_name <- as.factor(histoFH$sp1_name) #308 ASVs
unique(histoFH$sp1_name)
histoFHsP <-subset(histoFH, p_gt < 0.05) #p_lt is considered p-value
histoFHsP %>% count(sp1_name)
histoFHsP$sig <- c("Positive")
histoFHsP <-histoFHsP[,-c(12:16)]
histoFHsP <- merge(histoFHsP,ITStax.G, by.x = "sp1_name", by.y = 0)

histoFHsN <-subset(histoFH, p_lt < 0.05) #p_lt is considered p-value
histoFHsN %>% count(sp1_name)
histoFHsN$sig <- c("Negative")
histoFHsN <-histoFHsN[,-c(12:16)]
histoFHsN <- merge(histoFHsN,ITStax.G, by.x = "sp1_name", by.y = 0)

histoSigH <- rbind(histoFHsP,histoFHsN)
write.csv(histoSigH, "histoSigHealthy.csv")

mannF <- merge(BRD_bacfun, mann, by.x = "sp2_name", by.y = 0)

##negative associations
negBF.H <-subset(Healthy_bacfun, p_lt < 0.05) #negative 7887
randomH <- (281420 - 8734 - 7887)
negBF.H <-subset(negBF.H, obs_cooccur <1) 
negBF.H <-subset(negBF.H, prob_cooccur <0.05)
negBF.H <- merge(negBF.H, ITStax.G, by.x = "sp1_name", by.y = 0)
str(negBF.H)
colnames(negBF.H) <- c("sp1_name", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                     "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "sp2_name", "PhylumITS", "ClassITS","OrderITS", "FamilyITS","GenusITS", "SpeciesITS")
negBF.H <- merge(negBF.H,S16tax.G, by.x = "sp2_name", by.y = 0)
str(negBF.H)
colnames(negBF.H) <- c("sp2_name","sp1_name", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                     "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "PhylumITS", "ClassITS","OrderITS", "FamilyITS","GenusITS", "SpeciesITS",
                     "Phylum16S", "Class16S","Order16S", "Family16S","Genus16S", "Species16S")
write.csv(negBF.H, "negBFHealthy.csv")

##count the negative interactions 
negBFHealthy <- read.csv("negBFHealthy.csv", sep = ",", header = T)
str(negBFHealthy)
negBFHealthy$sp1_name <- as.factor(negBFHealthy$sp1_name)
negBFHealthy$GenusITS <- as.factor(negBFHealthy$GenusITS)
negBFHealthy$SpeciesITS <- as.factor(negBFHealthy$SpeciesITS)
negBFHealthy %>% tally()
negG <- negBFHealthy %>% count(GenusITS)
negS <- negBFHealthy %>% count(SpeciesITS)
negITS <- negBFHealthy %>% count(sp1_name)

negITS <- negITS %>% filter(n >= 30)
negITS <- merge(negITS, ITStax, by.x = "sp1_name", by.y = 0)
write.csv(negITS, "negITS.csv")

negBFBRD <- read.csv("negBFBRD.csv", sep = ",", header = T)
str(negBFBRD)
negBFBRD$sp1_name <- as.factor(negBFBRD$sp1_name)
negBFBRD$GenusITS <- as.factor(negBFBRD$GenusITS)
negBFBRD$SpeciesITS <- as.factor(negBFBRD$SpeciesITS)
negBFBRD %>% tally()
negGB <- negBFBRD %>% count(GenusITS)
negSB <- negBFBRD %>% count(SpeciesITS)
negITSB <- negBFBRD %>% count(sp1_name)

negITSB <- negITSB %>% filter(n >= 30)
negITSB <- merge(negITSB, ITStax, by.x = "sp1_name", by.y = 0)
write.csv(negITSB, "negITSB.csv")


