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
setwd("~/Desktop/eunice/PhD/PhDProject/FungalBRD/Qiime/QiimeFungalComplete/exported/")
metadata <- read.csv("qPCRmetadataCorrect06.07.csv", na.strings = c("","NA"), header=TRUE) #metadata containing the corect qPCR data for the IN beef samples (Eunice data)
metadata2 <- metadata
str(metadata2)
metadata$BRD <- as.factor(metadata$BRD)
#metadata$Date.Collection <- mdy(metadata$Date.Collection)
order_groups <- metadata$ID
rownames(metadata) <- metadata[,1]
metadata = metadata[,-c(1)]
str(metadata)

#Random forest for ITS
ITSs <- read_qza("tableFungal.qza") #6930 ITSs
ITS_s <- as.data.frame(ITSs$data) #7050
ITS_table <- as.data.frame(ITSs$data) 
ITS_table$ITSnos <- paste0("ITS", 1:nrow(ITS_table))
ITS_table$ITSstring <- rownames(ITS_table)
rownames(ITS_table) <- ITS_table$ITSnos ##We change the ITS name created in Qiime to ITSn
ITSkey <- ITS_table[, (ncol(ITS_table)-1):ncol(ITS_table)] #the key withe the names
ITS_table <- ITS_table[,-(ncol(ITS_table)-1):-ncol(ITS_table)]
ITS_table <- t(ITS_table)

#Taxonomy of each OTU
##Adding taxonomy
#Taxonomy of each OTU
tax <- read_qza("taxonomyFungal.qza")
tax <- as.data.frame(tax$data)
tax2 = separate(tax, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#This warning means that some cells are empty and that R is replacing empty cells with NA. Now there are other cells that are unclassified that  say, for example `s__` or `g__`. All these need to be removed and replaced with `NA`. 
#All this is OK except that in future use of the taxonomy table, these ASVs will be ignored because they are not classified. Why are ASVs not classified? Its because there is not a close enough match in the database. Just because there is not a good match in the database does not mean they don’t exist, so I wanted to make sure this data was not lost. So in my new code, from lines 300 – 334 I make it so that ASVs that are unclassified at any level are classified as the lowest taxonomic level for which there is a classification.

#All the strings that need to be removed and replaced with NA
na_strings <- c(" s__", " g__", " f__", " o__", " c__")

tax3 = replace_with_na_all(tax2, condition = ~.x %in% na_strings)

#This code is great because an ASV that is unclassified at a certain level are all listed as `NA`.
#Unfortunately this command changed ou Feature.ID names

#Next, all these `NA` classifications with the last level that was classified
tax3[] <- t(apply(tax3, 1, zoo::na.locf))
tax3 <- as.data.frame(tax3)
row.names(tax3) <- tax3[,1]
tax3 = tax3[,-c(1:2)]
tax.clean <- as.data.frame(tax3)
tax.clean$OTUs <- rownames(tax.clean)
#Would be good to check here to make sure the order of the two data frames was the same. You should do this on your own.

###Remove all the OTUs that don't occur in our OTU.clean data set
tax.final = tax.clean[row.names(tax.clean) %in% row.names(ITS_s),]

##Remove unneccesary information from the taxonomy names
tax.final$Phylum <- sub("k__*", "", tax.final[,1])
tax.final$Phylum <- sub("p__*", "", tax.final[,1])
tax.final$Class <- sub("p__*", "", tax.final[,2])
tax.final$Class <- sub("c__*", "", tax.final[,2])
tax.final$Class <- sub("k__*", "", tax.final[,2])
tax.final$Order <- sub("p__*", "", tax.final[,3])
tax.final$Order <- sub("c__*", "", tax.final[,3])
tax.final$Order <- sub("o__*", "", tax.final[,3])
tax.final$Order <- sub("k__*", "", tax.final[,3])
tax.final$Family <- sub("p__*", "", tax.final[,4])
tax.final$Family <- sub("c__*", "", tax.final[,4])
tax.final$Family <- sub("o__*", "", tax.final[,4])
tax.final$Family <- sub("f__*", "", tax.final[,4])
tax.final$Family <- sub("k__*", "", tax.final[,4])
tax.final$Genus <- sub("p__*", "", tax.final[,5])
tax.final$Genus <- sub("c__*", "", tax.final[,5])
tax.final$Genus <- sub("o__*", "", tax.final[,5])
tax.final$Genus <- sub("f__*", "", tax.final[,5])
tax.final$Genus <- sub("g__*", "", tax.final[,5])
tax.final$Genus <- sub("k__*", "", tax.final[,5])
tax.final$Species <- sub("p__*", "", tax.final[,6])
tax.final$Species <- sub("c__*", "", tax.final[,6])
tax.final$Species <- sub("o__*", "", tax.final[,6])
tax.final$Species <- sub("f__*", "", tax.final[,6])
tax.final$Species <- sub("g__*", "", tax.final[,6])
tax.final$Species <- sub("s__*", "", tax.final[,6])
tax.final$Species <- sub("k__*", "", tax.final[,6])

TaxITS <- merge(tax.final, ITSkey, by.x = 0, by.y = "ITSstring")
row.names(TaxITS) <- TaxITS[,10]
TaxITS = TaxITS[,-c(1,10)]

### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(ITS_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxITS))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq #[ 7050 taxa and 131 samples ] ITS data

colnames(tax_table(physeq_deseq))
## Filter any non-baxteria, chloroplast and mitochondria
# Step 1: subset the samples based on healthy or sick

#Pruning the data
# Set prunescale 
prunescale = 0.0001 #We will do this by eliminating OTUs with an average relative abundance below 0.0001
minlib = 40420 # this is the value used to rarified the data
# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean <- taxa_sums(physeq_deseq)/nsamples(physeq_deseq)
sites.prune <- prune_taxa(tax.mean > prunescale*minlib, physeq_deseq)
prunetable<- phyloseq_to_df(sites.prune, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance") ## the row names always have to start with letter and not with number

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]
#write.csv(NewTax, "NewTaxFungiCooccur.csv")


sites.prune # [ 1109 taxa and 131 samples ]
#[ 1236 taxa and 131 samples ]
# Make a dataframe of training data with OTUs as column and samples as rows
predictors <- (otu_table(sites.prune))
dim(predictors) #131 1109 

# Make one column for our outcome/response variable 
response <- as.factor(sample_data(sites.prune)$BRD) ##classification
response

# Combine them into 1 data frame
rf.data <- data.frame(response, predictors) ## complete dataframe with the sample ID, classification and OTU table
str(rf.data)
#write.table(rf.data,"OtutableRF.txt",sep=",", row.names = TRUE) 

### other way
# Setting the testing and training set
str(metadata)
set.seed(123)
sample <- sample.int(n = nrow(metadata), size = floor(.60*nrow(metadata)), replace = F)
train <- metadata[sample, ]
str(train)
train %>% tally()
train %>% count(BRD)
#BRD 33
#Healthy 46

test  <- metadata[-sample, ]
test %>% tally()
test %>% count(BRD)
#BRD 25
#Healthy 29
str(test)

## performing random forest for the 16S data
trainITS <- merge(train, rf.data, by.x = 0, by.y = 0)
rownames(trainITS) <- trainITS[,1]
trainITS = trainITS[,-c(1:6)]

rm(ITS)
set.seed(123)
ITS <- randomForest(response ~ .,data=trainITS, importance=TRUE)
print(ITS)
plot(ITS)
error <- ITS$err.rate
set.seed(123)
ITS.1 <- randomForest(response ~ .,data=trainITS, importance=TRUE, ntree=17) # with 17 trees, the lowest classification errof for BRD and healthy
print(ITS.1)
plot(ITS.1) ### 17 we the the lowest error instead ofn 500 trees

#rfcv：Random Forest Cross Validation
result = rfcv(trainITS, trainITS$response, cv.fold=10)
result$error.cv #corresponding vector of error rates or MSEs at each step
#the mean squared error (MSE) or mean squared deviation (MSD) of an estimator 
#(of a procedure for estimating an unobserved quantity) measures the average of the squares of the errors—that is, the average squared difference between the estimated values and the actual value.
#plot
#10% of the data is used to test the function
# Draw validation, results
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))

# Simple visualization
varImpPlot(ITS.1)
varImpPlot(ITS.1, main = "Top 17-Feature OTU importance",n.var = 17)
varImpPlot(ITS.1, main = "Top 17-Feature OTU importance",n.var = 17, bg = par("bg"),
           color = par("fg"), gcolor = par("fg"), lcolor = "gray")
imp <- varImpPlot(ITS.1) # let's save the varImp object
imp <- as.data.frame(imp)
imp = imp[order(imp[,1],decreasing = T),]
imp = imp[c(1:17),]

# this part just creates the data.frame for the plot part
library(dplyr)
imp$varnames <- rownames(imp) # row names to column

# imp species name decomposition
# Add Class level to remain in line
str(TaxITS)
imp <- merge(imp, TaxITS, by.x = 0, by.y = 0)
imp
imp = imp[order(imp[,2],decreasing = T),]
imp$Row.names <- as.factor(imp$Row.names)
imp$ASV <-imp$Row.names
levels(imp$Row.names)
levels(imp$Row.names) <- list("ASV106"="ITS106","ASV115"="ITS115","ASV119"="ITS119", "ASV124"="ITS124",
                                         "ASV142"="ITS142", "ASV17"="ITS17", "ASV29"="ITS29", "ASV3"="ITS3", 
                                         "ASV344"="ITS344","ASV41"="ITS41", "ASV43"="ITS43", "ASV460"="ITS460", "ASV53"="ITS53",
                                         "ASV540"="ITS540", "ASV62"="ITS62", "ASV969"="ITS969", "ASV14"="ITS14"
)
imp$SpeciesASV <- paste(imp$Species, imp$Row.names, sep="_")
write.csv(imp , "imporITS.csv")

#Plotting variables 
ggplot(data = imp, mapping = aes(x=reorder(SpeciesASV, MeanDecreaseAccuracy), y= MeanDecreaseAccuracy,fill=Phylum)) +
  geom_bar(stat="identity")+coord_flip() +theme_bw() +xlab("Species") +
  scale_fill_manual(values = my_colors) +
  theme(legend.text = element_text(size=12, face = "italic")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12, face="italic")) 


trainITS <- as.data.frame(ITS.1[["predicted"]])
trainITS$ID <- rownames(trainITS)
str(trainITS)
trainR <- merge(trainITS, metadata, by.x = "ID", by.y = 0)


#testing data
testITS <- merge(test, rf.data, by.x = 0, by.y = 0)
rownames(testITS) <- testITS[,1]
testITS = testITS[,-c(1:6)]

predITS.2 = predict(ITS.1, newdata=testITS)#predictions, gives the classification
predITS.2
pITS.2 <- as.data.frame(predITS.2)
pITS.2$ID <- rownames(pITS.2)
pITS.2 %>% count(predITS.2)
aITS <- table(testITS$response, pITS.2$predITS.2)
abITS <- as.data.frame(aITS)

#Confusion matrix
matrixITS = aITS
matrixITS

#Misclasification rate
testITS %>%
  mutate(lda.pred = (pITS.2$pred)) %>%
  summarise(lda.error = mean(response != lda.pred))


#with this predition, we get the probabily of either being classified as healthy or BRD
str(testITS)
predITS = predict(ITS.1, newdata=testITS,type='prob')
predITS
pITS <- as.data.frame(predITS)
pITS$ID <- rownames(pITS)
pITS
colnames(pITS) <- c("Prob.BRDITS", "Prob.HealthyITS", "ID")

## we need to identify the abundance of the 17-spare ASVs that the model used to classify BRD and healthy animals
classification <- read.csv("RandomForestClass.csv", na.strings = c("","NA"), header=TRUE) #metadata containing the corect qPCR data for the IN beef samples (Eunice data)
classASV <- merge(classification, testITS, by.x = "ID", by.y = 0)
classASV <- classASV[,-c(2:4)]
rownames(classASV) <- classASV$ID
classASV <- classASV[,-c(1)]
classASV <- as.data.frame(t(classASV))
#now we merge with the 17-spare ASVs
classASV <- merge(classASV, imp, by.x = 0, by.y = "ASV")
rownames(classASV) <- classASV$Row.names
classASV <- classASV[,-c(1, 30:42)]
classASV <- as.matrix(t(classASV))
otu.summary <- prop.table(classASV, 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples
a <- as.data.frame(otu_abund)
sum(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:17)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
str(metadata)
str(melt_otu)
meta_otu <- merge(classification, melt_otu, by.x = "ID", by.y = "Sample")
str(meta_otu)
meta_otu_tax <- merge(meta_otu,imp, by.x = "ASV", by.y = "ASV")
str(meta_otu_tax)
order_groups <- classification$ID
meta_otu_tax$Row.names <- factor(meta_otu_tax$ID, levels = order_groups)
summary(meta_otu_tax$Row.names) ###to check that all the samples have the same number of OTUs (6199 total, same value from the taxonomy file) 
meta_otu_tax$Family <- factor(meta_otu_tax$Family)
meta_otu_tax$Classification <- factor(meta_otu_tax$Classification)
meta_otu_tax$Genus <- factor(meta_otu_tax$Genus)
meta_otu_tax$Phylum <- factor(meta_otu_tax$Phylum)
meta_otu_tax$RFStatus <- factor(meta_otu_tax$RFStatus)
meta_otu_tax$ASV <- factor(meta_otu_tax$ASV)
meta_otu_tax$SpeciesASV <- factor(meta_otu_tax$SpeciesASV)
str(meta_otu_tax)

### Abundance at a phlyum level
BRD <- subset(meta_otu_tax, RFStatus=="BRD")
sum(BRD$Abundance)
BRD$ID <- as.factor(BRD$ID) 
### Abundance at a genus level
GenusAB <- BRD %>% 
  group_by(ID, Classification,  SpeciesASV) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Classification,  SpeciesASV) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(GenusAB)
sum(GenusAB$taxa.average)
GenusAB$SpeciesASV <- factor(GenusAB$SpeciesASV)
write.csv(GenusAB , "GenusAB.csv")

my_colors <- c(
  '#a6cee3','#1f78b4','#b3df8a','#33a03c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab3d6','#6a3d9a','#ffff99','#b15938', 
  "#CBD588", "#5F7FC7", "orange","#DA5734", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14385", "#653936", "#C84348", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

ggplot(GenusAB, aes(x = Classification, y = taxa.average, fill =SpeciesASV)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  #facet_grid(Status~.)+
  #facet_wrap(vars(BRD), scales = "free") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.text = element_text(size=11, face = "italic")) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=13, face="bold"), axis.title.y = element_text(color="black", size=13, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 10)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='RF Classification as BRD')

##relative abundance per sample
GenusAB2 <- BRD %>% 
  group_by(ID, Classification,  SpeciesASV) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(ID, Classification,  SpeciesASV) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(GenusAB)
sum(GenusAB$taxa.average)
GenusAB$SpeciesASV <- factor(GenusAB$SpeciesASV)
write.csv(GenusAB , "GenusAB.csv")

#Healthy
healthy <- subset(meta_otu_tax, RFStatus=="Healthy")
sum(healthy$Abundance)
healthy$ID <- as.factor(healthy$ID) 
### Abundance at a genus level
genusH <- healthy%>% 
  group_by(SpeciesASV) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/15)*100) ## the total number of samples (51)
attach(genusH)
sum(genusH$Ave_Abundance)
genusH <- genusH[order(-Ave_Abundance),]
write.csv(genusH,"RF.ITShealthy.csv")

top10H <- merge(healthy, imp, by.x = "SpeciesASV", by.y="SpeciesASV")
top10H$SpeciesASV <- as.factor(top10H$SpeciesASV)
levels(top10H$SpeciesASV)
str(top10H)

GenusAH <- top10H %>% 
  group_by(ID, Classification,  SpeciesASV) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Classification,  SpeciesASV) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(GenusAH)
sum(GenusAH$taxa.average)
GenusAH$SpeciesASV <- factor(GenusAH$SpeciesASV)
write.csv(GenusAH , "GenusAH.csv")

ggplot(GenusAH, aes(x = Classification, y = taxa.average, fill =SpeciesASV)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  #facet_grid(Status~.)+
  #facet_wrap(vars(BRD), scales = "free") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.text = element_text(size=11, face = "italic")) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=13, face="bold"), axis.title.y = element_text(color="black", size=13, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 10)) +
  ylab(paste0("Relative Abundance")) +  labs(x='RF Classification as Healthy')

###to identify the fungal abundance in the animals identified as BRD+ V
names <- c("Dairy107", "Dairy27", "Dairy28", "Dairy28", "Dairy29", "Dairy93")
BRDclass <- GenusAB2[GenusAB2$ID %in% c(names),]
write.csv(BRDclass, "BRDclass.csv")


###QPCR data
trainqPCR <- merge(train, rf.data, by.x = 0, by.y = 0)
rownames(trainqPCR) <- trainqPCR[,1]
trainqPCR = trainqPCR[,-c(1,7:1116)]

set.seed(123)
qpcr <- randomForest(BRD~ .,data=trainqPCR, importance=TRUE)
print(qpcr)
plot(qpcr)
erro.rate2 <- qpcr[["err.rate"]]
rm(qpcr.1)
set.seed(123)
qpcr.1 <- randomForest(BRD~ .,data=trainqPCR, importance=TRUE, ntree=35) #with 35 trees, the lowest classification errof for BRD and healthy
print(qpcr.1)
plot(qpcr.1)
trainqpcr <- as.data.frame(qpcr.1[["predicted"]])
trainqpcr$ID <- rownames(trainqpcr)
str(trainqpcr)
trainqpcr <- merge(trainqpcr, metadata, by.x = 0, by.y = 0)

#rfcv：Random Forest Cross Validation
result2 = rfcv(trainqPCR, trainqPCR$BRD, cv.fold=10)
result2$error.cv #corresponding vector of error rates or MSEs at each step
with(result2, plot(n.var, error.cv, log="x", type="o", lwd=2))

# Export feature importance
varImpPlot(qpcr.1)
imp2 <- varImpPlot(qpcr.1) # let's save the varImp object
imp2 <- as.data.frame(imp2)
imp2 = imp2[order(imp2[,1],decreasing = T),]

#testing set
testqPCR <- merge(test, rf.data, by.x = 0, by.y = 0)
rownames(testqPCR) <- testqPCR[,1]
testqPCR = testqPCR[,-c(1,7:1116)]

predQPCR = predict(qpcr.1, newdata=testqPCR) #predictions
predQPCR
predQPCR <- as.data.frame(predQPCR)
predQPCR$ID <- rownames(predQPCR)
predQPCR %>% count(predQPCR)
aPCR <- table(testqPCR$BRD, predQPCR$predQPCR)
Apqpcr <- as.data.frame(aPCR)


#Confusion matrix
matrixqpcr = Apqpcr
matrixqpcr

#Misclasification rate
testqPCR %>%
  mutate(lda.pred = (predQPCR$pred)) %>%
  summarise(lda.error = mean(BRD != lda.pred))

predqpcr = predict(qpcr.1, newdata=testqPCR, type='prob') #predictions
predqpcr
pqpcr <- as.data.frame(predqpcr)
pqpcr$ID <- rownames(pqpcr)
pqpcr
colnames(pqpcr) <- c("Prob.BRDqPCR", "Prob.HealthyqPCR", "ID")

##check the abundance of BRD-pathobionts in the animals classified as BRD and healthy by ITS
classification2 <- read.csv("RandomForestClassqPCR.csv", na.strings = c("","NA"), header=TRUE) #metadata containing the corect qPCR data for the IN beef samples (Eunice data)
classPCR <- merge(classification2, testqPCR, by.x = "ID", by.y = 0)
classPCR <- classPCR[,-c(2:4)]
rownames(classPCR) <- classPCR$ID
classPCR <- classPCR[,-c(1)]
classPCR <- as.data.frame(t(classPCR))
classPCR <- as.matrix(t(classPCR))
otu.summary <- prop.table(classPCR, 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples
a <- as.data.frame(otu_abund)
sum(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:4)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
str(melt_otu)
meta_otu <- merge(classification2, melt_otu, by.x = "ID", by.y = "Sample")
str(meta_otu)
order_groups <- classification2$ID
meta_otu$Row.names <- factor(meta_otu$ID, levels = order_groups)
summary(meta_otu$Row.names) ###to check that all the samples have the same number of OTUs (6199 total, same value from the taxonomy file) 
meta_otu$Classification <- factor(meta_otu$Classification)
meta_otu$ID <- factor(meta_otu$ID)
meta_otu$RFStatus <- factor(meta_otu$RFStatus)
str(meta_otu)

### Abundance at a phlyum level
BRD <- subset(meta_otu, RFStatus=="BRD")
BRD$ID <- factor(BRD$ID)
str(BRD)
sum(BRD$Abundance)
### Abundance at a genus level
genus <- BRD%>% 
  group_by(ASV) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/12)*100) ## the total number of samples (51)
attach(genus)
sum(genus$Ave_Abundance)
genus <- genus[order(-Ave_Abundance),]
write.csv(genus,"RF.qPCRBRD.csv")

GenusAB <- BRD %>% 
  group_by(ID, Classification,  ASV) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Classification,  ASV) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(GenusAB)
sum(GenusAB$taxa.average)
GenusAB$Classification <- factor(GenusAB$Classification)
write.csv(GenusAB , "GenusABqPCR.csv")

my_colors <- c(
  '#a6cee3','#1f78b4','#b3df8a','#33a03c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab3d6','#6a3d9a','#ffff99','#b15938', 
  "#CBD588", "#5F7FC7", "orange","#DA5734", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14385", "#653936", "#C84348", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

ggplot(GenusAB, aes(x = Classification, y = taxa.average, fill =ASV)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  #facet_grid(Status~.)+
  #facet_wrap(vars(BRD), scales = "free") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.text = element_text(size=11)) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=13, face="bold"), axis.title.y = element_text(color="black", size=13, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 10)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='RF Classification as BRD')

#Healthy
healthy <- subset(meta_otu, RFStatus=="Healthy")
healthy$ID <- as.factor(healthy$ID) 
sum(healthy$Abundance)
str(healthy)
### Abundance at a genus level
genusH <- healthy%>% 
  group_by(ASV) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/15)*100) ## the total number of samples (51)
attach(genusH)
sum(genusH$Ave_Abundance)
genusH <- genusH[order(-Ave_Abundance),]
write.csv(genusH,"RF.qPCRhealthy.csv")

GenusAH <- healthy %>% 
  group_by(ID, Classification,  ASV) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Classification,  ASV) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(GenusAH)
sum(GenusAH$taxa.average)
GenusAH$ASV <- factor(GenusAH$ASV)
write.csv(GenusAH , "GenusAHqPCR.csv")

ggplot(GenusAH, aes(x = Classification, y = taxa.average, fill =ASV)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  #facet_grid(Status~.)+
  #facet_wrap(vars(BRD), scales = "free") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.text = element_text(size=11)) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=13, face="bold"), axis.title.y = element_text(color="black", size=13, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 10)) +
  ylab(paste0("Relative Abundance")) +  labs(x='RF Classification as Healthy')

#now check the abundance of the animals identified as BRD-V but not qPCR
#to determine if the fungal played a better job than qPCR
names #BDR fungal
tqPCR <- testqPCR
tqPCR$ID <- rownames(tqPCR)
BRDclassqPCR <- tqPCR[tqPCR$ID %in% c(names),]
BRDclassqPCR <- BRDclassqPCR[,-c(1,6)]
BRDclassqPCR <- as.data.frame(t(BRDclassqPCR))
BRDclassqPCR <- as.matrix(BRDclassqPCR)
otu.summary <- prop.table(BRDclassqPCR, 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples
a <- as.data.frame(otu_abund)
sum(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:4)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)
write.csv(melt_otu, "BRDpathobiontsabundanceBRDITS.csv")



#plot decision Tree random forest
pqpcr
pred <- merge(pqpcr, metadata, by.x = "ID", by.y = 0)
pred <- merge(pred, predQPCR, by.x = "ID", by.y = "ID")
pred <- merge(pred, pITS, by.x = "ID", by.y = "ID")
pred <- merge(pred, pITS.2, by.x = "ID", by.y = "ID")
write.table(pred,"RandomForestpred2CorrectqPCREuniceData.txt",sep=",", row.names = TRUE) 

M3 <- read.csv("RandomForestpred2.csv", na.strings = c("","NA"), header=TRUE)
str(M3)
M5 <- read.csv("RFvisual.ITS.csv", na.strings = c("","NA"), header=TRUE)
str(M5)
M5$predicted <- as.factor(M5$predicted)
#Prediction for healthy animals
str(M3)
M3 %>% tally()
M3 %>% count(visual)
#BRD 24
#Healthy 28

M5 %>% tally()
M5 %>% count(predicted)

#other tables
M6 <- read.csv("RFvisual.qPCR.csv", na.strings = c("","NA"), header=TRUE)
str(M6)
M6$predicted <- as.factor(M6$predicted)
M6 %>% tally()
M6 %>% count(predicted)


str(M5)
A <- ggplot() + 
  geom_point(data=M3, aes(x=Prob.BRDITS, y=Prob.BRDqPCR, color=visual), size=3) + 
  theme_classic() + 
  geom_point(data=M5, aes(x=Prob.BRDITS, y=Prob.BRDqPCR, size=3, shape=predicted)) + 
  guides(size=FALSE) +
  geom_vline(xintercept=0.50, color = "red")+
  geom_hline(yintercept=0.50, color = "red")+
  guides(fill=FALSE) +
  #geom_point(data= wuaxis, aes(x=groupX, y=groupX, shape=wuaxis4, size=5, fill=wuaxis4)) +
  scale_shape_manual(values=c(2, 8, 3, 1))+
  labs(color= "Visual Diagnosis") +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(shape= "Agreement 2 methods") +
  #geom_segment(data= wuaxis, aes(x=wuaxis1, y=wuaxis2, xend=groupX, yend=groupY, color= wuaxis4), size = .05) +
  labs(x='Predicted BRD-ITS sequencing', y= 'Predicted BRD-bacterial pathobiont qPCR') +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) 
A

B <- ggplot() + 
  geom_point(data=M3, aes(x=Prob.BRDITS, y=Prob.BRDqPCR, color=visual), size=3) + 
  theme_classic() + 
  geom_point(data=M6, aes(x=Prob.BRDITS, y=Prob.BRDqPCR, size=3, shape=predicted)) + 
  geom_vline(xintercept=0.50, color = "red")+
  geom_hline(yintercept=0.50, color = "red")+
  guides(size=FALSE) +
  guides(fill=FALSE) +
  xlim(0, 1) +
  ylim(0, 1) +
  scale_shape_manual(values=c(2, 8, 3, 1))+
  labs(color= "Visual Diagnosis") +
  labs(shape= "Agreement 2 methods") +
  labs(x='Predicted BRD-ITS sequencing', y= 'Predicted BRD-bacterial pathobiont qPCR') +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) 
B

ggarrange(A, B, labels = c("A", "B"),
          ncol = 2)
##plot the average relative abundance of the top 5 discriminant taxa when classifying healthy and BRD animals
class <- read.csv("RFASV_classification.csv", na.strings = c("","NA"), header=TRUE) #metadata containing the corect qPCR data for the IN beef samples (Eunice data)
str(class )

my_colorsBG <- c(
  'wheat2', "#673770", "#5F7FC7", "orange","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385","khaki1" , "#CD9BCD","darksalmon", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

ggplot(class, aes(x = ASV, y = Ave.Abundance, fill =ASV)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colorsBG) +
  facet_grid(.~Classification) +
  coord_flip() +
  ylim(0, 0.1) + 
  theme(legend.position = "none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = 1, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=13, face="bold"), axis.title.y = element_text(color="black", size=13, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 9), axis.text.y = element_text(color = "black", size = 10, face = "italic")) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Average Relative Abundance")) +  labs(x='ASV')



#Functions
phyloseq_to_df <- function(physeq, addtax = T, addtot = F, addmaxrank = F, sorting = "abundance"){
  
  # require(phyloseq)
  
  ## Data validation
  if(any(addtax == TRUE || sorting == "taxonomy")){
    if(is.null(phyloseq::tax_table(physeq, errorIfNULL = F))){
      stop("Error: taxonomy table slot is empty in the input data.\n")
    }
  }
  
  ## Prepare data frame
  if(taxa_are_rows(physeq) == TRUE){
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), phyloseq::otu_table(physeq), stringsAsFactors = F)
  } else {
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), t(phyloseq::otu_table(physeq)), stringsAsFactors = F)
  }
  
  ## Check if the sample names were silently corrected in the data.frame
  if(any(!phyloseq::sample_names(physeq) %in% colnames(res)[-1])){
    if(addtax == FALSE){
      warning("Warning: Sample names were converted to the syntactically valid column names in data.frame. See 'make.names'.\n")
    }
    
    if(addtax == TRUE){
      stop("Error: Sample names in 'physeq' could not be automatically converted to the syntactically valid column names in data.frame (see 'make.names'). Consider renaming with 'sample_names'.\n")
    }
  }
  
  ## Add taxonomy
  if(addtax == TRUE){
    
    ## Extract taxonomy table
    taxx <- as.data.frame(phyloseq::tax_table(physeq), stringsAsFactors = F)
    
    ## Reorder taxonomy table
    taxx <- taxx[match(x = res$OTU, table = rownames(taxx)), ]
    
    ## Add taxonomy table to the data
    res <- cbind(res, taxx)
    
    ## Add max tax rank column
    if(addmaxrank == TRUE){
      
      ## Determine the lowest level of taxonomic classification
      res$LowestTaxRank <- get_max_taxonomic_rank(taxx, return_rank_only = TRUE)
      
      ## Reorder columns (OTU name - Taxonomy - Max Rank - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), "LowestTaxRank", phyloseq::sample_names(physeq))]
      
    } else {
      ## Reorder columns (OTU name - Taxonomy - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), phyloseq::sample_names(physeq))]
      
    } # end of addmaxrank
  }   # end of addtax
  
  ## Reorder OTUs
  if(!is.null(sorting)){
    
    ## Sort by OTU abundance
    if(sorting == "abundance"){
      otus <- res[, which(colnames(res) %in% phyloseq::sample_names(physeq))]
      res <- res[order(rowSums(otus, na.rm = T), decreasing = T), ]
    }
    
    ## Sort by OTU taxonomy
    if(sorting == "taxonomy"){
      taxtbl <- as.data.frame( phyloseq::tax_table(physeq), stringsAsFactors = F )
      
      ## Reorder by all columns
      taxtbl <- taxtbl[do.call(order, taxtbl), ]
      # taxtbl <- data.table::setorderv(taxtbl, cols = colnames(taxtbl), na.last = T)
      res <- res[match(x = rownames(taxtbl), table = res$OTU), ]
    }
  }
  
  ## Add OTU total abundance
  if(addtot == TRUE){
    res$Total <- rowSums(res[, which(colnames(res) %in% phyloseq::sample_names(physeq))])
  }
  
  rownames(res) <- NULL
  return(res)
}

