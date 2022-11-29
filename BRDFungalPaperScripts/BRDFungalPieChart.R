#                                 R-Script for BRD Fungi
#                     Author: Ruth Eunice Centeno Martinez, 2020
#          Written at Dr. Tim Johnson Lab, Dept. of Animal Sciences,       #
#                             Purdue University, 2020
#                               rcenteno@purdue.edu
library(ggplot2)
library(tidyr) #separate function
library(reshape2) #melt function
library(dplyr)
library(naniar) # for replace_with_na_all function
library(data.table)
library(phyloseq)
library(qiime2R)
library(ggpubr)
library(forcats)

setwd("~/Desktop/eunice/PhD/PhDProject/FungalBRD/Qiime/QiimeFungalComplete/exported/")
##

##### Taxonomy barplot
#The OTU table as exported from qiime has a pound sign before the header row. You need to delete that pound sign in a text editor.
# Read the experimental design, and species classified documents
ASVs <- read_qza("tableFungal.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #3860 ASVs
#str(ASV_table)
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #the key withe the names
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
ASV_table <- t(ASV_table)
ASV_table <- ASV_table[-c(124),]
colSums(ASV_table)

#Importing the metadata file
metadata <- read.csv("BRD.metadadaFungal.csv")
str(metadata)
metadata$BRD <- factor(metadata$BRD) 
metadata$PenCode <- factor(metadata$PenCode)
order_groups <- metadata$ID
row.names(metadata) = metadata[,1]
#metadata = metadata[,-1] # we remove sample 92 because no sequences were identified in the analysis
metadata = metadata[-125,]
metadata$Status <-metadata$BRD 
metadata$Status <- as.factor(metadata$Status)

#merging the abundance of each OTU with the metadata and the taxonomy file
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
tax.final = tax.clean[row.names(tax.clean) %in% row.names(ASV_s),]

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

#write.table(tax.final,"taxonomyNasal.txt",sep=",", row.names = FALSE) 
TaxASV <- merge(tax.final, ASVkey, by.x = 0, by.y = "ASVstring")
row.names(TaxASV) <- TaxASV[,10]
TaxASV = TaxASV[,-c(1,10)]
# write.table(TaxASV,"TaxASV.txt",sep="\t", row.names = T, col.names = T)

H_OTU <- ASV_table[rownames(ASV_table) %in% rownames(metadata),] #51, 3838
str(H_OTU) #129 7050
metadata2 <- metadata[rownames(metadata) %in% rownames(H_OTU),] #51, 3838

#### Separate by BRD and healthy
## -----------------------------------PHYLUM LEVEL, Healthy period
### CALCULATION OF THE ABUNDANCE OF EACH OTU  
## check the prop.table, it gives me NA in the sample Dairy 92
otu.summary <- prop.table(H_OTU, 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples
a <- as.data.frame(otu_abund)
sum(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:7050)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
str(metadata)
str(melt_otu)
meta_otu <- merge(metadata2, melt_otu, by.x = 0, by.y = "Sample")
str(meta_otu)
meta_otu_tax <- merge(meta_otu, TaxASV, by.x = "ASV", by.y = 0)
str(meta_otu_tax)
order_groups <- metadata$ID
meta_otu_tax$Row.names <- factor(meta_otu_tax$Row.names, levels = order_groups)
summary(meta_otu_tax$Row.names) ###to check that all the samples have the same number of OTUs (6199 total, same value from the taxonomy file) 
meta_otu_tax$Family <- factor(meta_otu_tax$Family)
meta_otu_tax$Status <- factor(meta_otu_tax$Status)
meta_otu_tax$Genus <- factor(meta_otu_tax$Genus)
meta_otu_tax$Phylum <- factor(meta_otu_tax$Phylum)
meta_otu_tax$ASV <- factor(meta_otu_tax$ASV)
str(meta_otu_tax)


####subset by healthy or BRD
healthy <- subset(meta_otu_tax, BRD=="Healthy")
BRD <- subset(meta_otu_tax, BRD=="BRD")

###calculate community abundance for each of the two groups
### Abundance at a phlyum level
str(healthy)
healthy$ID <- as.factor(healthy$ID) #73 samples
str(healthy)
phylumH <- healthy %>% 
  group_by(Phylum) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/73)) ## the total number of samples (129)
attach(phylumH)
sum(phylumH$Ave_Abundance)
phylumH <- phylumH[order(-Ave_Abundance),]
write.csv(phylumH,"phylumHealthy.csv")

my_colorsH <- c(
  '#a6cee3','#1f78b4','#b3df8a','#33a03c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab3d6','#6a3d9a','#ffff99','#b15938', 
  "#CBD588", "#5F7FC7", "orange","#DA5734", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14385", "#653936", "#C84348", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)
my_colorsB <- c(
  '#a6cee3','#1f78b4','#b3df8a','#e31a1c','#fb9a99','#33a03c',
  '#fdbf6f','#cab3d6','#ff7f00','#ffff99','#6a3d9a','#b15938', 
  "#CBD588", "#5F7FC7", "orange","#DA5734", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14385", "#653936", "#C84348", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)


library(ggrepel)
a <- ggplot(phylumH, aes(x = "", y = Ave_Abundance, fill = fct_inorder(Phylum))) +
  geom_bar(width = 2, stat = "identity") +
  coord_polar(theta = "y", start = 0) +
  theme_minimal() +
  
  scale_fill_manual(values = my_colorsH) +
  guides(fill = guide_legend(title = "Phylum", keywidth = 1.5, keyheight = 1, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.text = element_text(size=11, face = "italic")) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 10)) +
  ggtitle("Healthy") + 
  theme(plot.title = element_text(color = "black", size = 16, face= "bold", hjust = 0.5))

#BRD phyla
###calculate community abundance for each of the two groups
### Abundance at a phlyum level
str(BRD)
BRD$ID <- as.factor(BRD$ID) #73 samples
str(BRD)
phylumB <- BRD %>% 
  group_by(Phylum) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/56)) ## the total number of samples (129)
attach(phylumB)
sum(phylumB$Ave_Abundance)
phylumB <- phylumB[order(-Ave_Abundance),]
write.csv(phylumB,"phylumBRD.csv")

b <- ggplot(phylumB, aes(x = "", y = Ave_Abundance, fill = fct_inorder(Phylum))) +
  geom_bar(width = 2, stat = "identity") +
  coord_polar(theta = "y", start = 0) +
  theme_minimal() +
  scale_fill_manual(values = my_colorsB) +
  guides(fill = guide_legend(title = "Phylum", keywidth = 1.5, keyheight = 1, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.text = element_text(size=11, face = "italic")) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 10)) +
  #theme(text=element_text(family="Times New Roman"))+ 
  ggtitle("BRD") + 
  theme(plot.title = element_text(color = "black", size = 16, face= "bold", hjust = 0.5))


ggarrange(a,b, labels = c("a", "b"),vjust = 0.9,
          ncol=2)
##now genus
str(healthy)
genusH <- healthy %>% 
  group_by(Genus) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/73)) ## the total number of samples (129)
attach(genusH)
sum(genusH$Ave_Abundance)
genusH <- genusH[order(-Ave_Abundance),]
write.csv(genusH,"genusHealthy.csv")

top10genus <- genusH[c(1:15),] #to select the top 10 most abundant phylum
sum(top10genus$Ave_Abundance) ### the total of the community composed of the top 10
top10G <- merge(healthy, top10genus, by.x = "Genus", by.y="Genus")
top10G$Genus <- as.factor(top10G$Genus)
levels(top10G$Genus)
str(top10G)

GenusAB <- top10G %>% 
  group_by(Row.names, Date.Collection, Genus) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Date.Collection, Genus) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(GenusAB)
sum(GenusAB$taxa.average)
GenusAB$Genus <- factor(GenusAB$Genus)
GenusAB$Date.Collection <- factor(GenusAB$Date.Collection)
levels(GenusAB$Date.Collection)
levels(GenusAB$Date.Collection) <- list("12/2"="12/2/20", "11/18"="11/18/20","11/11"="11/11/20","11/4"="11/4/20","10/28"="10/28/20",
                                        "10/14"="10/14/20","9/30"="9/30/20","9/10"="9/10/20","8/26"="8/26/20","8/19"="8/19/20",
                                        "8/12"="8/12/20","7/21"="7/21/20","7/14"="7/14/20")
write.csv(GenusAB,"genusHealthyNasalStatus.csv")

#Plot the graph 
my_colorsHG <- c(
  'wheat2', "#673770", "#5F7FC7", "orange","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385", "#653936", "#CD9BCD", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)
my_colorsBG <- c(
  'wheat2', "#673770", "#5F7FC7", "orange","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385","khaki1" , "#CD9BCD","darksalmon", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)
c <- ggplot(GenusAB, aes(x = Date.Collection, y = taxa.average, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colorsHG) +
  coord_flip() +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = 1, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.text = element_text(size=11, face = "italic")) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=13, face="bold"), axis.title.y = element_text(color="black", size=13, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 10)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Date of Collection')

# BRD genus
str(healthy)
genusB <- BRD %>% 
  group_by(Genus) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/56)) ## the total number of samples (129)
attach(genusB)
sum(genusB$Ave_Abundance)
genusB <- genusB[order(-Ave_Abundance),]
write.csv(genusB,"genusBRD.csv")

top10genusB <- genusB[c(1:15),] #to select the top 10 most abundant phylum
sum(top10genusB$Ave_Abundance) ### the total of the community composed of the top 10
top10GB <- merge(BRD, top10genusB, by.x = "Genus", by.y="Genus")
top10GB$Genus <- as.factor(top10GB$Genus)
levels(top10GB$Genus)
str(top10GB)

GenusABB <- top10GB %>% 
  group_by(Row.names, Date.Collection, Genus) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Date.Collection, Genus) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(GenusABB)
sum(GenusABB$taxa.average)
GenusABB$Genus <- factor(GenusABB$Genus)
GenusABB$Date.Collection <- factor(GenusABB$Date.Collection)
levels(GenusABB$Date.Collection)
levels(GenusABB$Date.Collection) <- list("12/2"="12/2/20", "11/18"="11/18/20","11/11"="11/11/20","11/4"="11/4/20","10/28"="10/28/20",
                                        "10/14"="10/14/20","9/30"="9/30/20","9/10"="9/10/20","8/26"="8/26/20","8/19"="8/19/20",
                                        "8/12"="8/12/20","7/21"="7/21/20","7/14"="7/14/20")

write.csv(GenusABB,"genusBRDNasalStatus.csv")

#Plot the graph 
d <- ggplot(GenusABB, aes(x = Date.Collection, y = taxa.average, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colorsBG) +
  coord_flip() +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = 1, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.text = element_text(size=11, face = "italic")) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=13, face="bold"), axis.title.y = element_text(color="black", size=13, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 10)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Date of Collection')

ggarrange(c,d, labels = c("c", "d"),vjust = 0.9,
          ncol=2)
ggarrange(a, b,c,d, labels = c("a","b","c", "d"),vjust = 0.9,
          ncol=2, nrow=2)
