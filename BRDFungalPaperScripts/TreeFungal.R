###creating phylogenetic tree
library(qiime2R)
library(phyloseq)
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("ggtree")
library(ggtree)
library(ape)
#install.packages("ggjoy")
library(ggjoy)
library(dplyr)
library(tidyr) #for separate function
library(naniar)# Ffor replace all function
library(tidyr)
library(tidytree)


### importing my phylogenetic tree created in Qiime2
setwd("~/Desktop/eunice/PhD/PhDProject/FungalBRD/Qiime/QiimeFungalComplete/exported/")

# we need to make the phyloseq object
metadata <-read.csv("BRD.metadadaFungal.csv", na.strings = c("","NA"), header=TRUE)
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

ASVs <- read_qza("tableFungal.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #7050 ASVs
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #the key withe the names
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
ASV_table <- t(ASV_table)
class(ASV_table)
ASV_table = ASV_table[-124,]

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
TaxASV = TaxASV[,-c(1)]

#subset the taxonomy based on the family of Issatckenckia and trichosporon
##sequences
seq <- read_qza("rep-seqsFungal.qza")
seq <- as.data.frame(seq$data)
seq <- merge(seq, TaxASV, by.x = 0, by.y = "OTUs")
write.csv(seq, "seq.csv")

Alternaria <- subset(seq, Genus=="Alternaria") #for trichosporon outgroup
write.csv(Alternaria, "Alternaria.csv")

fusarium<- subset(seq, Genus=="Fusarium") #for trichosporon outgroup
write.csv(fusarium, "fusarium.csv")

#import the significant sequence from Deseq
DeSeq = read.csv("sigtabFungi.csv")
str(DeSeq)
DeSeq$High_low <- factor(DeSeq$High_low)
DeSeq$Genus <- factor(DeSeq$Genus)
Tris <- subset(DeSeq, Genus =="Trichosporon")
signames <- Tris$SpeciesASV
#"Trichosporon-ASV45"  "Trichosporon-ASV16"  "Trichosporon-ASV179" "Trichosporon-ASV101"
#"Trichosporon-ASV33"  "Trichosporon-ASV1" 

##visualizing the tree
library(tidyverse)
library(ggtree)

##Trichosporon treee
## we need to import a Newick tree
tree <- read.tree("TrichosporonBootStrap1000.tre")
tree
ggtree(tree)
treeinfo <- read.csv("TreeDataTris.csv", na.strings = c("","NA"), header=TRUE)

#to color specific species
cls <- list(HighBRD=c("Trichosporon-ASV16", "Trichosporon-ASV45", "Trichosporon-ASV179"),
            HighHealthy=c("Trichosporon-ASV33", "Trichosporon-ASV1", "Trichosporon-ASV101"),
            Outgroup=c("Vanrija-humicola-ASV6801", "Vanrija-humicola-ASV3781"),
            RefSequence=c("Trichosporon-ovoides-2NF903Bctg_85-NCBI", "Trichosporon-sp-Bio4-Blast", "Trichosporon-sp-mYJrh51-Blast"))

#to root the tree
id <- as.data.frame(tree$tip.label)
tree <- root(tree, 382:383)
tree2 <- groupOTU(tree, cls)
ggtree(tree2)

lenght <- as.data.frame(tree2[["edge.length"]])
my_colors <- c(
  "gray28","lightcoral", "turquoise3", "lightgreen", "purple")

#plot the tree
p <- ggtree(tree2, aes(color=group)) + 
  scale_color_manual(values=my_colors)
p <- ggtree(tree2, layout = "rectangular", aes(color=group)) + 
  scale_color_manual(values=my_colors)
nodes <- p[["data"]]

#to get the node number
p3 <- ggtree(tree2, aes(color=group), branch.length = 'non') +
  scale_color_manual(values=my_colors)  +
  geom_text2(aes(subset=!isTip, label=node))
p3
# expand node 570, 706, 647
#compress node 908,984,1014,1058,598,763, 61
p2 <- p %>% collapse(node=881) +
  geom_point2(aes(subset=(node==881)), shape=19, size=3) 
p2 <- collapse(p2, node=686) +
  geom_point2(aes(subset=(node==686)), shape=19, size=3) 
p2 <- p2 + #geom_cladelab(node=381, label="Trichosporon-mYJrh51",fontsize=3,offset.text=0.1) +
  #xlim(0, 0.9) +
  #geom_cladelab(node=465, label="Trichosporon-ovoides-2NF903Bctg_85", fontsize=3,offset.text=0.1) + 
  #geom_cladelab(node=522, label="Trichosporon-Bio4", fontsize=3,offset.text=0.1) +
  #geom_cladelab(node=382, label="Vanrija-humicola-ASV6801", fontsize=3,offset.text=0.1) +
  #geom_cladelab(node=383, label="Vanrija-humicola-ASV3781", fontsize=3,offset.text=0.1) +
  geom_cladelab(node=8, label="Trichosporon-ASV1", fontsize=3,offset.text=0.1) +
  geom_cladelab(node=53, label="Trichosporon-ASV101", fontsize=3,offset.text=0.1) +
  geom_cladelab(node=314, label="Trichosporon-ASV33", fontsize=3,offset.text=0.1) +
  geom_cladelab(node=111, label="Trichosporon-ASV16", fontsize=3,offset.text=0.1) +
  geom_cladelab(node=195, label="Trichosporon-ASV179", fontsize=3,offset.text=0.1) +
  geom_cladelab(node=317, label="Trichosporon-ASV45", fontsize=3,offset.text=0.1) +
  labs(color= "Classification") + theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.text = element_text(size=11)) + theme_tree2() 
p2

##---- Issatchenkia treee
## we need to import a Newick tree
tree <- read.tree("IssatchenkiaBootStrap1000.tre")
tree
Issat <- subset(DeSeq, Genus=="Issatchenkia")

treeinfo <- read.csv("TreeDataIssat.csv", na.strings = c("","NA"), header=TRUE)

#to color specific species
cls <- list(HighBRD=c("Issatchenkia-ASV129", "Issatchenkia-ASV222", "Issatchenkia-ASV37"),
            HighHealthy=c("Issatchenkia-ASV191", "Issatchenkia-ASV24", "Issatchenkia-ASV27",
                          "Issatchenkia-ASV36", "Issatchenkia-ASV56", "Issatchenkia-ASV61",
                          "Issatchenkia-ASV73", "Issatchenkia-ASV75"),
            Outgroup=c("Kazachstania-ASV1688", "Kazachstania-ASV6787"),
            RefSequence=c("Pichia-kudriavzevii-VN9Y-Blast", "Pichia-kudriavzevii-ATCC6258-ITS", "Pichia-kudriavzevii-3Y3-Blast"))

#to root the tree
id <- as.data.frame(tree$tip.label)
tree <- root(tree, 281:282)
tree2 <- groupOTU(tree, cls)

#plot the tree
p <- ggtree(tree2, aes(color=group)) + 
  scale_color_manual(values=my_colors)
p
p <- ggtree(tree2, layout = "rectangular", aes(color=group)) + 
  scale_color_manual(values=my_colors) + geom_tiplab(as_ylab=TRUE, color='firebrick')
p


#to get the node number
p3 <- ggtree(tree2, aes(color=group),branch.length = 'none') +
  scale_color_manual(values=my_colors)  +
  geom_text2(aes(subset=!isTip, label=node))
p3
# expand node 570, 706, 647
#compress node 908,984,1014,1058,598,763, 61
nodes <- p[["data"]]
p2 <- p + geom_cladelab(node=21, label="Issatchenkia-ASV24",fontsize=3,offset.text=0.1) +
  #xlim(0, 0.9) +
  geom_cladelab(node=23, label="Issatchenkia-ASV56", fontsize=3,offset.text=0.1) + 
  geom_cladelab(node=139, label="Issatchenkia-ASV61", fontsize=3,offset.text=0.1) +
  geom_cladelab(node=141, label="Issatchenkia-ASV27", fontsize=3,offset.text=0.1) +
  geom_cladelab(node=180, label="Issatchenkia-ASV75", fontsize=3,offset.text=0.1) +
  geom_cladelab(node=183, label="Issatchenkia-ASV191", fontsize=3,offset.text=0.1) +
  geom_cladelab(node=219, label="Issatchenkia-ASV73", fontsize=3,offset.text=0.1) +
  geom_cladelab(node=223, label="Issatchenkia-ASV36", fontsize=3,offset.text=0.1) +
  geom_cladelab(node=9, label="Issatchenkia-ASV37", fontsize=3,offset.text=0.1) +
  geom_cladelab(node=98, label="Issatchenkia-ASV222", fontsize=3,offset.text=0.1) +
  geom_cladelab(node=102, label="Issatchenkia-ASV129", fontsize=3,offset.text=0.1) +
  labs(color= "Classification") + theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.text = element_text(size=12)) + theme_tree2() 
p2

