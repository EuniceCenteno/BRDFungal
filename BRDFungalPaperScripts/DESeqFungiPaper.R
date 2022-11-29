###Differential Abundance with DESeq2

Adapted from https://joey711.github.io/phyloseq-extensions/DESeq2.html
setwd("~/Desktop/eunice/PhD/PhDProject/FungalBRD/Qiime/QiimeFungalComplete/exported/")

rm(list = ls ())

library("DESeq2")
library(dplyr)
library(tidyr)
library(ape)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(plotly)
library(tidyr)
library(naniar)
library(zoo)
library(lubridate)
library(qiime2R)

#OTU table (shared file)
#The OTU table as exported from qiime has a pound sign before the header row. You need to delete that pound sign in a text editor.
metadata <- read.csv("BRD.metadadaFungal.csv")
#metadata <- metadata2[-1,]
str(metadata)
metadata$BRD <- factor(metadata$BRD) 
metadata$PenCode <- factor(metadata$PenCode)
order_groups <- metadata$ID
row.names(metadata) = metadata[,1]
metadata = metadata[,-1] # we remove sample 92 because no sequences were identified in the analysis
metadata = metadata[-125,]
metadata$Status <-metadata$BRD 
metadata$Status <- as.factor(metadata$Status)
levels(metadata$Status) <- list("0"="Healthy", "1"="BRD")

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
row.names(TaxASV) <- TaxASV[,10]
TaxASV = TaxASV[,-c(1,10)]

#write.csv(TaxASV,"taxonomy.csv", row.names = TRUE)

### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq #[ 7050 taxa and 129 samples ]

colnames(tax_table(physeq_deseq))

## Running DESeq
## creating a new phyloseq object
#### DESEq differentially```
#To use DESeq, we need no zeros in our OTU table. So we will edit the table + 1

NewASVTable2 <- ASV_table + 1

### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(NewASVTable2), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq #[ 7050 taxa and 129 samples ]

str(metadata)
levels(metadata$Status) 

# establishing the model
diagdds = phyloseq_to_deseq2(physeq_deseq, ~ Status)
#PenCode needs to be factor
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
head(diagdds)
resultsNames(diagdds)

my_contrast = c("Status", "1", "0")
res = results(diagdds, contrast = my_contrast, cooksCutoff = FALSE, alpha=0.05)
summary(res)
res
res <- as.data.frame(res)
res2 <- merge(res, TaxASV, by.x = 0, by.y = 0)
alpha = 0.05
#sigtab = res ### No significant results
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_deseq)[rownames(sigtab), ], "matrix"))
sigtab

sigtab$High_low <- ifelse(
  sigtab$log2FoldChange < -2.00, 'High in Healthy',
  ifelse(sigtab$log2FoldChange > 2.00, 'High in BRD',
         'Mid Change'))
write.table(sigtab,"sigtab2.2.txt",sep=",", row.names = TRUE)

#To manke the figures
DeSeq = read.csv("sigtabFungi.csv")
str(DeSeq)
DeSeq$High_low <- factor(DeSeq$High_low)

ggplot(data = DeSeq,aes(x = reorder(SpeciesASV, log2FoldChange), y = log2FoldChange, group = factor(High_low))) + coord_flip() +
  geom_bar(stat = "identity", aes(fill = factor(High_low)), position = position_dodge(width = 0.9)) +
  labs(fill= "Diagnosis") +
  theme_bw()+
  ylab("Log2 Fold Change") +xlab ("Differentially abundant ASVs") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12, face="italic")) 
  
