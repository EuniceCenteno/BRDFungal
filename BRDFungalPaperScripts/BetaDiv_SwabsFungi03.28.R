BETA diversity Swabs
#install.packages("raster")
library(vegan) 
library(ggplot2)
library(ggpubr)
library(data.table)
library(phyloseq)
library(qiime2R)
library(tidyr)
library(naniar)
library(raster)
library(dplyr)
#transpose function


setwd("~/Desktop/eunice/PhD/PhDProject/FungalBRD/Qiime/QiimeFungalComplete/exported/") #sets new working directory for Windows systems (remember to replace … with your filepath)
metadata <-read.csv("BRD.metadadaFungal.csv", na.strings = c("","NA"), header=TRUE)
#OTU table (shared file)
#The OTU table as exported from qiime has a pound sign before the header row. You need to delete that pound sign in a text editor.
str(metadata)
metadata$BRD <- factor(metadata$BRD) 
metadata$PenCode <- factor(metadata$PenCode)
metadata$seasonGroup <- factor(metadata$seasonGroup, levels=c('Jul-Aug', 'Sept-Oct', 'Nov-Dec'))
levels(metadata$seasonGroup)
order_groups <- metadata$ID
row.names(metadata) = metadata[,1]
metadata %>% tally()
metadata %>% count(seasonGroup)
#seasonGroup  n
#Season1 54
#Season2 38
#Season3 39

# Converting date of collection to numeric values

#OTU table, we use the rarified table
ASVs <- read_qza("rarefied_tableFungi.qza") #6930 ASVs
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) 
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #the key withe the names
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
ASV_table <- t(ASV_table)
#write.csv(ASV_table, "ASV_table")

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
#write.table(TaxASV,"TaxASV.txt",sep=",", row.names = FALSE)

### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq #[ 6930 taxa and 126 samples ]

colnames(tax_table(physeq_deseq))
colnames(sample_data(metadata))

## Calculating distances based on Bray-curtis and Weighted Unifrac using the physeq object
#subsetting the metadata with the same number of samples as the ASV_Table
names <- rownames(ASV_table)
metadata = metadata[row.names(metadata) %in% c(names),]
#making the bray curtis distance
dist.bray <- phyloseq::distance(physeq_deseq, method = "bray")
dist.bray <- as.dist(dist.bray)

## PERMANOVA
#Bray curtis
BC <- adonis(dist.bray ~ metadata$BRD + metadata$seasonGroup, permutations = 999)
BC #not different by BRD by time is different 

TimeBC <- pairwise.adonis2(dist.bray ~ seasonGroup,  data=metadata)
TimeBC

#Bray_curtis dispersion
BRD_BC <- betadisper(dist.bray, type = c("centroid"), group = metadata$seasonGroup)
BRD_BC
boxplot(BRD_BC)
TukeyHSD(BRD_BC) 
plot(BRD_BC)

pBRD_BC<- permutest(BRD_BC, permutations = 999)
pBRD_BC# not significant 

#Bray_curtis dispersion
BRD_B <- betadisper(dist.bray, type = c("centroid"), group = metadata$BRD)
BRD_B
boxplot(BRD_B)
TukeyHSD(BRD_B) 
plot(BRD_B)

pBRD_B<- permutest(BRD_B, permutations = 999)
pBRD_B# not significant 


## Bray_curtis
ordu.bc <- ordinate(physeq_deseq, "PCoA", "bray")
Bray <- plot_ordination(physeq_deseq, 
                        ordu.bc, color="seasonGroup") 
Bray1 <- Bray[["data"]][["Axis.1"]]
Bray1 <- as.data.frame(Bray1)
Bray1$number <- rownames(Bray1)
Bray2 <- Bray[["data"]][["Axis.2"]]
Bray2 <- as.data.frame(Bray2)
Bray2$number <- rownames(Bray2)
Bray3 <- Bray[["data"]][["ID"]]
Bray3 <- as.data.frame(Bray3)
Bray3$number <- rownames(Bray3)
Bray4 <- Bray[["data"]][["seasonGroup"]]
Bray4 <- as.data.frame(Bray4)
Bray4$number <- rownames(Bray4)

bray <- merge(Bray3, Bray4, by.x = "number", by.y = "number")
colnames(bray) <- c("number", "ID", "seasonGroup")
bray <- merge(bray, Bray1, by.x = "number", by.y = "number")
bray <- merge(bray, Bray2, by.x = "number", by.y = "number")

my_colors <- c(
  "gray12","#fb9a99","steelblue2",'palegreen3','#e31a1c',
  '#ff7f00','#cab3d6','#6a3d9a','gray','#b15938', 
  "#CBD588", "#5F7FC7","#DA5734", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14385", "#653936", "#C84348", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

Bray <- Bray + ggtitle("Bray Curtis") + geom_point(size = 2)
Bray <- Bray + theme_classic() +
  ggtitle("Bray-Curtis") +
  scale_color_manual(values = c(my_colors)) +
  theme(legend.text = element_text(size=12)) +
  labs(color='SeasonGroup') +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) 
print(Bray + stat_ellipse(level=0.95))

##other plot
str(bray)

centroids <- aggregate(bray[,4:5], list(Group=bray$seasonGroup), mean)
colnames(centroids) <- c('seasonGroup','groupX', 'groupY')

bray <- merge(bray, centroids, by.x = "seasonGroup", by.y = "seasonGroup")

#distance betweeen centroids
pointDistance(c(0.08405844, 0.059075765), c(-0.07284415, -0.008914383), lonlat=FALSE) #Jul-Aug - Nov-Dec
#0.1710002
pointDistance(c(0.08405844, 0.059075765), c(-0.04105142, -0.072273446), lonlat=FALSE) #Jul-Aug - Sept-Oct
#0.1813976
pointDistance(c(-0.07284415, -0.008914383), c(-0.04105142, -0.072273446), lonlat=FALSE) #Sept-Oct - Nov-Dec
#0.07088828

str(bray)
ggplot(bray, aes(x=Bray1, y=Bray2, color=seasonGroup)) + 
  geom_point(size=2) + 
  scale_color_manual(values = c(my_colors)) +
  scale_fill_manual(values = c(my_colors)) +
  theme_classic() + stat_ellipse() +
  guides(size=FALSE) +
  guides(fill=FALSE) +
  geom_point(data= bray, aes(x=groupX, y=groupX, shape=seasonGroup, size=5, fill=seasonGroup)) +
  scale_shape_manual(values=c(25, 22, 23))+
  ggtitle("Bray-Curtis") +
  labs(color= "Season") +
  labs(shape= "Centroids") +
  #geom_segment(data= wuaxis, aes(x=wuaxis1, y=wuaxis2, xend=groupX, yend=groupY, color= wuaxis4), size = .05) +
  labs(x='Axis 1 (5.5%)', y= 'Axis 2  (4%)') +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) 


## Weighted Unifrac
# We need to create first a tree using OTU and taxa table-- we do this by creating a phyloseq object 

##Making the tree
rm(OTU_Unifrac)
library("ape")
random_tree = rtree(ntaxa(physeq_deseq), rooted=TRUE, tip.label=taxa_names(physeq_deseq))
##Merging the tree with the phyloseq object
physeq1 = merge_phyloseq(physeq_deseq, random_tree)
physeq1 #[ 6930 taxa and 126 samples ]

OTU_Unifrac <- UniFrac(physeq1, weighted = TRUE) ## Calculating the weighted unifrac distances

## Permanova 
M3 <- adonis(OTU_Unifrac~ metadata$BRD + metadata$seasonGroup, permutations = 999)
M3 #not significant by BRD and season

TimeWU <- pairwise.adonis2(OTU_Unifrac ~ seasonGroup,  data=metadata)
TimeWU

## Dispersion test 
BRD_Wu <- betadisper(OTU_Unifrac, type = c("centroid"), group = metadata$seasonGroup)
BRD_Wu
boxplot(BRD_Wu)
TukeyHSD(BRD_Wu) ## calculate the distance from the centroids to one group to another
plot(BRD_Wu) #no difference based on betadisper

pBRD_WU<- permutest(BRD_Wu, permutations = 999)
pBRD_WU# not significant 

## Dispersion testBRD
BRD_W <- betadisper(OTU_Unifrac, type = c("centroid"), group = metadata$BRD)
BRD_W
boxplot(BRD_W)
TukeyHSD(BRD_W) ## calculate the distance from the centroids to one group to another
plot(BRD_W) #no difference based on betadisper

pBRD_W<- permutest(BRD_W, permutations = 999)
pBRD_W# not significant
## Weighted unifrac 
ordu.wt.uni <- ordinate(physeq1 , "PCoA", "unifrac", weighted=T)
wt.unifrac <- plot_ordination(physeq1, 
                              ordu.wt.uni, color="season") 
wuaxis1 <- wt.unifrac[["data"]][["Axis.1"]]
wuaxis1 <- as.data.frame(wuaxis1)
wuaxis1$number <- rownames(wuaxis1)
wuaxis2 <- wt.unifrac[["data"]][["Axis.2"]]
wuaxis2 <- as.data.frame(wuaxis2)
wuaxis2$number <- rownames(wuaxis2)
wuaxis3 <- wt.unifrac[["data"]][["ID"]]
wuaxis3 <- as.data.frame(wuaxis3)
wuaxis3$number <- rownames(wuaxis3)
wuaxis4 <- wt.unifrac[["data"]][["season"]]
wuaxis4 <- as.data.frame(wuaxis4)
wuaxis4$number <- rownames(wuaxis4)

wuaxis <- merge(wuaxis3, wuaxis4, by.x = "number", by.y = "number")
colnames(wuaxis) <- c("number", "ID", "season")
wuaxis <- merge(wuaxis, wuaxis1, by.x = "number", by.y = "number")
wuaxis <- merge(wuaxis, wuaxis2, by.x = "number", by.y = "number")

wt.unifrac <- wt.unifrac + ggtitle("Weighted UniFrac") + geom_point(size = 2)
wt.unifrac <- wt.unifrac + theme_classic() +
  ggtitle("Weighted UniFrac") +
  scale_color_manual(values = c(my_colors)) +
  theme(legend.text = element_text(size=12)) +
  labs(color='Season') +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) 

print(wt.unifrac + stat_ellipse())#


# Pairwise Adonis funtion by edro Martinez Arbizu & Sylvain Monteux
#https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 






