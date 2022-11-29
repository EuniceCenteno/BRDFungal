library(afex)
library(lme4)
library(emmeans)
library(lubridate)
library(ggplot2)
library("cowplot")
theme_set(theme_grey())
#install.packages("sjstats")
library(jtools)
library(ggpubr)
library(sjstats)
library(dplyr)
library(devtools) 



rm(list = ls ())

setwd("~/Desktop/eunice/PhD/PhDProject/FungalBRD/Qiime/QiimeFungalComplete/exported/") 

metadata <- read.csv("BRD.metadadaFungal.csv", na.strings = c("","NA"), header=TRUE)
faith_pd <- read.table("faith_pd.tsv", header=TRUE, row.names=1, sep="\t")
evenness <-read.table("evenness.tsv", header=TRUE, row.names=1, sep="\t") 
chao <- read.table("chao1.tsv", header=TRUE, row.names=1, sep="\t")
observed_otus <-read.table("observed_otus.tsv", header=TRUE, row.names=1, sep="\t")
alpha_diversity <- merge(faith_pd, evenness, by.x = 0, by.y = 0)
alpha_diversity <- merge(alpha_diversity, observed_otus, by.x = "Row.names", by.y = 0)
alpha_diversity <- merge(alpha_diversity, chao, by.x = "Row.names", by.y = 0)
metadata2 <- merge(metadata, alpha_diversity, by.x = "ID", by.y = "Row.names")
row.names(metadata2) <- metadata2$ID
metadata2 <- metadata2[,-1]

#assign numerical values to factors
metadata2$BRD <- as.factor(metadata2$BRD)
levels(metadata2$BRD) <- list("Healthy"="Healthy", "BRD"="BRD")
metadata2$PenCode <- as.factor(metadata2$PenCode)

str(metadata)

# Converting date of collection to numeric values
metadata2$date <- metadata2$Date.Collection
metadata2$Date.Collection <- as.Date(metadata2$Date.Collection, "%m/%d/%y")
d<- as.Date('12/31/2020', "%m/%d/%y") #use to calculate the days
metadata2$Date.Collection <- as.Date(d) -as.Date(metadata2$Date.Collection) 
metadata2$Date.Collection <- as.numeric(metadata2$Date.Collection)
str(metadata2$Date.Collection)
#the highest day value is the date of the samples collected first 

plot(metadata2$Date.Collection, metadata2$Age)

## Dependent factors in the model
str(metadata2)
#1. Observed OTUs
#2. Chao1 (measures richness of the environment)
#3. Pielou_e (measures evenness)
#4. Faith_pd (phylogenetic diversity)

set_sum_contrasts() # important for afex

# full model
str(metadata2)
#For dependent variable Observed OTUs
#install.packages("piecewiseSEM")

M1 <- mixed(observed_otus ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1) #not significant

M2 <- mixed(pielou_e ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2) ## BRD and date of collection

M3 <- mixed(chao1 ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3) #not sig

M4 <- mixed(faith_pd ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4) #NS

#checking assumptions
# way 1:
plot(M1$full_model)
plot(M2$full_model)
plot(M3$full_model)
plot(M4$full_model)

# this is for testing the normality of the residuals
qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model))
qqnorm(residuals(M2$full_model))
qqline(residuals(M2$full_model))
qqnorm(residuals(M3$full_model))
qqline(residuals(M3$full_model))
qqnorm(residuals(M4$full_model))
qqline(residuals(M4$full_model))

shapiro.test(metadata2$observed_otus)
shapiro.test(metadata2$faith_pd)
shapiro.test(metadata2$chao1)
shapiro.test(metadata2$pielou_e)
## data not normal

metadata2<- mutate(metadata2, observed_otuslog = log(observed_otus + 1))
metadata2 <- mutate(metadata2, observed_otussqrt = sqrt(observed_otus + 0.5))
metadata2 <- mutate(metadata2, observed_otuscube = observed_otus^(1/3))

M1 <- mixed(observed_otuslog ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1)

M2.SF <- mixed(observed_otussqrt ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2.SF) 

M3 <- mixed(observed_otuscube ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3) 

qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model))
qqnorm(residuals(M2$full_model))
qqline(residuals(M2$full_model))

shapiro.test(metadata2$observed_otus)
shapiro.test(metadata2$observed_otuslog)
shapiro.test(metadata2$observed_otussqrt) #normal
anova(M2.SF)
shapiro.test(metadata2$observed_otuscube) 

hist(metadata2$observed_otuslog)
hist(metadata2$observed_otussqrt)
hist(metadata2$observed_otuscube)

# faith
metadata2<- mutate(metadata2, faith_pdlog = log(faith_pd + 1))
metadata2 <- mutate(metadata2, faith_pdsqrt = sqrt(faith_pd + 0.5))
metadata2 <- mutate(metadata2, faith_pdcube = faith_pd^(1/3))
metadata2 <- mutate(metadata2, faith_log10 = log10(faith_pd))

M1 <- mixed(faith_pdlog ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1) #not significant

M2 <- mixed(faith_pdsqrt ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2) ## BRD and date of collection

M3 <- mixed(faith_pdcube ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3) ## BRD and date of collection

M4 <- mixed(faith_log10 ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4) ## BRD and date of collection


qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model))
qqnorm(residuals(M2$full_model))
qqline(residuals(M2$full_model))

shapiro.test(metadata2$faith_pd) #not normal
shapiro.test(metadata2$faith_pdlog) #not normal
shapiro.test(metadata2$faith_pdsqrt) #not normal
shapiro.test(metadata2$faith_pdcube) #not normal
shapiro.test(metadata2$faith_log10)

hist(metadata2$faith_pd)
hist(metadata2$faith_pdlog)
hist(metadata2$faith_pdsqrt)
hist(metadata2$faith_pdcube)
hist(metadata2$faith_log10)

# chao1
metadata2<- mutate(metadata2, chao1log = log(chao1 + 1))
metadata2 <- mutate(metadata2, chao1sqrt = sqrt(chao1 + 0.5))
metadata2 <- mutate(metadata2, chao1cube = chao1^(1/3))

M1 <- mixed(chao1log ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1) #not significant

M2.S <- mixed(chao1sqrt ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2.S) ## BRD and date of collection

M3 <- mixed(chao1cube ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2) ## BRD and date of collection

qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model))
qqnorm(residuals(M2$full_model))
qqline(residuals(M2$full_model))
qqnorm(residuals(M3$full_model))
qqline(residuals(M3$full_model))

shapiro.test(metadata2$chao1)
shapiro.test(metadata2$chao1log) 
shapiro.test(metadata2$chao1sqrt) #normal
anova(M2)
shapiro.test(metadata2$chao1cube)

hist(metadata2$chao1)
hist(metadata2$chao1log)
hist(metadata2$chao1sqrt)
hist(metadata2$chao1cube)

# pielou_e
metadata2<- mutate(metadata2, pielou_elog = log(pielou_e + 1))
metadata2 <- mutate(metadata2, pielou_esqrt = sqrt(pielou_e + 0.5))
metadata2 <- mutate(metadata2, pielou_pdcube = pielou_e^(1/3))
metadata2 <- mutate(metadata2, pielou_log10 = log10(pielou_e))

M1 <- mixed(pielou_elog ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1) #not significant

M2 <- mixed(pielou_esqrt ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2) ## BRD and date of collection

M3 <- mixed(pielou_pdcube ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3) ## BRD and date of collection

M4 <- mixed(pielou_log10 ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4) ## BRD and date of collection

qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model))
qqnorm(residuals(M2$full_model))
qqline(residuals(M2$full_model))

shapiro.test(metadata2$pielou_e) # not normal
shapiro.test(metadata2$pielou_elog) # not normal
shapiro.test(metadata2$pielou_esqrt) # not normal
shapiro.test(metadata2$pielou_pdcube) # not normal
shapiro.test(metadata2$pielou_log10) # not norma

hist(metadata2$pielou_e)
hist(metadata2$pielou_elog)
hist(metadata2$pielou_esqrt)
hist(metadata2$pielou_pdcube)
hist(metadata2$pielou_log10)

#install.packages("moments")
library(moments)
#choose one of these 3 with the skewness value closer to 0 (= closer to the normal distribution)

#For pielou
skewness(metadata2$pielou_e)
skewness(metadata2$pielou_elog)
skewness(metadata2$pielou_esqrt)
skewness(metadata2$pielou_pdcube)

#For faith
skewness(metadata2$faith_pd) #closer to 0
skewness(metadata2$faith_pdlog)
skewness(metadata2$faith_pdsqrt)
skewness(metadata2$faith_pdcube)

#check the type of distribuiton for each of the measurements
#check for a log normal distribution
qqnorm(log(metadata2$faith_pd))
qqnorm(log(metadata2$pielou_e))

#testing if the data is log normal
install.packages("goft")
library(goft)
lnorm_test(metadata2$faith_pd)
lnorm_test(metadata2$faith_pdlog) # not log normal

lnorm_test(metadata2$pielou_e)
lnorm_test(metadata2$pielou_elog) # not log normal

#### non parametric test for faith and pielo
library(ggpubr)
group_by(metadata2,BRD) %>%
  summarise(
    count = n(),
    median = median(faith_pd, na.rm = TRUE),
    IQR = IQR(faith_pd, na.rm = TRUE))

# loading package for boxplot
ggboxplot(metadata2, x = "BRD", y = "faith_pd",
          color = "BRD", palette = "Dark2",
          ylab = "faith_pd", xlab = "BRD")

res <- wilcox.test(faith_pd~ BRD,
                   data = metadata2,
                   exact = FALSE)
res # no difference in the faith values between BRD and healthy

#now pielou
group_by(metadata2,BRD) %>%
  summarise(
    count = n(),
    median = median(pielou_e, na.rm = TRUE),
    IQR = IQR(pielou_e, na.rm = TRUE))

# loading package for boxplot
ggboxplot(metadata2, x = "BRD", y = "pielou_e",
          color = "BRD", palette = "Dark2",
          ylab = "pielou_e", xlab = "BRD")

res1 <- wilcox.test(pielou_e~ BRD,
                   data = metadata2,
                   exact = FALSE)
res1 # no difference in the pielou values between BRD and healthy


# interpreting results
# BRD plots, we will use Observed_outs and chao1 with the squared transformation
anova(M2.S) #chao1
anova(M2.SF)

library(ggbeeswarm)
afex_plot(M2, x = "BRD", id = "PenCode", dodge = 0.7, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") +
  labs(y = "Evenness (Pielou)", x = "Health Status") +
  theme(legend.position="none") 

### Date collection plots
ggplot(data = metadata2, aes(x = Date.Collection, y = pielou_e)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("Evenness (Pielou)") +xlab ("Date of Collection") +
  theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
  theme(text=element_text(family="Times New Roman"))

## Correlation between average temperature and date of collection

str(metadata2)
cordata = metadata2[,c(3,4)]
corr <- round(cor(cordata), 1)
corr

str(cordata)
cor(cordata$Date.Collection, cordata$Ave.temp)
cor.test(cordata$Date.Collection, cordata$Ave.temp)

library("ggpubr")
ggscatter(cordata, x = "Date.Collection", y = "Ave.temp", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Days relative to the end of the study", ylab = "Average temperature (°F)")


## testing temperature
str(meta)
set_sum_contrasts() 
M1 <- mixed(observed_otus ~ Ave.temp + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1) #NS

M2 <- mixed(pielou_e ~ Ave.temp + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2) #NS

M3 <- mixed(chao1 ~ Ave.temp + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3) #NS

M4 <- mixed(faith_pd ~ Ave.temp + (1|PenCode), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4) #NS

