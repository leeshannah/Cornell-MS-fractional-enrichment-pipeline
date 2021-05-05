#Set notebook figure print parameters
options(repr.plot.width=20, repr.plot.height=80, repr.plot.res = 400, repr.plot.quality = 500, repr.plot.pointsize = 12)

# Set working directory
#setwd('/Users/leesh/Documents/MS analysis pipeline/Nates results 2nd May/')

#install.packages("ggplot2",repos='http://cran.us.r-project.org',dependencies = TRUE)
#install.packages("dplyr",repos='http://cran.us.r-project.org',dependencies = TRUE)
#install.packages("reshape",repos='http://cran.us.r-project.org',dependencies = TRUE)
#install.packages("readxl",repos='http://cran.us.r-project.org',dependencies = TRUE)
#install.packages("tidyr",repos='http://cran.us.r-project.org',dependencies = TRUE)
#install.packages("data.table",repos='http://cran.us.r-project.org',dependencies = TRUE)
#install.packages("reshape2",repos='http://cran.us.r-project.org',dependencies = TRUE)
#install.packages("Polychrome",repos='http://cran.us.r-project.org',dependencies = TRUE)
#install.packages("here",repos='http://cran.us.r-project.org',dependencies = TRUE)
#install.packages("factoextra",repos='http://cran.us.r-project.org',dependencies = TRUE)


# Load packages

library ('ggplot2')
library("dplyr")
library("reshape")
library('readxl')
library('tidyr')
library ("data.table")
library('reshape2')
library("Polychrome")
library('factoextra')
library('here')


## Import data and reshape into tidy format (rows are samples, columns are metabolite intensities)

# Import data from WCM Ms excel spreadsheet
MS_Int <- read_excel('/Users/leesh/Documents/MS analysis pipeline/Nates results 2nd May/KRAS 1wk metabolomics.xlsx')

# Create a new first column: metabolite ID
MS_Int <- tibble::rowid_to_column(MS_Int, "Metabolite_ID")
colnames(MS_Int)[2] <- "Degree_of_confidence"


# Reshape data wide > long 
long_MS_Int <- melt(MS_Int, id.vars = c("Metabolite_ID","Degree_of_confidence", "HMDB", "Compound", "Mass"))                    
colnames(long_MS_Int) <- c("Metabolite_ID","Degree_of_confidence", "HMDB", "Compound", "Mass", "sample", "intensity")              


# Split sample column into 2 columns: condition and replicate number
long_MS_Int <-long_MS_Int %>%
  separate(sample, c("Condition", "replicate"), "_")

# Summary of number of metabolites detected at different confidence levels

MetabolitesSummary <- MS_Int %>%
    group_by(Degree_of_confidence) %>%
    summarize(No.ofCompounds = n())
MetabolitesSummary  <- rbind(MetabolitesSummary , c("Total", colSums(MetabolitesSummary[,2])))

MetabolitesSummary 
# write a csv of this summary
write.csv(MetabolitesSummary, file = "MetabolitesSummary.csv")


long_MS_Int_zeroEntries <- filter (long_MS_Int, intensity =="0")
long_MS_Int_zeroEntries <- long_MS_Int_zeroEntries[order(long_MS_Int_zeroEntries$Compound),]


long_MS_Int_zeroEntries

# write a csv of this summary
write.csv(long_MS_Int_zeroEntries, file = "Zero_Entries.csv")

# Change order of conditons for graphing

#MS_detected$Condition <- factor(MS_detected$Condition, c("CtrlIntact","CtrlSC", "RESTSC", "CDH1SC"))

# Rename condition labels

#MS_detected$Condition <- recode_factor(MS_detected$Condition, CtrlIntact  = "Intact", CtrlSC = "Single cells", RESTSC = "shREST SC", CDH1SC = "shCDH1 SC")



# Filter to just look at high confidence compound assignments

High_Confidence <- long_MS_Int %>%
  filter(Degree_of_confidence == "1")


High_Confidence

MS_Int_HC <- MS_Int %>%
  filter(Degree_of_confidence == "1")


# Reshape data long > wide so that now each row is a sample and each column is a metabolite
wide_HC <- dcast(High_Confidence,  Condition + replicate ~ Metabolite_ID, value.var="intensity")

# Create a new first column: Sample ID
wide_HC <- tibble::rowid_to_column(wide_HC, "SampleID")

# Make first column row names
dataNamed_HC <- wide_HC[,-1]
rownames(dataNamed_HC) <- wide_HC[,1]

dataNamed_HC

# Select just the columns of X data (not qualitative metadata)
Data_X_HC<- dataNamed_HC[3:181]
head (Data_X_HC)

# Compute PCA (scale = TRUE for UV and FALSE for MC)
res.pca_HC <- prcomp(Data_X_HC, scale = TRUE )

# Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
#fviz_eig(res.pca)

options(repr.plot.width=8, repr.plot.height=8, repr.plot.res = 200, repr.plot.quality = 300, repr.plot.pointsize = 12)


fviz_eig(res.pca_HC, addlabels = TRUE)+
labs(title=(expression(paste("R" ^2,": variance explained by each component")))) 


# PCA Scores plot coloured by qualitative variables + 95% confidence elipses 


#groups <- Data_X[9]
fviz_pca_ind(res.pca_HC,
             fill.ind = dataNamed_HC$Condition, col.ind = "black",
             pointshape = 21, pointsize = 7,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             ellipse.level = 0.95,
             legend.title = "Groups",
             mean.point = FALSE,
              ellipse.border.remove = TRUE,
             repel = TRUE
             )+
labs(title="PCA scores of high confidence assigned compounds")

# PCA loadings plot 

# Graph of variables. Positive correlated variables point to the same side of the plot. 
# Negative correlated variables point to opposite sides of the graph.

fviz_pca_var(res.pca_HC,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )+
labs(title="PCA loadings of high confidence assigned compounds")


dataNamed_HC
attach (dataNamed_HC)

# Function to perform a Welch T-test on every column of data (compound) you specify, with Condition used as the groups


Welch_ttest <- function(compound) {
    
    Wttest <- t.test(compound ~ Condition)
     pval = Wttest$p.value
    
    return(pval)
}


# fill in dataNamed[3:1055] such that all columns of metabolite intensity data are calculated (use dataNamed)
# - do not include metadata columns:

Welch_ttest  <- sapply(dataNamed_HC[3:181], 'Welch_ttest')

# Create dataframe of results of the Welch t-test
TtestDF <- as.data.frame(Welch_ttest)

# Make row name metabolite ID column
TtestDF <- setDT(TtestDF, keep.rownames = "Metabolite_ID1")

# Create a new dataframe with a column showing the Benjamini-Hochberg-corrected p-values
BH_testDF <-mutate(TtestDF,BH_correctedVal = p.adjust(TtestDF$Welch_ttest,method="BH") )

# Bind metadata to this dataframe (select which columns from MS_Int are metadata)
MS_Int_HC_meta <- MS_Int_HC %>% select (1:5)
Ttest_meta_HC<-cbind(MS_Int_HC_meta, BH_testDF)  

# Order the dataframe according to the adjusted pvalues
Ttest_meta_HC_orderedByPval <- Ttest_meta_HC[order(Ttest_meta_HC$Welch_ttest),]

Ttest_meta_HC_orderedByPval

# write a csv of this dataframe
write.csv(Ttest_meta_HC_orderedByPval, file = "Welch_BH_orderedResults_HighConfidenceCompoundsOnly.csv")

# Create a DF subselection of compounds with uncorrected Welch t-test pval of ≤ 0.05
Ttest_Signif_HC <- Ttest_meta_HC_orderedByPval %>%
filter (Welch_ttest <= 0.05)

# Select only the compounds with uncorrected Welch t-test pval of ≤ 0.05
Significant_Compounds_HC <- MS_Int_HC[MS_Int_HC$Metabolite_ID %in% Ttest_Signif_HC$Metabolite_ID,]

# Append the t-test results
Significant_Compounds_tResults_HC <- select (Ttest_meta_HC_orderedByPval,1,7,8)
Significant_Compounds_HC=merge(Significant_Compounds_HC,Significant_Compounds_tResults_HC,by="Metabolite_ID")

# Order the dataframe according to the adjusted pvalues
Significant_Compounds_HC <- Significant_Compounds_HC[order(Significant_Compounds_HC$Welch_ttest),]
Significant_Compounds_HC <- within(Significant_Compounds_HC, Compound <- reorder(Compound, Welch_ttest))# orders the facets in increasing order of Welch t-test


# write a csv of this dataframe
write.csv(Significant_Compounds_HC, file = "High_Confid_assign_Signif_Compounds_uncorrectedPval.csv")

# Reshape data wide > long (for ggplot graphs with facet)
Significant_Compounds_HC_long <- melt(Significant_Compounds_HC, id.vars = c("Metabolite_ID","Degree_of_confidence", "HMDB", "Compound", "Mass",  "Welch_ttest","BH_correctedVal"))                    

colnames(Significant_Compounds_HC_long) <- c("Metabolite_ID","Degree_of_confidence", "HMDB", "Compound", "Mass", "Welch_ttest","BH_correctedVal","sample", "intensity")              

# Split sample column into 2 columns: condition and replicate number
Significant_Compounds_HC_long <-Significant_Compounds_HC_long %>%
 separate(sample, c("Condition", "replicate"), "_")


# Pool sizes for compounds with uncorrected Welch t-test pval of ≤ 0.05

options(repr.plot.width=15, repr.plot.height=20, repr.plot.res = 400, repr.plot.quality = 500, repr.plot.pointsize = 12)

p<-ggplot(data=Significant_Compounds_HC_long, aes(x=Condition, y=intensity, fill = Condition)) +
  geom_bar(position = 'dodge', stat = 'summary', fun = 'mean')+
geom_errorbar(stat="summary", 
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))}, position=position_dodge(.9), width = 0.3, colour="black")+

geom_point(aes(fill=Condition), size=2, shape=21, colour="black", stroke = 1, alpha =1, 
           position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0)) +


labs(title="Pool sizes for compounds (assigned with high confidence)
with significant uncorrected Welch t-test pval (0.05)
Mean ± SEM",x="Condition", y = "intensity")+theme_classic()+
theme(axis.text.x=element_text(angle=60, hjust=1))+
theme(legend.position="none")+
facet_wrap(~Compound, ncol =6, scales="free",labeller = label_wrap_gen())+
theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0))

p

pdf(file = "HighConfidence_poolSize_UncorrSignif.pdf", width = 20, height = 30, family = "Helvetica") 

p<-ggplot(data=Significant_Compounds_HC_long, aes(x=Condition, y=intensity, fill = Condition)) +
  geom_bar(position = 'dodge', stat = 'summary', fun = 'mean')+
geom_errorbar(stat="summary", 
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))}, position=position_dodge(.9), width = 0.3, colour="black")+

geom_point(aes(fill=Condition), size=2, shape=21, colour="black", stroke = 1, alpha =1, 
           position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0)) +


labs(title="Pool sizes for compounds (assigned with high confidence)
with significant uncorrected Welch t-test pval (0.05)
Mean ± SEM",x="Condition", y = "intensity")+theme_classic()+
theme(axis.text.x=element_text(angle=60, hjust=1))+
theme(legend.position="none")+
facet_wrap(~Compound, ncol =6, scales="free",labeller = label_wrap_gen())+
theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0))

p
dev.off()

# Reshape data long > wide so that now each row is a sample and each column is a metabolite
wide_All <- dcast(long_MS_Int,  Condition + replicate ~ Metabolite_ID, value.var="intensity")

# Create a new first column: Sample ID
wide_All <- tibble::rowid_to_column(wide_All, "SampleID")

# Make first column row names
dataNamed <- wide_All[,-1]
rownames(dataNamed) <- wide_All[,1]

dataNamed

# Select just the columns of X data (not qualitative metadata)
Data_X<- dataNamed[3:1055]
head (Data_X)

# Compute PCA (scale = TRUE for UV and FALSE for MC)
res.pca <- prcomp(Data_X, scale = TRUE )

# Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
#fviz_eig(res.pca)

options(repr.plot.width=8, repr.plot.height=8, repr.plot.res = 200, repr.plot.quality = 300, repr.plot.pointsize = 12)


fviz_eig(res.pca, addlabels = TRUE)+
labs(title=(expression(paste("R" ^2,": variance explained by each component")))) 


# PCA Scores plot coloured by qualitative variables + 95% confidence elipses 


#groups <- Data_X[9]
fviz_pca_ind(res.pca,
             fill.ind = dataNamed$Condition, col.ind = "black",
             pointshape = 21, pointsize = 7,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             ellipse.level = 0.95,
             legend.title = "Groups",
             mean.point = FALSE,
              ellipse.border.remove = TRUE,
             repel = TRUE
             )

# PCA loadings plot 

# Graph of variables. Positive correlated variables point to the same side of the plot. 
# Negative correlated variables point to opposite sides of the graph.

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )


dataNamed
attach (dataNamed)

# Function to perform a Welch T-test on every column of data (compound) you specify, with Condition used as the groups


Welch_ttest <- function(compound) {
    
    Wttest <- t.test(compound ~ Condition)
     pval = Wttest$p.value
    
    return(pval)
}


# fill in dataNamed[3:1055] such that all columns of metabolite intensity data are calculated (use dataNamed)
# - do not include metadata columns:

Welch_ttest  <- sapply(dataNamed[3:1411], 'Welch_ttest')

# Create dataframe of results of the Welch t-test
TtestDF <- as.data.frame(Welch_ttest)

# Make row name metabolite ID column
TtestDF <- setDT(TtestDF, keep.rownames = "Metabolite_ID1")

# Create a new dataframe with a column showing the Benjamini-Hochberg-corrected p-values
BH_testDF <-mutate(TtestDF,BH_correctedVal = p.adjust(TtestDF$Welch_ttest,method="BH") )

# Bind metadata to this dataframe (select which columns from MS_Int are metadata)
MS_Int_meta <- MS_Int %>% select (1:5)
Ttest_meta<-cbind(MS_Int_meta, BH_testDF)  
 
# Order the dataframe according to the adjusted pvalues
Ttest_meta_orderedByPval <- Ttest_meta[order(Ttest_meta$Welch_ttest),]

Ttest_meta_orderedByPval

# write a csv of this dataframe
write.csv(Ttest_meta_orderedByPval, file = "Welch_BH_orderedResults.csv")

# Create a DF subselection of compounds with uncorrected Welch t-test pval of ≤ 0.05
Ttest_Signif <- Ttest_meta_orderedByPval %>%
filter (Welch_ttest <= 0.05)

# Select only the compounds with uncorrected Welch t-test pval of ≤ 0.05
Significant_Compounds <- MS_Int[MS_Int$Metabolite_ID %in% Ttest_Signif$Metabolite_ID,]

# Create a new column to label grpahs: metabolite Id with name appended 
Significant_Compounds$Name_ID <- paste(Significant_Compounds$Compound, "-", Significant_Compounds$Metabolite_ID)

# Append the t-test results
Significant_Compounds_tResults <- select (Ttest_meta_orderedByPval,1,7,8)
Significant_Compounds=merge(Significant_Compounds,Significant_Compounds_tResults,by="Metabolite_ID")

# create a new column for graphnames whereby compound name is used for high confidence metabolites
# metabolite ID is used for lower-confidence assignments

Significant_Compounds$GraphName <- Significant_Compounds$Compound
Significant_Compounds$GraphName[Significant_Compounds$Degree_of_confidence==2] <- Significant_Compounds$Metabolite_ID[Significant_Compounds$Degree_of_confidence==2]
Significant_Compounds$GraphName[Significant_Compounds$Degree_of_confidence==3] <- Significant_Compounds$Metabolite_ID[Significant_Compounds$Degree_of_confidence==3]


# Order the dataframe according to the adjusted pvalues
Significant_Compounds <- Significant_Compounds[order(Significant_Compounds$Welch_ttest),]
Significant_Compounds <- within(Significant_Compounds,GraphName	 <- reorder(GraphName, Welch_ttest))# orders the facets in increasing order of Welch t-test


# write a csv of this dataframe
write.csv(Significant_Compounds, file = "Significant_Compounds_uncorrectedPval.csv")

# Reshape data wide > long (for ggplot graphs with facet)
Significant_Compounds_long <- melt(Significant_Compounds, id.vars = c("Metabolite_ID","Degree_of_confidence", "HMDB", "Compound", "Mass", "Name_ID", "Welch_ttest","BH_correctedVal","GraphName"))                    

colnames(Significant_Compounds_long) <- c("Metabolite_ID","Degree_of_confidence", "HMDB", "Compound", "Mass","Name_ID", "Welch_ttest","BH_correctedVal","GraphName","sample", "intensity")              

# Split sample column into 2 columns: condition and replicate number
Significant_Compounds_long <-Significant_Compounds_long %>%
 separate(sample, c("Condition", "replicate"), "_")


# Pool sizes for compounds with uncorrected Welch t-test pval of ≤ 0.05

options(repr.plot.width=15, repr.plot.height=100, repr.plot.res = 400, repr.plot.quality = 500, repr.plot.pointsize = 12)

p<-ggplot(data=Significant_Compounds_long, aes(x=Condition, y=intensity, fill = Condition)) +
  geom_bar(position = 'dodge', stat = 'summary', fun = 'mean')+
geom_errorbar(stat="summary", 
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))}, position=position_dodge(.9), width = 0.3, colour="black")+

geom_point(aes(fill=Condition), size=2, shape=21, colour="black", stroke = 1, alpha =1, 
           position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0)) +


labs(title="Pool sizes for compounds with significant uncorrected Welch t-test pval (0.05)
Mean ± SEM",x="Condition", y = "intensity")+theme_classic()+
theme(axis.text.x=element_text(angle=60, hjust=1))+
theme(legend.position="none")+
facet_wrap(~GraphName, ncol =6, scales="free",labeller = label_wrap_gen())+
theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0))


p

pdf(file = "poolSize_UncorrSignif.pdf", width = 20, height = 70, family = "Helvetica") 

p<-ggplot(data=Significant_Compounds_long, aes(x=Condition, y=intensity, fill = Condition)) +
  geom_bar(position = 'dodge', stat = 'summary', fun = 'mean')+
geom_errorbar(stat="summary", 
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))}, position=position_dodge(.9), width = 0.3, colour="black")+

geom_point(aes(fill=Condition), size=2, shape=21, colour="black", stroke = 1, alpha =1, 
           position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0)) +


labs(title="Pool sizes for compounds with significant uncorrected Welch t-test pval (0.05)
Mean ± SEM",x="Condition", y = "intensity")+theme_classic()+
theme(axis.text.x=element_text(angle=60, hjust=1))+
theme(legend.position="none")+
facet_wrap(~GraphName, ncol =6, scales="free",labeller = label_wrap_gen())+
theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0))


p
dev.off()

















