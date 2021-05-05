#Set notebook figure print parameters
options(repr.plot.width=20, repr.plot.height=80, repr.plot.res = 400, repr.plot.quality = 500, repr.plot.pointsize = 12)

# Set working directory
# setwd('/Users/leesh/Documents/MS analysis pipeline/Rachels results/')

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

# Import data from WCM Ms excel spreadsheet
MS_Int <- read_excel('/Users/leesh/Documents/MS analysis pipeline/Rachels results/RC755_IntensityReport.xlsx')

# Transpose and correct row and column names
t_MS_Int<- transpose(MS_Int)
colnames(t_MS_Int) <- rownames(MS_Int)
rownames(t_MS_Int) <- colnames(MS_Int)

# Make first row column names
names(t_MS_Int) <- as.matrix(t_MS_Int[1, ])
t_MS_Int <- t_MS_Int[-1, ]
t_MS_Int[] <- lapply(t_MS_Int, function(x) type.convert(as.character(x)))

t_MS_Int                            

# Create column name for sample labels
setDT(t_MS_Int, keep.rownames = "Condition_Replicate")
                     
# Split first column into 2 columns: condition and replicate number
t_MS_Int <-t_MS_Int %>%
  separate(Condition_Replicate, c("Condition", "replicate"), "_")

# Create a new first column: Sample ID
t_MS_Int <- tibble::rowid_to_column(t_MS_Int, "SampleID")

# Reshape data wide > long to be able to separate compound and labelling M+
long_MS_Int <- melt(t_MS_Int, id.vars = c("SampleID", "Condition", "replicate"))                    
colnames(long_MS_Int) <- c("SampleID", "Condition", "replicate", "compound_labeling", "intensity")

# Split compound labelling column into compound and labelling M+
long_MS_Int <-long_MS_Int %>%
  separate(compound_labeling, c("compound", "labeling"), "_")

# Reshape data long > wide so that now each column is a compound 
wide_MS_Int <- dcast(long_MS_Int, SampleID + Condition +replicate +labeling ~ compound, value.var="intensity")
                                        
# Rename compounds to clean odd characters that can't be parsed e.g. L-Serine
colnames(wide_MS_Int) <- make.names(colnames(wide_MS_Int))

# Remove metabolites that have not been detected
MS_detected <- wide_MS_Int [, colSums(wide_MS_Int != 0, na.rm = TRUE) > 0]   

MS_detected

# Change order of conditons for graphing

#MS_detected$Condition <- factor(MS_detected$Condition, c("CtrlIntact","CtrlSC", "RESTSC", "CDH1SC"))

# Rename condition labels

#MS_detected$Condition <- recode_factor(MS_detected$Condition, CtrlIntact  = "Intact", CtrlSC = "Single cells", RESTSC = "shREST SC", CDH1SC = "shCDH1 SC")



# Sum each sample (sampleID) M+ = total pool size intensity per sample
# Adds a row of these totals at bottom of table

MS_detected_Totals <-MS_detected %>%
group_by(SampleID, Condition, replicate) %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum,na.rm = TRUE),
                      across(where(is.character), ~"Total")))

# reshape data wide > long to be able to compare labelling
long_All <- melt(MS_detected_Totals, id.vars = c("SampleID", "Condition", "replicate", "labeling"))                    
colnames(long_All) <- c("SampleID", "Condition", "replicate", "labeling", "metabolite", "intensity")

# Reshape data long > wide so that now each column is 0, 1, 2...., Total labeling

wide_All <- dcast(long_All, SampleID + Condition + replicate + metabolite ~ labeling, value.var="intensity")

# Have to rename "0" column as R gets confused with columns labelled as numbers

colnames(wide_All)[colnames(wide_All) == "0"] <- "zero"

# Select only rows where total > 0 labeling = metabolites which show fractional enrichment

OnlyFE_metabolites <- wide_All %>%
filter (Total > zero)

# Create a summary table showing which metabolites are enriched and in how many samples. Save as a csv

OnlyFE_metabolitesSummary <- OnlyFE_metabolites %>%
    group_by(metabolite) %>%
    summarize(No.SamplesWithEnrichment = n())
OnlyFE_metabolitesSummary 

write.csv(OnlyFE_metabolitesSummary, file = "OnlyFE_metabolitesSummary.csv")

# Create a summary table for how many metabolites were measured, detected and enriched
summary <- data.frame (Number_of_metabolites  = c("Measured", "Detected", "Enriched"),
                  Count = c(ncol(wide_MS_Int)-4, ncol(MS_detected)-4,nrow(OnlyFE_metabolitesSummary))
                  )
summary


# Select only the metabolites that are enriched from MS_detected

Selected_FE <- MS_detected %>% select(one_of(dput(as.character(OnlyFE_metabolitesSummary$metabolite))))

# Bind metadata (select which columns from MS_detected are metadata)

MS_detected_meta <- MS_detected %>% select (1:4)

MS_FE <-cbind(MS_detected_meta, Selected_FE)  


MS_FE

set.seed(723451) # for reproducibility
P30_2 <- createPalette(30, c("#00ffff", "#ff00ff", "#ffff00"), M=1000)
#swatch(P30_2)
names(P30_2) <- NULL
P30_2

StackPalette_M1 <- c( '#EE2EE8','#F0E716','#FF0D45','#78454B','#2E51FB','#3DFD66','#E87600','#D10D7E','#C7A7FE','#689F60','#FCDAAF','#224F63','#D4DDEB','#FFA3D6','#1CA9FB','#B61CFC','#7A267D','#47F8B5','#E2767A','#651CB2','#97D61C','#EE515A','#E9A54F','#4FA2C9','#FF91FB','#2EA797','#CFDD7D','#66633D','#2A4D95')
#swatch(StackPalette_M1)

# Long version of MS_FE for facet-graphs

Long_MS_FE <- melt(MS_FE, id.vars = c("SampleID", "Condition", "replicate", "labeling"))                    
colnames(Long_MS_FE) <- c("SampleID", "Condition", "replicate", "labeling", "metabolite", "intensity")

Long_MS_FE$labeling <- factor(Long_MS_FE$labeling,levels= c("0",'1',	'2',	'3',	'4',	'5',	'6',	'7',	'8',	'9',	'10',	'11',	'12',	'13',	'14',	'15',	'16',	'17',	'18',	'19',	'20',	'21',	'22',	'23',	'24',	'25',	'26',	'27'))


# Absolute values: Stacked bar plot with y scale fixed across metabolites
ggplot(Long_MS_FE, aes(fill= labeling, y=intensity, x=Condition)) + 
 geom_bar(aes(fill = labeling), position = position_stack(reverse = TRUE), stat="identity") +

labs(title="Absolute values of C13 labelling
(Mean)",subtitle = "y scale fixed across metabolites",x="", y = "Absolute intensity")+theme_classic()+

theme(axis.text.x=element_text(size=rel(3), angle=90))+
scale_fill_manual(values = P30_2)+
facet_wrap(~metabolite, ncol =5)+
theme(strip.text.x = element_text(size = 18, colour = "black", angle = 0))



# Absolute values: Stacked bar plot with y scale adjusted across metabolites

ggplot(Long_MS_FE, aes(fill= labeling, y=intensity, x=Condition)) + 
 geom_bar(aes(fill = labeling), position = position_stack(reverse = TRUE), stat="identity") +
labs(title="Absolute values of C13 labelling
(Mean)",subtitle = "y scale relative to metabolite",x="", y = "Absolute intensity")+theme_classic()+
theme(axis.text.x=element_text(size=rel(2), angle=0))+
scale_fill_manual(values = P30_2)+
facet_wrap(~metabolite, ncol =5, scales="free")+
theme(strip.text.x = element_text(size = 18, colour = "black", angle = 0))

# Normalized values: Stacked bar plot with y scale fixed across metabolites
ggplot(Long_MS_FE, aes(fill= labeling, y=intensity, x=Condition)) + 
 geom_bar(aes(fill = labeling), position = position_fill(reverse = TRUE), stat="identity") +
labs(title="Relative values of C13 labelling
(Mean)",x="", y = "Absolute intensity")+theme_classic()+
theme(axis.text.x=element_text(size=rel(3), angle=90))+
scale_fill_manual(values = P30_2)+
facet_wrap(~metabolite, ncol =5)+
theme(strip.text.x = element_text(size = 18, colour = "black", angle = 0))

pdf(file = "AllFE.pdf", width = 20, height = 50, family = "Helvetica") 

# Absolute values: Stacked bar plot with y scale fixed across metabolites
ggplot(Long_MS_FE, aes(fill= labeling, color = labeling, y=intensity, x=Condition)) + 
 geom_bar(aes(fill = labeling), position = position_stack(reverse = TRUE), stat="identity") +
labs(title="Absolute values of C13 labelling
(Mean)",subtitle = "y scale fixed across metabolites", x="", y = "Absolute intensity")+theme_classic()+
theme(axis.text.x=element_text(size=rel(3), angle=90))+
scale_fill_manual(values = P30_2)+
scale_color_manual(values = P30_2)+
facet_wrap(~metabolite, ncol =5)+
theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0))


# Absolute values: Stacked bar plot with y scale adjusted across metabolites

ggplot(Long_MS_FE, aes(fill= labeling, color = labeling, y=intensity, x=Condition)) + 
 geom_bar(aes(fill = labeling), position = position_stack(reverse = TRUE), stat="identity") +
labs(title="Absolute values of C13 labelling
(Mean)",subtitle = "y scale relative to metabolite",x="", y = "Absolute intensity")+theme_classic()+
theme(axis.text.x=element_text(size=rel(2), angle=0))+
scale_fill_manual(values = P30_2)+
scale_color_manual(values = P30_2)+
facet_wrap(~metabolite, ncol =5, scales="free")+
theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0))



# Normalized values: Stacked bar plot with y scale fixed across metabolites
ggplot(Long_MS_FE, aes(fill= labeling, color = labeling, y=intensity, x=Condition)) + 
 geom_bar(aes(fill = labeling), position = position_fill(reverse = TRUE), stat="identity") +
labs(title="Relative values of C13 labelling
(Mean)",x="", y = "Absolute intensity")+theme_classic()+
theme(axis.text.x=element_text(size=rel(3), angle=90))+
scale_fill_manual(values = P30_2)+
scale_color_manual(values = P30_2)+
facet_wrap(~metabolite, ncol =5)+
theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0))

dev.off()


MS_FE$labeling <- as.numeric(MS_FE$labeling)

MS_FE_M_Higher <- MS_FE %>%
filter(!((labeling == 0)))

MS_FE_M_Higher$labeling <- as.character(MS_FE_M_Higher$labeling)

# Long version of MS_FE for facet-graphs

Long_MS_FE_M_Higher <- melt(MS_FE_M_Higher, id.vars = c("SampleID", "Condition", "replicate", "labeling"))                    
colnames(Long_MS_FE_M_Higher) <- c("SampleID", "Condition", "replicate", "labeling", "metabolite", "intensity")

Long_MS_FE_M_Higher$labeling <- factor(Long_MS_FE_M_Higher$labeling,levels= c('1',	'2',	'3',	'4',	'5',	'6',	'7',	'8',	'9',	'10',	'11',	'12',	'13',	'14',	'15',	'16',	'17',	'18',	'19',	'20',	'21',	'22',	'23',	'24',	'25',	'26',	'27'))


# Absolute values: Stacked bar plot with y scale fixed across metabolites
ggplot(Long_MS_FE_M_Higher, aes(fill= labeling, y=intensity, x=Condition)) + 
 geom_bar(aes(fill = labeling), position = position_stack(reverse = TRUE), stat="identity") +

labs(title="Absolute values of C13 labelling: M+1 isotopologues and above
(Mean)",subtitle = "y scale fixed across metabolites",x="", y = "Absolute intensity")+theme_classic()+

theme(axis.text.x=element_text(size=rel(3), angle=90))+
scale_fill_manual(values = StackPalette_M1)+
facet_wrap(~metabolite, ncol =5)+
theme(strip.text.x = element_text(size = 18, colour = "black", angle = 0))



# Absolute values: Stacked bar plot with y scale adjusted across metabolites
# Attempt to remove white horizontal lines 

ggplot(Long_MS_FE_M_Higher, aes(fill= labeling, y=intensity, x=Condition)) + 
 geom_bar(aes(fill = labeling), position = position_stack(reverse = TRUE), stat="identity") +
labs(title="Absolute values of C13 labelling: M+1 isotopologues and above
(Mean)",subtitle = "y scale relative to metabolite",x="", y = "Absolute intensity")+theme_classic()+
theme(axis.text.x=element_text(size=rel(2), angle=0))+
scale_fill_manual(values =StackPalette_M1)+
facet_wrap(~metabolite, ncol =5, scales="free")+
theme(strip.text.x = element_text(size = 18, colour = "black", angle = 0))




# Normalized values: Stacked bar plot with y scale fixed across metabolites
ggplot(Long_MS_FE_M_Higher, aes(fill= labeling, y=intensity, x=Condition)) + 
 geom_bar(aes(fill = labeling), position = position_fill(reverse = TRUE), stat="identity") +
labs(title="Relative values of C13 labelling: M+1 isotopologues and above
(Mean)",x="", y = "Absolute intensity")+theme_classic()+
theme(axis.text.x=element_text(size=rel(3), angle=90))+
scale_fill_manual(values = StackPalette_M1)+
facet_wrap(~metabolite, ncol =5)+
theme(strip.text.x = element_text(size = 18, colour = "black", angle = 0))

pdf(file = "OnlyFE_isotopologues.pdf", width = 20, height = 50, family = "Helvetica") 

# Absolute values: Stacked bar plot with y scale fixed across metabolites
ggplot(Long_MS_FE_M_Higher, aes(fill= labeling, color = labeling, y=intensity, x=Condition)) + 
 geom_bar(aes(fill = labeling), position = position_stack(reverse = TRUE), stat="identity") +
labs(title="Absolute values of C13 labelling: M+1 isotopologues and above
(Mean)",subtitle = "y scale fixed across metabolites", x="", y = "Absolute intensity")+theme_classic()+
theme(axis.text.x=element_text(size=rel(3), angle=90))+
scale_fill_manual(values = StackPalette_M1)+
scale_color_manual(values = StackPalette_M1)+
facet_wrap(~metabolite, ncol =5)+
theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0))


# Absolute values: Stacked bar plot with y scale adjusted across metabolites

ggplot(Long_MS_FE_M_Higher, aes(fill= labeling, color = labeling, y=intensity, x=Condition)) + 
 geom_bar(aes(fill = labeling), position = position_stack(reverse = TRUE), stat="identity") +
labs(title="Absolute values of C13 labelling: M+1 isotopologues and above
(Mean)",subtitle = "y scale relative to metabolite",x="", y = "Absolute intensity")+theme_classic()+
theme(axis.text.x=element_text(size=rel(2), angle=0))+
scale_fill_manual(values = StackPalette_M1)+
scale_color_manual(values = StackPalette_M1)+
facet_wrap(~metabolite, ncol =5, scales="free")+
theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0))



# Normalized values: Stacked bar plot with y scale fixed across metabolites
ggplot(Long_MS_FE_M_Higher, aes(fill= labeling, color = labeling, y=intensity, x=Condition)) + 
 geom_bar(aes(fill = labeling), position = position_fill(reverse = TRUE), stat="identity") +
labs(title="Relative values of C13 labelling: M+1 isotopologues and above
(Mean)",x="", y = "Absolute intensity")+theme_classic()+
theme(axis.text.x=element_text(size=rel(3), angle=90))+
scale_fill_manual(values = StackPalette_M1)+
scale_color_manual(values = StackPalette_M1)+
facet_wrap(~metabolite, ncol =5)+
theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0))

dev.off()


# Create a new tabel ("poolsizeTable") pool sizes (summed M+ values) per sample
poolsizeTable <- wide_MS_Int %>% 
  group_by(Condition, replicate, SampleID) %>% 
summarise(across(where(is.numeric), sum, na.rm = TRUE))%>% ungroup() 

poolsizeTable$Condition_replicate <- paste(poolsizeTable$Condition,"-", poolsizeTable$replicate)
poolsizeTable <- poolsizeTable %>%
  select(Condition_replicate, everything())

# Remove metabolites that have not been detected
poolsizeTable_noZ <- poolsizeTable[, colSums(poolsizeTable != 0) > 0]

# Remove columns you don't need
poolsizeTable_noZ <-  poolsizeTable_noZ %>% ungroup() %>% select(!c(replicate, SampleID, Condition))

# Transpose and correct row and column names
t_poolsizeTable_noZ <- transpose(poolsizeTable_noZ )
colnames(t_poolsizeTable_noZ) <- rownames(poolsizeTable_noZ)
rownames(t_poolsizeTable_noZ) <- colnames(poolsizeTable_noZ)

# Make first row column names
names(t_poolsizeTable_noZ) <- as.matrix(t_poolsizeTable_noZ[1, ])
t_poolsizeTable_noZ <- t_poolsizeTable_noZ[-1, ]
t_poolsizeTable_noZ[] <- lapply(t_poolsizeTable_noZ, function(x) type.convert(as.character(x)))

# Make row name compound name
t_poolsizeTable_noZ <- setDT(t_poolsizeTable_noZ, keep.rownames = "Compound")
                                
# Reshape data wide > long 
long_poolsize <- melt(t_poolsizeTable_noZ , id.vars = c("Compound"))                    
colnames(long_poolsize) <- c("Compound", "sample", "intensity")              

# Split sample column into 2 columns: condition and replicate number
long_poolsize <-long_poolsize %>%
  separate(sample, c("Condition", "replicate"), "-")

         


# Reshape data long > wide so that now each row is a sample and each column is a metabolite
wide_PS <- dcast(long_poolsize,  Condition + replicate ~ Compound, value.var="intensity")

# Create a new first column: Sample ID
wide_PS <- tibble::rowid_to_column(wide_PS, "SampleID")

# Make first column row names
dataNamed_PS <- wide_PS[,-1]
rownames(dataNamed_PS) <- wide_PS[,1]

dataNamed_PS

# Select just the columns of X data (not qualitative metadata)
Data_X_PS<- dataNamed_PS[3:148]
head (Data_X_PS)

# Compute PCA (scale = TRUE for UV and FALSE for MC)
res.pca <- prcomp(Data_X_PS, scale = TRUE )

# Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
#fviz_eig(res.pca)

options(repr.plot.width=8, repr.plot.height=8, repr.plot.res = 200, repr.plot.quality = 300, repr.plot.pointsize = 12)


fviz_eig(res.pca, addlabels = TRUE)+
labs(title=(expression(paste("R" ^2,": variance explained by each component")))) 


# PCA Scores plot coloured by qualitative variables + 95% confidence elipses 


#groups <- Data_X[9]
fviz_pca_ind(res.pca,
             fill.ind = dataNamed_PS$Condition, col.ind = "black",
             pointshape = 21, pointsize = 7,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             ellipse.level = 0.95,
             legend.title = "Groups",
             mean.point = FALSE,
              ellipse.border.remove = TRUE,
             repel = TRUE
             )+
labs(title="PCA scores of pool sizes")

# PCA loadings plot 

# Graph of variables. Positive correlated variables point to the same side of the plot. 
# Negative correlated variables point to opposite sides of the graph.

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )+
labs(title="PCA loadings of pool sizes")


dataNamed_PS
attach (dataNamed_PS)

# Function to perform a Welch T-test on every column of data (compound) you specify, with Condition used as the groups

Welch_ttest <- function(compound) {
    
    Wttest <- t.test(compound ~ Condition)
     pval = Wttest$p.value
    
    return(pval)
}


# fill in dataNamed[3:1055] such that all columns of metabolite intensity data are calculated (use dataNamed)
# - do not include metadata columns:

Welch_ttest  <- sapply(dataNamed_PS[3:148], 'Welch_ttest')

# Create dataframe of results of the Welch t-test
TtestDF <- as.data.frame(Welch_ttest)

# Make row name metabolite ID column
TtestDF <- setDT(TtestDF, keep.rownames = "Compound")

# Create a new dataframe with a column showing the Benjamini-Hochberg-corrected p-values
BH_testDF <-mutate(TtestDF,BH_correctedVal = p.adjust(TtestDF$Welch_ttest,method="BH") )

# Order the dataframe according to the adjusted pvalues
PoolSize_BH_testDF_orderedByPval <- BH_testDF[order(BH_testDF$Welch_ttest),]

PoolSize_BH_testDF_orderedByPval 
# write a csv of this dataframe
write.csv(PoolSize_BH_testDF_orderedByPval , file = "Welch_BH_orderedResults_poolSizes.csv")

# Create a DF subselection of compounds with uncorrected Welch t-test pval of ≤ 0.05
Ttest_Signif_PS <- PoolSize_BH_testDF_orderedByPval  %>%
filter (Welch_ttest <= 0.05)

# Select only the compounds with uncorrected Welch t-test pval of ≤ 0.05
Significant_Compounds_PS <- t_poolsizeTable_noZ [t_poolsizeTable_noZ$Compound %in% Ttest_Signif_PS$Compound,]
Significant_Compounds_PS

# Append the t-test results
Significant_Compounds_tResults_PS <- select (PoolSize_BH_testDF_orderedByPval,1,2,3)
Significant_Compounds_PS=merge(Significant_Compounds_PS,Significant_Compounds_tResults_PS, by="Compound")

# Order the dataframe according to the adjusted pvalues
Significant_Compounds_PS <- Significant_Compounds_PS[order(Significant_Compounds_PS$Welch_ttest),]
Significant_Compounds_PS <- within(Significant_Compounds_PS, Compound <- reorder(Compound, Welch_ttest))# orders the facets in increasing order of Welch t-test

# write a csv of this dataframe
write.csv(Significant_Compounds_PS, file = "Poolsize_Signif_Compounds_uncorrectedPval.csv")

# Reshape data wide > long (for ggplot graphs with facet)
Significant_Compounds_PS_long <- melt(Significant_Compounds_PS, id.vars = c("Compound","Welch_ttest","BH_correctedVal"))                    

colnames(Significant_Compounds_PS_long) <- c("Compound", "Welch_ttest","BH_correctedVal","sample", "intensity")              

# Split sample column into 2 columns: condition and replicate number
Significant_Compounds_PS_long <-Significant_Compounds_PS_long %>%
 separate(sample, c("Condition", "replicate"), "-")


# Pool sizes for compounds with uncorrected Welch t-test pval of ≤ 0.05

options(repr.plot.width=15, repr.plot.height=15, repr.plot.res = 400, repr.plot.quality = 500, repr.plot.pointsize = 12)

p<-ggplot(data=Significant_Compounds_PS_long, aes(x=Condition, y=intensity, fill = Condition)) +
  geom_bar(position = 'dodge', stat = 'summary', fun = 'mean')+
geom_errorbar(stat="summary", 
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))}, position=position_dodge(.9), width = 0.3, colour="black")+

geom_point(aes(fill=Condition), size=2, shape=21, colour="black", stroke = 1, alpha =1, 
           position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0)) +


labs(title="Pool sizes with significant uncorrected Welch t-test pval (0.05)
Mean ± SEM",x="Condition", y = "intensity")+theme_classic()+
theme(axis.text.x=element_text(angle=60, hjust=1))+
theme(legend.position="none")+
facet_wrap(~Compound, ncol =6, scales="free",labeller = label_wrap_gen())+
theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0))

p

pdf(file = "poolSize_UncorrSignif.pdf", width = 18, height = 10, family = "Helvetica") 

p<-ggplot(data=Significant_Compounds_PS_long, aes(x=Condition, y=intensity, fill = Condition)) +
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
facet_wrap(~Compound, ncol =6, scales="free",labeller = label_wrap_gen())+
theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0))

p
dev.off()
