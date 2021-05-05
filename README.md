# Cornell-MS-fractional-enrichment-pipeline
Pipeline to visualise MS intensity data and perform PCA and univariate stats


There are 2 pipelines. One for unlabelled MS and one for isotope-labelled data.
The PCA and bar charts will work for 2+ group data, the Welch-test section only works for 2-group data.

1. Import data from WCM Ms excel spreadsheet - change the file name within quotations. For fractional enrichment data, this is the intensity report not the fraction report e.g.
`MS_Int <- read_excel('/Users/leesh/Documents/MS analysis pipeline/Nates results 2nd May/KRAS 1wk metabolomics.xlsx')`

2. There are some commands where you have to select columns for analysis - use the dataframe object preceding the [] to decide which columns of data you want to analyze.

e.g 
`Select just the columns of X data (not qualitative metadata)
Data_X_HC<- dataNamed_HC[3:181]
head (Data_X_HC)`

And also

```
fill in dataNamed[3:1055] such that all columns of metabolite intensity data are calculated (use dataNamed)
- do not include metadata columns:
`Welch_ttest  <- sapply(dataNamed_HC[3:181], 'Welch_ttest')
```

3. You can alter the size of the PDF plot space:

e.g.
`pdf(file = "HighConfidence_poolSize_UncorrSignif.pdf", width = 20, height = 30, family = "Helvetica")`
If you have many compounds to visualise, you will need to increase the height.

You can also alter the R markdown notebook plot size (to see figures better in the notebook) by specifying within the chunk, e.g.:
`{r, fig.width=15, fig.height=20}`


4. Optional: You can alter details for graphs: order of the X variables / names of the X variables:
Change order of conditons for graphing

```
#MS_detected$Condition <- factor(MS_detected$Condition, c("CtrlIntact","CtrlSC", "RESTSC", "CDH1SC"))

Rename condition labels

#MS_detected$Condition <- recode_factor(MS_detected$Condition, CtrlIntact  = "Intact", CtrlSC = "Single cells", RESTSC = "shREST SC", CDH1SC = "shCDH1 SC")
```

Working directory:
No need to # Set working directory
Library (here): All the files the script outputs (csv and PDFs) will be saved in the folder where you have the .rmd r script notebook

