# Lab week 8. Working with summary statistics

## Manhattan Plots, QQplots and genetic correlations  matrices
## Alternative plots using package manhattanly






 Import the summary statistics in R
```
library(tidyverse)
EAgwasResults<-read_csv("https://www.dropbox.com/scl/fi/mjqf8n3700awgged9za67/EA4_results.txt?rlkey=ihu2ll6tzfyh4g1aq3bg2ohfk&dl=1")

head(EAgwasResults)
dim(EAgwasResults)
names(EAgwasResults)
```

This is a very large file, we can speed up things by selecting ony SNPS with Pvalue <0.005
```
EAgwasResults_sub<-subset(EAgwasResults, P<0.0005)
dim(EAgwasResults_sub)

```



We can use a library created to plot manhattan plots
Load the manhattanly library
```
install.packages("manhattanly")
library(manhattanly)
```

## Create interactive Manhattan plot 
```
help(manhattanly)
manhattanly(EAgwasResults_sub, snp = "rsID" , bp="BP", p="P", chr="Chr")
```


## Load the library qqman
```
install.packages("qqman")
 library(qqman)

```

Save the figure into an external png file format
```
png(file="manhattan_without_highlights.png" , width = 1200, height = 600)


manhattan(EAgwasResults_sub, 
            chr="Chr",
	 					bp="BP",
						snp="rsID",
						p="P",
						suggestiveline=F)
						 
						 
dev.off()
```


QQplot
```
qq(EAgwasResults_sub$P)
```

Highlighting significant results

```
 hits<-EAgwasResults_sub[EAgwasResults_sub$P<5e-8,]

EAgwasResults_sub$highlight.snps<-0

for ( i in 1: dim(hits)[1]){

chr<-as.numeric(hits[i,4])
loc_min<- as.numeric(hits[i, 5]-5000)
loc_max<- as.numeric(hits[i, 5]+5000)

neighbours.snps<- EAgwasResults_sub$rsID[EAgwasResults_sub$Chr==chr & EAgwasResults_sub$BP>loc_min & EAgwasResults_sub$BP<loc_max]

EAgwasResults_sub$highlight.snps[EAgwasResults_sub$rsID %in% neighbours.snps] <- 1
}
```



```
png(file="manhattan_with_highlights.png" , width = 1200, height = 600)
# add highlight command to the Manhattan plot 

manhattan(EAgwasResults_sub, 
          chr="Chr",
          bp="BP",
          snp="rsID",
          p="P", 
          highlight=EAgwasResults_sub$rsID[EAgwasResults_sub$highlight.snps==1],
                            suggestiveline=F)
dev.off()
```











## Plotting genetic correlations  in R
import data on genetic correlation

```

data_rg<-read.table("http://nicolabarban.com/sociogenomics2023/week8/LD-Hub_genetic_correlation_example.txt",fill =T, sep="\t", header=T, quote="") 
```



draw heatmap
```
install.packages("ggplot2")
library(ggplot2)
ggplot(data = data_rg, aes(Trait1, Trait2, fill = rg))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = 
                                 "white",  midpoint = 0, limit = 
            c(-1.1,1.1), space = "Lab",
              name="Genetic\nCorrelation") +
      theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
          size = 8, hjust = 1))+
 coord_fixed()
```
## Locuszoom


* http://locuszoom.org



```
devtools::install_github("myles-lewis/locuszoomr")

library(locuszoomr)
data(SLE_gwas_sub)

library(BiocInstaller)
biocLite("EnsDb.Hsapiens.v75")

loc <- locus(gene = 'UBE2L3', SLE_gwas_sub, flank = 1e5)
summary(loc)
locus_plot(loc)

# Or FTP download the full summary statistics from
# https://www.ebi.ac.uk/gwas/studies/GCST003156
library(data.table)
SLE_gwas <- fread('../bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv')

loc <- locus(gene = 'UBE2L3', SLE_gwas, flank = 1e5)
locus_plot(loc)

```

