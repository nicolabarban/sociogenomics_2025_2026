


## PLotting genetic correlations  in R
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
