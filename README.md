# STAD
R implementation of STAD (Simplified Topological Approximation of Data)

## Description

STAD generates an abstract representation of high-dimensional data by giving each data point a location in a graph which preserves the distances in the original high-dimensional space. The STAD graph is built upon the Minimum Spanning Tree (MST) to which new edges are added until the correlation between the graph and the original dataset is maximized. Additionally, STAD supports the inclusion of filter functions to focus the exploration and allow the analysis of data from new perspectives, emphasizing traits in data which otherwise would remain hidden. 

## Install

[GitHub](https://github.com) is not directly supported by the basic
`install.packages` command. You could use the
[devtools](http://cran.r-project.org/web/packages/devtools/index.html) package
to install the development version of `stad`.

```r
install.packages("devtools")
library("devtools")
install_github("vda-lab/stad")
```

## Example
```r
library(stad)
library(ggplot2)

# Circles dataset
data(circles)

ggplot(circles, aes(x,y, color = lens)) +
  geom_point()

circles_distance <- dist(circles[,c("x", "y")])

## STAD without lens
set.seed(10)
circles_nolens <- stad(circles_distance)
plot_graph(circles_nolens, layout = igraph::layout_with_kk )

## STAD with lens
set.seed(10)
circles_lens <- stad(circles_distance, filter_values = circles$lens, num_intervals = 5)
plot_graph(circles_lens, layout = igraph::layout_with_kk )
```
