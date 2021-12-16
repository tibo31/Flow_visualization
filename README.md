# Investigating graphical representations of origin-destination spatial flow data

This document provides the **R** codes used to reproduce the results included in the paper *Investigating graphical representations of origin-destination spatial flow data*. 

Several packages and plenty of lines of code are required for plotting the different figures we present in this paper. To help potential users to represent these figures, we create a unique function called `plot_flows()`, available on Github, which simplifies considerably the syntax code. However, the function still depends on the following R packages: **arcdiagram**, **circlize**, **classInt**, **colorspace**, **ggalluvial**, **gplots**, **igraph**, **sf**, **tidyverse**, **viridis**. 

Interested users can install the required packages from CRAN or Github as follows: 

```{r, message = F, eval = F}
install.packages(c(
  "circlize", # circular flows charts
  "classInt", # discretization of numeric variable
  "colorspace", # color palettes
  "devtools", # necessary to install Github packages 
  "ggalluvial", # Sankey diagram
  "ggridges",  # ridges plot
  "gplots", # heat maps,
  "igraph", # useful to construct adjencency matrix
  "plot.matrix", # heat maps
  "sf", # spatial norm
  "tidyverse", # tidyverse universe
  "viridis" # color palette for heat map
  )) 
```

```{r, eval = F, message = F}
devtools::install_github('gastonstat/arcdiagram')   # Arc diagram plot
devtools::install_github("LukeCe/spflow")           # flows data
```

Users can now load our function from Github as follows:

