# Investigating graphical representations of origin-destination spatial flow data

This document provides the **R** codes used to reproduce the results included in the paper *Investigating graphical representations of origin-destination spatial flow data*. 

Several packages and plenty of lines of code are required for plotting the different figures we present in this paper. To help potential users to represent these figures, we create a unique function called `plot_flows()`, available on Github, which simplifies considerably the syntax code. However, the function still depends on the following R packages: **arcdiagram**, **circlize**, **classInt**, **colorspace**, **ggalluvial**, **gplots**, **igraph**, **sf**, **tidyverse**, **viridis**. 
