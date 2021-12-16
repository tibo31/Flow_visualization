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

```{r, message = F, echo = F, eval = T}
source("codes/plot_flows.R")
```

The mandatory input arguments of the function `plot_flows()` are:

* `y` is the flow vector of size $N$;
* `index_o`, respectively `index_d`, is the vector of size $N$ containing the identifiers of the origins, respectively destinations; 
* `xy_sf` is a spatial object **sf** which contains the coordinates of the $S$ sites;
* `type_plot` is a character giving the name of the methods among `"Sankey"`, `"circular"`, `"arc"`, `"flow_map"`, `"griffith"`, `"heatmap"`

## Illustration used in our article

Importing the data (the data were originally presented in Laurent *et al.*, 2020): 

```{r, message = F}
load("data/data_long.RData") 
```

For reproducing the different graphics we present with any source of flows data, user needs to define the following objects :


* `y` is the vector of flows of size $N$:  

```{r}
y <- data_long[, "y"]
N <- length(y)
```

* `index_o` resp. `index_d` is the vector of size $N$ with the indices of the origin resp. destination:

```{r}
index_o <- data_long[, "cont_o"]
index_d <- substr(data_long[, "cont_d"], 1, nchar(data_long[, "cont_d"]) - 1)
```


* Finally, we need to compute the geographical coordinates as POINT, for each site $s\in S$; it may be stored in a **sf** object (**centroid_sf** hereafter) so that it can be transformed easily in another CRS. The id of the observations may be stored in `S` variable and must coincide with the notation used for the flows. In our example, each geographical region is not defined as an administrative area; thus, we use the centroid of an arbitrary country which belongs to the considered zones: 

```{r, message = F, warning = F}
library(sf)
centroid_sf <- st_centroid(world[
  sapply(c("Canada", "Mexico", "Brazil", "Bahamas",
        "Switzerland", "France", "Poland", "Kazakhstan",
        "Jordan", "Algeria", "Cameroon",
        "India", "Thailand", "China", "Australia"), 
        function(x) grep(x, world$name_long)), ])
```

Then we define the abbreviate and full names of the geographical region:

```{r, message = F, warning = F}
levels_S <- c("N.America", "C.America", "S.America", "Caribbean", "Eu-oth",
                   "UE-bef-2004", "UE-aft-2004", "ex-URSS", "M.East", "N.Africa",
                   "Africa", "S.Asia", "SE.Asia", "E.Asia", "Pacific")
labels_S <- c("North America", "Central America",
             "South America", "Caribbean", "Other European countries",
             "European Union (before 2004 exp)", 
             "European Union (after 2004 exp)",
             "ex URSS", "Middle East", "North Africa",
             "Sub-Saharan Africa", "South Asia", "South East Asia", "East Asia", "Pacific")
```

We associate the names to the variable `S`:

```{r, message = F, warning = F}
centroid_sf$S <- levels_S
```

The CRS must be defined; in our case it is WGS 84:

```{r}
st_crs(centroid_sf) <- 4326
st_crs(world) <- 4326
```


* Find the most appropriated CRS: in our case, we will use two different representations: the Mollweide projection  ("+proj=moll") on one hand which preserves area relationships and "epsg:3035" on other hand which use a projection more convenient for plotting the flows.  User must choose this criteria with respect to the locations of the data (see for instance http://magrit.cnrs.fr/docs/projection_list_fr.html). We apply the transformation on the boundaries of the world countries (that does not correspond to our spatial unit but it covers the same territory and it is helpful to visualize them) and the centroid of the sites $s\in S$

```{r}
poly_sf <- st_transform(world[-160, ], "+proj=moll")
xy_sf <- st_transform(centroid_sf, "+proj=moll")
```

The colors have been chosen with package **colorspace**, one unique color for each $s\in S$ :

```{r}
library("colorspace")
q4 <- qualitative_hcl(length(levels_S), palette = "Dark 3")
names(q4) <- levels_S
```

Version of the input parameters including the full names:

```{r}
index_o_l <- factor(index_o, levels = levels_S, labels = labels_S)
index_d_l <- factor(index_d, levels = levels_S, labels = labels_S)
xy_sf_l <- xy_sf
xy_sf_l$S <- factor(xy_sf_l$S, levels = levels_S, labels = labels_S)
q4_l <- qualitative_hcl(length(labels_S), palette = "Dark 3")
names(q4_l) <- labels_S
```


