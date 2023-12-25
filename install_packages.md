## Install packages
### Seurat v5
Github: https://github.com/satijalab/seurat

```
install.packages('Seurat')
library(Seurat)
```
### Additional packages: 

**SeuratData** - access to Seurat tutorial datasets.

```
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
```
**SeuratWappers** - collections of community-provided methods and extensions for Seurat.
```
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)
```
**BPCells** - high performance and efficient storage.  
```
remotes::install_github("bnprks/BPCells")
```
### tidyverse
Intro: https://www.tidyverse.org/packages/

A set of R packages for data science, including **ggplot2, dplyr, tidyr, stringr, ...**
```
install.packages("tidyverse")
```

### Viridis
Intro: https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html

Makes pretty plots, and improves readability for color blindness and/or color vision deficiency readers.
```
install.packages("viridis")
```
