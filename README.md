# ChromosomeFixation
 
Data and analysis for preprint: doi:

## Data Availability
### Analyzed files available here:
1. *something.csv*: Additive peaks called using the cybr_lmpeaks() function
2. *somethingelse.csv*: Interaction peaks, which are the highest points on each chromosome above the 5% FDR
3. *somethingelse*:

### Sequencing Files
4. link to SRR

### Raw Data
5. output.table files in the `Input` folder
   x something
   x something
   x something

### cybrBSA Package
Please download the package by visiting github: https://github.com/cbuzby/cybrBSA
```
library(devtools)
install_github("cbuzby/cybrBSA")
```


## Analysis Code
### Sequencing
Sequencing scripts are located in the Sequencing folder; we used Nextflow, but each step can be done individually. Generally, steps are as follows:

```

```

### Bulk Segregant Analysis
All analysis is done using the `CuPAPER_1byRep_Process.Rmd` file, which utilizes the cybrBSA package. Please download this package from github to use.

### Visualization
Scripts for each figure of the manuscript are included in the visualizations `CuPAPER_Figures.Rmd`. Each 

