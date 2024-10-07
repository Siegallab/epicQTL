# ChromosomeFixation
 
Data and analysis for preprint: doi:

## Data Availability
### Analyzed files available here:
1. **CuSO4_glmer_output.csv**: Additive peaks called using the cybr_lmpeaks() function
2. **Interaction_Peaks.csv**: Interaction peaks, which are the highest points on each chromosome above the 5% FDR
   
### Sequencing Files
4. link to SRR

### Raw Data
5. output.table files in the `Input` folder  
   x **AllCuSO4.REF_.SortedCat.vcf.output.zip**: raw data, zipped via Windows  
   x **Oak_VCF.txt** and **Wine_VCF.txt**: text files for parent data  
   x **experiment_names_524.csv**: index of experiment names, to be loaded in to analysis scripts  

### cybrBSA Package
Please download the package by visiting github: https://github.com/cbuzby/cybrBSA
```
library(devtools)
install_github("cbuzby/cybrBSA")
```


## Analysis Code
### Sequencing
Sequencing scripts are located in the Sequencing folder; we used Nextflow, but each step can be done individually. Generally, steps are as follows (and can be found in the `Aug2024_CGSB_template.nf` Workflow segment):  
```
trim(files_ch) | set{ trimmed_ch }
align(trimmed_ch) | set{ aligned_ch }
sort(aligned_ch) | set{ sorted_ch }
MakeBQSRTable(sorted_ch) | set{ bqsrtables_ch }
ApplyBQSR(bqsrtables_ch) | set{ bqsr_bam_ch }
```

Once raw sequences were processed into bam files for each bulk, we combined all bulks together and ran variant calling on each chromosome (containing all bulks) in parallel:
```
CB_2.3_merge.split.all.q
CB_3.0_Index.q
CB_4.0_CallVariants.q
CB_5.0_zip.concat.sort.q
```

Please see READ.ME in the `Sequencing` folder for execution code.

### Bulk Segregant Analysis
All analysis is done using the `CuPAPER_1byRep_Process.Rmd` file, which utilizes the cybrBSA (https://github.com/cbuzby/cybrBSA) package. Please download this package from github to use. Analysis scripts for output table data can all be found in `CuPAPER_1byRep_Process.Rmd`. Additional examples (and a sample dataset) can be found in the documentation for cybrBSA.

### Visualization
Scripts for each figure of the manuscript are included in the visualizations `CuPAPER_Figures.Rmd`. 

