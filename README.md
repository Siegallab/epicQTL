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
   a. **AllCuSO4.REF_.SortedCat.vcf.output.zip**: raw data, zipped via Windows  
   b. **Oak_VCF.txt** and **Wine_VCF.txt**: text files for parent data  
   c. **experiment_names_524.csv**: index of experiment names, to be loaded in to analysis scripts  

### cybrBSA Package
Please download the package by visiting github: https://github.com/cbuzby/cybrBSA
```
library(devtools)
install_github("cbuzby/cybrBSA")
```


## Analysis Code
### Sequencing
Sequencing scripts are located in the `SequencingScripts` folder; we used Nextflow, but each step can be run individually. Workflow is as follows (and can be found in the `Aug2024_CGSB_template.nf` Workflow segment):  
```
#trim
java -jar /share/apps/trimmomatic/0.36/trimmomatic-0.36.jar SE -phred33 \
 $reads $trimmed_file ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#align
bwa mem -Y -K 100000000 -R \"${readGroup}\" \
 $REF $read_1 > ${id}.aln-se.sam

#sort
java -jar $PICARD_JAR SortSam \
                  INPUT=$read_1 \
                  OUTPUT=${id}.aln-se.bam \
                  SORT_ORDER=coordinate

#MakeBQSRTable
gatk BaseRecalibrator -R $REF -I $read --known-sites $KNOWN_SITES -O ${id}_recal_data.table

#ApplyBQSR
gatk ApplyBQSR \
	   -R $REF \
	   -I ${bamfile} \
	   --bqsr-recal-file ${table} \
	   -O ${id}.bqsr.output.bam
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
Files for analyzing output tables for epicQTL are found in the `BSA_Analysis` folder, with folders for Input and Output data. All analysis is done using the `CuPAPER_1byRep_Process.Rmd` file, which utilizes the cybrBSA (https://github.com/cbuzby/cybrBSA) package. Please download this package from github to use. Analysis scripts for output table data can all be found in `CuPAPER_1byRep_Process.Rmd`. Additional examples (and a sample dataset) can be found in the documentation for cybrBSA.

### Visualization
Scripts for each figure of the manuscript are included in the visualizations `CuPAPER_Figures.Rmd` within the `BSA_Analysis` folder. 

