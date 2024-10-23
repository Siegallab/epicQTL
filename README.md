# epic-QTL Data and Analysis Scripts
 
Here we have included the data and analysis for preprint: doi:TBD

All raw data is linked or included below, and analysis scripts are included in the two folders: `SequencingScripts` and `BSA_Analysis`. Supplemental figures and text can be found in PDF form in publication (linked above).

***  

## Data Availability
### Analyzed files available here:
1. **CuSO4_glmer_output.csv**: Additive peaks called using the cybr_lmpeaks() function
2. **Interaction_Peaks.csv**: Interaction peaks, which are the highest points on each chromosome above the 5% FDR
   
### Raw Sequencing Files
4. SRA Repository PRJNA1175662

### Input Data
5. output.table files in the `Input` folder  
   a. **AllCuSO4.REF_.SortedCat.vcf.output.zip**: raw data, zipped via Windows  
   b. **Oak_VCF.txt** and **Wine_VCF.txt**: text files for parent data  
   c. **experiment_names_524.csv**: index of experiment names, to be loaded in to analysis scripts  

### cybrBSA Package
Please download the package by visiting the [cybrBSA github](https://github.com/cbuzby/cybrBSA) or by running the following code in R:
```
library(devtools)
install_github("cbuzby/cybrBSA")
```

***  

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
############## CB_2.3_merge.split.all.q #######################################

samtools merge AllCuSO4 $1 $2... ${70}
samtools index AllCuSO4
bamtools split -in AllCuSO4 -reference

############## CB_3.0_Index.q #################################################

samtools index $1

############## CB_4.0_CallVariants.q ##########################################

gatk HaplotypeCaller -I $1 -R $REF -ploidy 1 -O ${1}.vcf

############## CB_5.0_zip.concat.sort.q #######################################

for i in ${1}*vcf; do bgzip -c $i > ${i}.gz; done
echo ${1}*vcf.gz |  xargs -n1 tabix -p vcf
bcftools concat -o unsortedcat.vcf -a -D ${1}*vcf.gz
bcftools sort -Oz -o ${1}.SortedCat.vcf unsortedcat.vcf

myfile=${1}.SortedCat.vcf

gatk VariantsToTable \
     -V ${myfile} \
     -F CHROM -F POS -F REF -F ALT \
     -GF AD -GF DP -GF GQ -GF PL \
     -O ${myfile}.output.table
```

Please see READ.ME in the `Sequencing` folder for execution code.

### Bulk Segregant Analysis
Files for analyzing output tables for epicQTL are found in the [BSA_Analysis](https://github.com/Siegallab/epicQTL/tree/main/BSA_Analysis) folder, with folders for `Input` and `Output` folders for data. All analysis is done using the `CuPAPER_1byRep_Process.Rmd` file, which utilizes the [cybrBSA](https://github.com/cbuzby/cybrBSA) package. Analysis scripts for output table data can all be found in [/BSA_Analysis/CuPAPER_1byRep_Process.Rmd](https://github.com/Siegallab/epicQTL/blob/main/BSA_Analysis/CuPAPER_1byRep_Process.Rmd). Additional examples (and a sample dataset) can be found in the documentation for cybrBSA.

The output table is read in through `cybrInputGATKTable()`, and then "called" for parent alleles using Wine.txt and Oak.txt (which are included in supplemental data).
```
mydatatotest = "AllCuSO4.REF_.SortedCat.vcf.output"

cybrInputGATKTable(mydatatotest) %>% mutate(Coverage = as.numeric(AD.REF) + as.numeric(AD.ALT)) %>%
  select(POS, CHROM, Dataset, GQ, AD.REF, AD.ALT, Coverage) -> rawdata

parentSNPids <- cybrConvertParentalAlleles(Truncate = TRUE)

rawdata %>% 
  merge(parentSNPids) %>% 
  filter(grepl("HNGLVDRXY", Dataset)) %>% #Filter for smaller dataset
  mutate(REFW = as.numeric(Type == "Wine"), REFO = as.numeric(Type == "Oak")) %>%
  group_by(Dataset, CHROM, POS) %>%
  mutate(Wine = max(REFW * as.numeric(AD.REF), REFO * as.numeric(AD.ALT)),
         Oak = max(REFW * as.numeric(AD.ALT), REFO * as.numeric(AD.REF))) %>%
  select(Dataset, POS, CHROM, Coverage, Wine, Oak) %>%
  pivot_longer(c(Oak, Wine), names_to = "Allele", values_to = "Reads") -> rawdata_called
```

Next, data is smoothed to reduce the error from sequencing at each position, using a weighted gaussian:
```
rawdata_called %>% group_by(Dataset, CHROM, Allele) %>% arrange(POS) %>%
  reframe(POS = POS, 
          SmoothCount = ceiling(frollapply(Reads, n = 200, FUN = cybr_weightedgauss, align = "center"))
          ) -> rawdata_smoothed
```

Finally, z-scores are calculated for EACH position using `glmer_cb2_short()`. The contrasts in our experiment are adjusted using the `contrasts()` function in R.
```
rawdata_glm_prep %>% na.omit() %>% 
  filter(CHROM %in% c("M", "I") == FALSE) %>%
  group_by(CHROM, POS) %>%
  mutate_if(is.character, as.factor) %>%
  summarize(summary = glmer_cb2_short(Allele = Allele,
                             Bulk = Bulk,
                             Parent = Parent,
                             Rep = Rep, #any parameters not used should NOT be included
                             W = SmoothCount,
                             formula = "Allele ~ Bulk * Parent + (1 | Rep)",
                             outputlength = 4),
            #MAKE SURE THIS IS THE SAME LENGTH AS OUTPUT LENGTH
            Factor = (c("Intercept", "Bulk", "Parent", "Interaction"))) -> processed_glm_all
```

Permutations for the false discovery rate of 5% are done as described in _Materials and Methods_.

### Visualization
Scripts for each figure of the manuscript are included in the visualizations [/BSA_Analysis/CuPAPER_Figures.Rmd](https://github.com/Siegallab/epicQTL/blob/main/BSA_Analysis/CuPAPER_Figures.Rmd), and utilize the `Ouput` data within the folder. 

