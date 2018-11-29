# FEA_simplification_comparison_scripts


## DESCRIPTION

These repository contains functions designed to perform not only a Functional Enrichment Analysis (FEA) but also to simplify and compare different sets of FEA outputs.
The FEA is done using **clusterProfiler** R package [G Yu, et al. OMICS, 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3339379/?report=classic). The simplification is done by applying **gene-co-occurrence** and **semantic similarity** -based simplification algorithms. Finally the repository also provides scripts based on **ggplot2** R package H. Wickham, Springer-Verlag, 2016 to produce aesthetic and meaningfull plots of the previous analyses.  

More info about [clusterProfiler R package](http://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html)
More info about [ggplot2  R package](https://ggplot2.tidyverse.org/index.html)   
More details about the pipeline described in this repository [M.L Garcia-Vaquero et al. Sci. Rep, 2018](https://rdcu.be/bb8Sw) (supplementary material)

## R SCRIPTS FOLDER

1. simplest_FEA.R
2. FEA_simplification.R
3. FEA_comparison.R
4. FEA_ggplots.R

### simplest_FEA.R

The simplest way to do a FEA using **clusterProfiler** R package [G Yu, et al. OMICS, 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3339379/?report=classic).  **input_data** folder contains two examples of accepted gene list inputs. By changing the argument named "ont", the user can specify the type of ontology desired. More info at **FEA documentation**. Biological Process (BP) GO terms are the most commonly used. For KEGG annotation enrichment use the function named **enrichKEGG**.
 
**enrichGO* function requires **ENTREZ_GENE** IDs thus, the script also includes a previous step to mapp the gene ids from **SYMBOL** to **ENTREZ** ID.
**mapIds** function accepts several types of IDs. However, gene or Protein ID mapping is a major issue and care must be taken to get the highest id covering.  
There are numerous R packages to perform id mapping -in this case **org.Hs.eg.db**- and each retrieves information from varied repositories that in turn are updated very differently,   
therefore each method might return different outputs. More info at **FEA documentation**.

Thus before perfoming the FEA, please check if **mapIds** function returns a sufficient number of IDs. If not, I recommend to directly map your genes at www.uniprot.org website --> "retireve/ID mapping" tab.

Additionally, **enrichGO* function requires an argument with GENE-GO TERM information to correctly perfom the FEA. For well-studied species we can directly use   
"annotation" R libraries i.e **org.Hs.eg.db** for human, or **org.mm.eg.db** for mouse. However, there are some species with scarce infomation and so is the user who must introduce it manually. For this case, I also include a last example of how you should load your own annotation data.frame. Likewise it requires a manually constructed "universe" argument that usually contains all the proteins described at annotation data.frame. More info at **FEA documentation**.


### FEA_simplification.R

In this script, the function **FEA_and_FILTERING** already includes the ID mapping, it further filters FEA results according to an user given gene background threshold (I recomend 0.7). It also computes the Fold Enrichment (gene ratio/background ratio) and dditionally, to help the result interpretation, it returns a last output with ENTREZ_IDs mapped again to original ID.  


#### **ARGUMENTS**  
ORIGINAL_gene_vec
ORIGINAL_id - OUTPUT_id
ann = "BP", "MF", "CC" or "KEGG" annotation type
bg_threshold =  backgroung threhold defined by user (default 0.7)

#### **OUTPUTS**
ENTREZ_LIST_A_DF = original and final id mapping data.frame
FEA_result = FEA with Background already filtered and original mapping names


Then, **single_funct_simplif3** function will simplifiy the FEA according to a) gene co-occurrence and b) semantic similarity


#### **ARGUMENTS**  

RAW_FUNC_ENRICH_DATAFRAME_A= input FEA 
gene_cooc_threshold: values between 0 (all functions are merged) to 101 (no function will be merged)
semantic_threshold: values between 0 (all functions are merged) to 1.1 (no function will be merged)
semantic_alg: one of "Wang", "Resnik", "Rel", "Jiang", and "Lin". More info go to **FEA documentation**
semData_info= d_BP or d_MF of d_KEGG  (retrieved in the same script).  
Note that since KEGG has not hierarchy, we cannot apply semantic similarity algortihms thus, FEA could only be simplified by using gene-co occurrence algorithm


#### **OUTPUTS**
"original_FEA"
"gene_coocc_simp" simplification output only for the first algorithm (gene co-occurrence) 
"sem_sim_simp" simplification output for the first and second algorithm (gene co-occurrence + semantic similarity algorithm) 
"both_simp" simplification output ( (gene co-occurrence + semantic similarity algorithm) )



### FEA_comparison.R

This script does 3 things, First it performes the simplest FEA for two independent gene lists.**input_data** folder contains two examples of accepted gene list inputs. Then the comparison of these two simple FEAs and finally, if user desires, it performes a simplification of the combined FEA. 
The comparison consists in merging GO TERMS when their associated genes are present in both sets. The function named **plural_funct_simplif3** returns a single data frame in which the user can distinguish between 3 different susbsets, a) GO TERMS only associated to gene list A, b)  GO TERMS only associated to gene list B and c) GO TERMS associated common to both gene lists. 

The function **plural_funct_simplif3** requires the FEA returned from **enrichGO** [as in simplest_FEA.R] and it also requires the arguments for the FEA simplification [as in FEA_simplification.R]

#### **OUTPUTS**

"original_FEA_A"
"original_FEA_B"
"merged_FEA_AB" = combination of A and B sets (previous to the simplification)
"gene_coocc_simp"= simplification output only for the first algorithm (gene co-occurrence) (for AB set)
"sem_sim_simp" = simplification output for the first and second algorithm (gene co-occurrence + semantic similarity algorithm) (for AB set)
"both_simp"=  simplification of the combined output (AB set)



### FEA_ggplots.R

For the moment it only contains one type of plot, an combination of histograms of the "functional classes" representative of the GO TERMS associated uniquely to SET A, to SET B and to both sets commonly. More info in supplementary material [M.L Garcia-Vaquero et al. Sci. Rep, 2018](https://rdcu.be/bb8Sw)

The first part of the script is identical to [FEA_comparison.R]. Then it applies different functions to edit the output according to **ggplot** requirements.
**EDITING_TO_PLOT_noS2B**, **rearrange_df_toplot_function**

At this point, user should look to the description of the GO Terms retrieved by the analysis ir order to construct its own text mining system to define the most relevant "functional classes". The number of "functional classes" can be changed too. The colors can be tested using **pie** function as well.

**subsetting_summres_plots_func** function will count the times a GO TERMS description matches with the given key terms and returns. Its output is a data frame that can be used as summary.

Finally, **total_wo_others_histograms_function** returns the 3 histograms representing the number of GO GROUPS assigned to each GO CLASS in the three independent sets. Use **ggsave** function to save the plot as pdf file.**g_legend** will return the legend of the GO CLASSES created by the users. 



![]["/Users/Marina/Documents/R_work/FEA_simplification_comparison_scripts/input_data/FEA_comparison_histogram.pdf" "Example of ggplot output"]


