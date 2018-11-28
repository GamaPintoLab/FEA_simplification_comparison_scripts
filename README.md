# FEA_simplification_comparison_scripts


## DESCRIPTION

These scripts and functions are designed to perform  not only a Functional Enrichment Analysis (FEA) using **ClusterProfiler** R package  
but also to simplify the FEA results by applying gene-co-occurrence and semantic similarity-based simplification algorithms.  

Addditonaly, we can also perfom pair-wise comparisons between two FEA results to compare for example two experimental gene lists functional differences and similarities.  

Finally the repository provides scripts based on **ggplot2** R package to produce aesthetic and meaningfull plots to summarize these comparisons.  

## R SCRIPTS

### simplest_FEA.R

The simplest FEA. There are two example input gene sets available at **input_data** folder  
(gene_symbol_vecA_df.csv and gene_symbol_vecB_df.csv)

The function is written to return and enrichment for Biological Process (BP) GO terms.   
Biological Process (BP) GO terms are the most frequently used type of GO terms but **enrichGO** also accepts Cellular Component (CC) or Molecular Function (MF) as well.
If the user is interested in KEGG, functional annotation use **enrichKEGG** function.  

User must be aware that **enrichGO* function requires ENTREZ_GENE IDs thus, the script also includes a previous step to mapp the gene ids  
from **SYMBOL** to **ENTREZ** ID. Gene or Protein ID mapping is a major issue and care must be taken to get the highest id covering.  
There are numerous R packages to perform id mapping and each retrieves information from varied repositories that in turn are updated very differently,   
therefore each method might return different outputs.

Thus before perfoming the FEA, please check if id mapping function returns a sufficient number of IDs. If not, I recommend to directly  
map your genes at www.uniprot.org website --> "retireve/ID mapping" tab.

Additionally, **enrichGO* function requires an argument with GENE-GO TERM infomation to correctly perfom the FEA. For well-studied species we can directly use   
"annotation" R libraries i.e **org.Hs.eg.db** for human, or **org.mm.eg.db** for mouse.

However, there are some species that lack of this infomation and so is the user who must introduce it manually. 
In this case, load a manually created annotation data.frame with two columns ("GO TERM", "GENE") and then use it as "TERM2GENE" argument at the function called **enricher** as in the last example at the script.
Likewise it requires a manually constructed "universe" argument that usually contains all the proteins described at annotation data.frame.



### FEA_simplification.R

In this script, the function **FEA_and_FILTERING** already includes the id mapping, it further filters FEa reuslts according to an user given gene background threshold.  
It also computes the Fold Enrichment (gene ratio/background ratio)
Additionally, to help the results interpretation, it returns a last output with ENTREZ_IDs mapped again to original ID.  


**ARGUMENTS**  
ORIGINAL_gene_vec
ORIGINAL_id - OUTPUT_id
ann = "BP", "MF", "CC" or "KEGG" annotation type
bg_threshold =  backgroung threhold defined by user (default 0.7)

**OUTPUTS**
ENTREZ_LIST_A_DF = original and final id mapping data.frame
FEA_result = FEA with Background already filtered and original mapping names


Then, **single_funct_simplif3** function will simplifiy the FEA according to a) gene co-occurrence and b) semantic similarity


**ARGUMENTS**  

RAW_FUNC_ENRICH_DATAFRAME_A= input FEA 
gene_cooc_threshold: values between 0 (all functions are merged) to 101 (no function will be merged)
semantic_threshold: values between 0 (all functions are merged) to 1.1 (no function will be merged)
semantic_alg: one of "Wang", "Resnik", "Rel", "Jiang", and "Lin". More info go to **FEA documentation**
semData_info= d_BP or d_MF of d_KEGG  (retrieved in the same script).  
Note that since KEGG has not hierarchy, we cannot apply semantic similarity algortihms thus, FEA could only be simplified by using gene-co occurrence algorithm


**OUTPUTS**
"original_FEA"
"gene_coocc_simp" simplification output only for the first algorithm (gene co-occurrence) 
"sem_sim_simp" simplification output for the first and second algorithm (gene co-occurrence + semantic similarity algorithm) 
"both_simp" simplification output ( (gene co-occurrence + semantic similarity algorithm) )






