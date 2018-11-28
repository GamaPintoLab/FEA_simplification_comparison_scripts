# FEA_simplification_comparison_scripts


## DESCRIPTION

These repository contains functions designed to perform not only a Functional Enrichment Analysis (FEA) but also to simplify and compare different sets of FEA outputs.
The FEA is done using **ClusterProfiler** R package. The simplification is done by applying **gene-co-occurrence** and **semantic similarity** -based simplification algorithms.
More info at **FEA documentation**. Finally the repository also provides scripts based on **ggplot2** R package to produce aesthetic and meaningfull plots of the previous analyses.  

## R SCRIPTS FOLDER

1. simplest_FEA.R
2. FEA_simplification.R
3. FEA_comparison.R

### simplest_FEA.R

The simplest way to do a FEA using **ClusterProfiler** R package.  **input_data** folder contains two examples of accepted gene list inputs.
By changing the argument named "ont", the user can specify the type of ontology desired. More info at **FEA documentation**. Biological Process (BP) GO terms are the most commonly used. For KEGG annotation enrichment use the function named **enrichKEGG**.
 
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






