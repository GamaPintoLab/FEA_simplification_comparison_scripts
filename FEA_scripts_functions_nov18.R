
## LOAD LIBRARIES

library(org.Hs.eg.db) 
library("clusterProfiler")


## LOAD IMPUT DATA (GENE VECTOR TO FUNCTIONALLY ENRICH)

setwd("~/Documents/R_work/FEA_simplification_comparison_scripts")


setwd("~/Documents/MARINA LUQUE/JPND project/TDP43 interactome/TDP43_mass_spec_analysis_nov18_interm")
BirA_05_df_variable_output=read.csv("BirA_05_df_variable_output.csv", header = T, stringsAsFactors = F)
gene_symbol_vec=BirA_05_df_variable_output$Protein.IDs
gene_symbol_vecA=gene_symbol_vec[1:500]
gene_symbol_vecB=gene_symbol_vec[501:1000]

gene_symbol_vecA_df=data.frame(gene_symbol_vecA, stringsAsFactors = F)
gene_symbol_vecB_df=data.frame(gene_symbol_vecB, stringsAsFactors = F)

setwd("~/Documents/R_work/FEA_simplification_comparison_scripts/input_data")
#write.csv(gene_symbol_vecA_df, file="gene_symbol_vecA_df.csv", row.names = F)
#write.csv(gene_symbol_vecB_df, file="gene_symbol_vecB_df.csv", row.names = F)

gene_symbol_vecA_df=read.csv("gene_symbol_vecA_df.csv", header=T, stringsAsFactors = F)
gene_symbol_vecB_df=read.csv("gene_symbol_vecB_df.csv", header=T, stringsAsFactors = F)


uniprotA=unique(sel_a1h3S2B_resdf[1:30,1])

## SIMPLEST WAY TO DO THE FEA  
# enrichGO function requires ENTREZ_GENE IDs

# If necessary, These 4 lines map UNIPROT IDs to ENTREZ IDs
ENTREZ_LIST_A=mapIds(org.Hs.eg.db, keys=(as.vector(gene_symbol_vecA)),column="ENTREZID", keytype = "SYMBOL", multivals="first")
ENTREZ_LIST_A=as.vector(as.numeric(ENTREZ_LIST_A))

ENTREZ_LIST_B=mapIds(org.Hs.eg.db, keys=(as.vector(gene_symbol_vecB)),column="ENTREZID", keytype = "SYMBOL", multivals="first")
ENTREZ_LIST_B=as.vector(as.numeric(ENTREZ_LIST_B))

# Functional Enrichment Analysis (FEA)
Aentrezgenes_enrichGOres=enrichGO(ENTREZ_LIST_A,'org.Hs.eg.db',ont="BP")
Aentrezgenes_enrichGOres=as.data.frame(Aentrezgenes_enrichGOres)

Bentrezgenes_enrichGOres=enrichGO(ENTREZ_LIST_B,'org.Hs.eg.db',ont="BP")
Bentrezgenes_enrichGOres=as.data.frame(Bentrezgenes_enrichGOres)




## Functional Enrichment Analysis (FEA) and Background filtering


uniprotA_BP07_res_A=FEA_and_FILTERING(gene_symbol_vecA, "SYMBOL","BP", 0.7, "UNIPROT")
uniprotA_MF07_res_A=FEA_and_FILTERING(gene_symbol_vecA,"SYMBOL","MF", 0.7, "UNIPROT")
uniprotA_KEGG07_res_A=FEA_and_FILTERING(gene_symbol_vecA,"SYMBOL","KEGG", 0.7, "UNIPROT")

uniprotA_BP07_res_B=FEA_and_FILTERING(gene_symbol_vecB, "SYMBOL","BP", 0.7, "UNIPROT")



#uniprotA_KEGG07_res$filtered_FEA

uniprotA_BP07_res_A_r=uniprotA_BP07_res_A$filtered_FEA
uniprotA_MF07_res_A_r=uniprotA_MF07_res_A$filtered_FEA
uniprotA_KEGG07_res_A_r=uniprotA_KEGG07_res_A$filtered_FEA



##### FEA simplification

# there is an initial filtering by gene background (0.70 by default)
# simplification by a) gene co-occurrence, b) semantic similarity or c) both
# if you want to avoid one of themm , use a 1 cutoff

library("GOSemSim")
library("corpcor")

d_BP <- godata('org.Hs.eg.db', ont="BP", computeIC=T) 
d_MF <- godata('org.Hs.eg.db', ont="MF", computeIC=T) 
d_KEGG=0
class(d_BP)=="GOSemSimDATA"
class(d_KEGG)!="GOSemSimDATA"



# if you want both simplifications define your thresholds
# if you want only simp by gene-co occurrence , define semantic_threshold=1.1
# if your want only simp by sem sim, define gene_cooc_threshold=1.1



RAW_FUNC_ENRICH_DATAFRAME_A=uniprotA_BP07_res_A_r

gene_cooc_threshold=1 # IF >1 MEANS DONT DO GENE CO-OCCURRENCE MERGE
semantic_alg="Lin"
semantic_threshold=0.9
semData_info=d_BP
simp_order="cooc_semsim"


single_funct_simplif2=function(RAW_FUNC_ENRICH_DATAFRAME_A,simp_order, gene_cooc_threshold, semantic_alg, semantic_threshold, semData_info){
  
  Auni_size=length(RAW_FUNC_ENRICH_DATAFRAME_A[,1])
  Buni_size=0
  ABcom_size=0
  
  total=ABcom_size+Auni_size+Buni_size 
  doubleset=c(rep(1,ABcom_size),rep(0,total-ABcom_size))  
  maxscore=as.numeric(pmax(RAW_FUNC_ENRICH_DATAFRAME_A[,5],RAW_FUNC_ENRICH_DATAFRAME_A[,6])) #by fold score
  reorder=order(doubleset,maxscore,decreasing=T)
  ordered_ABfun_enr=RAW_FUNC_ENRICH_DATAFRAME_A[reorder,]
  
  ordered_ABfun_enr=cbind(ordered_ABfun_enr, 1,1,1,1,1,1)
  colnames(ordered_ABfun_enr)=c("GO_Terms_ID","Description","Original_Genes_A","Original_Genes_B","Original_Fold_A","Original_Fold_B","Merged_GOT","Merged_desc","Merged_A_genes","Merged_B_genes","Merged_A_fold","Merged_B_fold")
  
  
  if (simp_order=="cooc_semsim"){
    geneco_AB_func_enrichment=GENECO_MERGING_GOT(ordered_ABfun_enr,0,length(ordered_ABfun_enr[,1]),0, gene_cooc_threshold, gene_cooc_threshold)
    if (class(semData_info)=="GOSemSimDATA") {
      sem_AB_func_enrichment=SEMANTIC_ANALYSIS3(semantic_alg,semantic_threshold,geneco_AB_func_enrichment,semData_info)
    }
    if (class(semData_info)!="GOSemSimDATA"){
      sem_AB_func_enrichment=geneco_AB_func_enrichment
      print( "KEGG functions do not have hierarchy")
    }
    final_simplified_res=sem_AB_func_enrichment
    print(paste("from an original size of",length(ordered_ABfun_enr[,1]) ,
                "functions, we get a simplified list of",length(final_simplified_res[,1]),"functional groups",sep=" ", collapse = " "  ))
    
    #list("original_FEA"=RAW_FUNC_ENRICH_DATAFRAME_A,"gene_coocc_simp"=geneco_AB_func_enrichment,"sem_sim_simp"=sem_AB_func_enrichment,"both_simp"=final_simplified_res)
  }
  
  if (simp_order=="semsim_cooc"){
    if (class(semData_info)=="GOSemSimDATA") {
      sem_AB_func_enrichment=SEMANTIC_ANALYSIS3(semantic_alg,semantic_threshold,ordered_ABfun_enr,semData_info)
    }
    if (class(semData_info)!="GOSemSimDATA"){
      sem_AB_func_enrichment=ordered_ABfun_enr
      print( "KEGG functions do not have hierarchy")
    }
    geneco_AB_func_enrichment=GENECO_MERGING_GOT(sem_AB_func_enrichment,0,length(sem_AB_func_enrichment[,1]),0, gene_cooc_threshold, gene_cooc_threshold)    
    final_simplified_res=geneco_AB_func_enrichment
    print(paste("from an original size of",length(ordered_ABfun_enr[,1]) ,
                "functions, we get a simplified list of",length(final_simplified_res[,1]),"functional groups",sep=" ", collapse = " "  ))
    
    
  }
  list("original_FEA"=ordered_ABfun_enr,"gene_coocc_simp"=geneco_AB_func_enrichment,"sem_sim_simp"=sem_AB_func_enrichment,"both_simp"=final_simplified_res)
}
#sem_AB_func_enrichment1=sem_AB_func_enrichment
#sem_AB_func_enrichment2=sem_AB_func_enrichment

length(colnames(sem_AB_func_enrichment1))
length(colnames(sem_AB_func_enrichment2))


single_g09s09simp_BP=single_funct_simplif2(uniprotA_BP07_res_A_r,"cooc_semsim",1.1,"Lin", 0.9, d_BP)
# outputs list("original_FEA"=ordered_ABfun_enr,"gene_coocc_simp"=geneco_AB_func_enrichment,"sem_sim_simp"=sem_AB_func_enrichment,"both_simp"=final_simplified_res)

simplifies_output_A=single_g09s09simp_BP[[4]]

#single_g09s09simp_MF=single_funct_simplif2(uniprotA_MF07_res_r,"cooc_semsim",0.9,"Lin", 0.9, d_MF)





### these are the suggested thresholds and semantic algorithm
bg_threshold=0.1
gene_cooc_threshold=0.7
semantic_alg="Lin"
semantic_threshold=0.7




##### FEA simplification AND COMPARISON FOR TWO GENE-FUNCTIONS LISTS


# before running the simplification function, go to line 69-end and run ALL THE FUNCTIONS
SIMPandCOMP_OUTPUT=pluralwithoutsubset_funct_simplif(Aentrezgenes_enrichGOres,Bentrezgenes_enrichGOres,0.1,0.7,0.7)
SIMPandCOMP_OUTPUT_clean=EDITING_TO_PLOT_noS2B(length(ENTREZ_LIST_A),length(ENTREZ_LIST_B) ,SIMPandCOMP_OUTPUT[[4]])
SIMPandCOMP_OUTPUT_clean_res=list("simplification_results_list"=SIMPandCOMP_OUTPUT,"simplified_genefreqandfold_df"=SIMPandCOMP_OUTPUT_clean)

#OUTPUT: 
# SIMPandCOMP_OUTPUT_clean_res[[1]][[1]] --> FEA FILTERED BY FOLD ENRICHMENT
# SIMPandCOMP_OUTPUT_clean_res[[1]][[2]] --> FEA MERGED BY GENE CO.OCCURRENCE
# SIMPandCOMP_OUTPUT_clean_res[[1]][[3]] --> FEA MERGED BY SEMANTIC SIMILARITY
# SIMPandCOMP_OUTPUT_clean_res[[1]][[4]] --> FEA MERGED BY GENE CO.OCCURRENCE AND SEMANTIC SIMILARITY (save each GOT fold)
# SIMPandCOMP_OUTPUT_clean_res[[2]] --> FEA MERGED BY GENE CO.OCCURRENCE AND SEMANTIC SIMILARITY WITH RESULTING GENEFREQ AND MEDIAN FOLD of the GOT group

# function groups (rows) with values ONLY in columns denoted with "A" will be functions only associated to gene list A, 
# function groups (rows) with values ONLY in columns denoted with "B" will be functions only associated to gene list B, 
# function groups (rows) with values IN BOTH columns denoted with "A" and "B" will be functions associated to BOTH gene lists



##### FEA simplification FOR SINGLE GENE-FUNCTIONS LIST

# before running the simplification function, go to line 69-end and run ALL THE FUNCTIONS
SIMPLIFICATION_OUTPUT=single_funct_simplif(Aentrezgenes_enrichGOres,0,0.1,0.7,0.7)
SIMPLIFICATION_OUTPUT_clean=EDITING_TO_PLOT_noS2B(length(ENTREZ_LIST_A),0 ,SIMPLIFICATION_OUTPUT[[4]])
SIMPLIFICATION_OUTPUT_res=list("simplification_results_list"=SIMPLIFICATION_OUTPUT,"simplified_genefreqandfold_df"=SIMPLIFICATION_OUTPUT_clean)

#OUTPUT: 
# SIMPLIFICATION_OUTPUT_res[[1]][[1]] --> FEA FILTERED BY FOLD ENRICHMENT
# SIMPLIFICATION_OUTPUT_res[[1]][[2]] --> FEA MERGED BY GENE CO.OCCURRENCE
# SIMPLIFICATION_OUTPUT_res[[1]][[3]] --> FEA MERGED BY SEMANTIC SIMILARITY
# SIMPLIFICATION_OUTPUT_res[[1]][[4]] --> FEA MERGED BY GENE CO.OCCURRENCE AND SEMANTIC SIMILARITY (save each GOT fold)
# SIMPLIFICATION_OUTPUT_res[[2]] --> FEA MERGED BY GENE CO.OCCURRENCE AND SEMANTIC SIMILARITY WITH RESULTING GENEFREQ AND MEDIAN FOLD of the GOT group
# this function was adapted from the comparison function thus, results will present some columns ("denoted by B") in blank because there is no comparing set


## STEP BY STEP

### FEA SIMPLIFICATION


d <- godata('org.Hs.eg.db', ont="BP", computeIC=T) 

### these are the suggested thresholds and semantic algorithm
bg_threshold=0.1
gene_cooc_threshold=0.7
semantic_alg="Lin"
semantic_threshold=0.7

filtering_A_enrichGO_res=FILTERING_BACKGROUND(NotAFF_enrichGOBP,0.7) 
withfold_filtering_A_enrichGO_res=fold_computing(filtering_A_enrichGO_res[[1]])
geneco_AB_func_enrichment_res=GENECO_MERGING_GOT(withfold_filtering_A_enrichGO_res,0,
                                                 length(withfold_filtering_A_enrichGO_res[,1]),0, gene_cooc_threshold, gene_cooc_threshold)
sem_AB_func_enrichment_res=SEMANTIC_ANALYSIS2(semantic_alg,semantic_threshold,geneco_AB_func_enrichment_res)
final_simplified_res_res=sem_AB_func_enrichment_res[[1]]
sem_merge_summary_res=sem_AB_func_enrichment_res[[2]]
SIMPLIFICATION_OUTPUT_res=list(withfold_filtering_A_enrichGO_res,geneco_AB_func_enrichment_res,sem_merge_summary_res,final_simplified_res_res)

SIMPLIFICATION_OUTPUT_clean_res=EDITING_TO_PLOT_noS2B(length(length(union(NotAFF_uniprot_all,NotAFF_entrez_all))),0 ,SIMPLIFICATION_OUTPUT_res[[4]])
SIMPLIFICATION_OUTPUT_res_res=list("simplification_results_list"=SIMPLIFICATION_OUTPUT_res,"simplified_genefreqandfold_df"=SIMPLIFICATION_OUTPUT_clean_res)


#####





####### PLOT RESULTS 


library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)


original_df=SIMPandCOMP_OUTPUT_clean_res$simplification_results_list[[1]][[1]]
resulting_df=SIMPandCOMP_OUTPUT_clean_res$simplified_genefreqandfold_df



###### re-arrange df to plot
rearrange_df_toplot_function=function(ori_df, res_df){
  Total_GOT=vector()
  RepDesc_vec=vector()
  for (i in 1:length(resulting_df[,1])){
    DESCs=strsplit(resulting_df[i,8],split="/") #merged Desc
    GOTs=strsplit(resulting_df[i,7],split="/") #merged GOTs
    t_GOTs=length(which(GOTs[[1]]!="1"))
    vec_gene_numA=vector()
    vec_gene_numB=vector()
    for (j in 1:length(GOTs[[1]]) ) {
      ind_GOT=which(original_df[,1]==GOTs[[1]][j])
      if (length(ind_GOT)!=0){
        merged_genesA=strsplit(original_df[ind_GOT,3], split="/")
        merged_genesB=strsplit(original_df[ind_GOT,4], split="/")
        vec_gene_numA=c(vec_gene_numA,length(merged_genesA[[1]]))
        vec_gene_numB=c(vec_gene_numB,length(merged_genesB[[1]]))
      }
    }
    vec_both=vector()
    for (k in 1:length(GOTs[[1]]) ) {
      if (GOTs[[1]][k]!="1"){
        vec_both[k]=c(max(vec_gene_numA[k],vec_gene_numB[k]))
      }
    }
    ordered_vec_both=order(vec_both, decreasing = T) # final_desc_ind[1] would be the got with highest gene freq in both sets (to select as repdesc)
    final_desc_ind=ordered_vec_both[1]
    final_desc=DESCs[[1]][final_desc_ind]
    RepDesc_vec=c(RepDesc_vec,final_desc)
    Total_GOT=c(Total_GOT,t_GOTs)
  }
  edited_resulting_df=data.frame("RepDesc"=RepDesc_vec,"Total_GOT"=Total_GOT,resulting_df[,8],resulting_df[,1:7],resulting_df[,9:16],stringsAsFactors = F)
  names=colnames(resulting_df)
  names_ed=names
  names_ed=c("RepDesc","Total_GOT",names[8],"original_GO_Terms_ID","original_Description",names[3:7],names[9:16])
  colnames(edited_resulting_df)=names_ed
  edited_resulting_df=edited_resulting_df[order(edited_resulting_df[,2], decreasing = F),] # order rows according to the # of merged GOTs
  #df_toplot=data.frame()
  df_toplot=data.frame(edited_resulting_df$RepDesc, edited_resulting_df$RepDesc, edited_resulting_df$A_medianfold, edited_resulting_df$B_medianfold, 
                       rep(0, length(edited_resulting_df[,1])), rep(0, length(edited_resulting_df[,1])), stringsAsFactors = F)
  df_toplot[,1]=c(1:length(edited_resulting_df[,1]))
  colnames(df_toplot)=c("GOT_index","RepDesc", "Fold_setA", "Fold_setB", "sumGeneFreq", "Function_class" )
  for (i in 1:length(df_toplot[,1])){
    if (is.na(df_toplot[i,3])){
      df_toplot[i,3]=0
    }
    if (is.na(df_toplot[i,4])){
      df_toplot[i,4]=0
    }
    df_toplot[i,5]=edited_resulting_df$A_genesfreq[i]+edited_resulting_df$B_genesfreq[i]
  }
  list("df_toplot"=df_toplot,"edited_resulting_df"=edited_resulting_df)
}





rearrange_fun_res=rearrange_df_toplot_function(original_df,resulting_df)


rearrange_fun_res$df_toplot



# look to  rearrange_fun_res$edited_resulting_df[,3] GOT descriptions in order to define the best key terms to represent your data

#### DEFINE THE KEY TERMS FOR FUNCTIONAL CLASSES ASSIGMENT BY TEXT MINING 


#group1 will be considered OTHERS (GOTs that do not match with any key term)
group2=c("neuron", "synaptic", "axon", "microglial", "glial", "neural", "neuromuscular","neurogenesis", "nervous","Alzheimer's","Glioma", "(ALS)")
group3=c("immune", "host", "pathogen", "interferon-beta", "cytokine", "fungus", "interleukin-2", "interleukin-1", "leukocyte", "Viral")
group4=c("muscle")
group5=c("stress", "heat", "oxidative", "UV", "X-ray", "superoxide")
group6=c("aggregation", "folding")
group7=c("apoptosis", "apoptotic", "autophagy", "Apoptosis")
group8=c("cytoskeleton", "microtubule", "actin")
group9=c("RNA", "processing", "mRNA", "spliceosomal", "splice", "Spliceosome")
group10=c("transcription", "chromatin", "histone")
group11=c("DNA", "repair")
group12=c("degradation", "proteolysis", "ubiquitination", "deubiquitination", "ERAD")
group13=c("cycle", "mitotic", "cytokinesis")
group14=c("localization", "transport", "import", "export", "targeting")
group15=c("transduction", "cascade", "signaling", "signal")
group16=c("development", "developmental", "differentiation", "embryo", "embryonic", "morphogenesis")


## CAUTION, If you change the number of classes, please delete (or add) them in the following line called "GROUPS"
GROUPS=list(group2,group3,group4,group5,group6,group7,group8,group9,group10,group11,group12,group13,group14, group15, group16)


## to be plot in the plots legend, you will be asked to choose a) to plot all key terms ("group_names")
## or b) simpler functional class names that you have to edit in the following lines ("group_names_plot")

group_names_plot=c("Other","Nervous system",
                   "Immune system", 
                   "Muscle",
                   "Stress", 
                   "Folding",
                   "Apoptosis", 
                   "Cytoskeleton", 
                   "RNA processing",
                   "Transcription", 
                   "DNA repair", 
                   "Protein degradation",
                   "Cell cycle", 
                   "Protein export/import",
                   "Signaling", 
                   "Development")


## colors. If you prefer other colors to be ploted, change here the codex
## the order of the colors will match with the order of classes assigned above

colors16=c( "#546269" ,"red" ,    "#FF6600", "#FFCC00" ,"#CCFF00", "green",   "#71c100" ,"#a37458" ,"#00FFCC", "#00CCFF" ,
            "#0066FF", "blue" ,   "#6600FF" ,"#CC00FF","#FF00CC", "#e4acf6")
pie(rep(1,16),col=colors16) 

##

subsetting_summres_plots_func=function(GROUPS,edited_resulting_df,df_toplot){
  functions_groups_df=as.data.frame(matrix(0, nrow=length(edited_resulting_df[,1]), ncol=length(GROUPS)+1)) # nrow=#functional groups to assign , # ncol= #functional classes defined above + 1
  for ( i in 1:length(edited_resulting_df[,1])){
    funcpast=paste(edited_resulting_df[i,1],edited_resulting_df[i,3],sep="/")
    functions_groups_df[i,1]=funcpast
  }
  for ( i in 1:length(functions_groups_df[,1])){
    fun_count=0
    splited_funct=strsplit(functions_groups_df[i,1], "/")
    for (j in 1:length(GROUPS)){
      for ( h in 1:length(splited_funct[[1]])){
        splited_funct_sent=strsplit(splited_funct[[1]][[h]], "/| ")
        fun_count=length(which(is.element(splited_funct_sent[[1]],GROUPS[[j]])))
        if (fun_count>=1){
          count=(functions_groups_df[i,j+1])+1
          functions_groups_df[i,j+1]=count
        }
      }
    }
  }
  
  
  
  for (i in 1:length(functions_groups_df[,1])){
    ordered=order(functions_groups_df[i,2:length(functions_groups_df[1,])],decreasing=T) 
    max_group=ordered[1]
    val=functions_groups_df[i,max_group+1]
    if (val!=0){
      edited_resulting_df[i,19]=max_group
    }
    if (val==0){
      edited_resulting_df[i,19]=0
    }
  }
  names_ed=colnames(edited_resulting_df)
  names_ed[19]="Function_class"
  colnames(edited_resulting_df)=names_ed
  df_toplot$Function_class=edited_resulting_df[,19]
  
  group_names=vector()
  for (i in 1:length(GROUPS)){
    gname=GROUPS[[i]]
    complete_gname=paste(gname, collapse="/")
    group_names=c(group_names,complete_gname)
  }
  
  
  group_names=c("others",group_names)
  sizeord=order(df_toplot[,5], decreasing = T) ## according to gene freq
  ordbysize_df_toplot=df_toplot[sizeord,]
  
  
  ordbysize_df_toplot2=data.frame(ordbysize_df_toplot, stringsAsFactors = F)
  for (i in 1:length(ordbysize_df_toplot2[,1])){
    ordbysize_df_toplot2[i,6]=(ordbysize_df_toplot2[i,6]+1)
  }
  
  ordbysize_df_toplot2$Function_class=as.factor(ordbysize_df_toplot2$Function_class)
  ordbysize_df_toplot2$Function_class=factor(ordbysize_df_toplot2$Function_class,levels=c(1:16))
  ordbysize_df_toplot2_nother=ordbysize_df_toplot2[-which(ordbysize_df_toplot2[,6]==1),] 
  
  ### SUBSET FOR COMMON, UNIA, UNIB FUNCTIONS
  
  UNIA_edited_resulting_df=edited_resulting_df[which(edited_resulting_df$B_genesfreq==0),]
  UNIB_edited_resulting_df=edited_resulting_df[which(edited_resulting_df$A_genesfreq==0),]
  COMM_edited_resulting_df=edited_resulting_df[which(edited_resulting_df$A_genesfreq!=0&edited_resulting_df$B_genesfreq!=0),]
  
  UNIAordbysize_df_toplot2=ordbysize_df_toplot2[which(ordbysize_df_toplot2$Fold_setB==0),]
  UNIBordbysize_df_toplot2=ordbysize_df_toplot2[which(ordbysize_df_toplot2$Fold_setA==0),]
  COMMordbysize_df_toplot2=ordbysize_df_toplot2[which(ordbysize_df_toplot2$Fold_setA!=0&ordbysize_df_toplot2$Fold_setB!=0),]
  
  
  SUBSETS_ordbysize_df_toplot2_LIST=list(UNIAordbysize_df_toplot2,UNIBordbysize_df_toplot2,COMMordbysize_df_toplot2)
  SUBSETS_edited_resulting_df_LIST=list(UNIA_edited_resulting_df,UNIB_edited_resulting_df,COMM_edited_resulting_df)
  
  
  SUBSETS_functions_hist_toplotdf_LIST=list()
  for ( w in 1:length(SUBSETS_edited_resulting_df_LIST)){
    EACH_edited_resulting_df=SUBSETS_edited_resulting_df_LIST[[w]]
    functions_groups_df=as.data.frame(matrix(0, nrow=length(EACH_edited_resulting_df[,1]), ncol=length(GROUPS)+1)) # nrow=#functional groups to assign , # ncol= #functional classes defined above + 1
    colnames(functions_groups_df)[1]="RepDesc"
    
    for ( i in 1:length(EACH_edited_resulting_df[,1])){
      funcpast=paste(EACH_edited_resulting_df[i,1],EACH_edited_resulting_df[i,3],sep="/")
      functions_groups_df[i,1]=funcpast
    }
    for ( i in 1:length(functions_groups_df[,1])){
      fun_count=0
      splited_funct=strsplit(functions_groups_df[i,1], "/")
      for (j in 1:length(GROUPS)){
        for ( h in 1:length(splited_funct[[1]])){
          splited_funct_sent=strsplit(splited_funct[[1]][[h]], "/| ")
          fun_count=length(which(is.element(splited_funct_sent[[1]],GROUPS[[j]])))
          if (fun_count>=1){
            count=(functions_groups_df[i,j+1])+1
            functions_groups_df[i,j+1]=count
          }
        }
      }
    }
    functions_hist_toplotdf=data.frame("Group_names"=group_names, "Group"=c(1:16), "Total"=rep(0, 16),"Represented"=rep(0, 16))
    binary_functions_groups_df=functions_groups_df
    for ( i in 2:length(functions_groups_df[1,])){
      ind_to1=which(functions_groups_df[,i]>=1)
      binary_functions_groups_df[ind_to1,i]=1
    }
    EACHordbysize_df_toplot2=SUBSETS_ordbysize_df_toplot2_LIST[[w]]
    for ( i in 2:length(binary_functions_groups_df[1,])){
      functions_hist_toplotdf[i,3]=sum(binary_functions_groups_df[,i])
      functions_hist_toplotdf[i,4]=length(which(EACHordbysize_df_toplot2[,6]==i)) 
    }
    
    functions_hist_toplotdf[1,3]=length(which(EACHordbysize_df_toplot2[,6]==1)) 
    functions_hist_toplotdf[1,4]=length(which(EACHordbysize_df_toplot2[,6]==1)) 
    functions_hist_toplotdf$Group=as.factor(functions_hist_toplotdf$Group)
    functions_hist_toplotdf$Ratio=paste(functions_hist_toplotdf$Total, "/", functions_hist_toplotdf$Represented, sep = "")
    #  totnum=paste("Total functional groups",sum(functions_hist_toplotdf$Represented), sep=": ")
    SUBSETS_functions_hist_toplotdf_LIST[[w]]=functions_hist_toplotdf
    
  }
  names(SUBSETS_functions_hist_toplotdf_LIST)=c("uniqueA_set", "uniqueB_set", "common_set")
  SUBSETS_functions_hist_toplotdf_LIST
}


subsetting_fun_res=subsetting_summres_plots_func(GROUPS,rearrange_fun_res$edited_resulting_df, rearrange_fun_res$df_toplot)



#graphics.off()

## FOR GOT

# barplots counts only the most represented function for each GOgroup
wo_others_histograms_function=function( colors16,group_names,functions_hist_toplotdf){
  val10=25
  colorsx=colors16[2:16]
  group_namesx=group_names[1:15]
  functions_hist_toplotdf_A=functions_hist_toplotdf[[1]]
  functions_hist_toplotdf_noother=functions_hist_toplotdf_A[2:16,]
  max_fun=max(functions_hist_toplotdf_noother$Represented)
  totnum=paste("Unique A GO groups",sum(functions_hist_toplotdf_noother$Represented), sep=": ")
  for ( i in 1:length(functions_hist_toplotdf_noother[,1])){
    functions_hist_toplotdf_noother[i,2]=as.numeric(functions_hist_toplotdf_noother[i,2])-1
  }
  functions_hist_toplotdf_noother$Group=as.factor(functions_hist_toplotdf_noother$Group)
  
  UNIAhist_plot=ggplot(data=functions_hist_toplotdf_noother) + 
    geom_histogram(aes(x=Group,y=Represented), stat = "identity", fill=colorsx, color="black")+
    scale_fill_manual(values=colorsx, 
                      name="Function classes" ,
                      labels=c(group_namesx), guide=F)+
    annotate("text",label=totnum, x=16/2, y=10,size=10)+
    xlab("GO classes")+
    ylab("#GO groups")+
    scale_y_continuous(limits=c(0,15))+
    theme_base()+
    theme(legend.position = "none")+  
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "white", fill=NA, size=1),
          axis.title = element_text(size = val10),
          axis.text = element_text(size=val10-5))
  
  
  
  functions_hist_toplotdf_B=functions_hist_toplotdf[[2]]
  functions_hist_toplotdf_noother=functions_hist_toplotdf_B[2:16,]
  totnum=paste("Unique B GO groups",sum(functions_hist_toplotdf_noother$Represented), sep=": ")
  for ( i in 1:length(functions_hist_toplotdf_noother[,1])){
    functions_hist_toplotdf_noother[i,2]=as.numeric(functions_hist_toplotdf_noother[i,2])-1
  }
  functions_hist_toplotdf_noother$Group=as.factor(functions_hist_toplotdf_noother$Group)
  
  max_fun=max(functions_hist_toplotdf_noother$Represented)
  
  
  UNIBhist_plot=ggplot(data=functions_hist_toplotdf_noother) + 
    geom_histogram(aes(x=Group,y=Represented), stat = "identity", fill=colorsx, color="black")+
    scale_fill_manual(values=colorsx, 
                      name="Function classes" ,
                      labels=c(group_namesx), guide=F)+
    annotate("text",label=totnum, x=16/2, y=10,size=10)+
    xlab("GO classes")+
    ylab("#GO groups")+
    scale_y_continuous(limits=c(0,15))+
    theme_base()+
    theme(legend.position = "none")+  
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "white", fill=NA, size=1),
          axis.title = element_text(size = val10),
          axis.text = element_text(size=val10-5))
  
  functions_hist_toplotdf_C=functions_hist_toplotdf[[3]]
  
  functions_hist_toplotdf_noother=functions_hist_toplotdf_C[2:16,]
  totnum=paste("Common GO groups",sum(functions_hist_toplotdf_noother$Represented), sep=": ")
  for ( i in 1:length(functions_hist_toplotdf_noother[,1])){
    functions_hist_toplotdf_noother[i,2]=as.numeric(functions_hist_toplotdf_noother[i,2])-1
  }
  functions_hist_toplotdf_noother$Group=as.factor(functions_hist_toplotdf_noother$Group)
  
  
  
  max_fun=max(functions_hist_toplotdf_noother$Represented)
  
  COMMhist_plot=ggplot(data=functions_hist_toplotdf_noother) + 
    geom_histogram(aes(x=Group,y=Represented), stat = "identity", fill=colorsx, color="black")+
    scale_fill_manual(values=colorsx, 
                      name="Function classes" ,
                      labels=c(group_namesx), guide=F)+
    annotate("text",label=totnum, x=16/2, y=10,size=10)+
    xlab("GO classes")+
    ylab("#GO groups")+
    scale_y_continuous(limits=c(0,15))+
    theme_base()+
    theme(legend.position = "none")+  
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "white", fill=NA, size=1),
          axis.title = element_text(size = val10),
          axis.text = element_text(size=val10-5))
  
  ALLhist_plot=grid.arrange(UNIAhist_plot,COMMhist_plot,UNIBhist_plot,nrow=3)
  ALLhist_plot
}


#bar plot counts all functions (not only the most represented one) and mantains the count for the #groups

total_wo_others_histograms_function=function(colors16,group_names,functions_hist_toplotdf){
  val10=25
  colorsx=colors16[2:16]
  group_namesx=group_names[1:15]
  functions_hist_toplotdf_A=functions_hist_toplotdf[[1]]
  functions_hist_toplotdf_noother=functions_hist_toplotdf_A[2:16,]
  
  max_fun=max(functions_hist_toplotdf_noother$Total)
  totnum=paste("Unique A GO groups",sum(functions_hist_toplotdf_noother$Represented), sep=": ")
  for ( i in 1:length(functions_hist_toplotdf_noother[,1])){
    functions_hist_toplotdf_noother[i,2]=as.numeric(functions_hist_toplotdf_noother[i,2])-1
  }
  functions_hist_toplotdf_noother$Group=factor(functions_hist_toplotdf_noother$Group, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))
  
  UNIAhist_plot=ggplot(data=functions_hist_toplotdf_noother) + 
    geom_histogram(aes(x=Group,y=Total), stat = "identity", fill=colorsx, color="black")+
    scale_fill_manual(values=colorsx, 
                      name="Function classes" ,
                      labels=c(group_namesx), guide=F)+
    annotate("text",label=totnum, x=16/2, y=10,size=10)+
    xlab("GO classes")+
    ylab("#GO groups")+
    scale_y_continuous(limits=c(0,max_fun))+
    theme_base()+
    theme(legend.position = "none")+  
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "white", fill=NA, size=1),
          axis.title = element_text(size = val10),
          axis.text = element_text(size=val10-5))
  
  
  
  functions_hist_toplotdf_B=functions_hist_toplotdf[[2]]
  functions_hist_toplotdf_noother=functions_hist_toplotdf_B[2:16,]
  totnum=paste("Unique B GO groups",sum(functions_hist_toplotdf_noother$Represented), sep=": ")
  for ( i in 1:length(functions_hist_toplotdf_noother[,1])){
    functions_hist_toplotdf_noother[i,2]=as.numeric(functions_hist_toplotdf_noother[i,2])-1
  }
  functions_hist_toplotdf_noother$Group=factor(functions_hist_toplotdf_noother$Group, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))
  
  max_fun=max(functions_hist_toplotdf_noother$Total)
  
  
  UNIBhist_plot=ggplot(data=functions_hist_toplotdf_noother) + 
    geom_histogram(aes(x=Group,y=Total), stat = "identity", fill=colorsx, color="black")+
    scale_fill_manual(values=colorsx, 
                      name="Function classes" ,
                      labels=c(group_namesx), guide=F)+
    annotate("text",label=totnum, x=16/2, y=10,size=10)+
    xlab("GO classes")+
    ylab("#GO groups")+
    scale_y_continuous(limits=c(0,max_fun))+
    theme_base()+
    theme(legend.position = "none")+  
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "white", fill=NA, size=1),
          axis.title = element_text(size = val10),
          axis.text = element_text(size=val10-5))
  
  functions_hist_toplotdf_C=functions_hist_toplotdf[[3]]
  
  functions_hist_toplotdf_noother=functions_hist_toplotdf_C[2:16,]
  totnum=paste("Common GO groups",sum(functions_hist_toplotdf_noother$Represented), sep=": ")
  for ( i in 1:length(functions_hist_toplotdf_noother[,1])){
    functions_hist_toplotdf_noother[i,2]=as.numeric(functions_hist_toplotdf_noother[i,2])-1
  }
  functions_hist_toplotdf_noother$Group=factor(functions_hist_toplotdf_noother$Group, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))
  
  
  
  max_fun=max(functions_hist_toplotdf_noother$Total)
  
  COMMhist_plot=ggplot(data=functions_hist_toplotdf_noother) + 
    geom_histogram(aes(x=Group,y=Total), stat = "identity", fill=colorsx, color="black")+
    scale_fill_manual(values=colorsx, 
                      name="Function classes" ,
                      labels=c(group_namesx), guide=F)+
    annotate("text",label=totnum, x=16/2, y=10,size=10)+
    xlab("GO classes")+
    ylab("#GO groups")+
    scale_y_continuous(limits=c(0,max_fun))+
    theme_base()+
    theme(legend.position = "none")+  
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "white", fill=NA, size=1),
          axis.title = element_text(size = val10),
          axis.text = element_text(size=val10-5))
  
  ALLhist_plot=grid.arrange(UNIAhist_plot,COMMhist_plot,UNIBhist_plot,nrow=3)
  ALLhist_plot
}





group_names=vector()
for (i in 1:length(GROUPS)){
  gname=GROUPS[[i]]
  complete_gname=paste(gname, collapse="/")
  group_names=c(group_names,complete_gname)
}

test_wo_others_histograms_function=wo_others_histograms_function(colors16,group_names,subsetting_fun_res)





#ggsave(filename="UNIAhist_plot.pdf",plot=UNIAhist_plot, device="pdf",limitsize=F, scale=1)
#ggsave(filename="UNIBhist_plot.pdf",plot=UNIBhist_plot, device="pdf",limitsize=F, scale=1)
#ggsave(filename="COMMhist_plot.pdf",plot=COMMhist_plot, device="pdf",limitsize=F, scale=1)

#ggsave(filename="ALLhist_plot.pdf",plot=ALLhist_plot, device="pdf",limitsize=F, scale=2)


functions_hist_toplotdf_A=functions_hist_toplotdf[[1]]
functions_hist_toplotdf_noother=functions_hist_toplotdf_A[2:16,]


COMMhist_legend=ggplot(data=functions_hist_toplotdf_noother) + 
  geom_histogram(aes(x=Group,y=Represented), stat = "identity", fill=colors15, color="black")+
  scale_fill_manual(values=colors15, 
                    name="Function classes" ,
                    labels=c(group_names15))+
  annotate("text",label=totnum, x=16/2, y=10,size=10)+
  xlab("GO classes")+
  ylab("#GO groups")+
  scale_y_continuous(limits=c(0,15))+
  theme_base()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "white", fill=NA, size=1),
        axis.title = element_text(size = val10),
        axis.text = element_text(size=val10-5))




g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(test_wo_others_histograms_function)

ggsave(filename="simp_res_plot2_LEGEND.pdf",plot=mylegend, device="pdf",limitsize=F, scale=1.5) ## plot only legend





#ggsave(filename="UNIBhist_plot.pdf",plot=UNIBhist_plot, device="pdf",limitsize=F, scale=0.5)



##### PLOT BOXES 

clean_res=SIMPLIFICATION_OUTPUT_clean_res[,c(2,8,9)]

for ( i in 1:length(clean_res[,1])){
  vec_prot=unlist(strsplit(clean_res[i,3],split="/"))
  vec_prot_entrez=mapIds(org.Hs.eg.db, keys=(as.vector(vec_prot)),column="SYMBOL", keytype = "ENTREZID", multivals="first")
  vec_prot=paste(vec_prot_entrez,sep="/",collapse="/")
  clean_res[i,3]=vec_prot
}


vec_desc_size=vector()
vec_gene_size=vector()
for ( i in 1:length(clean_res[,1])){
  desc_count=length(unlist(strsplit(clean_res[i,2],split="/")))
  vec_desc_size=c(vec_desc_size,desc_count)
  gene_count=length(unique(unlist(strsplit(clean_res[i,3],split="/"))))
  vec_gene_size=c(vec_gene_size,gene_count)
}  

vec_desc_size
vec_gene_size

for (i in 1:length(clean_res[,1])){
  vec=unlist(strsplit(clean_res[i,2], split="/"))
  vec2=vec[2:length(vec)]
  clean_res[i,2]=paste(vec2, sep="/", collapse = "/")
}

clean_res2=cbind(clean_res,rep(0))

for (i in 1:length(clean_res2[,1])){
  clean_res2[i,4]=paste("group", i, sep=" " )
}


poem_title2 <- data.frame(
  text = c(clean_res2[,1]),
  xmin = 0,
  xmax = 100,
  ymin = 85,
  ymax = 100,
  fit = c( clean_res2[,4])
)

poem_desc2 <- data.frame(
  text = c(clean_res2[,2]),
  xmin = 0,
  xmax = 100,
  ymin = 40,
  ymax = 80,
  fit = c( clean_res2[,4])
)


poem_gene2 <- data.frame(
  text = c(clean_res2[,3]),
  xmin = 0,
  xmax = 100,
  ymin = 0,
  ymax = 30,
  fit = c( clean_res2[,4])
)


mytitle2="FFLLentrez_enrichGOMFdf_010707"
res_desc=desc_plotAllLayers(clean_res2,poem_title2,poem_desc2,poem_gene2, mytitle2)


desc_plotAllLayers<-function(df0,df1,df2,df3,mytittle){
  p<-ggplot(df0, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                     label = text)) 
  for(h in 1:length(df0[,1])){ 
    p=p+geom_fit_text(data = subset(df1, fit == df0[h,4]),aes(fontface="bold",family="Times"), reflow = TRUE,
                      min.size=0, size=30) 
  }
  
  for(i in 1:length(df0[,1])){ 
    p=p+geom_fit_text(data = subset(df2, fit == df0[i,4]),aes(family="Times"), reflow = TRUE,
                      min.size = 0, size=20) 
  }
  for(j in 1:length(df0[,1])){ 
    p=p+geom_fit_text(data = subset(df3, fit == df0[j,4]),aes(family="Times"), reflow = TRUE,
                      min.size = 0) 
  }
  p=p+ lims(x = c(0, 100), y = c(0, 100)) 
  p=p+ labs(x = "", y = "") 
  p=p+ facet_wrap(~ fit, ncol=5)
  p=p+ theme_igray()
  p=p+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
             axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())
  
  p=p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p=p+theme(legend.position="none")
  p=p+labs(title=mytittle)
  return(p)
  
}













###################### FUNCTIONS REQUIRED FOR MAIN FUNCTIONS ABOVE

single_funct_simplif=function(RAW_FUNC_ENRICH_DATAFRAME_A,bg_threshold, gene_cooc_threshold, semantic_alg, semantic_threshold){
  filtering_A_enrichGO=FILTERING_BACKGROUND(RAW_FUNC_ENRICH_DATAFRAME_A,bg_threshold) 
  withfold_filtering_A_enrichGO=fold_computing(filtering_A_enrichGO[[1]])
  geneco_AB_func_enrichment=GENECO_MERGING_GOT(withfold_filtering_A_enrichGO,0,length(withfold_filtering_A_enrichGO[,1]),0, gene_cooc_threshold, gene_cooc_threshold)
  sem_AB_func_enrichment=SEMANTIC_ANALYSIS2(semantic_alg,semantic_threshold,geneco_AB_func_enrichment)
  final_simplified_res=sem_AB_func_enrichment[[1]]
  sem_merge_summary=sem_AB_func_enrichment[[2]]
  list(withfold_filtering_A_enrichGO,geneco_AB_func_enrichment,sem_merge_summary,final_simplified_res)
}


pluralwithoutsubset_funct_simplif=function(RAW_FUNC_ENRICH_DATAFRAME_A, RAW_FUNC_ENRICH_DATAFRAME_B,bg_threshold, gene_cooc_threshold, semantic_alg, semantic_threshold){
  filtering_A_enrichGO=FILTERING_BACKGROUND(RAW_FUNC_ENRICH_DATAFRAME_A,bg_threshold) 
  filtering_B_enrichGO=FILTERING_BACKGROUND(RAW_FUNC_ENRICH_DATAFRAME_B,bg_threshold) 
  AB_func_enrichment_res=JOINING_COMPARING_DATASETS(filtering_A_enrichGO[[1]],filtering_B_enrichGO[[1]])
  geneco_AB_func_enrichment=GENECO_MERGING_GOT(AB_func_enrichment_res[[1]],AB_func_enrichment_res[[2]],AB_func_enrichment_res[[3]],AB_func_enrichment_res[[4]], gene_cooc_threshold, gene_cooc_threshold)
  sem_AB_func_enrichment=SEMANTIC_ANALYSIS2(semantic_alg,semantic_threshold,geneco_AB_func_enrichment)
  final_simplified_res=sem_AB_func_enrichment[[1]]
  sem_merge_summary=sem_AB_func_enrichment[[2]]
  result_simp=list(AB_func_enrichment_res,geneco_AB_func_enrichment,sem_merge_summary,final_simplified_res)
}


#########

EDITING_TO_PLOT_noS2B=function(Ageneset_size, Bgeneset_size,merged_fussioning_result_df){
  all_merged_fussioning_result_df=data.frame()
  all_merged_fussioning_result_df=merged_fussioning_result_df
  for (i in 1:length(merged_fussioning_result_df[,1])){
    for (z in 1:12){
      if (z<=6){
        p_vec=paste(merged_fussioning_result_df[i,z],merged_fussioning_result_df[i,z+6],sep="/",collapse = "/")
        all_merged_fussioning_result_df[i,z+6]=p_vec
      }
    }
  }
  comp_merged_fussioning_result_df=data.frame()
  comp_merged_fussioning_result_df=all_merged_fussioning_result_df
  del=c("0","1")
  empty=" "
  for (i in 1:length(all_merged_fussioning_result_df[,9])){
    genes_a=strsplit(all_merged_fussioning_result_df[i,9],"/")
    genes_a=setdiff(genes_a[[1]],del)
    if (length(genes_a)!=0){
      genes_f=setdiff(genes_a[[1]],empty)
      f=strsplit(genes_f,"")
      if (length(f[[1]])==0){
        genes_a=vector()
      }
    }
    genes_a_freq=(length(genes_a)/Ageneset_size)*100
    comp_merged_fussioning_result_df[i,13]=genes_a_freq
    
    genes_b=strsplit(all_merged_fussioning_result_df[i,10],"/")
    genes_b=setdiff(genes_b[[1]],del)
    if (length(genes_b)!=0){
      genes_ff=setdiff(genes_b[[1]],empty)
      ff=strsplit(genes_ff,"")
      if (length(ff[[1]])==0){
        genes_b=vector()
      }
    }
    genes_b_freq=(length(genes_b)/Bgeneset_size)*100
    comp_merged_fussioning_result_df[i,14]=genes_b_freq
    
    fold_a=strsplit(all_merged_fussioning_result_df[i,11],"/")
    fold_a=setdiff(fold_a[[1]],del)
    if(length(fold_a)>=5){
      median_a5fold=median(as.numeric(fold_a[1:5]))
      comp_merged_fussioning_result_df[i,15]=median_a5fold
    }else{
      median_afold=median(as.numeric(fold_a))
      comp_merged_fussioning_result_df[i,15]=median_afold
    }
    
    fold_b=strsplit(all_merged_fussioning_result_df[i,12],"/")
    fold_b=setdiff(fold_b[[1]],del)
    if(length(fold_b)>=5){
      median_b5fold=median(as.numeric(fold_b[1:5]))
      comp_merged_fussioning_result_df[i,16]=median_b5fold
    }else{
      median_bfold=median(as.numeric(fold_b))
      comp_merged_fussioning_result_df[i,16]=median_bfold
    }
  }
  colnames(comp_merged_fussioning_result_df)=c("GO_Terms_ID","Description","Original_Genes_A","Original_Genes_B" , "Original_Fold_A" , "Original_Fold_B" ,  "Merged_GO_T" ,"Merged_desc_T","Merged_A_genes" ,  "Merged_B_genes"  ,  "Merged_A_fold"  ,  "Merged_B_fold","A_genesfreq","B_genesfreq","A_medianfold","B_medianfold")
  comp_merged_fussioning_result_df
}
########
FILTERING_BACKGROUND=function(enrichGOdf,bg_threshold) {
  bg_enrichGOdf=data.frame()
  for (i in 1:length(enrichGOdf[,1])){
    a=strsplit(as.character(enrichGOdf[i,4]),"/")
    v=as.numeric(a[[1]][1])/as.numeric(a[[1]][2])
    enrichGOdf[i,10]=v
  }
  bg_enrichGOdf=enrichGOdf
  colnames(bg_enrichGOdf)=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Background")
  x_bg_index=which(enrichGOdf[,10]<=bg_threshold)
  filteredbg_enrichGO=enrichGOdf[x_bg_index,]
  filtered_GOTs=as.vector(filteredbg_enrichGO[,1])
  list(filteredbg_enrichGO, filtered_GOTs,bg_enrichGOdf) # list with 3 outputs (ALREADY FILTERED COMPLETE RESULT, ONLY GOT, ORIGINAL FUNCTIONAL ENRICHMENT WIHT BACKGROUND COMPUTED )
}

######################

JOINING_COMPARING_DATASETS=function(funct_enrA, funct_enrB){
  AB_colnames_vector=c("GO_T", "desc" ,"Genes_A","Genes_B", "Fold_A", "Fold_B","ad.p-value_A", "ad.p-value_B" )
  A_unique_colnames_vector=c("GO_T","desc","Genes_A","Fold_A","ad.p-value_A")
  B_unique_colnames_vector=c("GO_T","desc","Genes_B","Fold_B","ad.p-value_B")
  
  AB_common_terms=common_terms(funct_enrA,funct_enrB,AB_colnames_vector)
  A_unique_terms_againstB=unique_terms(funct_enrA,funct_enrB,A_unique_colnames_vector)
  B_unique_terms_againstA=unique_terms(funct_enrB,funct_enrA,B_unique_colnames_vector)
  
  ed_uniA=cbind("Go"=A_unique_terms_againstB[,1],"Desc"=A_unique_terms_againstB[,2],
                "Genes_A"=as.character(A_unique_terms_againstB[,3]),"Genes_B"=rep(0,length(A_unique_terms_againstB[,1])),
                "Fold_A"=A_unique_terms_againstB[,4],
                "Fold_B"=rep(0,length(A_unique_terms_againstB[,1])))
  
  ed_uniB=cbind("Go"=B_unique_terms_againstA[,1],"Desc"=B_unique_terms_againstA[,2],
                "Genes_A"=rep(0,length(B_unique_terms_againstA[,3])),"Genes_B"=as.character(B_unique_terms_againstA[,3]),
                "Fold_A"=rep(0,length(B_unique_terms_againstA[,3])),
                "Fold_B"=B_unique_terms_againstA[,4])
  
  ed_AB=cbind("Go"=AB_common_terms[,1],"Desc"=AB_common_terms[,2],
              "Genes_A"=as.character(AB_common_terms[,3]),"Genes_B"=as.character(AB_common_terms[,4]),
              "Fold_A"=AB_common_terms[,5],
              "Fold_B"=AB_common_terms[,6])
  
  AB_joined_data=rbind(ed_AB,ed_uniA,ed_uniB)
  AB_joined_data=data.frame(AB_joined_data,stringsAsFactors = F)
  common_size=length(ed_AB[,1])
  uniA_size=length(ed_uniA[,1])
  uniB_size=length(ed_uniB[,1])
  list(AB_joined_data, common_size, uniA_size, uniB_size)
}

######################

GENECO_MERGING_GOT=function(ABfun_enr,ABcom_size, Auni_size, Buni_size, gene_co_threshA, gene_co_threshB){
  total=ABcom_size+Auni_size+Buni_size 
  doubleset=c(rep(1,ABcom_size),rep(0,total-ABcom_size))  
  maxscore=pmax(ABfun_enr[,5],ABfun_enr[,6]) #by fold score
  maxscore=as.numeric(maxscore)
  reorder=order(doubleset,maxscore,decreasing=T)
  ordered_ABfun_enr=data.frame()
  ordered_ABfun_enr=ABfun_enr[reorder,]
  
  gcmerged_ABfun_enr_res=GOT_allfussioning(ordered_ABfun_enr,gene_co_threshA,gene_co_threshB)
  gcmerged_ABfun_enr=gcmerged_ABfun_enr_res[[2]]
  gcmerged_ABfun_enr
}

######################

SEMANTIC_ANALYSIS2=function(alg_measure,sim_thresh,initial_df){
  AB_GOT=initial_df[,1]
  AB_GOT=c(as.character(AB_GOT))
  AB_testing=termSim(AB_GOT, AB_GOT,semData=d,method = "Lin")
  tst_values=sm2vec(AB_testing, diag = FALSE)
  tst_index=sm.index(AB_testing, diag = F)
  order_tst_values=order(tst_values,decreasing = T)
  sorted_tst_values=tst_values[order_tst_values]
  sorted_tst_index=tst_index[order_tst_values,]
  e_df=rbind(sorted_tst_index[sorted_tst_values>=0.7,])
  
  if (length(e_df[,1])!=0 & (all(is.na(e_df[,1]))&  all(is.na(e_df[,2]))==F)){
    e_df_values=c(unlist(e_df[,1:2]))
    e_df_summary=cbind(1:max(e_df_values),rep(0,max(e_df_values)),rep(0,max(e_df_values)))
    colnames(e_df_summary)=c("GOT_ind_list", "received_fused","deleted")
    minitial_df=initial_df
    for (i in 1:length(e_df[,1])){
      indexX=e_df[i,2]
      indexY=e_df[i,1]
      e_df_summary[indexX,2]=e_df[i,1]
      e_df_summary[indexY,3]=1
      if (e_df_summary[indexX,3]==0){
        minitial_df[indexX,7]=paste(initial_df[indexY,1],initial_df[indexX,7],sep="/")
        minitial_df[indexX,8]=paste(initial_df[indexY,2],initial_df[indexX,8],sep="/")
        minitial_df[indexX,3]=paste(initial_df[indexY,3],initial_df[indexX,9],sep="/")
        minitial_df[indexX,4]=paste(initial_df[indexY,4],initial_df[indexX,10],sep="/")
        minitial_df[indexX,5]=paste(initial_df[indexY,5],initial_df[indexX,11],sep="/")
        minitial_df[indexX,6]=paste(initial_df[indexY,6],initial_df[indexX,12],sep="/")
      }
    }
    del_ind=which(e_df_summary[,3]==1) #####
    minitial_df=minitial_df[-del_ind,] #####
    list(minitial_df,e_df_summary)
  }
  else {
    minitial_df=initial_df
    e_df_summary="there is not enought semantic similarity"
    list(minitial_df,e_df_summary)
  }
}

######################


common_terms=function(dfA,dfB,df_colnames_vector){
  res_common_termsdf=data.frame(stringsAsFactors = F)
  paralell_index=0
  for (i in 1:length(dfA[,1])) {
    if (is.element(dfA[i,1],dfB[,1])) {
      bbb_index=which(dfB[,1]==dfA[i,1])
      paralell_index=paralell_index+1
      res_common_termsdf[paralell_index,1]=dfA[i,1]
      res_common_termsdf[paralell_index,2]=dfA[i,2]
      res_common_termsdf[paralell_index,3]=dfB[bbb_index,8]
      res_common_termsdf[paralell_index,4]=dfA[i,8]
      
      a_ftss=strsplit(as.character(dfB[bbb_index,3]),"/") #generatio
      v_ftss=as.numeric(a_ftss[[1]][1])/as.numeric(a_ftss[[1]][2]) #generatio divided
      b_ftss=strsplit(as.character(dfB[bbb_index,4]),"/")
      w_ftss=as.numeric(b_ftss[[1]][1])/as.numeric(b_ftss[[1]][2])
      ab_ftss=v_ftss/w_ftss
      res_common_termsdf[paralell_index,5]=ab_ftss
      
      a_s2b=strsplit(as.character(dfA[i,3]),"/") #generatio
      v_s2b=as.numeric(a_s2b[[1]][1])/as.numeric(a_s2b[[1]][2]) #generatio divided
      b_s2b=strsplit(as.character(dfA[i,4]),"/")
      w_s2b=as.numeric(b_s2b[[1]][1])/as.numeric(b_s2b[[1]][2])
      ab_s2b=v_s2b/w_s2b
      res_common_termsdf[paralell_index,6]=ab_s2b
      
      res_common_termsdf[paralell_index,7]=dfB[bbb_index,5]
      res_common_termsdf[paralell_index,8]=dfA[i,5]
    }
  }
  colnames(res_common_termsdf)=df_colnames_vector
  res_common_termsdf
}


######################

unique_terms=function(dfA,dfB,df_colnames_vector){
  res_unique_termsdf=data.frame(stringsAsFactors = F)
  paralell_index=0
  for (i in 1:length(dfA[,1])) {
    if (!is.element(dfA[i,1],dfB[,1])) {
      paralell_index=paralell_index+1
      res_unique_termsdf[paralell_index,1]=dfA[i,1]
      res_unique_termsdf[paralell_index,2]=dfA[i,2]
      res_unique_termsdf[paralell_index,3]=dfA[i,8]
      a_s2b=strsplit(as.character(dfA[i,3]),"/") #generatio
      v_s2b=as.numeric(a_s2b[[1]][1])/as.numeric(a_s2b[[1]][2]) #generatio divided
      b_s2b=strsplit(as.character(dfA[i,4]),"/")
      w_s2b=as.numeric(b_s2b[[1]][1])/as.numeric(b_s2b[[1]][2])
      ab_s2b=v_s2b/w_s2b
      res_unique_termsdf[paralell_index,4]=ab_s2b
      res_unique_termsdf[paralell_index,5]=dfA[i,5]
    }
  }
  colnames(res_unique_termsdf)=df_colnames_vector
  res_unique_termsdf
}


######################

GOT_allfussioning=function(df,ftss_coin_min,S2B_coin_min) {
  # order dataframe according to the best parameter
  fussioned_GOT_df=as.data.frame(df)
  fussioned_GOT_df[,7]=c(1)
  fussioned_GOT_df[,8]=c(1)
  fussioned_GOT_df[,9]=c(1)
  fussioned_GOT_df[,10]=c(1)
  fussioned_GOT_df[,11]=c(1)#
  fussioned_GOT_df[,12]=c(1)#
  n=length(df[,1])
  remove=vector()
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if ((!is.element(i,remove))&(!is.element(j,remove))){
        genes_ftss1=as.vector(df[i,3]) #(vec1)
        genes_ftss2=as.vector(df[j,3]) #(vec2)
        ftss_check=GOTcoincidence_analysischeck(genes_ftss1,genes_ftss2) 
        GOTfun_ftss=GOTcoincidence_analysis(genes_ftss1,genes_ftss2) #maxsize
        
        genes_S2B1=as.vector(df[i,4]) #(vec1)
        genes_S2B2=as.vector(df[j,4]) #(vec2)
        S2B_check=GOTcoincidence_analysischeck(genes_S2B1,genes_S2B2)
        GOTfun_S2B=GOTcoincidence_analysis(genes_S2B1,genes_S2B2)
        
        a1=(ftss_check==0)#
        a2=(S2B_check==0)#
        
        c1=(GOTfun_ftss>=ftss_coin_min)
        c2=(GOTfun_S2B>=S2B_coin_min)
        
        allcheck=((c1&c2)|(a1&c2)|(a2&c1))
        
        if (allcheck){   ######### for coincidences
          fussioned_GO_vec=vector()
          fussioned_desc_vec=vector()
          genes=vector()
          if (fussioned_GOT_df[i,7]==1){ ########################
            fussioned_GOT_df[i,7]=df[j,1]
            fussioned_GOT_df[i,8]=df[j,2]
            fussioned_GOT_df[i,9]=df[j,3]
            fussioned_GOT_df[i,10]=df[j,4]
            #genes=paste(df[j,3],df[j,4],sep="/")
            fussioned_GOT_df[i,11]=df[j,5]
            fussioned_GOT_df[i,12]=df[j,6]#
            
          } else {
            fussioned_GO_vec=paste(fussioned_GOT_df[i,7],df[j,1],sep="/")
            fussioned_desc_vec=paste(fussioned_GOT_df[i,8],df[j,2],sep="/")
            fussioned_gene_ftss_vec=paste(fussioned_GOT_df[i,9],df[j,3],sep="/")
            fussioned_gene_S2B_vec=paste(fussioned_GOT_df[i,10],df[j,4],sep="/")
            fussioned_fold_ftss_vec=paste(fussioned_GOT_df[i,11],df[j,5],sep="/")
            fussioned_fold_S2B_vec=paste(fussioned_GOT_df[i,12],df[j,6],sep="/")
            #genes=paste(df[j,3],df[j,4],sep="/")
            fussioned_GOT_df[i,7]=fussioned_GO_vec
            fussioned_GOT_df[i,8]=fussioned_desc_vec
            fussioned_GOT_df[i,9]=fussioned_gene_ftss_vec
            fussioned_GOT_df[i,10]=fussioned_gene_S2B_vec
            fussioned_GOT_df[i,11]=fussioned_fold_ftss_vec
            fussioned_GOT_df[i,12]=fussioned_fold_S2B_vec
            
            
          }
          remove=c(remove, j)
        }
      }
    }
    
  }
  real_indexi=setdiff((1:n), remove) 
  colnames(fussioned_GOT_df)=c("GO_Terms_ID","Description","Original_Genes_A","Original_Genes_B","Original_Fold_A","Original_Fold_B","Merged_GOT","Merged_desc","Merged_A_genes","Merged_B_genes","Merged_A_fold","Merged_B_fold")
  depurated_fussioned_GOT_df=fussioned_GOT_df[real_indexi,]
  for (i in 1:length(depurated_fussioned_GOT_df[,1])){
    depurated_fussioned_GOT_df[i,9]=GOTrepeat_gene_elim(depurated_fussioned_GOT_df[i,9])
    depurated_fussioned_GOT_df[i,10]=GOTrepeat_gene_elim(depurated_fussioned_GOT_df[i,10])
  }
  list(fussioned_GOT_df,depurated_fussioned_GOT_df)
  
}

######################

GOTcoincidence_analysischeck=function(vec1,vec2){
  checkval=1
  if ((vec1=="0")|(vec2=="0")){
    checkval=0
  } 
  checkval
}

######################

GOTcoincidence_analysis=function(vec1,vec2){
  genes1=strsplit(as.character(vec1),"/")
  genes2=strsplit(as.character(vec2),"/")
  inter=intersect(genes1[[1]],genes2[[1]])
  sizes=c(length(genes1[[1]]),length(genes2[[1]]))
  if (min(sizes)==0){
    coincid=0
  } else {
    coincid=length(inter)*100/min(sizes)
  }
}

######################

GOTrepeat_gene_elim=function(vec1){
  genes=strsplit(as.character(vec1),"/")
  unigenes=unique(genes[[1]])
  unigenes=setdiff(unigenes,c("0"))
  vec=paste(unigenes,sep="/",collapse="/")
}

######################





FILTERING_BACKGROUND2=function(enrichGOdf,bg_threshold) {
  bg_enrichGOdf=cbind(enrichGOdf, "Background"=0)
  for (i in 1:length(bg_enrichGOdf[,1])){
    a=strsplit(as.character(bg_enrichGOdf[i,4]),"/")
    v=as.numeric(a[[1]][1])/as.numeric(a[[1]][2])
    bg_enrichGOdf[i,10]=v
  }
  x_bg_index=which(bg_enrichGOdf[,10]<=bg_threshold)
  bg_enrichGOdf=bg_enrichGOdf[x_bg_index,]
  bg_enrichGOdf
}



######################
fold_computing2=function(tes_B_func_enrichment){
  tes2_B_func_enrichment=data.frame(matrix(ncol=6))
  colnames(tes2_B_func_enrichment)=c("GO_T", "desc" ,"Genes_Aa","Genes_Ab", "Fold_Aa", "Fold_Ab")
  for (i in 1:length(tes_B_func_enrichment[,1])){
    a_ftss=unlist(strsplit(as.character(tes_B_func_enrichment[i,3]),"/")) #generatio
    b_ftss=unlist(strsplit(as.character(tes_B_func_enrichment[i,4]),"/"))
    ab_ftss=as.numeric(a_ftss[1])/as.numeric(a_ftss[2])/as.numeric(b_ftss[1])/as.numeric(b_ftss[2])     #generatio divided 
    tes2_B_func_enrichment[i,]=c(tes_B_func_enrichment[i,1],tes_B_func_enrichment[i,2],tes_B_func_enrichment[i,8], "0", ab_ftss, "0")
  }
  tes2_B_func_enrichment
}
################

FEA_and_FILTERING=function(ORIGINAL_gene_vec, ORIGINAL_id, ann, bg_threshold, OUTPUT_id){
  ENTREZ_LIST_A=mapIds(org.Hs.eg.db, keys=(ORIGINAL_gene_vec),column="ENTREZID", keytype = ORIGINAL_id, multivals="first")
  ENTREZ_LIST_A_DF=data.frame(ORIGINAL_gene_vec, ENTREZ_LIST_A, stringsAsFactors = F)
  if (ann!="KEGG"){
    Aentrezgenes_enrichres=as.data.frame(enrichGO(ENTREZ_LIST_A_DF[,2],'org.Hs.eg.db',ont=ann))
    filtering_A_enrichGO=FILTERING_BACKGROUND2(Aentrezgenes_enrichres,0.7) 
    filtering_A_enrichGO_res=fold_computing2(filtering_A_enrichGO)
    num_out=length(Aentrezgenes_enrichres[,1])-length(filtering_A_enrichGO_res[,1])
    num_notmapp=length(which(is.na(ENTREZ_LIST_A_DF[,2])))
    filtering_A_enrichGO_res2=cbind(filtering_A_enrichGO_res,"genes"=0)
    for ( i in 1:length(filtering_A_enrichGO_res2[,1])){
      vec_prot=unlist(strsplit(filtering_A_enrichGO_res2[i,3],split="/"))
      vec_prot_entrez=mapIds(org.Hs.eg.db, keys=(as.vector(vec_prot)),column=OUTPUT_id, keytype = "ENTREZID", multivals="first")
      vec_prot=paste(vec_prot_entrez,sep="/",collapse="/")
      filtering_A_enrichGO_res2[i,7]=vec_prot
    }
  }

  if (ann=="KEGG"){
    Aentrezgenes_enrichres=as.data.frame(enrichKEGG(ENTREZ_LIST_A, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH"))
    filtering_A_enrichGO=FILTERING_BACKGROUND2(Aentrezgenes_enrichres,0.7) 
    filtering_A_enrichGO_res=fold_computing2(filtering_A_enrichGO)
    num_out=length(Aentrezgenes_enrichres[,1])-length(filtering_A_enrichGO_res[,1])
    num_notmapp=length(which(is.na(ENTREZ_LIST_A_DF[,2])))
    filtering_A_enrichGO_res2=cbind(filtering_A_enrichGO_res,"genes"=0)
    for ( i in 1:length(filtering_A_enrichGO_res2[,1])){
      vec_prot=unlist(strsplit(filtering_A_enrichGO_res2[i,3],split="/"))
      vec_prot_entrez=mapIds(org.Hs.eg.db, keys=(as.vector(vec_prot)),column=OUTPUT_id, keytype = "ENTREZID", multivals="first")
      vec_prot=paste(vec_prot_entrez,sep="/",collapse="/")
      filtering_A_enrichGO_res2[i,7]=vec_prot
    }
  }
  print(paste(num_notmapp, "/",length(ORIGINAL_gene_vec),  ORIGINAL_id , "ids not mapped to ENTREZ ID" , sep=" ", collapse = " "))
  print(paste(num_out, "/",length(Aentrezgenes_enrichres[,1]), ann,  "functions with gene background higher or equal to", bg_threshold, sep=" ", collapse = " "))
  list("entrez_id_df"=ENTREZ_LIST_A_DF ,     "filtered_FEA"=filtering_A_enrichGO_res, "FEA_result"=filtering_A_enrichGO_res2)
}



################


GENECO_MERGING_GOT=function(ABfun_enr,ABcom_size, Auni_size, Buni_size, gene_co_threshA, gene_co_threshB){
  gene_co_threshA=gene_co_threshA*100
  gene_co_threshB=gene_co_threshB*100
  gcmerged_ABfun_enr_res=GOT_allfussioning(ABfun_enr,gene_co_threshA,gene_co_threshB)
  gcmerged_ABfun_enr=gcmerged_ABfun_enr_res[[2]]
  
  num_cooc_merged=length(ABfun_enr[,1])-length(gcmerged_ABfun_enr[,1])
  print(paste(num_cooc_merged, "/",length(ABfun_enr[,1]), "functions merged by gene co-occurrence with a threshold of",gene_co_threshA/100  , sep=" ", collapse = " "))
  gcmerged_ABfun_enr
}

#df=ordered_ABfun_enr
#ftss_coin_min=gene_co_threshA
#S2B_coin_min=gene_co_threshB

GOT_allfussioning=function(df,ftss_coin_min,S2B_coin_min) {
  # order dataframe according to the best parameter
  #fussioned_GOT_df=as.data.frame(df)
  #fussioned_GOT_df[,7]=c(1)
  #fussioned_GOT_df[,8]=c(1)
  #fussioned_GOT_df[,9]=c(1)
  #fussioned_GOT_df[,10]=c(1)
  #fussioned_GOT_df[,11]=c(1)#
  #fussioned_GOT_df[,12]=c(1)#
  
  fussioned_GOT_df=df
  n=length(df[,1])
  remove=vector()
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if ((!is.element(i,remove))&(!is.element(j,remove))){
        genes_ftss1=as.vector(df[i,3]) #(vec1)
        genes_ftss2=as.vector(df[j,3]) #(vec2)
        ftss_check=GOTcoincidence_analysischeck(genes_ftss1,genes_ftss2) 
        GOTfun_ftss=GOTcoincidence_analysis(genes_ftss1,genes_ftss2) #maxsize
        
        genes_S2B1=as.vector(df[i,4]) #(vec1)
        genes_S2B2=as.vector(df[j,4]) #(vec2)
        S2B_check=GOTcoincidence_analysischeck(genes_S2B1,genes_S2B2)
        GOTfun_S2B=GOTcoincidence_analysis(genes_S2B1,genes_S2B2)
        
        a1=(ftss_check==0)#
        a2=(S2B_check==0)#
        
        c1=(GOTfun_ftss>=ftss_coin_min)
        c2=(GOTfun_S2B>=S2B_coin_min)
        
        allcheck=((c1&c2)|(a1&c2)|(a2&c1))
        
        if (allcheck){   ######### for coincidences
          fussioned_GO_vec=vector()
          fussioned_desc_vec=vector()
          genes=vector()
          if (fussioned_GOT_df[i,7]==1){ ########################
            fussioned_GOT_df[i,7]=df[j,1]
            fussioned_GOT_df[i,8]=df[j,2]
            fussioned_GOT_df[i,9]=df[j,3]
            fussioned_GOT_df[i,10]=df[j,4]
            #genes=paste(df[j,3],df[j,4],sep="/")
            fussioned_GOT_df[i,11]=df[j,5]
            fussioned_GOT_df[i,12]=df[j,6]#
            
          } else {
            fussioned_GO_vec=paste(fussioned_GOT_df[i,7],df[j,1],sep="/")
            fussioned_desc_vec=paste(fussioned_GOT_df[i,8],df[j,2],sep="/")
            fussioned_gene_ftss_vec=paste(fussioned_GOT_df[i,9],df[j,3],sep="/")
            fussioned_gene_S2B_vec=paste(fussioned_GOT_df[i,10],df[j,4],sep="/")
            fussioned_fold_ftss_vec=paste(fussioned_GOT_df[i,11],df[j,5],sep="/")
            fussioned_fold_S2B_vec=paste(fussioned_GOT_df[i,12],df[j,6],sep="/")
            #genes=paste(df[j,3],df[j,4],sep="/")
            fussioned_GOT_df[i,7]=fussioned_GO_vec
            fussioned_GOT_df[i,8]=fussioned_desc_vec
            fussioned_GOT_df[i,9]=fussioned_gene_ftss_vec
            fussioned_GOT_df[i,10]=fussioned_gene_S2B_vec
            fussioned_GOT_df[i,11]=fussioned_fold_ftss_vec
            fussioned_GOT_df[i,12]=fussioned_fold_S2B_vec
            
            
          }
          remove=c(remove, j)
        }
      }
    }
    
  }
  real_indexi=setdiff((1:n), remove) 
  colnames(fussioned_GOT_df)=c("GO_Terms_ID","Description","Original_Genes_A","Original_Genes_B","Original_Fold_A","Original_Fold_B","Merged_GOT","Merged_desc","Merged_A_genes","Merged_B_genes","Merged_A_fold","Merged_B_fold")
  depurated_fussioned_GOT_df=fussioned_GOT_df[real_indexi,]
  for (i in 1:length(depurated_fussioned_GOT_df[,1])){
    depurated_fussioned_GOT_df[i,9]=GOTrepeat_gene_elim(depurated_fussioned_GOT_df[i,9])
    depurated_fussioned_GOT_df[i,10]=GOTrepeat_gene_elim(depurated_fussioned_GOT_df[i,10])
  }
  list(fussioned_GOT_df,depurated_fussioned_GOT_df)
  
}



GOT_NOTfussioning=function(df,ftss_coin_min,S2B_coin_min) {
  # order dataframe according to the best parameter
  fussioned_GOT_df=as.data.frame(df)
  fussioned_GOT_df[,7]=c(1)
  fussioned_GOT_df[,8]=c(1)
  fussioned_GOT_df[,9]=c(1)
  fussioned_GOT_df[,10]=c(1)
  fussioned_GOT_df[,11]=c(1)#
  fussioned_GOT_df[,12]=c(1)#
  n=length(df[,1])
  remove=vector()
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if ((!is.element(i,remove))&(!is.element(j,remove))){
        genes_ftss1=as.vector(df[i,3]) #(vec1)
        genes_ftss2=as.vector(df[j,3]) #(vec2)
        ftss_check=GOTcoincidence_analysischeck(genes_ftss1,genes_ftss2) 
        GOTfun_ftss=GOTcoincidence_analysis(genes_ftss1,genes_ftss2) #maxsize
        
        genes_S2B1=as.vector(df[i,4]) #(vec1)
        genes_S2B2=as.vector(df[j,4]) #(vec2)
        S2B_check=GOTcoincidence_analysischeck(genes_S2B1,genes_S2B2)
        GOTfun_S2B=GOTcoincidence_analysis(genes_S2B1,genes_S2B2)
        
        a1=(ftss_check==0)#
        a2=(S2B_check==0)#
        
        c1=(GOTfun_ftss>=ftss_coin_min)
        c2=(GOTfun_S2B>=S2B_coin_min)
        
        allcheck=((c1&c2)|(a1&c2)|(a2&c1))
        
        if (allcheck){   ######### for coincidences
          fussioned_GO_vec=vector()
          fussioned_desc_vec=vector()
          genes=vector()
          if (fussioned_GOT_df[i,7]==1){ ########################
            fussioned_GOT_df[i,7]=df[j,1]
            fussioned_GOT_df[i,8]=df[j,2]
            fussioned_GOT_df[i,9]=df[j,3]
            fussioned_GOT_df[i,10]=df[j,4]
            #genes=paste(df[j,3],df[j,4],sep="/")
            fussioned_GOT_df[i,11]=df[j,5]
            fussioned_GOT_df[i,12]=df[j,6]#
            
          } else {
            fussioned_GO_vec=paste(fussioned_GOT_df[i,7],df[j,1],sep="/")
            fussioned_desc_vec=paste(fussioned_GOT_df[i,8],df[j,2],sep="/")
            fussioned_gene_ftss_vec=paste(fussioned_GOT_df[i,9],df[j,3],sep="/")
            fussioned_gene_S2B_vec=paste(fussioned_GOT_df[i,10],df[j,4],sep="/")
            fussioned_fold_ftss_vec=paste(fussioned_GOT_df[i,11],df[j,5],sep="/")
            fussioned_fold_S2B_vec=paste(fussioned_GOT_df[i,12],df[j,6],sep="/")
            #genes=paste(df[j,3],df[j,4],sep="/")
            fussioned_GOT_df[i,7]=fussioned_GO_vec
            fussioned_GOT_df[i,8]=fussioned_desc_vec
            fussioned_GOT_df[i,9]=fussioned_gene_ftss_vec
            fussioned_GOT_df[i,10]=fussioned_gene_S2B_vec
            fussioned_GOT_df[i,11]=fussioned_fold_ftss_vec
            fussioned_GOT_df[i,12]=fussioned_fold_S2B_vec
            
            
          }
          remove=c(remove, j)
        }
      }
    }
    
  }
  real_indexi=setdiff((1:n), remove) 
  colnames(fussioned_GOT_df)=c("GO_Terms_ID","Description","Original_Genes_A","Original_Genes_B","Original_Fold_A","Original_Fold_B","Merged_GOT","Merged_desc","Merged_A_genes","Merged_B_genes","Merged_A_fold","Merged_B_fold")
  depurated_fussioned_GOT_df=fussioned_GOT_df[real_indexi,]
  for (i in 1:length(depurated_fussioned_GOT_df[,1])){
    depurated_fussioned_GOT_df[i,9]=GOTrepeat_gene_elim(depurated_fussioned_GOT_df[i,9])
    depurated_fussioned_GOT_df[i,10]=GOTrepeat_gene_elim(depurated_fussioned_GOT_df[i,10])
  }
  list(fussioned_GOT_df,depurated_fussioned_GOT_df)
  
}





GOTcoincidence_analysis=function(vec1,vec2){
  genes1=strsplit(as.character(vec1),"/")
  genes2=strsplit(as.character(vec2),"/")
  inter=intersect(genes1[[1]],genes2[[1]])
  sizes=c(length(genes1[[1]]),length(genes2[[1]]))
  if (min(sizes)==0){
    coincid=0
  } else {
    coincid=length(inter)*100/min(sizes)
  }
}

GOTcoincidence_analysischeck=function(vec1,vec2){
  checkval=1
  if ((vec1=="0")|(vec2=="0")){
    checkval=0
  } 
  checkval
}



GOTrepeat_gene_elim=function(vec1){
  genes=strsplit(as.character(vec1),"/")
  unigenes=unique(genes[[1]])
  unigenes=setdiff(unigenes,c("0"))
  vec=paste(unigenes,sep="/",collapse="/")
}



uniprotA_BP07_res_r=uniprotA_BP07_res$filtered_FEA
################
alg_measure="Lin"
sim_thresh=0.7
initial_df=uniprotA_BP07_res_r
semData_info=d_BP


SEMANTIC_ANALYSIS3=function(alg_measure,sim_thresh,initial_df, semData_info){
  AB_GOT=initial_df[,1]
  AB_GOT=c(as.character(AB_GOT))
  AB_testing=termSim(AB_GOT, AB_GOT,semData=semData_info,method = alg_measure)
  tst_values=sm2vec(AB_testing, diag = FALSE)
  tst_index=sm.index(AB_testing, diag = F)
  order_tst_values=order(tst_values,decreasing = T)
  sorted_tst_values=tst_values[order_tst_values]
  sorted_tst_index=tst_index[order_tst_values,]
  e_df=rbind(sorted_tst_index[sorted_tst_values>=sim_thresh,])
  
  if (length(e_df[,1])!=0 & (all(is.na(e_df[,1])) &  all(is.na(e_df[,2]))) ==F) {
    e_df[which(is.na(e_df[,1])),1]=0
    e_df[which(is.na(e_df[,2])),2]=0
    e_df_values=c(unlist(e_df[,1:2]))

    e_df_summary=cbind(1:max(e_df_values),rep(0,max(e_df_values)),rep(0,max(e_df_values)))
    colnames(e_df_summary)=c("GOT_ind_list", "received_fused","deleted")
    minitial_df=initial_df
    for (i in 1:length(e_df[,1])){
      indexX=e_df[i,2]
      indexY=e_df[i,1]
      if (indexX!=0 && indexY!=0){
        e_df_summary[indexX,2]=e_df[i,1]
        e_df_summary[indexY,3]=1
        if (e_df_summary[indexX,3]==0){
          minitial_df[indexX,7]=paste(initial_df[indexY,1],initial_df[indexX,7],sep="/")
          minitial_df[indexX,8]=paste(initial_df[indexY,2],initial_df[indexX,8],sep="/")
          minitial_df[indexX,3]=paste(initial_df[indexY,3],initial_df[indexX,9],sep="/")
          minitial_df[indexX,4]=paste(initial_df[indexY,4],initial_df[indexX,10],sep="/")
          minitial_df[indexX,5]=paste(initial_df[indexY,5],initial_df[indexX,11],sep="/")
          minitial_df[indexX,6]=paste(initial_df[indexY,6],initial_df[indexX,12],sep="/")
        }
      }
    }
    del_ind=which(e_df_summary[,3]==1) #####
    minitial_df=minitial_df[-del_ind,] #####
    num_semsim_merged=length(initial_df[,1])-length(minitial_df[,1])
    print(paste(num_semsim_merged, "/",length(initial_df[,1]), "functions merged by semantic similarity with a threshold of",sim_thresh, "and", alg_measure, "algorithm" , sep=" ", collapse = " "))
    minitial_df
  }
  else {
    minitial_df=initial_df
    print(paste( "0 functions merged by semantic similarity with a threshold of",sim_thresh, "and", alg_measure, "algorithm" , sep=" ", collapse = " "))
    minitial_df
  }
}



#########

EDITING_TO_PLOT_noS2B=function(Ageneset_size, Bgeneset_size,merged_fussioning_result_df){
  all_merged_fussioning_result_df=data.frame()
  all_merged_fussioning_result_df=merged_fussioning_result_df
  for (i in 1:length(merged_fussioning_result_df[,1])){
    for (z in 1:12){
      if (z<=6){
        p_vec=paste(merged_fussioning_result_df[i,z],merged_fussioning_result_df[i,z+6],sep="/",collapse = "/")
        all_merged_fussioning_result_df[i,z+6]=p_vec
      }
    }
  }
  comp_merged_fussioning_result_df=data.frame()
  comp_merged_fussioning_result_df=all_merged_fussioning_result_df
  del=c("0","1")
  empty=" "
  for (i in 1:length(all_merged_fussioning_result_df[,9])){
    genes_a=strsplit(all_merged_fussioning_result_df[i,9],"/")
    genes_a=setdiff(genes_a[[1]],del)
    if (length(genes_a)!=0){
      genes_f=setdiff(genes_a[[1]],empty)
      f=strsplit(genes_f,"")
      if (length(f[[1]])==0){
        genes_a=vector()
      }
    }
    genes_a_freq=(length(genes_a)/Ageneset_size)*100
    comp_merged_fussioning_result_df[i,13]=genes_a_freq
    
    genes_b=strsplit(all_merged_fussioning_result_df[i,10],"/")
    genes_b=setdiff(genes_b[[1]],del)
    if (length(genes_b)!=0){
      genes_ff=setdiff(genes_b[[1]],empty)
      ff=strsplit(genes_ff,"")
      if (length(ff[[1]])==0){
        genes_b=vector()
      }
    }
    genes_b_freq=(length(genes_b)/Bgeneset_size)*100
    comp_merged_fussioning_result_df[i,14]=genes_b_freq
    
    fold_a=strsplit(all_merged_fussioning_result_df[i,11],"/")
    fold_a=setdiff(fold_a[[1]],del)
    if(length(fold_a)>=5){
      median_a5fold=median(as.numeric(fold_a[1:5]))
      comp_merged_fussioning_result_df[i,15]=median_a5fold
    }else{
      median_afold=median(as.numeric(fold_a))
      comp_merged_fussioning_result_df[i,15]=median_afold
    }
    
    fold_b=strsplit(all_merged_fussioning_result_df[i,12],"/")
    fold_b=setdiff(fold_b[[1]],del)
    if(length(fold_b)>=5){
      median_b5fold=median(as.numeric(fold_b[1:5]))
      comp_merged_fussioning_result_df[i,16]=median_b5fold
    }else{
      median_bfold=median(as.numeric(fold_b))
      comp_merged_fussioning_result_df[i,16]=median_bfold
    }
  }
  colnames(comp_merged_fussioning_result_df)=c("GO_Terms_ID","Description","Original_Genes_A","Original_Genes_B" , "Original_Fold_A" , "Original_Fold_B" ,  "Merged_GO_T" ,"Merged_desc_T","Merged_A_genes" ,  "Merged_B_genes"  ,  "Merged_A_fold"  ,  "Merged_B_fold","A_genesfreq","B_genesfreq","A_medianfold","B_medianfold")
  comp_merged_fussioning_result_df
}
########
FILTERING_BACKGROUND=function(enrichGOdf,bg_threshold) {
  bg_enrichGOdf=data.frame()
  for (i in 1:length(enrichGOdf[,1])){
    a=strsplit(as.character(enrichGOdf[i,4]),"/")
    v=as.numeric(a[[1]][1])/as.numeric(a[[1]][2])
    enrichGOdf[i,10]=v
  }
  bg_enrichGOdf=enrichGOdf
  colnames(bg_enrichGOdf)=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Background")
  x_bg_index=which(enrichGOdf[,10]<=bg_threshold)
  filteredbg_enrichGO=enrichGOdf[x_bg_index,]
  filtered_GOTs=as.vector(filteredbg_enrichGO[,1])
  list(filteredbg_enrichGO, filtered_GOTs,bg_enrichGOdf) # list with 3 outputs (ALREADY FILTERED COMPLETE RESULT, ONLY GOT, ORIGINAL FUNCTIONAL ENRICHMENT WIHT BACKGROUND COMPUTED )
}

######################

JOINING_COMPARING_DATASETS=function(funct_enrA, funct_enrB){
  AB_colnames_vector=c("GO_T", "desc" ,"Genes_A","Genes_B", "Fold_A", "Fold_B","ad.p-value_A", "ad.p-value_B" )
  A_unique_colnames_vector=c("GO_T","desc","Genes_A","Fold_A","ad.p-value_A")
  B_unique_colnames_vector=c("GO_T","desc","Genes_B","Fold_B","ad.p-value_B")
  
  AB_common_terms=common_terms(funct_enrA,funct_enrB,AB_colnames_vector)
  A_unique_terms_againstB=unique_terms(funct_enrA,funct_enrB,A_unique_colnames_vector)
  B_unique_terms_againstA=unique_terms(funct_enrB,funct_enrA,B_unique_colnames_vector)
  
  ed_uniA=cbind("Go"=A_unique_terms_againstB[,1],"Desc"=A_unique_terms_againstB[,2],
                "Genes_A"=as.character(A_unique_terms_againstB[,3]),"Genes_B"=rep(0,length(A_unique_terms_againstB[,1])),
                "Fold_A"=A_unique_terms_againstB[,4],
                "Fold_B"=rep(0,length(A_unique_terms_againstB[,1])))
  
  ed_uniB=cbind("Go"=B_unique_terms_againstA[,1],"Desc"=B_unique_terms_againstA[,2],
                "Genes_A"=rep(0,length(B_unique_terms_againstA[,3])),"Genes_B"=as.character(B_unique_terms_againstA[,3]),
                "Fold_A"=rep(0,length(B_unique_terms_againstA[,3])),
                "Fold_B"=B_unique_terms_againstA[,4])
  
  ed_AB=cbind("Go"=AB_common_terms[,1],"Desc"=AB_common_terms[,2],
              "Genes_A"=as.character(AB_common_terms[,3]),"Genes_B"=as.character(AB_common_terms[,4]),
              "Fold_A"=AB_common_terms[,5],
              "Fold_B"=AB_common_terms[,6])
  
  AB_joined_data=rbind(ed_AB,ed_uniA,ed_uniB)
  AB_joined_data=data.frame(AB_joined_data,stringsAsFactors = F)
  common_size=length(ed_AB[,1])
  uniA_size=length(ed_uniA[,1])
  uniB_size=length(ed_uniB[,1])
  list(AB_joined_data, common_size, uniA_size, uniB_size)
}

######################

GENECO_MERGING_GOT=function(ABfun_enr,ABcom_size, Auni_size, Buni_size, gene_co_threshA, gene_co_threshB){
  total=ABcom_size+Auni_size+Buni_size 
  doubleset=c(rep(1,ABcom_size),rep(0,total-ABcom_size))  
  maxscore=pmax(ABfun_enr[,5],ABfun_enr[,6]) #by fold score
  maxscore=as.numeric(maxscore)
  reorder=order(doubleset,maxscore,decreasing=T)
  ordered_ABfun_enr=data.frame()
  ordered_ABfun_enr=ABfun_enr[reorder,]
  
  gcmerged_ABfun_enr_res=GOT_allfussioning(ordered_ABfun_enr,gene_co_threshA,gene_co_threshB)
  gcmerged_ABfun_enr=gcmerged_ABfun_enr_res[[2]]
  gcmerged_ABfun_enr
}

######################

SEMANTIC_ANALYSIS2=function(alg_measure,sim_thresh,initial_df){
  AB_GOT=initial_df[,1]
  AB_GOT=c(as.character(AB_GOT))
  AB_testing=termSim(AB_GOT, AB_GOT,semData=d,method = "Lin")
  tst_values=sm2vec(AB_testing, diag = FALSE)
  tst_index=sm.index(AB_testing, diag = F)
  order_tst_values=order(tst_values,decreasing = T)
  sorted_tst_values=tst_values[order_tst_values]
  sorted_tst_index=tst_index[order_tst_values,]
  e_df=rbind(sorted_tst_index[sorted_tst_values>=0.7,])
  
  if (length(e_df[,1])!=0 & (all(is.na(e_df[,1]))&  all(is.na(e_df[,2]))==F)){
    e_df_values=c(unlist(e_df[,1:2]))
    e_df_summary=cbind(1:max(e_df_values),rep(0,max(e_df_values)),rep(0,max(e_df_values)))
    colnames(e_df_summary)=c("GOT_ind_list", "received_fused","deleted")
    minitial_df=initial_df
    for (i in 1:length(e_df[,1])){
      indexX=e_df[i,2]
      indexY=e_df[i,1]
      e_df_summary[indexX,2]=e_df[i,1]
      e_df_summary[indexY,3]=1
      if (e_df_summary[indexX,3]==0){
        minitial_df[indexX,7]=paste(initial_df[indexY,1],initial_df[indexX,7],sep="/")
        minitial_df[indexX,8]=paste(initial_df[indexY,2],initial_df[indexX,8],sep="/")
        minitial_df[indexX,3]=paste(initial_df[indexY,3],initial_df[indexX,9],sep="/")
        minitial_df[indexX,4]=paste(initial_df[indexY,4],initial_df[indexX,10],sep="/")
        minitial_df[indexX,5]=paste(initial_df[indexY,5],initial_df[indexX,11],sep="/")
        minitial_df[indexX,6]=paste(initial_df[indexY,6],initial_df[indexX,12],sep="/")
      }
    }
    del_ind=which(e_df_summary[,3]==1) #####
    minitial_df=minitial_df[-del_ind,] #####
    list(minitial_df,e_df_summary)
  }
  else {
    minitial_df=initial_df
    e_df_summary="there is not enought semantic similarity"
    list(minitial_df,e_df_summary)
  }
}

######################


common_terms=function(dfA,dfB,df_colnames_vector){
  res_common_termsdf=data.frame(stringsAsFactors = F)
  paralell_index=0
  for (i in 1:length(dfA[,1])) {
    if (is.element(dfA[i,1],dfB[,1])) {
      bbb_index=which(dfB[,1]==dfA[i,1])
      paralell_index=paralell_index+1
      res_common_termsdf[paralell_index,1]=dfA[i,1]
      res_common_termsdf[paralell_index,2]=dfA[i,2]
      res_common_termsdf[paralell_index,3]=dfB[bbb_index,8]
      res_common_termsdf[paralell_index,4]=dfA[i,8]
      
      a_ftss=strsplit(as.character(dfB[bbb_index,3]),"/") #generatio
      v_ftss=as.numeric(a_ftss[[1]][1])/as.numeric(a_ftss[[1]][2]) #generatio divided
      b_ftss=strsplit(as.character(dfB[bbb_index,4]),"/")
      w_ftss=as.numeric(b_ftss[[1]][1])/as.numeric(b_ftss[[1]][2])
      ab_ftss=v_ftss/w_ftss
      res_common_termsdf[paralell_index,5]=ab_ftss
      
      a_s2b=strsplit(as.character(dfA[i,3]),"/") #generatio
      v_s2b=as.numeric(a_s2b[[1]][1])/as.numeric(a_s2b[[1]][2]) #generatio divided
      b_s2b=strsplit(as.character(dfA[i,4]),"/")
      w_s2b=as.numeric(b_s2b[[1]][1])/as.numeric(b_s2b[[1]][2])
      ab_s2b=v_s2b/w_s2b
      res_common_termsdf[paralell_index,6]=ab_s2b
      
      res_common_termsdf[paralell_index,7]=dfB[bbb_index,5]
      res_common_termsdf[paralell_index,8]=dfA[i,5]
    }
  }
  colnames(res_common_termsdf)=df_colnames_vector
  res_common_termsdf
}


######################

unique_terms=function(dfA,dfB,df_colnames_vector){
  res_unique_termsdf=data.frame(stringsAsFactors = F)
  paralell_index=0
  for (i in 1:length(dfA[,1])) {
    if (!is.element(dfA[i,1],dfB[,1])) {
      paralell_index=paralell_index+1
      res_unique_termsdf[paralell_index,1]=dfA[i,1]
      res_unique_termsdf[paralell_index,2]=dfA[i,2]
      res_unique_termsdf[paralell_index,3]=dfA[i,8]
      a_s2b=strsplit(as.character(dfA[i,3]),"/") #generatio
      v_s2b=as.numeric(a_s2b[[1]][1])/as.numeric(a_s2b[[1]][2]) #generatio divided
      b_s2b=strsplit(as.character(dfA[i,4]),"/")
      w_s2b=as.numeric(b_s2b[[1]][1])/as.numeric(b_s2b[[1]][2])
      ab_s2b=v_s2b/w_s2b
      res_unique_termsdf[paralell_index,4]=ab_s2b
      res_unique_termsdf[paralell_index,5]=dfA[i,5]
    }
  }
  colnames(res_unique_termsdf)=df_colnames_vector
  res_unique_termsdf
}


######################

GOT_allfussioning=function(df,ftss_coin_min,S2B_coin_min) {
  # order dataframe according to the best parameter
  fussioned_GOT_df=as.data.frame(df)
  fussioned_GOT_df[,7]=c(1)
  fussioned_GOT_df[,8]=c(1)
  fussioned_GOT_df[,9]=c(1)
  fussioned_GOT_df[,10]=c(1)
  fussioned_GOT_df[,11]=c(1)#
  fussioned_GOT_df[,12]=c(1)#
  n=length(df[,1])
  remove=vector()
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if ((!is.element(i,remove))&(!is.element(j,remove))){
        genes_ftss1=as.vector(df[i,3]) #(vec1)
        genes_ftss2=as.vector(df[j,3]) #(vec2)
        ftss_check=GOTcoincidence_analysischeck(genes_ftss1,genes_ftss2) 
        GOTfun_ftss=GOTcoincidence_analysis(genes_ftss1,genes_ftss2) #maxsize
        
        genes_S2B1=as.vector(df[i,4]) #(vec1)
        genes_S2B2=as.vector(df[j,4]) #(vec2)
        S2B_check=GOTcoincidence_analysischeck(genes_S2B1,genes_S2B2)
        GOTfun_S2B=GOTcoincidence_analysis(genes_S2B1,genes_S2B2)
        
        a1=(ftss_check==0)#
        a2=(S2B_check==0)#
        
        c1=(GOTfun_ftss>=ftss_coin_min)
        c2=(GOTfun_S2B>=S2B_coin_min)
        
        allcheck=((c1&c2)|(a1&c2)|(a2&c1))
        
        if (allcheck){   ######### for coincidences
          fussioned_GO_vec=vector()
          fussioned_desc_vec=vector()
          genes=vector()
          if (fussioned_GOT_df[i,7]==1){ ########################
            fussioned_GOT_df[i,7]=df[j,1]
            fussioned_GOT_df[i,8]=df[j,2]
            fussioned_GOT_df[i,9]=df[j,3]
            fussioned_GOT_df[i,10]=df[j,4]
            #genes=paste(df[j,3],df[j,4],sep="/")
            fussioned_GOT_df[i,11]=df[j,5]
            fussioned_GOT_df[i,12]=df[j,6]#
            
          } else {
            fussioned_GO_vec=paste(fussioned_GOT_df[i,7],df[j,1],sep="/")
            fussioned_desc_vec=paste(fussioned_GOT_df[i,8],df[j,2],sep="/")
            fussioned_gene_ftss_vec=paste(fussioned_GOT_df[i,9],df[j,3],sep="/")
            fussioned_gene_S2B_vec=paste(fussioned_GOT_df[i,10],df[j,4],sep="/")
            fussioned_fold_ftss_vec=paste(fussioned_GOT_df[i,11],df[j,5],sep="/")
            fussioned_fold_S2B_vec=paste(fussioned_GOT_df[i,12],df[j,6],sep="/")
            #genes=paste(df[j,3],df[j,4],sep="/")
            fussioned_GOT_df[i,7]=fussioned_GO_vec
            fussioned_GOT_df[i,8]=fussioned_desc_vec
            fussioned_GOT_df[i,9]=fussioned_gene_ftss_vec
            fussioned_GOT_df[i,10]=fussioned_gene_S2B_vec
            fussioned_GOT_df[i,11]=fussioned_fold_ftss_vec
            fussioned_GOT_df[i,12]=fussioned_fold_S2B_vec
            
            
          }
          remove=c(remove, j)
        }
      }
    }
    
  }
  real_indexi=setdiff((1:n), remove) 
  colnames(fussioned_GOT_df)=c("GO_Terms_ID","Description","Original_Genes_A","Original_Genes_B","Original_Fold_A","Original_Fold_B","Merged_GOT","Merged_desc","Merged_A_genes","Merged_B_genes","Merged_A_fold","Merged_B_fold")
  depurated_fussioned_GOT_df=fussioned_GOT_df[real_indexi,]
  for (i in 1:length(depurated_fussioned_GOT_df[,1])){
    depurated_fussioned_GOT_df[i,9]=GOTrepeat_gene_elim(depurated_fussioned_GOT_df[i,9])
    depurated_fussioned_GOT_df[i,10]=GOTrepeat_gene_elim(depurated_fussioned_GOT_df[i,10])
  }
  list(fussioned_GOT_df,depurated_fussioned_GOT_df)
  
}

######################

GOTcoincidence_analysischeck=function(vec1,vec2){
  checkval=1
  if ((vec1=="0")|(vec2=="0")){
    checkval=0
  } 
  checkval
}

######################

GOTcoincidence_analysis=function(vec1,vec2){
  genes1=strsplit(as.character(vec1),"/")
  genes2=strsplit(as.character(vec2),"/")
  inter=intersect(genes1[[1]],genes2[[1]])
  sizes=c(length(genes1[[1]]),length(genes2[[1]]))
  if (min(sizes)==0){
    coincid=0
  } else {
    coincid=length(inter)*100/min(sizes)
  }
}

######################

GOTrepeat_gene_elim=function(vec1){
  genes=strsplit(as.character(vec1),"/")
  unigenes=unique(genes[[1]])
  unigenes=setdiff(unigenes,c("0"))
  vec=paste(unigenes,sep="/",collapse="/")
}

######################
fold_computing=function(tes_B_func_enrichment){
  tes2_B_func_enrichment=data.frame()
  for (i in 1:length(tes_B_func_enrichment[,1])){
    a_ftss=unlist(strsplit(as.character(tes_B_func_enrichment[i,3]),"/")) #generatio
    v_ftss=as.numeric(a_ftss[1])/as.numeric(a_ftss[2]) #generatio divided
    b_ftss=unlist(strsplit(as.character(tes_B_func_enrichment[i,4]),"/"))
    w_ftss=as.numeric(b_ftss[1])/as.numeric(b_ftss[2])
    ab_ftss=v_ftss/w_ftss
    tes2_B_func_enrichment[i,1]=tes_B_func_enrichment[i,1]
    tes2_B_func_enrichment[i,2]=tes_B_func_enrichment[i,2]
    tes2_B_func_enrichment[i,3]=tes_B_func_enrichment[i,8]
    tes2_B_func_enrichment[i,4]="0"
    tes2_B_func_enrichment[i,5]=ab_ftss
    tes2_B_func_enrichment[i,6]="0"
    colnames_vector=c( "GO_T", "desc" ,"Genes_Aa","Genes_Ab", "Fold_Aa", "Fold_Ab")
    colnames(tes2_B_func_enrichment)=colnames_vector
  }
  tes2_B_func_enrichment
}
################
single_funct_simplif2=function(RAW_FUNC_ENRICH_DATAFRAME_A,simp_order, gene_cooc_threshold, semantic_alg, semantic_threshold, semData_info){
  
  Auni_size=length(RAW_FUNC_ENRICH_DATAFRAME_A[,1])
  Buni_size=0
  ABcom_size=0
  
  total=ABcom_size+Auni_size+Buni_size 
  doubleset=c(rep(1,ABcom_size),rep(0,total-ABcom_size))  
  maxscore=as.numeric(pmax(RAW_FUNC_ENRICH_DATAFRAME_A[,5],RAW_FUNC_ENRICH_DATAFRAME_A[,6])) #by fold score
  reorder=order(doubleset,maxscore,decreasing=T)
  ordered_ABfun_enr=RAW_FUNC_ENRICH_DATAFRAME_A[reorder,]
  
  ordered_ABfun_enr=cbind(ordered_ABfun_enr, 1,1,1,1,1,1)
  colnames(ordered_ABfun_enr)=c("GO_Terms_ID","Description","Original_Genes_A","Original_Genes_B","Original_Fold_A","Original_Fold_B","Merged_GOT","Merged_desc","Merged_A_genes","Merged_B_genes","Merged_A_fold","Merged_B_fold")
  
  
  if (simp_order=="cooc_semsim"){
    geneco_AB_func_enrichment=GENECO_MERGING_GOT(ordered_ABfun_enr,0,length(ordered_ABfun_enr[,1]),0, gene_cooc_threshold, gene_cooc_threshold)
    if (class(semData_info)=="GOSemSimDATA") {
      sem_AB_func_enrichment=SEMANTIC_ANALYSIS3(semantic_alg,semantic_threshold,geneco_AB_func_enrichment,semData_info)
    }
    if (class(semData_info)!="GOSemSimDATA"){
      sem_AB_func_enrichment=geneco_AB_func_enrichment
      print( "KEGG functions do not have hierarchy")
    }
    final_simplified_res=sem_AB_func_enrichment
    print(paste("from an original size of",length(ordered_ABfun_enr[,1]) ,
                "functions, we get a simplified list of",length(final_simplified_res[,1]),"functional groups",sep=" ", collapse = " "  ))
    
    #list("original_FEA"=RAW_FUNC_ENRICH_DATAFRAME_A,"gene_coocc_simp"=geneco_AB_func_enrichment,"sem_sim_simp"=sem_AB_func_enrichment,"both_simp"=final_simplified_res)
  }
  
  if (simp_order=="semsim_cooc"){
    if (class(semData_info)=="GOSemSimDATA") {
      sem_AB_func_enrichment=SEMANTIC_ANALYSIS3(semantic_alg,semantic_threshold,ordered_ABfun_enr,semData_info)
    }
    if (class(semData_info)!="GOSemSimDATA"){
      sem_AB_func_enrichment=ordered_ABfun_enr
      print( "KEGG functions do not have hierarchy")
    }
    geneco_AB_func_enrichment=GENECO_MERGING_GOT(sem_AB_func_enrichment,0,length(sem_AB_func_enrichment[,1]),0, gene_cooc_threshold, gene_cooc_threshold)    
    final_simplified_res=geneco_AB_func_enrichment
    print(paste("from an original size of",length(ordered_ABfun_enr[,1]) ,
                "functions, we get a simplified list of",length(final_simplified_res[,1]),"functional groups",sep=" ", collapse = " "  ))
    
    
  }
  list("original_FEA"=ordered_ABfun_enr,"gene_coocc_simp"=geneco_AB_func_enrichment,"sem_sim_simp"=sem_AB_func_enrichment,"both_simp"=final_simplified_res)
}
