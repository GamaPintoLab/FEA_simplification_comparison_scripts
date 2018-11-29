

## LOAD LIBRARIES

library(org.Hs.eg.db) 
library("clusterProfiler")


## LOAD IMPUT DATA (GENE VECTOR TO FUNCTIONALLY ENRICH)

gene_symbol_vecA_df=read.csv("gene_symbol_vecA_df.csv", header=T, stringsAsFactors = F)
gene_symbol_vecB_df=read.csv("gene_symbol_vecB_df.csv", header=T, stringsAsFactors = F)



# BEFORE doing any analysis, go to line titled "# ASSOCIATED FUNCTIONS " and run all the lines until the end of the file

# mapping and simplest FEA for the two sets independently

# SET A

# If necessary, These 4 lines map UNIPROT IDs to ENTREZ IDs
ENTREZ_LIST_A=mapIds(org.Hs.eg.db, keys=(as.vector(gene_symbol_vecA_df[,1])),column="ENTREZID", keytype = "SYMBOL", multivals="first")
ENTREZ_LIST_A=as.vector(as.numeric(ENTREZ_LIST_A))

# Functional Enrichment Analysis (FEA) - BIOLOGICAL PROCESS
Aentrezgenes_enrichGOres=as.data.frame(enrichGO(ENTREZ_LIST_A,'org.Hs.eg.db',ont="BP"))

# SET B

# If necessary, These 4 lines map UNIPROT IDs to ENTREZ IDs
ENTREZ_LIST_B=mapIds(org.Hs.eg.db, keys=(as.vector(gene_symbol_vecB_df[,1])),column="ENTREZID", keytype = "SYMBOL", multivals="first")
ENTREZ_LIST_B=as.vector(as.numeric(ENTREZ_LIST_B))

# Functional Enrichment Analysis (FEA) - BIOLOGICAL PROCESS
Bentrezgenes_enrichGOres=as.data.frame(enrichGO(ENTREZ_LIST_B,'org.Hs.eg.db',ont="BP"))




##### FEA simplification

#  by a) gene co-occurrence, b) semantic similarity or c) both
# if you want to avoid one of themm , use a 1000 threshold

library("GOSemSim")
library("corpcor")

d_BP <- godata('org.Hs.eg.db', ont="BP", computeIC=T) 
#d_MF <- godata('org.Hs.eg.db', ont="MF", computeIC=T) 
#d_KEGG=0


# if you want both simplifications define your thresholds
# if you want only simp by gene-co occurrence , define semantic_threshold= 1.1
# if your want only simp by sem sim, define gene_cooc_threshold= 101


#####

##### FEA simplification AND COMPARISON FOR TWO GENE-FUNCTIONS LISTS

# before running the simplification function, go to line 69-end and run ALL THE FUNCTIONS



FEAs_comparison_outputs=plural_funct_simplif3(Aentrezgenes_enrichGOres,Bentrezgenes_enrichGOres,0.7,101,"Lin",0.7, d_BP)
class(FEAs_comparison_outputs$merged_FEA_AB)

  
## IF YOU WANT ONLY THE COMPARISON PLOTS
FEAs_simplified_comparison_editing=EDITING_TO_PLOT_noS2B(length(ENTREZ_LIST_A),length(ENTREZ_LIST_B), FEAs_comparison_outputs$merged_FEA_AB) ###  without simplification!!  
FEAs_simplified_comparison_editing=list("FEAs_comparison_outputs_list"=FEAs_comparison_outputs,"edited_output"=FEAs_comparison_editing)

## IF YOU WANT THE COMPARISON AND SIMPLIFICATION PLOTS
FEAs_simplified_comparison_editing=EDITING_TO_PLOT_noS2B(length(ENTREZ_LIST_A),length(ENTREZ_LIST_B) ,FEAs_comparison_outputs$both_simp) ###  simplified
FEAs_simplified_comparison_edited_to_plot=list("FEAs_comparison_outputs_list"=FEAs_comparison_outputs,"edited_output"=FEAs_simplified_comparison_editing)


#OUTPUT: 

# FEAs_comparison_edited_to_plot[[1]][[1]] --> original_FEA_A
# FEAs_comparison_edited_to_plot[[1]][[2]] --> original_FEA_B
# FEAs_comparison_edited_to_plot[[1]][[3]] --> merged_FEA_AB
# FEAs_comparison_edited_to_plot[[1]][[4]] --> geneco_AB_func_enrichment
# FEAs_comparison_edited_to_plot[[1]][[5]] --> sem_AB_func_enrichment
# FEAs_comparison_edited_to_plot[[1]][[6]] --> final_AB_simplified_res

# FEAs_comparison_edited_to_plot[[2]] --> final_AB_simplified_res edited to be ploted



####### PLOT RESULTS 


library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)


#original_df=FEAs_comparison_edited_to_plot$FEAs_comparison_outputs_list$merged_FEA_AB 
#resulting_df=FEAs_comparison_edited_to_plot$edited_output


original_df=FEAs_simplified_comparison_edited_to_plot$FEAs_comparison_outputs_list$merged_FEA_AB
resulting_df=FEAs_simplified_comparison_edited_to_plot$edited_output


rearrange_fun_res=rearrange_df_toplot_function(original_df,resulting_df)

#### DEFINE THE KEY TERMS FOR FUNCTIONAL CLASSES ASSIGMENT BY TEXT MINING 
# ALL THE GO TERMS Description. Change the keyterms bellow that best adjust to your results and interests
# The number of groups can be edited as well but might affect to plot rendering
# only group1 cannot be changed, it will contain "other" GO TERMS that do not match with your specifications

rearrange_fun_res$edited_resulting_df[,3] 

## KEY TERMS ARGUMENT: 

#group1 will be considered OTHERS (GOTs that do not match with any key term) ALWAYS KEEP IT THAT WAY

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



## CAUTION: DON'T FORGUET TO CHANGE THE GROUPS HERE ASWELL
GROUPS=list(group2,group3,group4,group5,group6,group7,group8,group9,group10,group11,group12,group13,group14, group15, group16 )



## GROUP NAMES ARGUMENT: to be plot in the plots legend, you will be asked to choose a) to plot all key terms ("group_names")
## or b) simpler functional class names that you have to edit in the following lines ("group_names_plot")

group_names_plot=c("0. Other",
                   "1. Nervous system",
                   "2. Immune system", 
                   "3. Muscle",
                   "4. Stress", 
                   "5. Folding",
                   "6. Apoptosis", 
                   "7. Cytoskeleton", 
                   "8. RNA processing",
                   "9. Transcription", 
                   "10. DNA repair", 
                   "11. Protein degradation",
                   "12. Cell cycle", 
                   "13. Protein export/import",
                   "14. Signaling", 
                   "15. Development")
                   
                   
## COLORS ARGUMENT: If you prefer other colors to be ploted, change here the codex
## the order of the colors will match with the order of classes assigned above
# use pie function to check the colors

colors_def=c( "#546269" ,"red" ,    "#FF6600", "#FFCC00" ,"#CCFF00", "green",   "#71c100" ,"#a37458" ,"#00FFCC", "#00CCFF" ,
            "#0066FF", "blue" ,   "#6600FF" ,"#CC00FF","#FF00CC", "#e4acf6")

pie(rep(1,length(colors_def)),col=colors_def) 


subsetting_fun_res=subsetting_summres_plots_func(GROUPS,rearrange_fun_res$edited_resulting_df, rearrange_fun_res$df_toplot,group_names_plot)

wo_others_histograms_plots=total_wo_others_histograms_function(colors_def,group_names_plot,subsetting_fun_res)

graphics.off()

#ggsave(filename="FEA_comparison_histogram.pdf",plot=wo_others_histograms_plots, device="pdf",limitsize=F, scale=2)


## PLOT THE LEGEND 

hist_example_df=subsetting_fun_res[[1]]
hist_example_df=hist_example_df[2:length(hist_example_df[,1]),]

max_fun=max(hist_example_df$Total)
colorsx=colors_def[2:length(colors_def)]
group_namesx=group_names_plot[2:length(group_names_plot)]
GROUPS_user=GROUPS

hist_example_plot=ggplot(data=hist_example_df) + 
  geom_histogram(aes(x=Group,y=Total,fill=colorsx), stat = "identity", color="black")+
  scale_fill_manual(values=colorsx, 
                    name="Function classes" ,
                    labels=c(group_namesx))+
  theme_light()+
  theme(legend.title = element_text(colour="black", size=16, face="bold"))


mylegend<-g_legend(hist_example_plot)

ggsave(filename="FEA_comparison_histogram_legend.pdf",plot=mylegend, device="pdf",limitsize=F, scale=1.5) 



# ASSOCIATED FUNCTIONS - run these lines BEFORE doing any analysis


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


######

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

#######


subsetting_summres_plots_func=function(GROUPS_user,edited_resulting_df, df_toplot, group_names_user){
  functions_groups_df=as.data.frame(matrix(0, nrow=length(edited_resulting_df[,1]), ncol=length(GROUPS_user)+1)) # nrow=#functional groups to assign , # ncol= #functional classes defined above + 1
  for ( i in 1:length(edited_resulting_df[,1])){
    funcpast=paste(edited_resulting_df[i,1],edited_resulting_df[i,3],sep="/")
    functions_groups_df[i,1]=funcpast
  }
  for ( i in 1:length(functions_groups_df[,1])){
    fun_count=0
    splited_funct=strsplit(functions_groups_df[i,1], "/")
    for (j in 1:length(GROUPS_user)){
      for ( h in 1:length(splited_funct[[1]])){
        splited_funct_sent=strsplit(splited_funct[[1]][[h]], "/| ")
        fun_count=length(which(is.element(splited_funct_sent[[1]],GROUPS_user[[j]])))
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
  
  
  sizeord=order(df_toplot[,5], decreasing = T) ## according to gene freq
  ordbysize_df_toplot=df_toplot[sizeord,]
  
  
  ordbysize_df_toplot2=data.frame(ordbysize_df_toplot, stringsAsFactors = F)
  for (i in 1:length(ordbysize_df_toplot2[,1])){
    ordbysize_df_toplot2[i,6]=(ordbysize_df_toplot2[i,6]+1)
  }
  
  ordbysize_df_toplot2$Function_class=as.factor(ordbysize_df_toplot2$Function_class)
  ordbysize_df_toplot2$Function_class=factor(ordbysize_df_toplot2$Function_class,levels=c(1:length(GROUPS_user)))
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
    functions_groups_df=as.data.frame(matrix(0, nrow=length(EACH_edited_resulting_df[,1]), ncol=length(GROUPS_user)+1)) # nrow=#functional groups to assign , # ncol= #functional classes defined above + 1
    colnames(functions_groups_df)[1]="RepDesc"
    
    for ( i in 1:length(EACH_edited_resulting_df[,1])){
      funcpast=paste(EACH_edited_resulting_df[i,1],EACH_edited_resulting_df[i,3],sep="/")
      functions_groups_df[i,1]=funcpast
    }
    for ( i in 1:length(functions_groups_df[,1])){
      fun_count=0
      splited_funct=strsplit(functions_groups_df[i,1], "/")
      for (j in 1:length(GROUPS_user)){
        for ( h in 1:length(splited_funct[[1]])){
          splited_funct_sent=strsplit(splited_funct[[1]][[h]], "/| ")
          fun_count=length(which(is.element(splited_funct_sent[[1]],GROUPS_user[[j]])))
          if (fun_count>=1){
            count=(functions_groups_df[i,j+1])+1
            functions_groups_df[i,j+1]=count
          }
        }
      }
    }
    functions_hist_toplotdf=data.frame("Group_names"=group_names_user, "Group"=c(1:length(group_names_user)), "Total"=rep(0, length(group_names_user)),"Represented"=rep(0, length(group_names_user)))
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
    #totnum=paste("Total functional groups",sum(functions_hist_toplotdf$Represented), sep=": ")
    SUBSETS_functions_hist_toplotdf_LIST[[w]]=functions_hist_toplotdf
    
  }
  names(SUBSETS_functions_hist_toplotdf_LIST)=c("uniqueA_set", "uniqueB_set", "common_set")
  SUBSETS_functions_hist_toplotdf_LIST
}

#######

total_wo_others_histograms_function=function(colors_user,group_names_user,functions_hist_toplotdf){
  val10=25
  colorsx=colors_user[2:length(colors_user)]
  group_namesx=group_names_user[2:length(group_names_user)]
  functions_hist_toplotdf_A=functions_hist_toplotdf[[1]]
  functions_hist_toplotdf_noother=functions_hist_toplotdf_A[2:length(colors_user),]
  
  max_fun=max(functions_hist_toplotdf_noother$Total)
  totnum=paste("Unique A GO groups",sum(functions_hist_toplotdf_noother$Represented), sep=": ")
  for ( i in 1:length(functions_hist_toplotdf_noother[,1])){
    functions_hist_toplotdf_noother[i,2]=as.numeric(functions_hist_toplotdf_noother[i,2])-1
  }
  functions_hist_toplotdf_noother$Group=factor(functions_hist_toplotdf_noother$Group, levels=c(as.character(1:length(GROUPS_user))))
  
  UNIAhist_plot=ggplot(data=functions_hist_toplotdf_noother) + 
    geom_histogram(aes(x=Group,y=Total), stat = "identity", fill=colorsx, color="black")+
    scale_fill_manual(values=colorsx, 
                      name="Function classes" ,
                      labels=c(group_namesx), guide=F)+
    annotate("text",label=totnum, x=length(GROUPS_user)/2, y=max_fun+2,size=val10-18)+
    xlab("GO classes")+
    ylab("#GO groups")+
    scale_y_continuous(limits=c(0,max_fun+4))+
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
  functions_hist_toplotdf_noother=functions_hist_toplotdf_B[2:length(colors_user),]
  totnum=paste("Unique B GO groups",sum(functions_hist_toplotdf_noother$Represented), sep=": ")
  for ( i in 1:length(functions_hist_toplotdf_noother[,1])){
    functions_hist_toplotdf_noother[i,2]=as.numeric(functions_hist_toplotdf_noother[i,2])-1
  }
  functions_hist_toplotdf_noother$Group=factor(functions_hist_toplotdf_noother$Group, levels=c(as.character(1:length(GROUPS_user))))
  
  max_fun=max(functions_hist_toplotdf_noother$Total)
  
  
  UNIBhist_plot=ggplot(data=functions_hist_toplotdf_noother) + 
    geom_histogram(aes(x=Group,y=Total), stat = "identity", fill=colorsx, color="black")+
    scale_fill_manual(values=colorsx, 
                      name="Function classes" ,
                      labels=c(group_namesx), guide=F)+
    annotate("text",label=totnum, x=length(GROUPS_user)/2, y=max_fun+2, size=val10-18)+
    xlab("GO classes")+
    ylab("#GO groups")+
    scale_y_continuous(limits=c(0,max_fun+4))+
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
  
  functions_hist_toplotdf_noother=functions_hist_toplotdf_C[2:length(colors_user),]
  totnum=paste("Common GO groups",sum(functions_hist_toplotdf_noother$Represented), sep=": ")
  for ( i in 1:length(functions_hist_toplotdf_noother[,1])){
    functions_hist_toplotdf_noother[i,2]=as.numeric(functions_hist_toplotdf_noother[i,2])-1
  }
  functions_hist_toplotdf_noother$Group=factor(functions_hist_toplotdf_noother$Group, levels=c(as.character(1:length(GROUPS_user))))
  
  
  
  max_fun=max(functions_hist_toplotdf_noother$Total)
  
  COMMhist_plot=ggplot(data=functions_hist_toplotdf_noother) + 
    geom_histogram(aes(x=Group,y=Total), stat = "identity", fill=colorsx, color="black")+
    scale_fill_manual(values=colorsx, 
                      name="Function classes" ,
                      labels=c(group_namesx), guide=F)+
    annotate("text",label=totnum, x=length(GROUPS_user)/2, y=max_fun+2,size=val10-18)+
    xlab("GO classes")+
    ylab("#GO groups")+
    scale_y_continuous(limits=c(0,max_fun+4))+
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

########

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
