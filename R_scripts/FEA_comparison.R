
## LOAD LIBRARIES

library(org.Hs.eg.db) 
library("clusterProfiler")


## LOAD IMPUT DATA (GENE VECTOR TO FUNCTIONALLY ENRICH)
setwd("~/Documents/R_work/FEA_simplification_comparison_scripts/input_data")
gene_symbol_vecA_df=read.csv("gene_symbol_vecA_df.csv", header=T, stringsAsFactors = F)
gene_symbol_vecB_df=read.csv("gene_symbol_vecB_df.csv", header=T, stringsAsFactors = F)


# BEFORE doing any analysis, go to line titled "# ASSOCIATED FUNCTIONS " and run all the lines until the end of the file

# mapping and simplest FEA for the two sets independently

# SET A

# If necessary, These 4 lines map UNIPROT IDs to ENTREZ IDs
ENTREZ_LIST_A=mapIds(org.Hs.eg.db, keys=(as.vector(gene_symbol_vecA_df[,1])),column="ENTREZID", keytype = "SYMBOL", multivals="first")
ENTREZ_LIST_A=as.vector(as.numeric(ENTREZ_LIST_A))
# Functional Enrichment Analysis (FEA) - BIOLOGICAL PROCESS
uniprotA_BP07_res_A=as.data.frame(enrichGO(ENTREZ_LIST_A,'org.Hs.eg.db',ont="BP"))
uniprotA_KEGG07_res_A=as.data.frame(enrichKEGG(ENTREZ_LIST_A,organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH"))
# SET B

# If necessary, These 4 lines map UNIPROT IDs to ENTREZ IDs
ENTREZ_LIST_B=mapIds(org.Hs.eg.db, keys=(as.vector(gene_symbol_vecB_df[,1])),column="ENTREZID", keytype = "SYMBOL", multivals="first")
ENTREZ_LIST_B=as.vector(as.numeric(ENTREZ_LIST_B))
# Functional Enrichment Analysis (FEA) - BIOLOGICAL PROCESS
uniprotA_BP07_res_B=as.data.frame(enrichGO(ENTREZ_LIST_B,'org.Hs.eg.db',ont="BP"))
uniprotA_KEGG07_res_B=as.data.frame(enrichKEGG(ENTREZ_LIST_B,organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH"))




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

FEAs_comparison_outputs_KEGG=plural_funct_simplif3(uniprotA_KEGG07_res_A,uniprotA_KEGG07_res_B,80,"Lin",0.7, d_KEGG)
FINAL_OUTPUT=FEAs_comparison_outputs_KEGG$both_simp

FEAs_comparison_outputs=plural_funct_simplif3(uniprotA_BP07_res_A,uniprotA_BP07_res_B,101,"Lin",0.7, d_BP)
FINAL_OUTPUT=FEAs_comparison_outputs$both_simp

# function groups (rows) with values ONLY in columns denoted with "A" will be functions only associated to gene list A, 
# function groups (rows) with values ONLY in columns denoted with "B" will be functions only associated to gene list B, 
# function groups (rows) with values IN BOTH columns denoted with "A" and "B" will be functions associated to BOTH gene lists




plural_funct_simplif3=function(SET_DATAFRAME_A, SET_DATAFRAME_B,gene_cooc_threshold, semantic_alg, semantic_threshold, semData_info){

  AB_func_enrichment_res=JOINING_COMPARING_DATASETS(SET_DATAFRAME_A,SET_DATAFRAME_B)
  geneco_AB_func_enrichment=GENECO_MERGING_GOT(AB_func_enrichment_res[[1]],AB_func_enrichment_res[[2]],AB_func_enrichment_res[[3]],AB_func_enrichment_res[[4]], gene_cooc_threshold, gene_cooc_threshold)
  co_oc_message=paste((length(AB_func_enrichment_res[[1]][,1])-length(geneco_AB_func_enrichment[,1])),"/", length(AB_func_enrichment_res[[1]][,1]), "functions merged by gene co-occurrence with a threshold of",gene_cooc_threshold, sep=" ", collapse = " ")
  
  if (class(semData_info)=="GOSemSimDATA") {
    sem_AB_func_enrichment=SEMANTIC_ANALYSIS5(semantic_alg,semantic_threshold,geneco_AB_func_enrichment,semData_info)
    final_simplified_res=sem_AB_func_enrichment[[1]]
    print(paste("the FEAs comparison returns", (AB_func_enrichment_res[[3]]),"/", length(SET_DATAFRAME_A[,1]), "unique functions associated to SET A, ",
                (AB_func_enrichment_res[[4]]), "/", length(SET_DATAFRAME_B[,1]), " to SET B, and ", AB_func_enrichment_res[[2]], "common to both AB . This gives a total of ", length(AB_func_enrichment_res[[1]][,1]),
                "not simplified functions. Then, after simplification, we get ",length(final_simplified_res[,1]),"functional groups",sep=" ", collapse = " "  ))
    print(co_oc_message)
    print(sem_AB_func_enrichment[[2]])
    
  }
  if (class(semData_info)!="GOSemSimDATA"){
    sem_AB_func_enrichment=0
    final_simplified_res=geneco_AB_func_enrichment
    print(paste("the FEAs comparison returns", (AB_func_enrichment_res[[3]]),"/", length(SET_DATAFRAME_A[,1]), "unique functions associated to SET A, ",
                (AB_func_enrichment_res[[4]]), "/", length(SET_DATAFRAME_B[,1]), " to SET B, and ", AB_func_enrichment_res[[2]], "common to both AB . This gives a total of ", length(AB_func_enrichment_res[[1]][,1]),
                "not simplified functions. Then, after simplification, we get ",length(final_simplified_res[,1]),"functional groups",sep=" ", collapse = " "  ))
    
    sem_message=paste("KEGG functions do not have hierarchy thus,", 0, "functions merged by semantic similarity with a threshold of",semantic_threshold, "and", semantic_alg, "algorithm" , sep=" ", collapse = " ")
    print(co_oc_message)
    print(sem_message)
  
  }
  
  list("original_FEA_A"=SET_DATAFRAME_A,"original_FEA_B"=SET_DATAFRAME_B,"merged_FEA_AB"=AB_func_enrichment_res[[1]], "gene_coocc_simp"=geneco_AB_func_enrichment,"sem_sim_simp"=sem_AB_func_enrichment,"both_simp"=final_simplified_res)

}



#SIMPandCOMP_OUTPUT_clean=EDITING_TO_PLOT_noS2B(length(ENTREZ_LIST_A),length(ENTREZ_LIST_B) ,SIMPandCOMP_OUTPUT[[4]])
#SIMPandCOMP_OUTPUT_clean_res=list("simplification_results_list"=SIMPandCOMP_OUTPUT,"simplified_genefreqandfold_df"=SIMPandCOMP_OUTPUT_clean)






# ASSOCIATED FUNCTIONS - run these lines BEFORE doing any analysis


##############

DATAFORMAT_CHANGE=function(Aentrezgenes_enrichres, bg_threshold){
  filtering_A_enrichGO=FILTERING_BACKGROUND2(Aentrezgenes_enrichres,bg_threshold) 
  filtering_A_enrichGO_res=fold_computing2(filtering_A_enrichGO)
  filtering_A_enrichGO_res2=cbind(filtering_A_enrichGO_res,"genes"=0)
  list("entrez_id_df"=0 , "filtered_FEA"=filtering_A_enrichGO_res, "FEA_result"=filtering_A_enrichGO_res2)
}

###

SEMANTIC_ANALYSIS5=function(alg_measure,sim_thresh,initial_df, semData_info){
  AB_GOT=initial_df[,1]
  AB_GOT=c(as.character(AB_GOT))
  AB_testing=termSim(AB_GOT, AB_GOT,semData=d_BP,method = alg_measure)
  tst_values=sm2vec(AB_testing, diag = FALSE)
  tst_index=sm.index(AB_testing, diag = F)
  order_tst_values=order(tst_values,decreasing = T)
  sorted_tst_values=tst_values[order_tst_values]
  sorted_tst_index=tst_index[order_tst_values,]
  e_df=rbind(sorted_tst_index[sorted_tst_values>=0.7,])
  
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
          #minitial_df[indexX,3]=paste(initial_df[indexY,3],initial_df[indexX,9],sep="/")
          #minitial_df[indexX,4]=paste(initial_df[indexY,4],initial_df[indexX,10],sep="/")
          minitial_df[indexX,9]=paste(initial_df[indexY,3],initial_df[indexX,9],sep="/")
          minitial_df[indexX,10]=paste(initial_df[indexY,4],initial_df[indexX,10],sep="/")
          #minitial_df[indexX,5]=paste(initial_df[indexY,5],initial_df[indexX,11],sep="/")
          #minitial_df[indexX,6]=paste(initial_df[indexY,6],initial_df[indexX,12],sep="/")
          minitial_df[indexX,11]=paste(initial_df[indexY,5],initial_df[indexX,11],sep="/")
          minitial_df[indexX,12]=paste(initial_df[indexY,6],initial_df[indexX,12],sep="/")
        }
      }
    }
    del_ind=which(e_df_summary[,3]==1) #####
    minitial_df=minitial_df[-del_ind,] #####
    num_semsim_merged=length(initial_df[,1])-length(minitial_df[,1])
    sem_message=paste(num_semsim_merged, "/",length(initial_df[,1]), "functions merged by semantic similarity with a threshold of",sim_thresh, "and", alg_measure, "algorithm" , sep=" ", collapse = " ")
    minitial_df
  }
  else {
    minitial_df=initial_df
    sem_message=paste( "0 functions merged by semantic similarity with a threshold of",sim_thresh, "and", alg_measure, "algorithm" , sep=" ", collapse = " ")
    minitial_df
  }
  list(minitial_df, sem_message)
}
###########
SEMANTIC_ANALYSIS4=function(alg_measure,sim_thresh,initial_df, semData_info){
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
    sem_message=paste(num_semsim_merged, "/",length(initial_df[,1]), "functions merged by semantic similarity with a threshold of",sim_thresh, "and", alg_measure, "algorithm" , sep=" ", collapse = " ")
    minitial_df
  }
  else {
    minitial_df=initial_df
    sem_message=paste( "0 functions merged by semantic similarity with a threshold of",sim_thresh, "and", alg_measure, "algorithm" , sep=" ", collapse = " ")
    minitial_df
  }
  list(minitial_df, sem_message)
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


#####


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

#########

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


#######

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
