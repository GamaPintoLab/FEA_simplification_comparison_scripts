
## LOAD LIBRARIES

library(org.Hs.eg.db) 
library("clusterProfiler")

# BEFORE doing any analysis, go to line titled "# ASSOCIATED FUNCTIONS " and run all the lines until the end of the file


## LOAD IMPUT DATA (GENE VECTOR TO FUNCTIONALLY ENRICH)

gene_symbol_vecA_df=read.csv("gene_symbol_vecA_df.csv", header=T, stringsAsFactors = F)

# filtering by gene background (0.70 by default)

# Functional Enrichment Analysis (FEA) - BIOLOGICAL PROCESS
uniprotA_BP07_res_A=FEA_and_FILTERING(gene_symbol_vecA_df[,1],"SYMBOL","BP", 0.7, "UNIPROT")

# Functional Enrichment Analysis (FEA) - MOLECULAR FUNCTION
uniprotA_MF07_res_A=FEA_and_FILTERING(gene_symbol_vecA_df[,1],"SYMBOL","MF", 0.7, "UNIPROT")

# Functional Enrichment Analysis (FEA) - KEGG ANNOTATION
uniprotA_KEGG07_res_A=FEA_and_FILTERING(gene_symbol_vecA_df[,1],"SYMBOL","KEGG", 0.7, "UNIPROT")

uniprotA_BP07_res_A_r=uniprotA_BP07_res_A$filtered_FEA
uniprotA_MF07_res_A_r=uniprotA_MF07_res_A$filtered_FEA
uniprotA_KEGG07_res_A_r=uniprotA_KEGG07_res_A$filtered_FEA


##### FEA simplification

#  by a) gene co-occurrence, b) semantic similarity or c) both
# if you want to avoid one of themm , use a 1000 threshold

library("GOSemSim")
library("corpcor")

d_BP <- godata('org.Hs.eg.db', ont="BP", computeIC=T) 
d_MF <- godata('org.Hs.eg.db', ont="MF", computeIC=T) 
d_KEGG=0


# if you want both simplifications define your thresholds
# if you want only simp by gene-co occurrence , define semantic_threshold= 1.1
# if your want only simp by sem sim, define gene_cooc_threshold= 101


semantic_alg="Lin"
#semantic_threshold (from 0 to 1.1)
#gene_cooc_threshold (from 0 to 101)
semData_info=d_BP

# Functional Enrichment Analysis (FEA) SIMPLIFICATION - BIOLOGICAL PROCESS
single_g101s07simp_BP=single_funct_simplif3(uniprotA_BP07_res_A_r,101,"Lin", 0.7, d_BP)
simplifies_output_A=single_g101s07simp_BP$both_simp

# Functional Enrichment Analysis (FEA) SIMPLIFICATION - MOLECULAR FUNCTION
single_g101s07simp_MF=single_funct_simplif3(uniprotA_MF07_res_A_r,101,"Lin", 0.7, d_MF)
simplifies_output_A=single_g101s07simp_MF$both_simp

# Functional Enrichment Analysis (FEA) SIMPLIFICATION - KEGG ANNOTATION
single_g101s07simp_KEGG=single_funct_simplif3(uniprotA_KEGG07_res_A_r,101,"Lin", 0.7, d_KEGG)
simplifies_output_A=single_g101s07simp_KEGG$both_simp



# ASSOCIATED FUNCTIONS - run these lines BEFORE doing any analysis

single_funct_simplif3=function(RAW_FUNC_ENRICH_DATAFRAME_A, gene_cooc_threshold, semantic_alg, semantic_threshold, semData_info){
  
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
  
  list("original_FEA"=ordered_ABfun_enr,"gene_coocc_simp"=geneco_AB_func_enrichment,"sem_sim_simp"=sem_AB_func_enrichment,"both_simp"=final_simplified_res)
}



###

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