
## LOAD LIBRARIES

library(org.Hs.eg.db) 
library("clusterProfiler")


## LOAD IMPUT DATA (GENE VECTOR TO FUNCTIONALLY ENRICH)

gene_symbol_vecA_df=read.csv("gene_symbol_vecA_df.csv", header=T, stringsAsFactors = F)


## SIMPLEST WAY TO DO THE FEA  
# enrichGO function requires ENTREZ_GENE IDs


## GENE OR PROTEIN IDS ACCEPTED
columns(org.Hs.eg.db)

# If necessary, These 4 lines map UNIPROT IDs to ENTREZ IDs
ENTREZ_LIST_A=mapIds(org.Hs.eg.db, keys=(as.vector(gene_symbol_vecA_df[,1])),column="ENTREZID", keytype = "SYMBOL", multivals="first")
ENTREZ_LIST_A=as.vector(as.numeric(ENTREZ_LIST_A))


# Functional Enrichment Analysis (FEA) - BIOLOGICAL PROCESS
Aentrezgenes_enrichGOres=as.data.frame(enrichGO(ENTREZ_LIST_A,'org.Hs.eg.db',ont="BP"))

# Functional Enrichment Analysis (FEA) - CELLULAR COMPONENT
Aentrezgenes_enrichGOres=as.data.frame(enrichGO(ENTREZ_LIST_A,'org.Hs.eg.db',ont="CC"))

# Functional Enrichment Analysis (FEA) - MOLECULAR FUNCTION
Aentrezgenes_enrichGOres=as.data.frame(enrichGO(ENTREZ_LIST_A,'org.Hs.eg.db',ont="MF"))

# Functional Enrichment Analysis (FEA) - KEGG ANNOTATION
Aentrezgenes_enrichGOres=as.data.frame(enrichKEGG(ENTREZ_LIST_A, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH"))

### FEA OTHER SPECIES

gene_symbol_vecA_df=read.csv("gene_symbol_vecA_df.csv", header=T, stringsAsFactors = F)
DATA_FRAME_UNIVERSE_TERM_GENES=read.csv("DATA_FRAME_UNIVERSE_TERM_GENES.csv", header=T, stringsAsFactors = F)
VECTOR_UNIVERSE_GENES=as.vector(DATA_FRAME_UNIVERSE_TERM_GENES[,2])

FEA_other_specie=as.data.frame(enricher(as.vector(c(gene_symbol_vecA_df[,1])), universe =VECTOR_UNIVERSE_GENES, TERM2GENE = DATA_FRAME_UNIVERSE_TERM_GENES))




