###GSE116899
DF <- read.table("Data/GSE116899_ra-neut-counts-EnsembIDs-GRCh37.p10.txt", header=T )

Annotation = data.frame(Sample = colnames(DF)[7:36])

Annotation$Donor = Annotation$Sample
  Annotation$Donor <- gsub("X|VDBR|b|sf", "", gsub("_C2CUGACXX_[AGCT]*_L00[12].accepted_hits.bam", "", Annotation$Donor))
  # Annotation$Tissue <- gsub("X|VDBR|[0-9]+", "", gsub("_C2CUGACXX_[AGCT]*_L00[12].accepted_hits.bam", "", Annotation$Tissue))
  # Annotation$Tissue <- gsub("b", "Blood", Annotation$Tissue)
  # Annotation$Tissue <- gsub("sf", "SynovialFluid", Annotation$Tissue)
  
  ##################
  #Create Expression set
  ##################
  gset <- new("ExpressionSet", exprs = ((as.matrix(DF[,7:ncol(DF)]))))
  fData(gset) <- DF[,1:6]
  #Add max value per gene
  fData(gset)$maxExpr <- apply(exprs(gset), 1, max)
  pData(gset) <- data.frame(Annotation)
    gset$Disease <- "rheumatoidarthritis"
    gset$Disease <- ifelse(gset$Donor %in% c(387,396,10,40,43), "Healthy", "RheumatoidArthritis")
    gset$Cell <- "Nuetrohils"
  library(biomaRt)
  mart<- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))
  Ann <-   getBM(filters="ensembl_gene_id",
                 attributes=c("ensembl_gene_id",   "hgnc_symbol"),
                 values=fData(gset)$Geneid[1:10],
                 mart=mart)
  
  #Test Row names and column names match  
  table(row.names(pData(gset)) == colnames(exprs(gset)))


  Ann <-   read.table("/app/Shiny/RAExpressionDashboard/Data/GeneSymbol.txt", sep="\t", header=T)
  names(Ann) <- c("Geneid", "Gene.name")
  library(plyr)
  jj <- join(fData(gset)[,c(1,6,7)], Ann, by="Geneid")
  fData(gset) <- jj
  
  
  save (gset, file="Data/GSE116899.rdata")
  
  