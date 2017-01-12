#' A DraAnno function for adding Deinococcus radiodurans annotation
#'
#'DraAnno.edgeR takes 2 argument: the name of the edgeR processed csv file(containing "logFC","logCPM","LR","PValue","FDR" variables), and in second argument, select "up" and "down" to get either up regulation or down regulation (cutoff parameters are fixed at FDR < 0.05, logFC <= -1 or >=1). Default of write is TRUE, so it automatically writes an annotion term added csv file to your current directory. You can defined the csv file name in "...".
#'
#' @import biomartr
#' @import dplyr
#' @import KEGGREST
#'
#'
#'
#'
#'@examples
#'\dontrun{DraAnno.edgeR(drad_24_diff_edgeR, "up", "up24hrannotation")}
#'
#'
#'
#'
#'@export
DraAnno.deseq <- function(dataset, select, ..., write = T ){



        Sig <- read.csv(dataset)
        names(Sig)[1] <- "refID"

        up <-  Sig %>%
                filter(padj <= 0.05 & log2FoldChange >= 1) %>%
                select(refID, log2FoldChange)

        down <-  Sig %>%
                filter(padj <= 0.05 & log2FoldChange <= -1) %>%
                select(refID, log2FoldChange)

        ##=============================================================
        ## Construct annotation database from NCBI


        dra.proteome <- getProteome(db = "refseq", organism = "Deinococcus radiodurans", path = file.path("_ncbi_downloads", "proteome"))

        pro <- read_proteome(dra.proteome, format = "fasta")
        pro.anno <- strsplit(pro@ranges@NAMES, ".1 ", fixed = T)

        pro.index <- list(method = "vector")
        for(i in 1:length(pro.anno)){
                pro.index[i] <- pro.anno[[i]][1]
        }

        pro.index <- unlist(pro.index)

        pro.annotation <- list(method = "vector")
        for(i in 1:length(pro.anno)){
                pro.annotation[i] <- pro.anno[[i]][2]
        }

        pro.annotation <- unlist(pro.annotation)

        # set strings as characters
        annotation_df <- data.frame(proteinID = pro.index, annotation = pro.annotation, stringsAsFactors = F)

        ##===========================================================
        ## Construct a dataset that can bridge proteinID and refID
        kg <- keggConv("dra", "ncbi-proteinid")
        kg_split <- strsplit(names(kg), ":", fixed = TRUE)


        ##protein ID
        kg_proteinID <- list(method = "vector")
        for(i in 1:length(kg)){
                kg_proteinID[i] <- kg_split[[i]][2]
        }
        kg_proteinID <- unlist(kg_proteinID)

        ##refID
        kg_refID <- list(method = "vector")
        for(i in 1:length(kg)){
                split <- strsplit(kg, ":", fixed = TRUE)
                kg_refID[i] <- split[[i]][2]
        }

        kg_refID <- unlist(kg_refID)

        # set strings as characters
        kg.index <- data.frame(refID = kg_refID, proteinID = kg_proteinID, stringsAsFactors = F)
        ##========================================================

        ## gene ID
        geneid <- keggConv("dra","ncbi-geneid")


        #x has to be a keggConv object
        get_geneID <- function(x = geneid){
                #x has to be a keggConv object
                genex <- vector()
                refx <- vector()

                for ( i in 1:length(x)){
                        genex[i] <- sapply(strsplit(names(x[i]), ":"), "[", 2)
                        refx[i] <- sapply(strsplit(x[i], ":"), "[", 2)
                }

                data.frame(refID = refx, geneID = genex, stringsAsFactors = F)
        }

        kg.geneid <- get_geneID(geneid)


        ### annotation function

        dra_anno <- function(x){
                # our dataset has to set first column to refID
                # has to set as characters or it returns NA
                sample_ref <- x$refID
                lfc <- x$logFC
                kg_ref <- kg.index$refID
                gene_ref <- kg.geneid$refID # geneid
                anno_ref <- annotation_df$proteinID

                logfc <- vector()
                proID <- vector()
                geneid <- vector()
                feature <- vector()

                for(i in 1:length(sample_ref)){

                        logfc[i] <- lfc[i]

                        if(sample_ref[i] %in% kg_ref){


                                proID[i] <-  kg.index[kg.index$refID == sample_ref[i],]$proteinID


                        } else {

                                proID[i] <- "Not in db"
                        }

                        if(sample_ref[i] %in% gene_ref ){

                                geneid[i] <- kg.geneid[kg.geneid$refID == sample_ref[i],]$geneID

                        } else {

                                geneid[i] <- "Not in db"

                        }

                        if(proID[i] %in% anno_ref){

                                feature[i] <- annotation_df[annotation_df$proteinID == proID[i],]$annotation

                        } else {

                                feature[i] <- "Not in db"
                        }
                }

                data.frame(refseq_locus = sample_ref, geneID = geneid, proteinID = proID, logFC = logfc, annotation = feature)

        }


        if(select == "up" & write == TRUE){

                annotated <- dra_anno(up)
                write.csv(annotated, ...)


        } else if( select == "down" & write == TRUE){

                annotated <- dra_anno(down)
                write.csv(annotated, ...)


        } else if( select == "up" & write == FALSE){

                annotated <- dra_anno(up)


                warning("Return 'up' regulated genes but without creating a csv file")

        } else if( select == "down" & write == FALSE){

                annotated <- dra_anno(down)


                warning("Return 'down' regulated genes but without creating a csv file")

        } else {

                warning("Define data set of 'up' or 'down' genes . If Write = F: not to write a .csv file of the annotated data.")
        }

        annotated

}


