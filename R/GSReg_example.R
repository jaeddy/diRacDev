library(GSBenchMark)
library(GSReg)
library(dplyr)

data(diracpathways)
data(GSBenchMarkDatasets)
print(GSBenchMark.Dataset.names)

dataSetName <- GSBenchMark.Dataset.names[[4]]
data(list = dataSetName)

nanGenes <- rowSums(is.nan(exprsdata)) > 0
exprsdata <- exprsdata[!nanGenes, ]

geneNames <- rownames(exprsdata)
head(geneNames)

?GSReg.GeneSets.DIRAC

diracResult <- GSReg.GeneSets.DIRAC(exprsdata, diracpathways, phenotypes, 10)
help(GSReg:::GSReg.DIRAC.Pathways)

sigPathways <- which(diracResult$pvalues < 0.05)
dysregulatedPathways <- rbind(diracResult$mu1[sigPathways],
                              diracResult$mu2[sigPathways], 
                              diracResult$pvalues[sigPathways]);
rownames(dysregulatedPathways) <- c("mu1", "mu2", "pvalues");
print(dysregulatedPathways[,1:5])

# convert results matrix to data frame
topPathways <- data.frame(t(dysregulatedPathways))
topPathways <- cbind(pathway = row.names(topPathways), topPathways)

gs <- as.character(topPathways$pathway[1])

# function to map gene expression values to an individual gene set
map_gene_to_gs <- function(gs, gene_mat) {
    gene_mat %>%
        data.frame() %>%
        mutate(gene = as.factor(row.names(exprsdata))) %>%
        filter(gene %in% diracpathways[[gs]]) %>%
        group_by(gene) %>%
        summarise_each_(funs(max), list(quote(-gene)))
}

gs_df <- map_gene_to_gs(gs, exprsdata)

label_samples <- function(gs_df, phenotypes) {
    classes <- as.numeric(phenotypes)
    labels <- paste0(classes, 
                     rep("_", length(classes)),
                     as.character(c(1:sum(classes == 1), 
                                    1:sum(classes == 2))))
    labels <- c("gene", labels)
    names(gs_df) <- labels    
    row.names(gs_df) <- gs_df$gene
    gs_df %>%
        select(-gene)
}

gs_df <- label_samples(gs_df, phenotypes)
