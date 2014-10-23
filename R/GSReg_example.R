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

length(diracpathways[[gs]])
sum(row.names(exprsdata) %in% diracpathways[[gs]])

test <- exprsdata %>%
    data.frame() %>%
    mutate(gene = as.factor(row.names(exprsdata))) %>%
    filter(gene %in% diracpathways[[gs]]) %>%
    group_by(gene) %>%
    summarise_each_(funs(max), list(quote(-gene)))

test2 <- test
names(test) <- as.character(phenotypes)

# function to map gene expression values to an individual gene set
map_gene_to_gs <- function(gs, gene_df) {
    gene_mat %>%
        mutate
        filter(gene %in% gs$genes)
}


