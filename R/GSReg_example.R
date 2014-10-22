#library(GSBenchMark)
#library(GSReg)

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



# function to convert gene set collection to data frame
gsc_to_df <- function(gsc) {
    df <- data.frame()
    for (gs in gsc) {
        df <- rbind(df, data.frame(gs = setName(gs), genes = geneIds(gs)))
    }
    df
}

gsc_df <- gsc_to_df(bcgsc)

# convert expression set to data frame
es_df <- as.data.frame(sample.ExpressionSet)
names(es_df) <- gsub("^X", "", names(es_df))
es_df$sample <- row.names(es_df)

# use melt to put data in 'long' form
es_melt <- es_df %>%
    melt(id.vars = c("sex", "type", "score", "sample"),
         variable.name = "probe", value.name = "expression")

# function to convert probe IDs to gene symbols
get_symbol <- function(x, data) {
    x <- as.character(x)
    getSYMBOL(x, data)
}

# get max expression across probes for each sample and gene
gene_df <- es_melt %>%
    mutate(gene = get_symbol(probe, "hgu95av2")) %>%
    group_by(sample, sex, type, score, gene) %>%
    summarise(expression = max(expression)) %>%
    as.data.frame()

# gene set for testing
for (bcgs in bcgsc) {
    gs_name <- setName(bcgs)
    gs <- gsc_df %>%
        filter(gs == gs_name)
    if (sum(gs$genes %in% gene_df$gene) > 2) {
        print(gs_name)
        print(sum(gs$genes %in% gene_df$gene))
    }
}
gs <- gsc_df %>%
    filter(gs == "BIOCARTA_TEL_PATHWAY")
gs

# function to count genes mapped to each gene set (for testing)
count_mapped_genes <- function(gs, gene_df) {
    gs_df <- gene_df %>%
        filter(gene %in% gs$genes) %>%
        select(gene) %>%
        distinct()
    nrow(gs_df)
}

count_mapped_genes(gs, gene_df)

# function to map gene expression values to an individual gene set
map_gene_to_gs <- function(gs, gene_df) {
    gene_df %>%
        filter(gene %in% gs$genes)
}

gs_df <- map_gene_to_gs(gs, gene_df)