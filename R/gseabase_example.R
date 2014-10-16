
library(GEOquery)
library(GSEABase)
library(hgu95av2)
library(dplyr)
library(reshape2)

# Expression data for 500 features and 26 samples
data(sample.ExpressionSet)

# Extract 50 features
exprs <- sample.ExpressionSet[201:250, ]

# Create gene set collection
biocartaFile <- "c2.cp.biocarta.v4.0.symbols.gmt"
# fileAddress <- paste0("http://www.broadinstitute.org/gsea/msigdb/download_file",
#                      ".jsp?filePath=/resources/msigdb/4.0/", biocartaFile)
# download.file(fileAddress, paste0("data/", biocartaFile), method = "curl")

bcgsc <- getGmt(paste0("data/", biocartaFile),
                collectionType = BroadCollection(category = "c2"),
                geneIdType = SymbolIdentifier())
bcgsc[2]

fData(sample.ExpressionSet) <- featureData(sample.ExpressionSet) %>%
    row.names() %>%
    getSYMBOL("hgu95av2") %>%
    as.data.frame() %>%
    select(symbol = starts_with("featureData"))

for (i in 1:length(bcgsc)) {
    gs <- bcgsc[i]
    if (sum(geneIds(gs) %in% f) > 1) {
        print(i)
    }
}

f <- fData(sample.ExpressionSet)
gs <- bcgsc[["BIOCARTA_GLYCOLYSIS_PATHWAY"]]


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

gs_df <- map_gene_to_gs(gs, gene_df) # this could potentially be done with a join...


gs_mat <- gs_df %>%
    acast(gene ~ sample)

gs_ranks <- gs_mat %>%
    apply(2, rank)

gs_kt <- cor(gs_ranks, method = "kendall")


