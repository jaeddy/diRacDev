
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
                      ".jsp?filePath=/resources/msigdb/4.0/", biocartaFile)
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

mapIdx <- f$symbol %in% geneIds(gs)
sum(mapIdx)
# probably replacing this with data frame functions
map_eset_to_geneset <- function(gs, es) {
    map_idx <- fData(es)$symbol %in% geneIds(gs)
    df <- data.frame(gene = fData(es)$symbol[map_idx],
                     exprs(es[map_idx])) %>%
        group_by(gene) %>%
        summarise_each(funs(max))
    df
}

test <- map_eset_to_geneset(gs, sample.ExpressionSet)

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
test <- es_melt %>%
    mutate(gene = get_symbol(probe, "hgu95av2")) %>%
    group_by(sample, gene, sex, type, score) %>%
    summarise(expression = max(expression))
