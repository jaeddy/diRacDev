
library(GEOquery)
library(GSEABase)
library(hgu95av2)
library(dplyr)

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

f <- fData(sample.ExpressionSet)

allGenes <- unique(as.character(unlist(geneIds(bcgsc))))
bcGenes <- geneSymbols[geneSymbols %in% allGenes]

map_expression_to_set <- function(gs, es) {
    mapIdx <- geneIds(gs) %in% fData(es)$symbol
}

for (gs in bcgsc) {
    if (sum(geneIds(gs) %in% geneSymbols) > 1) {
        print(gs)
    }
}


gs <- bcgsc[["BIOCARTA_GLYCOLYSIS_PATHWAY"]]
geneIds(gs) %in% geneSymbols
sum(geneSymbols %in% geneIds(gs))
mapIdx <- geneSymbols %in% geneIds(gs)

test <- data.frame(gene = geneSymbols[mapIdx], 
                   exprs(sample.ExpressionSet[mapIdx]))

test2 <- test %>%
    group_by(gene) %>%
    summarise_each(funs(max))
