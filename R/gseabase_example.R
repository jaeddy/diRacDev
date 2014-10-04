<<<<<<< HEAD
<<<<<<< HEAD
library(GEOquery)
library(GSEABase)
library(hgu95av2)
library(dplyr)
=======
library(GSEABase)
library(hgu95av2)
>>>>>>> a26c92e3306cb55aca4f989abeb28b209846232e
=======
library(GSEABase)
library(hgu95av2)
>>>>>>> b9ac84f46928c2accee71eefc597608cdcaa1fad

# Expression data for 500 features and 26 samples
data(sample.ExpressionSet)

# Extract 50 features
exprs <- sample.ExpressionSet[201:250, ]

# Create gene set collection
<<<<<<< HEAD
<<<<<<< HEAD
biocartaFile <- "c2.cp.biocarta.v4.0.symbols.gmt"
fileAddress <- paste0("http://www.broadinstitute.org/gsea/msigdb/download_file",
                      ".jsp?filePath=/resources/msigdb/4.0/", biocartaFile)
download.file(fileAddress, paste0("data/", biocartaFile), method = "curl")

bcgsc <- getGmt(paste0("data/", biocartaFile),
                collectionType = BroadCollection(category = "c2"),
                geneIdType = SymbolIdentifier())
bcgsc[2]

features <- fData(exprs)

egsc <- GeneSetCollection(exprs, idType = "missing", setType = BroadCollection())

library(KEGG.db)
lst <- head(as.list(KEGGEXTID2PATHID))
gsc <- GeneSetCollection(mapply(function(geneIds, keggId) {
    GeneSet(geneIds, geneIdType=EntrezIdentifier(),
            collectionType=KEGGCollection(keggId),
            setName=keggId)
}, lst, names(lst)))

geneSymbols <- featureData(exprs) %>%
    row.names() %>%
    getSYMBOL("hgu95av2")
=======
gsc <- GeneSetCollection(exprs, setType = BroadCollection())
source("http://bioconductor.org/biocLite.R")
>>>>>>> a26c92e3306cb55aca4f989abeb28b209846232e
=======
gsc <- GeneSetCollection(exprs, setType = BroadCollection())
source("http://bioconductor.org/biocLite.R")
>>>>>>> b9ac84f46928c2accee71eefc597608cdcaa1fad
