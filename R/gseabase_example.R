library(GSEABase)
library(hgu95av2)

# Expression data for 500 features and 26 samples
data(sample.ExpressionSet)

# Extract 50 features
exprs <- sample.ExpressionSet[201:250, ]

# Create gene set collection
gsc <- GeneSetCollection(exprs, setType = BroadCollection())
source("http://bioconductor.org/biocLite.R")
