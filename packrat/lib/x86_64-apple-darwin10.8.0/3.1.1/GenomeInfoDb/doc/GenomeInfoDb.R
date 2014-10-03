
## ----style, eval=TRUE, echo=FALSE, results="asis"------------------------
BiocStyle::latex()


## ----preliminaries, echo=FALSE, message=FALSE----------------------------
library(GenomeInfoDb)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)


## ----genomeStyles1-------------------------------------------------------
seqmap <- genomeStyles()
head(seqmap,n=2)


## ----name----------------------------------------------------------------
names(genomeStyles())


## ----genomeStyles2-------------------------------------------------------
head(genomeStyles("Homo_sapiens"),5)


## ----style-present-------------------------------------------------------
"UCSC" %in% names(genomeStyles("Homo_sapiens"))


## ----extractSeqlevels----------------------------------------------------
extractSeqlevels(species="Arabidopsis_thaliana", style="NCBI")


## ----extractSeqlevelsgroup-----------------------------------------------
extractSeqlevelsByGroup(species="Arabidopsis_thaliana", style="NCBI",
                         group="auto")


## ----seqlevelsStyle------------------------------------------------------
seqlevelsStyle(paste0("chr",c(1:30)))
seqlevelsStyle(c("2L","2R","X","Xhet"))


## ----keepChr-txdb--------------------------------------------------------
newchr <- paste0("chr",c(1:22,"X","Y","M","1_gl000192_random","4_ctg9_hap1"))
seqlevelsInGroup(newchr, group="sex")
seqlevelsInGroup(newchr, group="auto")
seqlevelsInGroup(newchr, group="circular")
seqlevelsInGroup(newchr, group="sex","Homo_sapiens","UCSC")


## ----check2--------------------------------------------------------------
seqnames <- c("chr1", "chr9", "chr2", "chr3", "chr10")
all(seqnames %in% extractSeqlevels("Homo_sapiens", "UCSC"))


## ----orderSeqlevels------------------------------------------------------
seqnames <- c("chr1","chr9", "chr2", "chr3", "chr10")
orderSeqlevels(seqnames)


## ----rankSeqlevels-------------------------------------------------------
seqnames <- c("chr1","chr9", "chr2", "chr3", "chr10")
rankSeqlevels(seqnames)


## ----find----------------------------------------------------------------
mapSeqlevels(c("chrII", "chrIII", "chrM"), "NCBI")


## ----quick-style---------------------------------------------------------
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
seqlevels(txdb)
genomeStyles("Drosophila melanogaster")
mapSeqlevels(seqlevels(txdb), "NCBI")


## ----sequence, eval=FALSE------------------------------------------------
## sequence <- seqlevels(x)
## 
## ## sequence is in UCSC format and we want NCBI style
## newStyle <- mapSeqlevels(sequence,"NCBI")
## newStyle <- newStyle[complete.cases(newStyle)] # removing NA cases.
## 
## ## rename the seqlevels
## x <- renameSeqlevels(x,newStyle)
## 
## ## keep only the seqlevels you want (say autosomes)
## auto <- extractSeqlevelsByGroup(species="Homo sapiens", style="NCBI",
##                                 group="auto")
## x <- keepSeqlevels(x,auto)


