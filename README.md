# GeneNetworkBuilder

[![platforms](http://bioconductor.org/shields/availability/devel/GeneNetworkBuilder.svg)](http://bioconductor.org/packages/devel/bioc/html/GeneNetworkBuilder.html)
[![build](http://bioconductor.org/shields/build/devel/bioc/GeneNetworkBuilder.svg)](http://bioconductor.org/packages/devel/bioc/html/GeneNetworkBuilder.html)

Build Regulatory Network from ChIP-chip/ChIP-seq and Expression Data

Appliation for discovering direct or indirect targets of transcription factors using ChIP-chip or ChIP-seq, and microarray or RNA-seq gene expression data. Inputting a list of genes of potential targets of one TF from ChIP-chip or ChIP-seq, and the gene expression results, GeneNetworkBuilder generates a regulatory network of the TF.

## Installation

To install this package, start R and enter:

```r
library(BiocManager)
BiocManager::install("GeneNetworkBuilder")
```

## Documentation

To view documentation of GeneNetworkBuilder, start R and enter:
```r
browseVignettes("GeneNetworkBuilder")
```

