# scRNA-seq data sets

## `library(scRNAseq)`

### Mouse retina data

- Consists of data sets from:
    1. Macosko et al.
    2. Shekhar et al.
- Two batches: Different laboratories both using Drop-seq
- Description found in Tran et al. (2020): Data set 7
- Details: [hemberg-lab](https://hemberg-lab.github.io/scRNA.seq.datasets/mouse/retina/)

```{r}
library(scRNAseq)

macosko <- MacoskoRetinaData()
shekhar <- ShekharRetinaData()
```

## Human pancreas data set

- Consists of data sets from:
    1. Baron et al.
    2. Muraro et al.
    3. Segerstolpe et al.
    4. Wang et al.
    5. Xin et al.
    6. Grun et al.
    7. Lawlor et al.
- Details: [hemberg-lab](https://hemberg-lab.github.io/scRNA.seq.datasets/human/pancreas/)
- Description found in Tran et al. (2020): Data set 4 

```{r}
library(scRNAseq)

baron <- BaronPancreasData()
muraro <- MuraroPancreasData()
seger <- SegerstolpePancreasData()
wang <- readRDS("processed/wang.rds")
xin <- XinPancreasData()
grun <- GrunPancreasData()
lawlor <- LawlorPancreasData()
```

## `library(SeuratData)`

### panc8

- Consists of eight pancreas data sets
- Five batches: Different sequencing platforms
- Description found in Tran et al. (2020): Data set 4 
- 13 cell types
- Dim: (34363, 14890)

```{r}
library(SeuratData)

pbmc3k <- InstallData("pbmc3k")
pbmc3k <- LoadData("pbmc3k")
```
## Gene Expression Omnibus / Others

### Human dendritic cells

- Villani et al. (2017)
- Description found in Tran et al. (2020): Data set 1
- Only data from the discovery set is used
- Sequencing platform: SmartSeq2 (no UMI tag) - TPM values reported
- Metadata obtained from cell IDs
- Two batches: Different runs
- Batch 1: Plates 7, 8, 9, 10; Batch 2: Plates 3, 4, 13, 14
- Dim: (26593, 768)

## Other sources

- `library(ExperimentHub)`
