# scRNA-seq data sets

## Sources

- `library(SeuratData)`
- `library(scRNAseq)`
- `library(ExperimentHub)`

## Data sets

## Human dendritic cells

- Villani et al. (2017) 
- Tran et al. (2020): Data set 1

## Human pancreas data set

- [Details: hemberg-lab](https://hemberg-lab.github.io/scRNA.seq.datasets/human/pancreas/)
- Tran et al. (2020): Data set 5

### References

1. Baron
2. Muraro
3. Segerstolpe
4. Wang
5. Xin
6. Grun
7. Lawlor

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

## Mouse retina data set

- [Details: hemberg-lab](https://hemberg-lab.github.io/scRNA.seq.datasets/mouse/retina/)
- Tran et al. (2020): Data set 7

### References

```{r}
library(scRNAseq)

macosko <- MacoskoPancreasData()
shekhar <- ShekharPancreasData()
```

## PBMC (pbmc3k)

```{r}
library(SeuratData)

pbmc3k <- InstallData("pbmc3k")
```

## Pancreas data set (panc8)

## IFNB-stimulated and control PBMC (ifnb)
