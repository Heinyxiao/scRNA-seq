### Load single-cell RNA-seq dataset

**10X data**
```
scRNA.counts <- Read10X(data.dir = "filtered_feature_bc_matrix")
```

**Expression Matrix**
```
Patient_scRNA.counts <- read.delim("GSE130000/GSM3729170_P1_dge.txt.gz", row.names = 1)
```

### Create Seurat Object
```
Patient_OC <- CreateSeuratObject(counts = Patient_scRNA.counts, project = "T59")
```

### Merge Seurat Objects
```
Patient_OC <- merge(ASCAF_scRNA, 
                   y = CAF23_scRNA, 
                   add.cell.ids = c("ASCAF", "CAF23"), 
                   project = "Patient_OC",
                   merge.data = TRUE)
#Check if the name was added
head(colnames(Patient_OC))
#Check if two origin was added
table(Patient_OC$orig.ident) 
```

### Start Quality Control
```
dir.create("QC")

# Mitochondria genes % 
Patient_OC[["percent.mt"]] <- PercentageFeatureSet(Patient_OC, pattern = "^MT-")

# Red blood cell genes % 
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(Patient_OC@assays$RNA)) 
HB.genes <- rownames(Patient_OC@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
Patient_OC[["percent.HB"]]<-PercentageFeatureSet(Patient_OC, features=HB.genes) 
head(Patient_OC@meta.data)
col.num <- length(levels(Patient_OC@active.ident))

# Draw violin graph, before QC 
violin <- VlnPlot(Patient_OC,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, 
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("QC/vlnplot_before_qc.png", plot = violin, width = 12, height = 6)  
plot1 <- FeatureScatter(Patient_OC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Patient_OC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(Patient_OC, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
ggsave("QC/pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5) 
ggsave("QC/pearplot_before_qc.png", plot = pearplot, width = 12, height = 5)

# Setup threshold for QC
minGene=500
maxGene=8000
pctMT=10

# Start QC 
Patient_OC <- subset(Patient_OC, 
                    subset = 
                      nFeature_RNA > minGene
                    & nFeature_RNA < maxGene 
                    & percent.mt < pctMT)
col.num <- length(levels(Patient_OC@active.ident))
violin <-VlnPlot(Patient_OC,
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.1, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

# Draw violin graph, after QC
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("QC/vlnplot_after_qc.png", plot = violin, width = 12, height = 6)

# Save RDS file
saveRDS(Patient_OC, file="Patient_OC.rds")

```