## Import package
```
library("AtaCNV")
```

## Example data
The example data is collected from a public scATAC-seq dataset (GEO: GSE129785, sample SU008_tumor_post and SU008_tcell_post). Here, `count.rds` and `count_paired.rds` are the single-cell read count matrices for tumor and non-tumor samples, respectively, while `cell_info.rds` provides the annotations for cell types.
```
cell_info <- readRDS("./example/cell_info.rds")
count <- readRDS("./example/count.rds")
count_paired <- readRDS("./example/count_paired.rds")
```

## Input generation
The main input for AtaCNV is a bin-by-cell count matrix. The matrix is derived from the preprocessing results of 10X Cell Ranger ATAC, and rely on the fragment file and cell barcodes. In the paper, we further perform quality control using R package ArchR (https://www.archrproject.com/). This includes removing barcodes with low TSS (transcription start site) enrichment and potential doublets. \
For the example data, you can directly download the inputs from https://github.com/Aelita-Stone/AtaCNV/tree/main/example.
```
fragment_file <- "your_fragment_file.tsv.gz"
cell_barcodes <- read.csv("your_cell_barcodes.tsv", header=FALSE)
count <- generate_input(fragment_file = fragment_file, 
                        cell_barcodes = cell_barcodes,
                        genome = "hg19",
                        output_dir = "./")
```

## Normalization
Using count matrix as input, AtaCNV removes confounding factors like GC content to deconvolute CNA signals. Based on the available information, we provide four different modes to estimate the baseline read counts for non-tumor cells.
### Option1: mode="matched normal sample"
Matched non-tumor sample is available
```
norm_re1 <- AtaCNV::normalize(count,
                              genome="hg19", mode="matched normal sample",
                              cell_cluster=cell_info$cluster,
                              count_paired=count_paired,
                              output_dir="./",
                              output_name="norm_re2.rds"
)
```

### Option2: mode="normal cells"
The normal cells in the sample are known.
```
norm_re2 <- AtaCNV::normalize(count,
                              genome="hg19", mode="normal cells",
                              normal_cells=(cell_info$cell_type=="normal"),
                              output_dir="./",
                              output_name="norm_re2.rds"
)
```

### Option3: mode="all cells"
Construct the non-tumor baseline using all cells. This mode is only applicable when the tumor purity is low.
```
norm_re3 <- AtaCNV::normalize(count,
                              genome="hg19", mode="all cells",
                              cell_cluster=cell_info$cluster,
                              output_dir="./",
                              output_name="norm_re3.rds"
)
```

### Option4: mode="none"
No additional information. AtaCNV will automatically identify the most likely normal cells in the sample.
```
norm_re4 <- AtaCNV::normalize(count,
                              genome="hg19", mode="none",
                              cell_cluster=cell_info$cluster,
                              output_dir="./",
                              output_name="norm_re4.rds"
)
```

**Note**: Mode 1 is the most reliable in ensuring the accuracy. A feasible approach is to use software like ArchR to estimate marker gene expression from scATAC-seq data for a preliminary identification of normal cells.\
**Note**: Make sure that the `genome` parameter is consistent with the reference genome used during the sequence alignment. It can be set to hg19, hg38, or mm10.\
**Optional**: If there is no pre-existing clustering result, you can ommit the `cell_cluster` parameter, and AtaCNV will automatically perform clustering based on the count matrix.

## Segmentation
At this step, AtaCNV uses the BICseq algorithm to infer CNA breakpoints and estimate the copy number ratio for each segment. The input is the `norm_re` object obtained from the previous step.
```
seg_re <- calculate_CNV(norm_count=norm_re1$norm_count,
                        baseline=norm_re1$baseline,
                        output_dir="./",
                        output_name="seg_re.rds")
```



## Copy state inference
At this step, AtaCNV infers shared CNA states (loss, neutral or gain) for all tumor cells. This step requires providing normal and tumor cell annotation to the `cell_type` parameter. If there is no pre-defined annotation, you can classify normal and tumor cells after clustering cells based on copy ratios. In the resulting matrix, 0.5, 1, and 1.5 correspond to copy number loss, neutral, and gain, respectively.\
Additionally, if the heatmap indicates the presence of subclones within tumor cells, clustering can also be employed to obtain subclone labels, which can then be assigned to `cell_type`. Then the function will infer the shared CNA states for each subclone class.
```
CN_state <- estimate_cnv_state_cluster(count = count,
                                       genome = "hg19",
                                       copy_ratio = norm_re$copy_ratio,
                                       bkp = seg_re$bkp,
                                       label = cell_info$cell_type)
```

## Heatmap visualization
After estimating copy ratio, you can plot it using heatmap. The `cell_cluster` parameter is optional. If it is removed, the cell order will be determined by hierarchical clustering.
```
plot_heatmap(copy_ratio=seg_re$copy_ratio, 
             cell_cluster=cell_info$cluster,
             output_dir="./",
             output_name="copy_ratio.png")
```
Similarly, you can plot heatmap for copy state result.
```
plot_heatmap(copy_ratio=CN_state$cns, 
             cell_cluster=cell_info$cluster,
             output_dir="./",
             output_name="copy_state.png")
```

