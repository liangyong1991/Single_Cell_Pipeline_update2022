本流程是从 https://github.com/johnstonmj/single-cell-rna-seq-template 上下载的，修改了其中的部分代码、完善每个rule的分析环境，以及添加了安装、使用说明。

另外，本流程的输入是段细胞转录组数据的每个细胞中各个基因表达文件、metadata文件，以及marker文件，来对细胞进行分型和差异分析的流程。

## 所需环境
- snakemake(version=6.15.5)
- singularity version 3.8.7-1.el7
- docker
- git


## 安装
#### 从github上下载源代码
```{shell}
cd /你的工作目录

#下载源代码
git clone git@github.com:liangyong1991/Single_Cell_Pipeline_update2022.git

#进入到下载的repository
cd Single_Cell_Pipeline_update2022
```

#### 构建单个rule需要的分析环境
1. 直接使用别人构建好的singularity环境（.sif文件）
可以联系未知君梁勇（liangyong@xbiome.com）拷贝.sif文件（4个环境大概5G左右, 在未知君服务器位置为：/share/work/HPC/work_tmp/liangyong/github/Single_Cell_Pipeline_update2022/singularity_test/recipe）

2. 从yaml文件构建
- 环境eval
```{shell}
cd singularity_test/recipe
sudo singularity build eval.sif eval.def
```

- 环境cellassign(需要先封装成docker image；然后再用docker image转换成singularity image；当然也鼓励大家直接使用singularity build直接构建环境)
```{shell}
cd singularity_test/recipe/cellassign
docker build -t cellassign:t

SINGULARITY_NOHTTPS=1 

#docker image转换 singularity image(这里需要上传到远端的docker仓库)
sudo singularity build cellassign.sif docker:cellassign:t
```

- 环境heatmap
```{shell}
cd singularity_test/recipe
sudo singularity build heatmap.sif heatmap.def
```

- 环境edger
```{shell}
cd singularity_test/recipe
sudo singularity build edger.sif edger.def
```

#### 将每个rule(*smk)中需要的 singularity替换成新构建的image
这些*.smk文件的路径为:
```{shell}
ls -l  rules/*smk
-rw-rw-r-- 1 liangyong liangyong 1215 11月 14 15:18 rules/cell-cycle.smk
-rw-rw-r-- 1 liangyong liangyong 3081 11月 16 17:23 rules/cell-type.smk
-rw-rw-r-- 1 liangyong liangyong 1325 11月 13 14:39 rules/common.smk
-rw-rw-r-- 1 liangyong liangyong  728 11月 14 15:05 rules/counts.smk
-rw-rw-r-- 1 liangyong liangyong 1992 11月 16 17:52 rules/diffexp.smk
-rw-rw-r-- 1 liangyong liangyong  730 11月 14 15:10 rules/filtration.smk
-rw-rw-r-- 1 liangyong liangyong 1319 11月 14 15:14 rules/normalization.smk
-rw-rw-r-- 1 liangyong liangyong 3259 11月 14 15:12 rules/qc.smk
-rw-rw-r-- 1 liangyong liangyong 3556 11月 14 15:36 rules/variance.smk
```
将这些文件中的image替换成新构建的。

#### 配置分析的参数文件 -- config.yaml
###### 设置gene count table
```{shell}
counts:
  # specify count table (rows: genes/transcripts/spikes, cols: cells)
  path: test_realdata/GSE131907_Lung_Cancer_raw_UMI_matrix.txt
  # define which kind of features are described in each row (must be a term understood by Ensembl biomart, e.g. ensembl_gene_id, hgnc_symbol)
  feature_ids: ensembl_gene_id
  # BioMart host to use (e.g. www.ensembl.org, useast.ensembl.org, ...)
  biomart: www.ensembl.org

```

###### 设置cell分组文件
```{shell}
cells: test_realdata/GSE131907_Lung_Cancer_cell_annotation.txt
#注意：这里的cell文件的rownames要与 counts文件的colnames一致，包括顺序
```

###### 设置marker文件
```{shell}
celltype:
  # Table describing markers for assignment of cell types.
  # Columns: name (cell type name), parent (parent cell type name),
  #          genes (comma-separated list of gene names/ids, as listed in the
  #          count matrix)
  # Thereby, parent is usually empty. If not, it means that assignment for that
  # type happens recursively only on those cells that have been assigned to the
  # parent type.
  markers: test_realdata/markers.tsv
  # Minimum gamma score for assigned cell type (resembles a posterior) to be
  # considered as correctly assigned. Cells where the certainty of cellassign
  # does not pass this threshold will show as celltype=NA.
  min_gamma: 0.9
  # Genes to create expression plots stratified by celltype for.
  # This can be used to find the right selection of marker genes for cellassign.
  expression-plot-genes:
    - CD2
    - ENG
```

###### 设置差异分析参数
```{shell}
# Comment out to not do differential expression analysis.
diffexp:
  # Add one entry per comparison here. The key below can be an arbitrary name.
  #a-vs-b:
  LUNG_N01-vs-LUNG_N06:
    # EdgeR design formula.
    # Refer to any colData from SingleCellExperiment here.
    # In addition, you can use celltype and detection_rate
    # (number of expressed genes in cell divided by total
    # number of genes in experiment).
    #design: "~ test.condition"
    design: "~ Cell_type"
    # Which coefficients of the model should be tested equal to zero.
    # E.g., 2 to test the first coefficient after the implicit intercept
    # (i.e., celltype in the example above).
    coef: 2
    # False discovery rate to control for.
    fdr: 0.05
    # Optional: constrain to cell types (comment out to use all cell types).
    constrain-celltypes:
      celltypes:
        # - Endothelial-cell
        #- "IFNy+"
        #- "IFNy-"
        - "T-cell"
      # Optional: constrain cells to those with the given covariate occurring in all celltypes 
      # (comment in if needed).
      # This can be used to avoid confounding of an important batch variable.
      # E.g., if you want the differential expression across cell types, and sample is a 
      # batch variable to control for, you need to ensure that each sample contains all
      # considered cell types.
      # common: sample
    # Genes to plot
    genes_of_interest:
      - Mitf
      - Mycn
```

###### 设置基因与基因比较参数
```{shell}
gene-vs-gene-plots:
  all-malignant:
    # uncomment below to perform a correlation of given type (pearson, spearman, ...)
    correlation: spearman
    # uncomment below to perform a regression with given formula
    # regression: "y ~ x"
    # constrain to cells of the following types (comment out if not needed)
    constrain-celltypes:
      # - Malignant
      - "Endothelial-cell"
      - "T-cell"
      - "B-cell"
    pairs:
      x: CD2
      y:
        - CD3D
        - CD3E
```

#### 实际运行(真实数据)
数据的地址为：`/share/work/runtime/PipelineTestdata_SIAT/SingleCell_RNA`，里面包含4个文件，如下：
```{shell}
$ls -l /share/work/runtime/PipelineTestdata_SIAT/SingleCell_RNA
总用量 67489
-rw-r--r-- 1 root root     5428 11月 20 16:32 config.yaml
-rw-r--r-- 1 root root    30559 11月 20 16:33 GSE131907_Lung_Cancer_cell_annotation.txt.1000cols.part.sort
-rw-r--r-- 1 root root 34507758 11月 20 16:33 GSE131907_Lung_Cancer_raw_UMI_matrix.txt.1000cols.ensemble_id
-rw-r--r-- 1 root root      174 11月 20 16:33 markers.tsv.1
```

- config.yaml 为配置文件，用户可以调整其中的分析参数、设置差异分析组，以及想研究的基因表达情况
- GSE131907_Lung_Cancer_cell_annotation.txt.1000cols.part.sort 为细胞分组文件
- GSE131907_Lung_Cancer_raw_UMI_matrix.txt.1000cols.ensemble_id 为基因reads count文件
- markers.tsv.1 为用于区分细胞类型的基因名字

```{shell}
snakemake --cores 100 --use-singularity
```


#### 注意事项
1. 细胞分组文件中的细胞顺序，应该与gene count文件中的细胞顺序一致
2. gene count文件中的id应该是ensemble_gene_id
3. marker gene表格的内容应该是基因名字（如：CD8），而不是ensemble gene id






