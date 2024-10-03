# 1.monocle安装
```
BiocManager::install("monocle")
library(monocle)
```
# 2.加载数据
这里使用的数据是示例数据，**来自SeuratData**，数据已经经过了细胞分群
```
library(pbmc3k.SeuratData)
data("pbmc3k") 
sce=pbmc3k.final 
```
```
#####提取要做monocle的细胞群子集
sce1=subset(sce, seurat_annotations %in% c('CD14+ Mono','FCGR3A+ Mono')) 
sce1 <- NormalizeData(sce1) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)  
  ```
# 3.转化为 monocle CDS数据格式
 ## 3.1直接转换
 一步转换，很轻松
 ```
 ### 2.1 Seurat转为cds数据
HSMM <- as.CellDataSet(sce1)
HSMM
```
![](2024-10-02-15-37-56.png)
 ## 3.2自行构建cds
 ```
 ## 3.2 如果是其他格式的数据，需要自行构建cds格式
# 3.2.1 表型数据
sample_ann <- sce1@meta.data  
# 3.2.2 基因信息
gene_ann <- data.frame(
  gene_short_name = rownames(sce1@assays$RNA), 
  row.names = rownames(sce1@assays$RNA) 
)
# 3.2.3 表达矩阵
pd <- new("AnnotatedDataFrame",
          data=sample_ann)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
ct=as.data.frame(sce1@assays$RNA@counts)
# 3.2.4 构建cds对象
HSMM2 <- newCellDataSet(
  as.matrix(ct), 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
HSMM2
```
# 4.运行monocle
## 4.1计算size factor
```
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
```
<font color="red">
<font size=2>
Warning Messages:

1. <code>group_by_()</code> was deprecated in <code>dplyr</code> 0.7.0.  
   <i>Please use <code>group_by()</code> instead.</i>  
   <i>See vignette('programming') for more help</i>
   
   The deprecated feature was likely used in the <code>monocle</code> package.  
   Please report the issue to the authors.  
   
   This warning is displayed once every 8 hours. Call <code>lifecycle::last_lifecycle_warnings()</code> to see where this warning was generated.

2. <code>select_()</code> was deprecated in <code>dplyr</code> 0.7.0.  
   <i>Please use <code>select()</code> instead.</i>  
   <i>See vignette('programming') for more help</i>
   
   The deprecated feature was likely used in the <code>monocle</code> package.  
   Please report the issue to the authors.  
   
   This warning is displayed once every 8 hours. Call <code>lifecycle::last_lifecycle_warnings()</code> to see where this warning was generated.
</font>
</font>
报错了，问题出在monocle包的代码
## 修改monocle代码
monocle代码有点问题，现在我直接修改monocle2.30.1的代码
```
trace(monocle:::estimateDispersionsForCellDataSet, edit = T)
```
修改对应的代码即可，这种方法修改代码只能用一次
## 4.2数据质控
```
HSMM <- detectGenes(HSMM, min_expr = 1)
print(head(fData(HSMM)))
```                 
| gene_short_name      | num_cells_expressed |
|----------------------|---------------------|
| AL627309.1           | 0                   |
| AP006222.2           | 0                   |
| RP11-206L10.2        | 0                   |
| RP11-206L10.9        | 0                   |
| LINC00115            | 0                   |
| NOC2L                | 7                   |
```
#########
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

head(pData(HSMM))
```
| orig.ident    | nCount_RNA | nFeature_RNA | seurat_annotations | percent.mt |
|---------------|------------|--------------|--------------------|------------|
| AAACCGTGCTTCCG | 2639       | 960          | CD14+ Mono          | 1.743085   |
| AAACGCTGTTTCTG | 1103       | 550          | FCGR3A+ Mono         | 2.901179   |
| AAAGAGACGCGAGA | 3033       | 1058         | CD14+ Mono          | 1.417738   |
| AAAGCAGATATCGG | 4584       | 1422         | CD14+ Mono          | 1.396161   |
| AAAGTTTGTAGCGT | 2683       | 877          | CD14+ Mono          | 2.497205   |
| AAATCAACCCTATT | 5676       | 1541         | FCGR3A+ Mono         | 2.431290   |
```
length(expressed_genes)
```
1803
## 4.3选择输入的基因用于排序细胞
第一步是决定用哪些基因来进行聚类（特征选择）。我们可以使用所有的基因，但是我们将包括很多没有表达到足够高水平来提供有意义的信号的基因，它们只会给系统增加噪音。因此，并不是所有的基因都有作用，所以先进行挑选合适的基因用来进行聚类。

这里选择基因的方式有很多，说明文档中建议以下4种选择基因的方式
（1）选择离散程度高的基因（例如Seurat的高变基因）；
（2）选择clusters差异表达基因（算法预测）；
（3）选择发育差异表达基因（需要结合背景知识）；
（4）自定义发育marker基因（需要结合背景知识）。
### 4.3.1 根据高变基因(最常用)
```
HSMM <- setOrderingFilter(HSMM, VariableFeatures(sce1))
plot_ordering_genes(HSMM)
```
![](2024-10-03-11-38-30.png)
### 4.3.2根据cluster
```
## 4.3.2 选择clusters差异表达基因
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_annotations",
                                      cores = 10)
```
```
# 挑选差异最显著的基因进行可视化
ordering_genes <- subset(diff_test_res, qval < 0.05)
# 按 p 值对挑选出的基因进行排序
ordering_genes = ordering_genes[order(ordering_genes$pval),]
# 显示排序后的基因的前几行，包括基因名、p 值和 q 值
head(ordering_genes[,c("gene_short_name", "pval", "qval")] )
# 提取前几个基因的名称作为字符向量
cg = as.character(head(ordering_genes$gene_short_name))
cg
```
"S100A9" "S100A8" "FCGR3A" "LYZ"    "IFITM2" "LGALS2"
```
# 从差异测试结果中挑选校正后的p值（qval）小于0.05的基因
# 并将这些基因的行名（通常是基因名）保存在ordering_genes变量中
ordering_genes <- row.names(subset(diff_test_res, qval < 0.05))

# 计算挑选出的基因数量
gene_count <- length(ordering_genes)
# 打印基因数量
print(gene_count)

# 使用setOrderingFilter函数将挑选出的基因设置为过滤条件
# 这里只取了前3000个基因进行后续分析
HSMM <- setOrderingFilter(HSMM, ordering_genes[1:3000])

# 绘制挑选出的基因的排序图
plot_ordering_genes(HSMM)
```
![](2024-10-03-11-30-43.png)
### 4.3根据表达量进行注释
```
# 使用 dispersionTable 函数计算 HSMM 对象中基因的离散度表
disp_table <- dispersionTable(HSMM)

# 从离散度表中筛选出平均表达量大于等于 0.1 的基因
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)

# 使用 setOrderingFilter 函数设置 HSMM 对象的排序过滤条件
# 这里使用筛选出的基因的 ID 列表作为过滤条件
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
```
![](2024-10-03-11-41-55.png)
## 4.4 降维 & 排序
```
# 使用 reduceDimension 函数对 HSMM 对象进行降维分析
# max_components 设置为 2，表示结果中只保留前两个主成分
# num_dim 设置为 20，表示要计算的主成分数量
# method 设置为 'DDRTree'，使用 DDRTree 方法进行降维
HSMM <- reduceDimension(HSMM,
                        max_components = 2,
                        num_dim = 20,
                        # 如果数据中存在批次效应，可以取消注释下面的行来指定批次
                        # residualModelFormulaStr = "~SampleID", # 如果存在批次则指定批次
                        method = 'DDRTree') 

# 使用 orderCells 函数对细胞进行排序，以便于后续的可视化
HSMM <- orderCells(HSMM)

pData(HSMM) %>% head()
```
| orig.ident        | nCount_RNA | nFeature_RNA | seurat_annotations | percent.mt   |
|-------------------|------------|--------------|--------------------|--------------|
| AAACCGTGCTTCCG     | 2639       | 960          | CD14+ Mono          | 1.743085     |
| AAACGCTGTTTCTG     | 1103       | 550          | FCGR3A+ Mono         | 2.901179     |
| AAAGAGACGCGAGA     | 3033       | 1058         | CD14+ Mono          | 1.417738     |
| AAAGCAGATATCGG     | 4584       | 1422         | CD14+ Mono          | 1.396161     |
| AAAGTTTGTAGCGT     | 2683       | 877          | CD14+ Mono          | 2.497205     |
| AAATCAACCCTATT     | 5676       | 1541         | FCGR3A+ Mono         | 2.431290     |
```
# 定义一个函数，用于根据细胞发展状态（cds）和起始点（starting_point）来确定细胞的排序
GM_state <- function(cds, starting_point, cluster){
    # 如果细胞状态（cds$State）不唯一，即存在多个状态，那么执行if语句
    if (length(unique(cds$State)) > 1){
        # 创建一个表格，统计每个状态在特定聚类（cluster）中的出现次数
        T0_counts <- table(cds$State, cds@phenoData@data[,cluster])[,starting_point]
        # 找出出现次数最多的状态的名称
        return(as.numeric(names(T0_counts)[which
                                           (T0_counts == max(T0_counts))]))
    } else {
        # 如果所有细胞状态相同，即只有一个状态，返回状态1
        return (1)
    }
}
```
```
root_start = GM_state(cds = HSMM,starting_point = "CD14+ Mono",cluster = "celltype")
HSMM <- monocle::orderCells(HSMM, root_state = root_start)
```
# 5.可视化
```
#配色
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

#主题
if(T){
  text.size = 12
  text.angle = 45
  text.hjust = 1
  legend.position = "right"
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"), 
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                   panel.grid=element_blank(), # 去网格线
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size)
  )
}
```
### 5.1 基于各种“类型”可视化
```
# 分别为根据 seurat cluster ，State ，Pseudotime 和 singleR注释后的cell type 着色。
a1 <- plot_cell_trajectory(HSMM, color_by = "seurat_annotations") + scale_color_manual(values = colour)
a2 <- plot_cell_trajectory(HSMM, color_by = "State") + scale_color_manual(values = colour)
a3 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime") + ggsci::scale_color_gsea()
options(repr.plot.width = 14, repr.plot.height = 4.5)
wrap_plots(a1, a2, a3)
```
![](2024-10-03-16-36-41.png)
## 5.2 分面展示
```
options(repr.plot.width = 4, repr.plot.height = 6)
plot_cell_trajectory(HSMM, color_by = "Pseudotime") +
    facet_wrap(~celltype, nrow = 2) + ggsci::scale_color_gsea()
```
![](2024-10-03-16-47-58.png)
### 5.3 添加“树形图”
```
# plot_complex_cell_trajectory函数添加“树形图”
p1 <- plot_cell_trajectory(HSMM, x = 1, y = 2, color_by = "celltype") + 
  theme(legend.position='none',panel.border = element_blank()) + 
  scale_color_manual(values = colour) 
p2 <- plot_complex_cell_trajectory(HSMM, x = 1, y = 2,
                                   color_by = "celltype")+
  scale_color_manual(values = colour) +
  theme(legend.title = element_blank()) 

options(repr.plot.width = 4.5, repr.plot.height = 9)
wrap_plots(p1, p2, ncol = 1, heights = c(2,1.5))
```
![](2024-10-03-16-49-14.png)
### 5.4 箱式图
```
input.data = data.frame(celltype = HSMM$celltype,
                        Pseudotime = HSMM$Pseudotime)

options(repr.plot.width = 4.5, repr.plot.height = 4)
ggboxplot(data = input.data,
          fill = "celltype", 
          x = "celltype",  
          y = "Pseudotime")+mytheme
```
![](2024-10-03-16-53-06.png)
### 5.5 gene表达
```
input.gene <- c("CD14","S100A8", "FCGR3A")
cds_subset <- HSMM[input.gene,]

options(repr.plot.width = 12, repr.plot.height = 3)
plot_genes_in_pseudotime(cds_subset, color_by = "celltype", ncol = 3)
```
![](2024-10-03-16-59-28.png)
```
# 提取名为 FCGR3A 的基因表达数据，并将其添加到 HSMM 对象的元数据中
pData(HSMM)$FCGR3A <- as.numeric(GetAssayData(object = sce1, assay = "RNA", slot = "data")["FCGR3A",])

# 提取名为 CD14 的基因表达数据，并将其添加到 HSMM 对象的元数据中
pData(HSMM)$CD14 <- as.numeric(GetAssayData(object = sce1, assay = "RNA", slot = "data")["CD14",])

# 使用 monocle 包的 plot_cell_trajectory 函数绘制细胞发展轨迹
# 并使用 ggsci 包的 scale_color_gsea 函数为轨迹图添加颜色
p1 <- plot_cell_trajectory(HSMM, color_by = "FCGR3A") + ggsci::scale_color_gsea()

# 同样绘制细胞发展轨迹，但这次根据 CD14 基因表达数据进行颜色编码
p2 <- plot_cell_trajectory(HSMM, color_by = "CD14") + ggsci::scale_color_gsea()

options(repr.plot.width = 12, repr.plot.height = 4.5)
wrap_plots(p1, p2, ncol = 2)
```
![](2024-10-03-17-23-28.png)
### 5.6 gene/通路的相关性折线图
```
input.data = data.frame(celltype = HSMM$celltype,
                        Pseudotime = HSMM$Pseudotime)
## 基因
input.data$FCGR3A = as.numeric(GetAssayData(object = sce, assay = "RNA",slot = "data")["FCGR3A",])
options(repr.plot.width = 6, repr.plot.height = 4)
## Vis
ggplot(input.data, aes(x = Pseudotime, y = FCGR3A)) +
  labs(x="Pseudotime",y = "Expression level (log2)")+
  ggpubr::stat_cor(label.sep = "\n",
                   label.y = 3,
                   label.y.npc = "top",
                   size = 4,
                   method = "sp")  +
  geom_smooth(method='loess',size=0.8, color = "black") + theme_bw() +
  ggfun::facet_set(label = "FCGR3A")+
  mytheme + theme(legend.position = "none")
```
![](2024-10-03-17-39-33.png)
```
## 通路
sce = AddModuleScore(object = sce, features = 
                       list(test = c("LGALS2", "FCGR3A", "S100A9", "CCL3", "CD14", "LYZ", "FOLR3", "ECHDC1", "S100A8")),
                     name = "test_pathway")
input.data$test_pathway = sce$test_pathway1
ggplot(input.data, aes(x = Pseudotime, y = test_pathway)) +
  labs(x="Pseudotime",y = "Pathway activity")+
  ggpubr::stat_cor(label.sep = "\n",
                   label.y = 0,
                   label.y.npc = "top",
                   size = 4,
                   method = "sp")  +
  geom_smooth(method='loess',size=0.8, color = "black") + theme_bw() +
  ggfun::facet_set(label = "Test pathway")+
  mytheme + theme(legend.position = "none")
  ```
  ![](2024-10-03-17-40-55.png)
# 6. 感兴趣基因的热图
```
input.gene <- c("LGALS2", "FCGR3A", "S100A9", "CCL3", "CD14", "LYZ", "FOLR3", "ECHDC1", "S100A8")
options(repr.plot.width = 6, repr.plot.height = 4)
plot_pseudotime_heatmap(HSMM[input.gene,], 
                             num_cluster = 4, 
                             show_rownames = T, 
                             return_heatmap = T)
```
![](2024-10-03-17-47-09.png)
# 6.BEAM分析
确定了分化起点后，Monocle可以模拟出每个细胞所处的分化时间，并寻找随着分化时间逐渐升高或降低的基因，即Beam分析。
```
# 使用 BEAM 函数对 HSMM 对象进行分支轨迹分析
BEAM_res <- BEAM(HSMM, branch_point = 1, 
                 cores = 4, 
                 progenitor_method = "duplicate")
# 按照 qval 对 BEAM 分析结果进行排序
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
# 选择包含基因名称、p 值和 q 值的列
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
# 查看 BEAM 分析结果的维度，即结果的行数和列数
dim(BEAM_res)
# 将 BEAM 分析结果保存为 CSV 文件
write.csv(BEAM_res, file = "./Outdata/Step1.BEAM_res.csv")
# 从 BEAM 分析结果中提取 q 值小于 0.05 的基因
input.gene = row.names(subset(BEAM_res, qval < 0.05))
# 打印提取的基因列表
input.gene
# 计算并打印提取的基因数量
length(input.gene)
# 设置图形输出选项，定义图形的宽度和高度
options(repr.plot.width = 4.5, repr.plot.height = 6)
# 使用 plot_genes_branched_heatmap 函数绘制分支热图
# 函数参数指定了分支点、聚类数、使用的核心数以及是否显示基因名称和行名
plot_genes_branched_heatmap(HSMM[input.gene,],
                            branch_point = 1,
                            num_clusters = 6,
                            cores = 2,
                            use_gene_short_name = T,
                            show_rownames = T)
```
![alt text]({AFDA0733-9DE8-410F-A588-DFA859355EF9}.png)
上图将Top差异基因聚成6类，展示了pre-branch向cell fate1和cell fate2状态分化相关的命运决定基因。                   