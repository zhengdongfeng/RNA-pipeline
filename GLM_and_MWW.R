####contents:GLM和MWW法寻找细胞型差异表达基因和伪时间差异表达基因

###GLM：
##celltype DEG: fit_model(model_formula_str = "~celltype")---coefficient_table---获得表格id/term/estimate/q_value列
##pseudotime DEG: order_cells---coldata加入pseudotime列---fit_model(model_formula_str = "~pseudotime")---coefficient_table---获得表格id/term/estimate/q_value列
#pericycle: XPP vs LRP, XPP vs MAT, XPP vs MAT pseudotime, XPP vs LRP pseudotime, xpp_lrp_mat_pseudotime(plot_cells画出4组基因表达umap：两个轨迹均差异表达，两个轨迹均不差异表达，只在LRP轨迹差异表达，只在MAT轨迹差异表达)
#stele: Ambiguous Stele Cells vs. all other stele cells(coldata列加入is_ambiguous列，fit_model(model_formula_str = "~is_ambiguous")) 
#ce: endodermis vs cortex, lre vs cortex, endodermis_main_branch_pseudotime, endodermis_lre_branch_pseudotime, endodermis_pseudotime(plot_cells画出4组基因表达umap：两个分支轨迹均差异表达，两个分支轨迹均不差异表达，只在main_branch轨迹差异表达，只在lre_branch轨迹差异表达)

###MWW:
##先标准化表达矩阵：norm_mat <- t(t(counts(cds))/size_factors(cds))
##wilcox_p_values(wilcox.test)--wilcox_p_values_adj(p.adjust)---分组stats(means, var, n)---table("id", "xpp_mean", "xpp_stand_err", "lrp_mean", "lrp_stand_err", "xpp_lrp_log2FoldChange", "wilcox_p_values_adj")
#pericycle: XPP vs LRP, XPP vs MP
#ce: cortex vs endodermis, cortex vs lre



#!/usr/bin/env Rscript
### Part 0: Load Packages
library(monocle3)
library(Matrix)
library(gplots)
library(ggplot2)
library(pheatmap)
library(viridis)
library(dplyr)
options(stringsAsFactors=F)


### GLM 
# load pericycle data
cds = readRDS("../../data/GSE158761_pericycle_cells_labeled.rds")#18868*2243

## XPP vs. Lateral Root Primordia
lm = fit_models(cds[, colData(cds)$celltype == "Xylem Pole Pericycle" | colData(cds)$celltype == "Lateral Root Primordia" ],  model_formula_str = "~celltype", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)## view(lm),其中model_formula_str就是要比较的分组对象，如果想获得不同的cluster或者partition的差异基因，就用model_formula_str = "~cluster"或者model_formula_str = "~partition"；另外还支持添加多个变量，比如考虑到批次效应 model_formula_str = "~embryo.time + batch"，embryo.time为colData中的某列名;大多数研究使用负二项分布，这种分布往往是符合测序reads或UMI count值的，这个方法也被应用在RNA-seq分析（如DESeq2）中,多种分布，默认是：quasipoisson ，它和负二项分布相似，只不过比负二项分布准确性低一点、速度会更快，对于细胞数量多更适用。
saveRDS(lm, "../../data/stele_XPPvLRP_lm.rds") # save linear model
ct <- coefficient_table(lm)# view(ct),coefficient_table()默认使用 Benjamini and Hochberg（BH）方法进行了p值的校正，得到了q值
#View(ct):term一列intercept作为参照组，LRP作为实验组，estimate值大于0且q值足够小说明基因在LRP中上调
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv  
write.table(ot, "../../data/stele_XPPvLRP_glm.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

## XPP vs. Mature Pericycle 
lm = fit_models(cds[, colData(cds)$celltype == "Xylem Pole Pericycle" | colData(cds)$celltype == "Mature Pericycle" ],  model_formula_str = "~celltype", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "../../data/stele_XPPvMAT_lm.rds") # save linear model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "../../data/stele_XPPvMATUREPERICYCLE_glm.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

## XPP and Mature Pericycle pseudotime
cds_xpp_mat = cds[, colData(cds)$celltype == "Xylem Pole Pericycle" | colData(cds)$celltype == "Mature Pericycle" ]#plot_cells(cds_xpp_mat,color_cells_by="celltype"),figure: 5xpp_mp_umap.png
cds_xpp_mat <- order_cells(cds_xpp_mat) # this is interactive!!!! Select XPP cells as root
colData(cds_xpp_mat)$pseudotime <- cds_xpp_mat@principal_graph_aux$UMAP$pseudotime # save pseudotime info; #plot_cells(cds_xpp_mat,color_cells_by="pseudotime",label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE),figure: 6xpp_mp_pseudotime.png
lm <- fit_models(cds_xpp_mat,  model_formula_str = "~pseudotime", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "../../data/pseudotime_lm_xpp_mat_pseudotime.rds") # save model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "../../data/stele_XPPtoMAT_pseudotime.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

## XPP and Lateral Root pseudotime
cds_xpp_lrp = cds[, colData(cds)$celltype == "Xylem Pole Pericycle" | colData(cds)$celltype == "Lateral Root Primordia" ]
cds_xpp_lrp <- order_cells(cds_xpp_lrp) # this is interactive!!!! Select XPP cells as root
colData(cds_xpp_lrp)$pseudotime <- cds_xpp_lrp@principal_graph_aux$UMAP$pseudotime # save pseudotime info
lm <- fit_models(cds_xpp_lrp,  model_formula_str = "~pseudotime", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "../../data/pseudotime_lm_xpp_lrp_pseudotime.rds") # save model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "../../data/stele_XPPtoLRP_pseudotime.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

# plotting
xpp_tab = read.table("./Pseudotime_lists/xpp_lrp_mat_pseudotime.txt", sep="\t", header=T) # slim version of supplemental table 1
plot_cells(cds, genes = unique(xpp_tab), show_trajectory_graph=T, label_branch_points=F, label_leaves=F, cell_size=1.5, label_cell_groups=F) + xlim(-4.5, 2) + ylim(-6, 0.5) + coord_equal() + scale_color_viridis(na.value=rgb(1, 1, 1, 0.01))#xlim(-4.5, 2)+ylim(-6, 0.5)截取对应的横纵坐标范围, +theme( title =element_text (size =14))修改文字大小。figure: 7xpp_lrp_mat_pseudotime.png


## Ambiguous Stele Cells vs. all other stele cells
# load stele data
cds = readRDS("../../data/GSE158761_stele_cells_labeled.rds")
# create feature column to run linear model
colData(cds)$is_ambiguous = colData(cds)$celltype == "Ambiguous Stele Cells"
lm = fit_models(cds,  model_formula_str = "~is_ambiguous", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "../../data/stele_ASCvall_lm.rds") # save model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "../../data/stele_ASCvall_glm.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

## Cortex/Enodermis/Lateral Root Endodermis
# load cortex / endodermis data
cds = readRDS("../../data/GSE158761_ce_cells_labeled.rds")#plot_cells(cds,color_cells_by = 'Round2_clusters',label_cell_groups = T,label_branch_points = F,label_leaves = F,label_roots = F), figure: 9ce_cluster.png; plot_cells(cds,color_cells_by = 'celltype',label_cell_groups = T,label_branch_points = F,label_leaves = F,label_roots = F), figure: 8ce_celltype_umap.png;
# Endodermis vs. Cotex
lm = fit_models(cds[, colData(cds)$celltype == "Cortex" | colData(cds)$celltype == "Endodermis" ],  model_formula_str = "~celltype", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "../../data/stele_ENDODERMISvCORTEX_lm.rds") # save linear model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "../../data/stele_ENDODERMISvCORTEX_glm.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

# Lateral Root Endodermis vs. Cortex
lm = fit_models(cds[, colData(cds)$celltype == "Cortex" | colData(cds)$celltype == "Lateral Root Endodermis" ],  model_formula_str = "~celltype", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "../../data/stele_LREvCORTEX_lm.rds") # save linear model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "../../data/stele_LREvCORTEX_glm.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

# Endodermis Pseudotime
cds_branch1 <- cds[, !( colData(cds)$Round2_clusters %in% c(1, 2, 5, 7, 8, 10, 13, 14) ) ] # plot_cells(cds_branch1,color_cells_by = 'Round2_clusters',label_cell_groups = T,label_branch_points = F,label_leaves = F,label_roots = F)+xlim(-5, 6)+ylim(-2.5, 3.5), figure: 10endo_mainbranch_cluster.png; plot_cells(cds_branch1,color_cells_by = 'celltype',label_cell_groups = T,label_branch_points = F,label_leaves = F,label_roots = F)+xlim(-5, 6)+ylim(-2.5, 3.5), figure: 11endo_mainbranch_celltype_umap.png
cds_branch2 <- cds[, !( colData(cds)$Round2_clusters %in% c(2, 3, 10, 15) ) ]# plot_cells(cds_branch2,color_cells_by = 'Round2_clusters',label_cell_groups = T,label_branch_points = F,label_leaves = F,label_roots = F)+xlim(-5, 6)+ylim(-2.5, 3.5), figure: 12endo_lre_branch_cluster.png; plot_cells(cds_branch2,color_cells_by = 'celltype',label_cell_groups = T,label_branch_points = F,label_leaves = F,label_roots = F)+xlim(-5, 6)+ylim(-2.5, 3.5), figure: 13endo_lreb_branch_celltype_umap.png

cds_branch1 <- order_cells(cds_branch1) # this is interactive!!!! Select Endodermis cells as root
colData(cds_branch1)$pseudotime <- cds_branch1@principal_graph_aux$UMAP$pseudotime # saveRDS(cei_r2_branch1, 'cei_branch1.rds')#plot_cells(cds_branch1,color_cells_by="pseudotime",label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+xlim(-5, 6)+ylim(-2.5, 3.5),figure: 14endodermis_mainbranch_pseudotime.png
cds_branch2 <- order_cells(cds_branch2) # this is interactive!!!! Select Endodermis cells as root
colData(cds_branch2)$pseudotime <- cds_branch2@principal_graph_aux$UMAP$pseudotime # saveRDS(cei_r2_branch2, 'cei_branch2.rds')#plot_cells(cds_branch2,color_cells_by="pseudotime",label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+xlim(-5, 6)+ylim(-2.5, 3.5),figure: 15endodermis_lre_branch_pseudotime.png

lm <- fit_models(cds_branch1,  model_formula_str = "~pseudotime", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "../../data/xlm_endodermis_mainbranch_pseudotime.rds") # save model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "../../data/endodermis_mainbranch_pseudotime.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

lm <- fit_models(cds_branch2,  model_formula_str = "~pseudotime", expression_family="negbinomial", cores=1, clean_model=T, verbose=T)
saveRDS(lm, "../../data/lm_endodermis_lrebranch_pseudotime.rds") # save model
ct <- coefficient_table(lm)
ot <- ct %>% filter(term != '(Intercept)') %>% select(id, term, estimate, q_value) # extract columns to save as tsv
write.table(ot, "../../data/endodermis_lrebranch_pseudotime.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)

# plotting
endo_tab = read.table("./Pseudotime_lists/endodermis_pseudotime.txt", sep="\t", header=T) # slim version of supplemental table 3
plot_cells(cds, genes=endo_tab, cell_size=1.5, label_branch_points=F, label_leaves=F, label_cell_groups=F)+xlim(-5, 6)+ylim(-2.5, 3.5) + coord_equal()  + scale_color_viridis(na.value=rgb(1, 1, 1, 0.01))#figure: 16endodermis_pseudotime.png


### MWW
# load pericycle data
cds = readRDS("../../data/GSE158761_pericycle_cells_labeled.rds")

# size factor normalized count matrix
norm_mat <- t(t(counts(cds))/size_factors(cds))

## XPP vs. LRP MWW
wilcox_p_values <- apply(norm_mat, 1, function(x) return(wilcox.test(x[colData(cds)$celltype == "Xylem Pole Pericycle"], x[colData(cds)$celltype == "Lateral Root Primordia"])[['p.value']]))#两个独立样本的 t-test 是检验两个样本的均值是否相等。而相对的，两个独立样本的 wilcoxon test 则是检验两个样本的中位数。返回每个基因的p值
wilcox_p_values_adj <- p.adjust(wilcox_p_values, method='BH')

# calculate LRP stats
lrp_mean <- rowMeans(norm_mat[, colData(cds)$celltype == "Lateral Root Primordia"])
lrp_var <- apply(norm_mat[, colData(cds)$celltype == "Lateral Root Primordia"], 1, var)
lrp_n <- sum(colData(cds)$celltype == "Lateral Root Primordia")

# calculate XPP stats
xpp_mean <- rowMeans(norm_mat[, colData(cds)$celltype == "Xylem Pole Pericycle"])
xpp_var <- apply(norm_mat[, colData(cds)$celltype == "Xylem Pole Pericycle"], 1, var)
xpp_n <- sum(colData(cds)$celltype == "Xylem Pole Pericycle")

# create output table
xpp_lrp_log2FoldChange <- log2(xpp_mean/lrp_mean)
out <- cbind(rownames(cds), xpp_mean, sqrt(xpp_var/xpp_n), lrp_mean, sqrt(lrp_var/lrp_n), xpp_lrp_log2FoldChange, wilcox_p_values_adj)
colnames(out) <- c("id", "xpp_mean", "xpp_stand_err", "lrp_mean", "lrp_stand_err", "xpp_lrp_log2FoldChange", "wilcox_p_values_adj")
write.table(out, "../../data/XPPvLRP_normalized_counts_wilcoxon_adj.tsv", append=F, quote=F, sep="\t", col.names=T, row.names=F)

## XPP vs. Mature Pericycle MWW
wilcox_p_values <- apply(norm_mat, 1, function(x) return(wilcox.test(x[colData(cds)$celltype == "Xylem Pole Pericycle"], x[colData(cds)$celltype == "Mature Pericycle"])[['p.value']]))
wilcox_p_values_adj <- p.adjust(wilcox_p_values, method='BH')

# calculate Mature Pericycle stats
mat_mean <- rowMeans(norm_mat[, colData(cds)$celltype == "Mature Pericycle"])
mat_var <- apply(norm_mat[, colData(cds)$celltype == "Mature Pericycle"], 1, var)
mat_n <- sum(colData(cds)$celltype == "Mature Pericycle")

# create output table
xpp_mat_log2FoldChange <- log2(xpp_mean/mat_mean)
out <- cbind(rownames(cds), xpp_mean, sqrt(xpp_var/xpp_n), mat_mean, sqrt(mat_var/mat_n), xpp_mat_log2FoldChange, wilcox_p_values_adj)
colnames(out) <- c("id", "xpp_mean", "xpp_stand_err", "mat_mean", "mat_stand_err", "xpp_mat_log2FoldChange", "wilcox_p_values_adj")
write.table(out, "../../data/XPPvMATUREPERICYCLE_normalized_counts_wilcoxon_adj.tsv", append=F, quote=F, sep="\t", col.names=T, row.names=F)

# load cortex / endodermis data
cds = readRDS("../../data/ce_cells_labeled.rds")

# size factor normalized count matrix
norm_mat <- t(t(counts(cds))/size_factors(cds))

## cortex vs. endodermis
endodermis_wilcox_p_values <- apply(norm_mat, 1, function(x) return(wilcox.test(x[colData(cds)$celltype == "Cortex"], x[colData(cds)$celltype == "Endodermis"])[['p.value']]))
endodermis_wilcox_p_values_adj <- p.adjust(endodermis_wilcox_p_values, method='BH')

## cortex vs. lre
lre_wilcox_p_values <- apply(norm_mat, 1, function(x) return(wilcox.test(x[colData(cds)$celltype == "Cortex"], x[colData(cds)$celltype == "Lateral Root Endodermis"])[['p.value']]))
lre_wilcox_p_values_adj <- p.adjust(lre_wilcox_p_values, method='BH')

# calculate cortex stats
cortex_mean <- rowMeans(norm_mat[, colData(cds)$celltype == "Cortex"])
cortex_var <- apply(norm_mat[, colData(cds)$celltype == "Cortex"], 1, var)
cortex_n <- sum(colData(cds)$celltype == "Cortex")

# calculate endodermis stats
endodermis_mean <- rowMeans(norm_mat[, colData(cds)$celltype == "Endodermis"])
endodermis_var <- apply(norm_mat[, colData(cds)$celltype == "Endodermis"], 1, var)
endodermis_n <- sum(colData(cds)$celltype == "Endodermis")
write.table(cbind(rownames(cds), endodermis_mean, cortex_mean, endodermis_wilcox_p_values_adj), "ENDODERMISvCORTEX_normalized_counts_wilcoxon_adj.tsv", append=F, quote=F, sep="\t", col.names=T, row.names=F)

# calculate lre stats
lre_mean <- rowMeans(norm_mat[, colData(cds)$celltype == "Lateral Root Endodermis"])
lre_var <- apply(norm_mat[, colData(cds)$celltype == "Lateral Root Endodermis"], 1, var)
lre_n <- sum(colData(cds)$celltype == "Lateral Root Endodermis")
write.table(cbind(rownames(cds), lre_mean, cortex_mean, lre_wilcox_p_values_adj), "LREvCORTEX_normalized_counts_wilcoxon_adj.tsv", append=F, quote=F, sep="\t", col.names=T, row.names=F)

#zdfcode: 韦恩图绘制GLM，MWW，Vision三种方法得到的DEG的交集情况
library(VennDiagram)
vennlist<-list(glm=glm_genenames,mww=mww_genenames,vision=vision_genenames)
venn.diagram(vennlist,fill=c("red","green","blue"),filename = '../../data/venn.png')
##取交集元素
inter <- get.venn.partitions(vennlist)






