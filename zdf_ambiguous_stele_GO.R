library(ggplot2)
library(stringr)
library(clusterProfiler)
library(org.At.tair.db)
stele_ASCvall_glm<-read.table("../../data/processed_data/stele_ASCvall_glm.tsv", header = T)

#寻找差异基因
rownames(stele_ASCvall_glm)<-stele_ASCvall_glm$id
diff_genelist<-rownames(stele_ASCvall_glm[stele_ASCvall_glm$q_value<0.1,])
up_genelist<-rownames(stele_ASCvall_glm[stele_ASCvall_glm$q_value<0.1&stele_ASCvall_glm$estimate>0,])
down_genelist<-rownames(stele_ASCvall_glm[stele_ASCvall_glm$q_value<0.1&stele_ASCvall_glm$estimate<0,])
length(diff_genelist)#2134
length(up_genelist)#776

#火山图验证
stele_ASCvall_glm$change=as.factor(ifelse(stele_ASCvall_glm$q_value<0.1&stele_ASCvall_glm$estimate!=0,ifelse(stele_ASCvall_glm$estimate>0,'up','down'),'not'))
this_tile <- paste0('stele_isambiguous DGE',
                    '\nCutoff for estimate is 0',
                    '\nThe number of up gene is ',nrow(stele_ASCvall_glm[stele_ASCvall_glm$change =='up',]) ,
                    '\nThe number of down gene is ',nrow(stele_ASCvall_glm[stele_ASCvall_glm$change =='down',])
)
g_volcano = ggplot(data=stele_ASCvall_glm, aes(x=estimate, y=-log10(q_value), color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+#theme_bw设置白背景灰网格线，base_size控制字体大小
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile  ) + 
  theme(plot.title = element_text(size=15,hjust = 0.5),#hjust:字体水平距离，0-1之间
        axis.text=element_text(size=14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=10))+
  scale_colour_manual(values = c('blue','black','red'))  ## corresponding to the levels(stele_ASCvall_glm$change)
print(g_volcano)#figure: 17stele_isambiguous_volcano.png

#基因ID转换symbol，entrezid
ls('package:org.At.tair.db')
id2symbol=toTable(org.At.tairSYMBOL);colnames(id2symbol)[1]<-'id'
id2entrezid=toTable(org.At.tairENTREZID);colnames(id2entrezid)<-c('id','entrezid')
stele_ASCvall_glm<-merge(stele_ASCvall_glm,id2entrezid,by='id')
stele_ASCvall_glm<-merge(stele_ASCvall_glm,id2symbol,by='id')

#GO注释
up_genelist<-stele_ASCvall_glm[stele_ASCvall_glm$q_value<0.1&stele_ASCvall_glm$estimate>0,]$entrezid
egoCC <- enrichGO(gene         = up_genelist,
                OrgDb         = org.At.tair.db, 
                ont           ='CC' ,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.2,
                qvalueCutoff  = 0.9)
egox <- setReadable(egoCC, 'org.At.tair.db', 'ENTREZID')#在go分析中，输入gene entrezid,输出的结果中，egoCC@result的geneID列为entrezid，使用setreadable函数可使geneid列的entrezid变为symbol

egoBP <- enrichGO(gene         = up_genelist,
                  OrgDb         = org.At.tair.db, 
                  ont           ='BP' ,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.2,
                  qvalueCutoff  = 0.9)
egoMF <- enrichGO(gene         = up_genelist,
                  OrgDb         = org.At.tair.db, 
                  ont           ='MF' ,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.2,
                  qvalueCutoff  = 0.9)
egoCC@result$annotation<-'cellular component'
egoBP@result$annotation<-'biological process'
egoMF@result$annotation<-'molecular function'
ego<-rbind(egoBP@result[,1:10],egoCC@result[,1:10],egoMF@result[,1:10])
write.csv(ego,'../../data/zdf_ego_stele_isambiguous.csv',col.names = T,row.names = F,quote = F)

#Go可视化
dotplot(egoBP, showCategory=10,title="UP biological process in stele_isambiguous GO enrichment")#figure: 20stele_isambiguous_upBP_go_entichment_dotplot.png

barplot(egoBP,showCategory = 10) + ggtitle("barplot for Biological process")

ggplot(data=egoBP[1:10,], aes(x=reorder(Description,-log10(qvalue)),y=-log10(qvalue))) + #y也可以是count
  geom_bar(stat="identity", width=0.8,fill='salmon1') + 
  coord_flip() +  xlab("GO term") + ylab("-log_qvalue") +
  theme_bw()#去除黑背景色

#可视化基因与通路间的网络联系:与一般功能富集柱状图和气泡图相比较的话，其优势在于可以清晰的展示出基因与功能的关系，哪些基因作用于什么功能。
foldchange<-stele_ASCvall_glm[stele_ASCvall_glm$q_value<0.1&stele_ASCvall_glm$estimate>0,]$estimate
names(foldchange)<-stele_ASCvall_glm[stele_ASCvall_glm$q_value<0.1&stele_ASCvall_glm$estimate>0,]$symbol#必须加上名字，否则没颜色
cnetplot(egox, showCategory = 5, categorySize="pvalue",foldChange=foldchange)
y <- c("Fatty acid biosynthesis","Fatty acid metabolism","PPAR signaling pathway","Insulin signaling pathway","Fatty acid degradation","Glucagon signaling pathway","Insulin resistance")
cnetplot(egox, showCategory = y, categorySize="pvalue",foldChange=foldchange)

cnetplot(egox, categorySize="pvalue", circular = TRUE, colorEdge = TRUE, foldChange=foldchange)

p1 <- cnetplot(egox, node_label="category") 
p2 <- cnetplot(egox, node_label="gene") 
p3 <- cnetplot(egox, node_label="all") 
p4 <- cnetplot(egox, node_label="none") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

heatplot(egox, foldChange=stele_ASCvall_glm[stele_ASCvall_glm$q_value<0.1&stele_ASCvall_glm$estimate>0,]$estimate)

#KEGG注释
up.kk <- enrichKEGG(gene         = up_genelist,
                    organism     = 'ath',#https://www.genome.jp/kegg/catalog/org_list.html中查看
                    pvalueCutoff = 0.99,
                    qvalueCutoff=0.99
)
up.kk@result#ID: pathway通路ID, GeneRatio: 分子是富集到该条目的gene数，分母为所有输入的做富集分析的gene数即差异表达分析得到的gene, BgRatio:分母是拟南芥的所有表达基因中有GO注释的gene数，分子是拟南芥所有表达基因中注释到该条目的gene数

#KEGG可视化
dotplot(up.kk, showCategory=10,title="KEGG_Up_biological_stele_is_ambiguous")#默认x为generatio, color为p.adjust,orderBy='x',figure: 18stele_isambiguous_upkegg_dotplot

barplot(up.kk, showCategory=10)

ggplot(up.kk@result[1:10,], aes(x=reorder(Description,-log10(qvalue)), y=-log10(qvalue), fill=-log10(qvalue))) + #reorder用于ggplot图中，aes (x=reorder(age,money),y=money)，即按 money 对 age 排序
  geom_bar(stat="identity") + #条形图高度表示y值，当stat为count时，高度等于每组中的数据的个数
  #scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
  scale_x_discrete(name ="Pathway names") +
  scale_y_continuous(name ="-log10P-value") +
  coord_flip() +
  ggtitle("KEGG_Up_biological_stele_is_ambiguous")#figure: 19stele_isambiguous_upkegg_geombar





