library(ggplot2)
library(Seurat)

#############################GSE133486数据 iWAT##############################
###################nuclei的RT数据
##读入数据
setwd('/Users/guo/Desktop/GEO_data/GSE133486/ko_efficiency_experiment/')
nuc <- CreateSeuratObject(counts = Read10X(data.dir = '/Users/guo/Desktop/GEO_data/GSE133486/iWAT_nuclei'), min.cells = 5, min.features = 200, project = 'GSE133486')
nuc[["percent.mt"]] <- PercentageFeatureSet(nuc, pattern = "^mt-")     #17987*52773
nuc <- subset(nuc, subset = nFeature_RNA > 200 & percent.mt < 10)     #按照文章的条件，筛选之后是17987*52765
FeatureScatter(nuc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 1.5)   #查看线性，0.97

##筛选contral_RT组的数据，在geo上传的metadata中给出Ad6N是negative RT组
nuc_adn_rt <- subset(x = nuc,idents='Ad6N')   #17987*4573
##统计Adipoq基因在adn_RT的表达分布，确定count值大于1为Adipoq有表达
Adipoq <- nuc_adn_rt@assays$RNA@counts['Adipoq',]
ggplot2::qplot(Adipoq)

###################计算所有基因在adipoq-cre条件下的theoretical ko_efficiency
result <- c()
data <- nuc_adn_rt@assays$RNA@counts
nuc_adn_rt_adi <- subset(x=nuc_adn_rt, Adipoq>0, slot = 'counts')     #Adipoq在所有细胞中有表达的数据，386个cells
data_adi <- nuc_adn_rt_adi@assays$RNA@counts
genes <- row.names(nuc_adn_rt_adi@assays$RNA@counts)
for (i in genes[1:5]) {
  n_all <- sum(data[i,]>0)       #统计基因在所有的adn中表达值大于0的个数，既在多少个细胞中有表达
  n_adi <- sum(data_adi[i,]>0)   #统计基因在Adipoq基因有表达的细胞中也有表达，有多少个细胞是这个基因也有表达的
  n <- sum(data_adi[i,])         #统计基因，在与Adipoq都有表达的细胞中的所有表达值
  m <- sum(data[i,])             #统计基因，在所有细胞中的所有表达值
  z = n/m
  res = cbind(i,n_all,n_adi,n,m,z)   #cbind依据列对数据进行合并
  result <- rbind(result,res)        #rbind依据行对数据进行合并
}
result <- data.frame(result)
colnames(result) <- c('gene','gene_exp_in_all','gene_exp_in_adi',
                      'count_in_adi','count_in_all','ko_efficiency')
write.csv(result,file = 'all_gene_ko_efficiency.csv',
          quote = FALSE,row.names = FALSE,sep = '\t')

