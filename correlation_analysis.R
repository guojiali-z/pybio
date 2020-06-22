library(ggplot2)
###################相关性分析
setwd('/Users/guo/Desktop/GEO_data/GSE133486/ko_efficiency_experiment/20200616')
ko_data <- read.csv(file = 'scwat_pcr.csv')
##########pearson corr coef
pearson_cor_coef <- cor(ko_data$koe_sc,ko_data$koe_pcr,method = 'pearson')        
pearson_p_value <- cor.test(ko_data$koe_sc,ko_data$koe_pcr,method = 'pearson')   
##########spearman corr coef
spearman_cor_coef <- cor(ko_data$koe_sc,ko_data$koe_pcr,method = 'spearman')
spearman_p_value <- cor.test(ko_data$koe_sc,ko_data$koe_pcr,method = 'spearman')
##########散点图
ggplot(ko_data,mapping = aes(x=koe_sc,y=koe_pcr)) + 
  geom_point(colour = 'pink',size=3) + 
  geom_text(label=paste(ko_data$gene),check_overlap = TRUE)

##########lm()
x <- ko_data$koe_sc
y <- ko_data$koe_pcr
lm.model <- lm(y~x)                     #lm(y~x)默认是计算截距，同lm(y~x+1),而lm(y~x-1)表示没有截距项
summary(lm.model)                       
coefficients(lm.model)                  #coefficients(lm.model)表示输出模型的参数估计值
#####对story genes做预测
koe_sc <- c(0.09859,0.15493,0.10811,0.13253,0.16260,0.06005,0.15094,0.06593,
            0.11864,0.10927,0.10739,0.12593,0.07246,0.12162,0.2,0.13226,0.07715,0.11189)
point <- data.frame(x=koe_sc)
rownames(point) <- c('Gsdmd','Kat8','Pcif1','Ctcf','Nsun2','Cyfip1','Fah','Map2k1',
                     'Bnip3l','Dnajc5','Nktr','Acss2','Eif2a','Prmt1','Spi1','Ldb1','Rbpj','Tardbp')
pred <- predict(lm(y~x),point,interval = "prediction",level = 0.95)
pred <- data.frame(pred)
numerator <- c(14,11,12,55,20,23,16,6,21,33,321,17,15,9,3,82,26,160)
denominator <- c(142,71,111,415,123,383,106,91,177,302,2989,135,207,74,15,620,337,1430)
pred <- cbind(numerator,denominator,koe_sc,pred[,1:ncol(pred)])
write.csv(pred,file = 'story_gene_koe_pred.csv',quote = FALSE,sep = '\t')




