#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 10:31:26 2020

@author: guo
"""

import pandas as pd
from collections import Counter
import os

path = '/Users/guo/Desktop/GEO_data/GSE133486/ko_efficiency_experiment/20200616'
os.chdir(path)
#os.getcwd()

############################################6月16号师姐整理的所有数据
#########################read curation data
df = pd.read_excel('KO_Efficiency_20200616.xlsx', index_col='symbol')
df = df[df.index.notnull()]  
df1 = pd.DataFrame(df.iloc[:,5:].fillna(value = '0'))

#########################找出符合记录要求的皮下wat的pcr数据
gene = []
evalue = [] 
num = len(df1.index)
for i in range(num):
    col = sorted(df1.iloc[i,:],reverse = True)   #排序的目的是首先让iwat出现
    for j in col:
        data = str.lower(j)                      #为了方便，全部改成小写
        if 'pcr:iwat:' in data:
            if ':~' in data:
                evl = data.split(':~')[-1]
            else:
                evl = data.split(':')[-1]
            gene.append(df1.index[i])
            evalue.append(evl)
            break
        else:
            if 'pcr:wat:' in data or 'pcr:scwat:' in data or 'pcr:iwat(chow):' in data:
                if ':~' in data:
                    evl = data.split(':~')[1]
                else:
                    evl = data.split(':')[-1]
                gene.append(df1.index[i])
                evalue.append(evl)

df_out = pd.DataFrame(evalue,index=gene,columns=['koe_pcr'])     #将选出的gene和pcr生成df数据
df_out = df_out.astype('float')                                  #koe_pcr改成float
df_out = df_out[df_out['koe_pcr']>0]                             #删除<=0之后有67个基因
frequency = Counter(df_out.index.to_list())                      #统计每个基因出现的频次,有56个基因,主要是查看有多篇文献记录的基因
print(Counter(list(frequency.values())))                         #一共有48个基因是有唯一的记录

##########输出唯一记录的基因名
select_unique_data = []
for key, value in frequency.items():
    if value == 1:
        select_unique_data.append(key)
    else:
        print(key)                   #Nr3c1,Kdm1a,Mtor,Insr,Stk11,Gps2,Klb,Ghr这8个基因有多个结果
df_out1 = df_out.loc[select_unique_data,:]

#########################将基因结果与koe_sc结果合并
df_all_ko = pd.read_csv('all_gene_ko_efficiency.tsv',index_col='gene',header = 0)
gene_select = []
for i in df_out1.index.to_list():
    if i in df_all_ko.index:
        gene_select.append(i)                         #筛选出scrna有数据的基因,有47个基因
    else:
        print(i)                                      #Tp53没有scRNA数据

df_gene = df_all_ko.loc[gene_select,'count_in_adi':]
df_out1 = df_out1.loc[gene_select,:]
df_out1 = pd.concat([df_gene,df_out1],axis=1)         #将列合并，总列数增加
df_out1.columns = ['count_adi','count_all','koe_sc','koe_pcr']
df_out1 = df_out1.fillna(value=0)                     #将缺失值填充为0
df_out1 = df_out1[-((df_out1['count_adi']==0)&(df_out1['count_all']==0)&(df_out1['koe_sc']==0))]  #删除都为0的行
df_out1.to_csv('scwat_pcr.csv')

#########################给陈一龙老师的数据，有pcr的记录，没有的记录基因
gene_nopcr = []
for i in df.index.to_list():
    if i not in df_out1.index and i in df_all_ko.index:       #判断在不在有pcr的结果中，且在有all_ko的基因中
        gene_nopcr.append(i)
df_nopcr = df_all_ko.loc[gene_nopcr,'count_in_adi':]
df_nopcr.columns = ['count_adi','count_all','koe_sc']
df_nopcr = pd.DataFrame(df_nopcr.fillna(value=0))             #填充缺失值，主要目的是找出全为0的基因进行删除
df_nopcr = df_nopcr[-((df_nopcr['count_adi']==0)&(df_nopcr['count_all']==0)&(df_nopcr['koe_sc']==0))]
df_out_all = pd.concat([df_out1,df_nopcr],axis=0)             #行合并，行数增加
df_out_all1 = df_out_all.drop_duplicates(keep='first')        #去重
df_out_all1.to_csv('curation_gene.csv')








