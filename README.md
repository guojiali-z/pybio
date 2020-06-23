# pybio
1. 先用ko_efficiency_analysis_code.R处理GSE133486RT条件下的单细胞数据，统计每个基因在adipoq-cre系统中的理论敲除效率，得到  all_gene_ko_efficiency.csv文件。
2. 再用scwat_koe_pcr.py脚本，处理依据实验数据整理的pcr敲除效率，得到皮下白色脂肪敲除实验的数据，并与单细胞数据计算的敲除效率合并。
3. 最后使用correlation_analysis.R对koe_sc和koe_pcr做相关性分析，lm模型分析（包括预测story genes的敲除效率）。
4. 在整理的代码里面，运行script.sh就可以按照顺序实现上面三个脚本的程序。
