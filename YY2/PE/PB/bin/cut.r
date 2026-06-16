### vars_stat_SKII5264.txt ###
result1 = read.table("vars_stat.txt", sep="\t", header=T)
sample_names = read.table("sampleB.tsv", header=F)
result2 = result1[result1$sampleID %in% sample_names[,1],]
write.table(result2, "vars_stat_PB.txt", sep="\t", col.names=T, row.names=F, quote=F)

### qc_stat_SKII5264.txt ###
result1 = read.table("qc_stat_SKII5264.txt", sep="\t", header=F)
sample_names = read.table("sample_name.tsv", header=F)
result2 = result1[result1$V3 %in% sample_names[,1],]
write.table(result2, "qc_stat_SKII5264_PA.txt", sep="\t", col.names=F, row.names=F, quote=F)

### result.txt ###
result1 = read.csv("result.txt", sep="\t", header=F)
sample_names = read.table("sampleB.tsv", header=F)
result2 = result1[result1$V3 %in% sample_names[,1],]
write.table(result2, "result_PB.txt", sep="\t", col.names=F, row.names=F, quote=F)

### stat ###
library(dplyr)
library(stringr)
# 1. 读取数据
df <- read.table("result.txt", sep="\t", stringsAsFactors = FALSE, header = F)
# 2. 定义关键字
keywords <- c("SKII15277", "SKII15275", "SKII5264", "SKII15266", "SKII15265")
# 筛选
keywords_regex <- paste(keywords, collapse="|")
sub_df <- df[grep(keywords_regex, df$V1), ]
# 提取分组标签
sub_df$group <- regmatches(sub_df$V1, regexpr(keywords_regex, sub_df$V1))
# 清洗数值列 (以第6列为例)
cols_to_fix <- c(6, 17, 18)
for(i in cols_to_fix) {
  sub_df[,i] <- as.numeric(gsub("%", "", as.character(sub_df[,i])))
}
# 分组求平均
aggregate(cbind(V6, V12, V17, V18) ~ group, data = sub_df, FUN = mean, na.rm = TRUE)