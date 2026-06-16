library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(reshape2)

# 1. 加载数据
sample_data <- read.csv("summary_by_sample.csv")
target_data <- read.csv("summary_by_target.csv")

pdf("R_analysis_report.pdf", width=10, height=8)

# Plot 1: 样本层级 Ratio 分布
p1 <- ggplot(sample_data, aes(x=Sample, y=Ratio)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_hline(yintercept=1, linetype="dashed", color="red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title="Sample-level SE1/SE2 Depth Ratio", y="Ratio (SE1/SE2)")
print(p1)

# Plot 2: Target 层级散点图 (Avg SE1 vs SE2)
p2 <- ggplot(target_data, aes(x=Avg_SE1, y=Avg_SE2)) +
  geom_point(alpha=0.5, color="darkgreen") +
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed") +
  scale_x_log10() + scale_y_log10() +
  theme_bw() +
  labs(title="Target-level Avg Depth: SE1 vs SE2 (Log Scale)", 
       x="Average SE1 Depth", y="Average SE2 Depth")
print(p2)

# Plot 3: 筛选比值偏差最大的前 50 个 Target 进行热图展示
# 重新加载原始宽表以计算每个样本在每个 Target 的 Ratio
raw_df <- read.table("../SKII15275_PA_depth_processed.txt", sep="\t", header=TRUE, check.names=FALSE)
se1_idx <- grep("_Depth_SE1$", colnames(raw_df))
se2_idx <- grep("_Depth_SE2$", colnames(raw_df))

# 计算 Ratio 矩阵
ratio_mat <- as.matrix(raw_df[, se1_idx]) / as.matrix(raw_df[, se2_idx])
rownames(ratio_mat) <- raw_df$Target
colnames(ratio_mat) <- gsub("_Depth_SE1", "", colnames(raw_df)[se1_idx])

# 选取偏差大的 Target
target_data$abs_log_ratio <- abs(log2(target_data$Ratio))
top_targets <- target_data %>% arrange(desc(abs_log_ratio)) %>% head(50) %>% pull(Target)

heatmap_data <- log2(ratio_mat[top_targets, ])
heatmap_data[is.infinite(heatmap_data)] <- NA

pheatmap(heatmap_data, 
         main="Log2(SE1/SE2 Ratio) Heatmap - Top 50 Divergent Targets",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols = TRUE, cluster_rows = TRUE,
         na_col = "grey")

dev.off()
print("R 可视化报告生成完成。")
