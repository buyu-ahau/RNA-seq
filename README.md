# =======================
# DESeq2 差异表达分析流程脚本
# Author: buyu-ahau
# =======================

# 1. 加载必要的 R 包 -----------------------------------
library(DESeq2)
library(tidyverse)

# 2. 读取 Counts 文件 -----------------------------------
# 设置文件路径
file_path <- "/disk192/users_dir/buyu/范书豪RNA-seq/gene_counts.txt"

# 读取数据
count_data <- read.table(
  file_path,
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  row.names = 1
)

# 3. 清理数据（修正版） --------------------------------
# 只保留 count 矩阵部分（假设前5列为注释，第1列已做行名）
counts_matrix <- count_data[, 6:ncol(count_data)]

# 清理列名，去除冗余部分
original_colnames <- colnames(counts_matrix)
clean_colnames <- gsub("..06_bam_files.", "", original_colnames, fixed = TRUE)
clean_colnames <- gsub("_Aligned.sortedByCoord.out.bam", "", clean_colnames, fixed = TRUE)
colnames(counts_matrix) <- clean_colnames

print("【修正后】的样本名:")
print(colnames(counts_matrix))

# 4. 创建样本信息表 -------------------------------------
sample_names <- colnames(counts_matrix)
sample_conditions <- gsub("_\\d+", "", sample_names)
colData <- data.frame(
  row.names = sample_names,
  condition = as.factor(sample_conditions)
)

print("【修正后】生成的样本信息表:")
print(colData)

# 检查样本名是否匹配
if (!all(rownames(colData) == colnames(counts_matrix))) {
  stop("错误：样本名不匹配！请检查代码。")
} else {
  print("样本名匹配成功，可以继续分析。")
}

# 5. 构建 DESeqDataSet 对象 -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = colData,
  design = ~ condition
)

# 6. 预过滤低表达基因 -----------------------------------
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
print("预过滤完成。")

# 7. 运行 DESeq2 核心分析 -------------------------------
print("正在运行 DESeq2，这可能需要一两分钟...")
dds <- DESeq(dds)
print("DESeq2 分析完成！")

# 8. 获取差异表达结果 ------------------------------------
# 指定对比组，condition 的两个水平为：DON 和 CK
res <- results(dds, contrast = c("condition", "DON", "CK"))
print("差异表达结果已成功生成。")

# 9. 整理并查看结果 --------------------------------------
res_ordered <- as.data.frame(res[order(res$pvalue), ])
print("查看排名最靠前的差异基因：")
print(head(res_ordered))

# 10. 导出完整和显著差异基因结果 -------------------------
write.csv(res_ordered, file = "DESeq2_DON_vs_CK_all_genes.csv")

# 筛选显著差异基因
significant_genes <- subset(res_ordered, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(significant_genes, file = "DESeq2_DON_vs_CK_significant_genes.csv")

cat("\n分析流程全部完成！\n")
cat("共找到", nrow(significant_genes), "个显著差异基因 (padj < 0.05 & |log2FC| > 1)。\n")
cat("结果已保存到 'DESeq2_DON_vs_CK_all_genes.csv' 和 'DESeq2_DON_vs_CK_significant_genes.csv' 文件中。\n")
