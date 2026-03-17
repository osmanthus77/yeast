#!/usr/bin/env Rscript

library(getopt)
library(sqldf)
library(dplyr)
library(tidyr)
library(gsubfn)
library(proto)
library(RSQLite)


# 参数配置
spec = matrix(
    c(
        "help",      "h", 0, "logical",   "brief help message",
        "suffix",    "s", 1, "character", "input file suffix, like DMS|PARS",
        "input_dir", "i", 1, "character", "input base path, like /scratch/wangq/wsn/T2T_pars/result",
        "name",      "n", 1, "character", "input name",
        "region",    "r", 1, "character", "region need to run"
    ),
    byrow = TRUE,
    ncol = 5
)
opt = getopt(spec)

# 设置默认值与路径
suffix   <- opt$suffix
base_path <- opt$input_dir
name <- opt$name
result_path <- paste0(base_path, "/", name)
region <- opt$region

setwd(result_path)

# 数据读取与预处理
# 使用第一阶段 count_position_gene.pl 生成的文件
file_SNPs <- paste0(result_path, "/data_SNPs_", suffix, "_", region, "_pos.tsv")
data_SNPs <- read.csv(file_SNPs, header = TRUE, sep = "\t")

data_SNPs <- data_SNPs %>%
  mutate(
    sample_size = nchar(as.character(occured)),
    rel_freq = freq / sample_size
  )

# 初始化用于存储所有长度分组结果的列表
categories <- c("SNPs", "gene", "stem", "loop",
                "stem_AT_GC", "loop_AT_GC", "stem_GC_AT", "loop_GC_AT")
master_data <- list()
for(cat in categories) {
    # 预先建立频率区间列
    master_data[[cat]] <- data.frame(name = c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%",
                                             "50-60%", "60-70%", "70-80%", "80-90%", "90-100%"))
}

# 统计
group_names <- paste0(region, "_", 1:15)
bin_bounds <- seq(0, 1, by = 0.1)

for (g in group_names) {
    # 筛选对应 island_length 的数据
    len_val <- as.numeric(sub(paste0(region, "_"), "", g))
    t <- subset(data_SNPs, island_length == len_val)
    
    if(nrow(t) > 0) {
        
        # 临时存储当前长度 (mRNA_i) 的 10 个频率区间统计
        tmp_stats <- list()
        for(cat in categories) tmp_stats[[cat]] <- numeric(10)
        
        # 统计每个区间的SNP和gene
        for(i in 1:10) {
            # 频率分箱
            lower <- bin_bounds[i]
            upper <- bin_bounds[i+1]
            
            # 第一个区间包含 0，后续区间为左开右闭 (lower, upper]
            if(i == 1) {
                n <- subset(t, rel_freq >= lower & rel_freq <= upper)
            } else {
                n <- subset(t, rel_freq > lower & rel_freq <= upper)
            }
                        
            if(nrow(n) == 0) next
            
            # 统计总的SNP和gene
            tmp_stats[["SNPs"]][i] <- nrow(n)
            tmp_stats[["gene"]][i] <- n_distinct(n$gene)
            
            # Stem/Loop 区分
            data_stem <- subset(n, structure == "stem")
            data_loop <- subset(n, structure == "loop")
            
            tmp_stats[["stem"]][i] <- nrow(data_stem)
            tmp_stats[["loop"]][i] <- nrow(data_loop)
            
            # 使用 SQL 进行突变方向筛选
            tmp_stats[["stem_AT_GC"]][i] <- nrow(sqldf('SELECT * FROM data_stem WHERE mutant_to IN ("A->G", "A->C", "T->G", "T->C")'))
            tmp_stats[["stem_GC_AT"]][i] <- nrow(sqldf('SELECT * FROM data_stem WHERE mutant_to IN ("G->A", "C->A", "G->T", "C->T")'))
            tmp_stats[["loop_AT_GC"]][i] <- nrow(sqldf('SELECT * FROM data_loop WHERE mutant_to IN ("A->G", "A->C", "T->G", "T->C")'))
            tmp_stats[["loop_GC_AT"]][i] <- nrow(sqldf('SELECT * FROM data_loop WHERE mutant_to IN ("G->A", "C->A", "G->T", "C->T")'))
        }
        
        # 将当前长度的结果合并到 master_data 中
        for(cat in categories) {
            new_col <- data.frame(tmp_stats[[cat]])
            colnames(new_col) <- g
            master_data[[cat]] <- cbind(master_data[[cat]], new_col)
        }
    } else {
        # 如果该长度无数据，填充 0
        for(cat in categories) {
            master_data[[cat]][[g]] <- 0
        }
    }
}

# 输出
output_dir <- paste0(result_path, "/freq_10/stem_length")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for(cat in categories) {
    output_file <- paste0(output_dir, "/", suffix, "_", region, "_stat_", cat, "_freq_10.csv")
    write.table(master_data[[cat]], file = output_file, sep = ",", row.names = FALSE, quote = FALSE)
}
