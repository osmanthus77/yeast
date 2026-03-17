#!/usr/bin/env Rscript

library(getopt)
library(dplyr)
library(readr)
library(tidyr)

# 参数配置
spec = matrix(
    c(
        "help",      "h", 0, "logical",   "help message",
        "suffix",    "s", 1, "character", "input file suffix, like DMS|PARS",
        "input_dir", "i", 1, "character", "input base path, like /scratch/wangq/wsn/T2T_pars/result",
        "name",      "n", 1, "character", "input name, like Scer_spar"
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

setwd(result_path)

# 创建输出目录
dir.create("freq_each", showWarnings = FALSE, recursive = TRUE)
dir.create("freq_10", showWarnings = FALSE, recursive = TRUE)

# 突变分类函数
classify_mut <- function(df) {
  df %>% mutate(
    mut_cat = case_when(
      mutant_to %in% c("A->G", "A->C", "T->G", "T->C") ~ "AT_GC",
      mutant_to %in% c("G->A", "C->A", "G->T", "C->T") ~ "GC_AT",
      TRUE ~ "Other"
    )
  )
}

# 卡方检验函数
run_stats <- function(df_sub) {
  d_st <- filter(df_sub, structure == "stem")
  d_lp <- filter(df_sub, structure == "loop")
  
  sa <- sum(d_st$mut_cat == "AT_GC"); sg <- sum(d_st$mut_cat == "GC_AT")
  la <- sum(d_lp$mut_cat == "AT_GC"); lg <- sum(d_lp$mut_cat == "GC_AT")
  
  # 计算 Ratio，防止除以 0
  ratio_s <- if ((sa + sg) > 0) sa / (sa + sg) else 0
  ratio_l <- if ((la + lg) > 0) la / (la + lg) else 0
  
  # 执行卡方
  mat <- matrix(c(sa, sg, la, lg), nrow = 2, byrow = TRUE)
  c_val <- NA; p_val <- NA
  if (all(rowSums(mat) > 0) && all(colSums(mat) > 0)) {
    test <- tryCatch(chisq.test(mat, correct = FALSE), error = function(e) list(statistic=NA, p.value=NA))
    c_val <- test$statistic; p_val <- test$p.value
  }
  
  # 返回结构化数据
  data.frame(
    structure = c("stem", "loop"),
    AT_GC = c(sa, la),
    GC_AT = c(sg, lg),
    AT_GC_ratio = c(ratio_s, ratio_l),
    X2 = c(NA, c_val),
    p_value = c(NA, p_val)
  )
}

groups <- c("cds", "utr", "syn", "nsy", "mRNA")

for (g in groups) {
  # 动态构建输入文件名：例如 data_SNPs_DMS_cds.tsv
  file_name <- paste0("data_SNPs_", suffix, "_", g, ".tsv")
  if (!file.exists(file_name)) next
  
  raw_data <- read_tsv(file_name, show_col_types = FALSE) %>%
      classify_mut() %>%
      mutate(
        sample_size = nchar(occured),
        rel_freq = freq / sample_size
      )
  max_f <- max(raw_data$freq)
  
  # --- 第一部分：按单一 Frequency 统计 (freq_each) ---
  
  stat_list_each <- list()
  snp_res <- list(); gene_res <- list(); stem_res <- list(); loop_res <- list()
  s_atgc_res <- list(); s_gcat_res <- list(); l_atgc_res <- list(); l_gcat_res <- list()
  each_chi <- list()

  for (i in 1:max_f) {
    n <- filter(raw_data, freq == i)
    d_stem <- filter(n, structure == "stem")
    d_loop <- filter(n, structure == "loop")
    
    # 基础统计量计算
    snp_res[[i]]  <- data.frame(name = i, SNPs = nrow(n))
    gene_res[[i]] <- data.frame(name = i, gene = n_distinct(n$gene))
    stem_res[[i]] <- data.frame(name = i, SNPs = nrow(d_stem))
    loop_res[[i]] <- data.frame(name = i, SNPs = nrow(d_loop))
    
    # AT <-> GC 计数
    s_atgc <- sum(d_stem$mut_cat == "AT_GC"); s_gcat <- sum(d_stem$mut_cat == "GC_AT")
    l_atgc <- sum(d_loop$mut_cat == "AT_GC"); l_gcat <- sum(d_loop$mut_cat == "GC_AT")
    
    s_atgc_res[[i]] <- data.frame(name = i, SNPs = s_atgc)
    s_gcat_res[[i]] <- data.frame(name = i, SNPs = s_gcat)
    l_atgc_res[[i]] <- data.frame(name = i, SNPs = l_atgc)
    l_gcat_res[[i]] <- data.frame(name = i, SNPs = l_gcat)
    
    stat_list_each[[i]] <- data.frame(
      structure = c("stem", "loop"),
      AT_GC = c(s_atgc, l_atgc),
      GC_AT = c(s_gcat, l_gcat)
    )
    
    # 卡方检验
    res <- run_stats(n)
    each_chi[[i]] <- res %>% mutate(freq_point = i)
  }
  
  # 批量输出各统计文件
  write_csv(bind_rows(stat_list_each),  paste0('freq_each/', suffix, '_', g, '_stat.csv'))
  write_csv(bind_rows(snp_res),    paste0('freq_each/', suffix, '_', g, '_stat_SNPs.csv'))
  write_csv(bind_rows(gene_res),   paste0('freq_each/', suffix, '_', g, '_stat_gene.csv'))
  write_csv(bind_rows(stem_res),   paste0('freq_each/', suffix, '_', g, '_stat_stem.csv'))
  write_csv(bind_rows(loop_res),   paste0('freq_each/', suffix, '_', g, '_stat_loop.csv'))
  write_csv(bind_rows(s_atgc_res), paste0('freq_each/', suffix, '_', g, '_stat_stem_AT_GC.csv'))
  write_csv(bind_rows(l_atgc_res), paste0('freq_each/', suffix, '_', g, '_stat_loop_AT_GC.csv'))
  write_csv(bind_rows(s_gcat_res), paste0('freq_each/', suffix, '_', g, '_stat_stem_GC_AT.csv'))
  write_csv(bind_rows(l_gcat_res), paste0('freq_each/', suffix, '_', g, '_stat_loop_GC_AT.csv'))
  write_csv(bind_rows(each_chi),   paste0('freq_each/', suffix, '_', g, '_stat_chi_square.csv'))

  # --- 第二部分：10% 分段统计 (freq_10)  ---
  
  stat_list_10 <- list()
  s10 <- list(); g10 <- list(); st10 <- list(); lp10 <- list()
  st_atgc10 <- list(); st_gcat10 <- list(); lp_atgc10 <- list(); lp_gcat10 <- list()
  freq10_chi <- list()
  
  for (i in 1:10) {
        lower <- (i - 1) / 10
        upper <- i / 10
        bin_label <- paste0((i-1)*10, "-", i*10, "%")
        
        # 根据每行计算出的 rel_freq 进行筛选
        # 第一个分箱包含 0，其余为左开右闭 (lower, upper]
        if (i == 1) {
            n_bin <- filter(raw_data, rel_freq >= lower & rel_freq <= upper)
        } else {
            n_bin <- filter(raw_data, rel_freq > lower & rel_freq <= upper)
        }
                
                
        d_st <- filter(n_bin, structure == "stem")
        d_lp <- filter(n_bin, structure == "loop")
        
        # 10% 分段基础统计
        s10[[i]]  <- data.frame(name = bin_label, SNPs = nrow(n_bin))
        g10[[i]]  <- data.frame(name = bin_label, gene = n_distinct(n_bin$gene))
        st10[[i]] <- data.frame(name = bin_label, SNPs = nrow(d_st))
        lp10[[i]] <- data.frame(name = bin_label, SNPs = nrow(d_lp))
        
        # 突变分类计数
        sa10 <- sum(d_st$mut_cat == "AT_GC"); sg10 <- sum(d_st$mut_cat == "GC_AT")
        la10 <- sum(d_lp$mut_cat == "AT_GC"); lg10 <- sum(d_lp$mut_cat == "GC_AT")
        
        st_atgc10[[i]] <- data.frame(name = bin_label, SNPs = sa10)
        st_gcat10[[i]] <- data.frame(name = bin_label, SNPs = sg10)
        lp_atgc10[[i]] <- data.frame(name = bin_label, SNPs = la10)
        lp_gcat10[[i]] <- data.frame(name = bin_label, SNPs = lg10)
        
        stat_list_10[[i]] <- data.frame(
          structure = c("stem", "loop"),
          AT_GC = c(sa10, la10),
          GC_AT = c(sg10, lg10)
        )
        
        # 卡方检验
        res10 <- run_stats(n_bin)
        freq10_chi[[i]] <- res10 %>% mutate(bin = bin_label)
  }

      # 输出 freq_10 目录下所有统计文件
      write_csv(bind_rows(stat_list_10), paste0('freq_10/', suffix, '_', g, '_stat_freq_10.csv'))
      write_csv(bind_rows(s10),       paste0('freq_10/', suffix, '_', g, '_stat_SNPs_freq_10.csv'))
      write_csv(bind_rows(g10),       paste0('freq_10/', suffix, '_', g, '_stat_gene_freq_10.csv'))
      write_csv(bind_rows(st10),      paste0('freq_10/', suffix, '_', g, '_stat_stem_freq_10.csv'))
      write_csv(bind_rows(lp10),      paste0('freq_10/', suffix, '_', g, '_stat_loop_freq_10.csv'))
      write_csv(bind_rows(st_atgc10), paste0('freq_10/', suffix, '_', g, '_stat_stem_AT_GC_freq_10.csv'))
      write_csv(bind_rows(lp_atgc10), paste0('freq_10/', suffix, '_', g, '_stat_loop_AT_GC_freq_10.csv'))
      write_csv(bind_rows(st_gcat10), paste0('freq_10/', suffix, '_', g, '_stat_stem_GC_AT_freq_10.csv'))
      write_csv(bind_rows(lp_gcat10), paste0('freq_10/', suffix, '_', g, '_stat_loop_GC_AT_freq_10.csv'))
      write_csv(bind_rows(freq10_chi),paste0('freq_10/', suffix, '_', g, '_stat_freq_10_chi_square.csv'))
      
}
