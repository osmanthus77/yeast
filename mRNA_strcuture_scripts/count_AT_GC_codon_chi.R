#!/usr/bin/env Rscript

library(getopt)
library(dplyr)
library(readr)
library(tidyr)

# --- 参数配置 ---
spec = matrix(
    c(
        "help",      "h", 0, "logical",   "brief help message",
        "suffix",    "s", 1, "character", "input file suffix, like DMS|PARS",
        "name",      "n", 1, "character", "input name (e.g., Scer_Spar)",
        "input_dir", "i", 1, "character", "base data path, slike /scratch/wangq/wsn/T2T_pars/result"
    ),
    byrow = TRUE,
    ncol = 5
)
opt = getopt(spec)

# 设置默认路径与工作空间
suffix   <- opt$suffix
name <- opt$name
base_path <- opt$input_dir
result_path <- paste0(base_path, "/", name)
setwd(result_path)

# 创建输出目录
dir.create("freq_each", showWarnings = FALSE, recursive = TRUE)
dir.create("freq_10", showWarnings = FALSE, recursive = TRUE)


# 定义 4D 和 tRNA 密码子集合
codons_4D <- c(
    "CCT->CCC", "CCT->CCA", "CCT->CCG", "CCC->CCT", "CCC->CCA", "CCC->CCG", "CCA->CCC", "CCA->CCT", "CCA->CCG", "CCG->CCC", "CCG->CCT", "CCG->CCA",
    "GTT->GTC", "GTT->GTA", "GTT->GTG", "GTC->GTT", "GTC->GTA", "GTC->GTG", "GTA->GTC", "GTA->GTT", "GTA->GTG", "GTG->GTC", "GTG->GTT", "GTG->GTA",
    "ACT->ACC", "ACT->ACA", "ACT->ACG", "ACC->ACT", "ACC->ACA", "ACC->ACG", "ACA->ACC", "ACA->ACT", "ACA->ACG", "ACG->ACC", "ACG->ACT", "ACG->ACA",
    "GCT->GCC", "GCT->GCA", "GCT->GCG", "GCC->GCT", "GCC->GCA", "GCC->GCG", "GCA->GCC", "GCA->GCT", "GCA->GCG", "GCG->GCC", "GCG->GCT", "GCG->GCA",
    "GGT->GGC", "GGT->GGA", "GGT->GGG", "GGC->GGT", "GGC->GGA", "GGC->GGG", "GGA->GGC", "GGA->GGT", "GGA->GGG", "GGG->GGC", "GGG->GGT", "GGG->GGA"
)

codons_tRNA <- c(
    "CTA->CTG", "CTG->CTA", "GTT->GTC", "GTT->GTA", "GTC->GTT", "GTC->GTA", "GTA->GTT", "GTA->GTC",
    "GGT->GGC", "GGC->GGT", "GCT->GCC", "GCT->GCA", "GCC->GCT", "GCC->GCA", "GCA->GCT", "GCA->GCC",
    "AGA->AGG", "AGG->AGA", "CGT->CGC", "CGT->CGA", "CGC->CGT", "CGC->CGA", "CGA->CGT", "CGA->CGC",
    "AAA->AAG", "AAG->AAA", "GAA->GAG", "GAG->GAA", "GAT->GAC", "GAC->GAT", "ACT->ACC", "ACT->ACA",
    "ACC->ACT", "ACC->ACA", "ACA->ACT", "ACA->ACC", "TAT->TAC", "TAC->TAT", "TCT->TCC", "TCT->TCA",
    "TCC->TCT", "TCC->TCA", "TCA->TCT", "TCA->TCC", "TCA->TCG", "TCG->TCA", "CAT->CAC", "CAC->CAT",
    "TTT->TTC", "TTC->TTT", "TGT->TGC", "TGC->TGT"
)

# 突变方向分类函数
classify_mut <- function(df) {
  df %>% mutate(
    mut_cat = case_when(
      mutant_to %in% c("A->G", "A->C", "T->G", "T->C") ~ "AT_GC",
      mutant_to %in% c("G->A", "C->A", "G->T", "C->T") ~ "GC_AT",
      TRUE ~ "Other"
    )
  )
}

# 统计与卡方检验函数
run_stats_with_chi <- function(df_sub) {
  d_st <- filter(df_sub, structure == "stem")
  d_lp <- filter(df_sub, structure == "loop")
  
  sa <- sum(d_st$mut_cat == "AT_GC"); sg <- sum(d_st$mut_cat == "GC_AT")
  la <- sum(d_lp$mut_cat == "AT_GC"); lg <- sum(d_lp$mut_cat == "GC_AT")
  
  # 执行卡方检验
  mat <- matrix(c(sa, sg, la, lg), nrow = 2, byrow = TRUE)
  c_val <- NA; p_val <- NA
  if (all(rowSums(mat) > 0) && all(colSums(mat) > 0)) {
    test <- chisq.test(mat, correct = FALSE)
    c_val <- test$statistic
    p_val <- test$p.value
  }
  
  data.frame(
    structure = c("stem", "loop"),
    AT_GC = c(sa, la),
    GC_AT = c(sg, lg),
    AT_GC_ratio = c(sa/(sa+sg), la/(la+lg)),
    X2 = c(NA, c_val),
    p_value = c(NA, p_val)
  )
}

# --- 数据读入与预处理 ---
raw_file <- paste0(result_path, "/data_SNPs_", suffix, "_syn_codon.tsv")
data_SNPs <- read_tsv(raw_file, show_col_types = FALSE) %>%
  classify_mut() %>%
  mutate(
    sample_size = nchar(as.character(occured)),
    rel_freq = freq / sample_size
  )


# --- 分组处理 (4D 和 tRNA) ---
groups <- list("4D" = codons_4D, "tRNA" = codons_tRNA)

for (g_name in names(groups)) {
  # 初始过滤密码子
  t_data <- data_SNPs %>% filter(Codon_to %in% groups[[g_name]])
  max_f <- max(t_data$freq, na.rm = TRUE)
  
  # --- 按单一 Frequency 统计 (freq_each) ---
  res_each <- list()
  snp_list <- list(); gene_list <- list(); stem_list <- list(); loop_list <- list()
  s_atgc_list <- list(); s_gcat_list <- list(); l_atgc_list <- list(); l_gcat_list <- list()

  for (i in 1:max_f) {
    n <- filter(t_data, freq == i)
    if(nrow(n) == 0) next
    
    # 基础统计
    snp_list[[i]]  <- data.frame(name = i, SNPs = nrow(n))
    gene_list[[i]] <- data.frame(name = i, gene = n_distinct(n$gene))
    stem_list[[i]] <- data.frame(name = i, SNPs = sum(n$structure == "stem"))
    loop_list[[i]] <- data.frame(name = i, SNPs = sum(n$structure == "loop"))
    
    # 统计分类计数用于输出
    stats_res <- run_stats_with_chi(n)
    res_each[[i]] <- stats_res %>% mutate(freq_point = i)
    
    s_atgc_list[[i]] <- data.frame(name = i, SNPs = stats_res$AT_GC[1])
    l_atgc_list[[i]] <- data.frame(name = i, SNPs = stats_res$AT_GC[2])
    s_gcat_list[[i]] <- data.frame(name = i, SNPs = stats_res$GC_AT[1])
    l_gcat_list[[i]] <- data.frame(name = i, SNPs = stats_res$GC_AT[2])
  }
  
  # 批量保存 freq_each 文件
  write_csv(bind_rows(res_each) %>% select(structure, AT_GC, GC_AT), file.path('freq_each', paste0(suffix, "_", g_name, '_stat.csv')))
  write_csv(bind_rows(res_each), file.path('freq_each', paste0(suffix, "_", g_name, '_stat_chi_square.csv')))
  write_csv(bind_rows(snp_list), file.path('freq_each', paste0(suffix, "_", g_name, '_stat_SNPs.csv')))
  write_csv(bind_rows(gene_list), file.path('freq_each', paste0(suffix, "_", g_name, '_stat_gene.csv')))
  write_csv(bind_rows(stem_list), file.path('freq_each', paste0(suffix, "_", g_name, '_stat_stem.csv')))
  write_csv(bind_rows(loop_list), file.path('freq_each', paste0(suffix, "_", g_name, '_stat_loop.csv')))
  write_csv(bind_rows(s_atgc_list), file.path('freq_each', paste0(suffix, "_", g_name, '_stat_stem_AT_GC.csv')))
  write_csv(bind_rows(l_atgc_list), file.path('freq_each', paste0(suffix, "_", g_name, '_stat_loop_AT_GC.csv')))
  write_csv(bind_rows(s_gcat_list), file.path('freq_each', paste0(suffix, "_", g_name, '_stat_stem_GC_AT.csv')))
  write_csv(bind_rows(l_gcat_list), file.path('freq_each', paste0(suffix, "_", g_name, '_stat_loop_GC_AT.csv')))

  # --- 10% 分段统计 (freq_10) ---
  if (nrow(t_data) > 0) {
    res_10 <- list()
    s10 <- list(); g10 <- list(); st10 <- list(); lp10 <- list()
    st_atgc10 <- list(); lp_atgc10 <- list(); st_gcat10 <- list(); lp_gcat10 <- list()
    
    bin_bounds <- seq(0, 1, by = 0.1)
    
    for (i in 1:10) {
      bin_label <- paste0((i-1)*10, "-", i*10, "%")
      lower <- bin_bounds[i]
      upper <- bin_bounds[i+1]
      
      if(i == 1) {
          n_bin <- filter(t_data, rel_freq >= lower & rel_freq <= upper)
      } else {
          n_bin <- filter(t_data, rel_freq > lower & rel_freq <= upper)
      }

      if(nrow(n_bin) == 0) next
      
      stats_10 <- run_stats_with_chi(n_bin)
      res_10[[i]] <- stats_10 %>% mutate(bin = bin_label)
      
      s10[[i]]  <- data.frame(name = bin_label, SNPs = nrow(n_bin))
      g10[[i]]  <- data.frame(name = bin_label, gene = n_distinct(n_bin$gene))
      st10[[i]] <- data.frame(name = bin_label, SNPs = sum(n_bin$structure == "stem"))
      lp10[[i]] <- data.frame(name = bin_label, SNPs = sum(n_bin$structure == "loop"))
      
      st_atgc10[[i]] <- data.frame(name = bin_label, SNPs = stats_10$AT_GC[1])
      lp_atgc10[[i]] <- data.frame(name = bin_label, SNPs = stats_10$AT_GC[2])
      st_gcat10[[i]] <- data.frame(name = bin_label, SNPs = stats_10$GC_AT[1])
      lp_gcat10[[i]] <- data.frame(name = bin_label, SNPs = stats_10$GC_AT[2])
    }
          
    # 批量保存 freq_10 文件
    write_csv(bind_rows(res_10) %>% select(structure, AT_GC, GC_AT), file.path('freq_10', paste0(suffix, "_", g_name, '_stat_freq_10.csv')))
    write_csv(bind_rows(res_10), file.path('freq_10', paste0(suffix, "_", g_name, '_stat_freq_10_chi_square.csv')))
    write_csv(bind_rows(s10),  file.path('freq_10', paste0(suffix, "_", g_name, '_stat_SNPs_freq_10.csv')))
    write_csv(bind_rows(g10),  file.path('freq_10', paste0(suffix, "_", g_name, '_stat_gene_freq_10.csv')))
    write_csv(bind_rows(st10), file.path('freq_10', paste0(suffix, "_", g_name, '_stat_stem_freq_10.csv')))
    write_csv(bind_rows(lp10), file.path('freq_10', paste0(suffix, "_", g_name, '_stat_loop_freq_10.csv')))
    write_csv(bind_rows(st_atgc10), file.path('freq_10', paste0(suffix, "_", g_name, '_stat_stem_AT_GC_freq_10.csv')))
    write_csv(bind_rows(lp_atgc10), file.path('freq_10', paste0(suffix, "_", g_name, '_stat_loop_AT_GC_freq_10.csv')))
    write_csv(bind_rows(st_gcat10), file.path('freq_10', paste0(suffix, "_", g_name, '_stat_stem_GC_AT_freq_10.csv')))
    write_csv(bind_rows(lp_gcat10), file.path('freq_10', paste0(suffix, "_", g_name, '_stat_loop_GC_AT_freq_10.csv')))
  }
}
