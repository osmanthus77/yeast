import pandas as pd
import numpy as np
import os
import argparse
from scipy.special import hyp1f1 # 新增引入合流超几何函数

def preprocess_data(file_path, col_name, manual_label=None):
    """
    预处理：提取标签、计算SNP比例并转置
    """
    # 提取标签
    if manual_label:
        label = manual_label
    else:
        file_name = os.path.basename(file_path)
        parts = file_name.split('_')
        label = parts[1] if len(parts) > 1 else "unknown"
    
    # 读取 CSV 并计算比例
    df = pd.read_csv(file_path)
    
    if col_name not in df.columns:
        available_cols = ", ".join(df.columns)
        raise ValueError(f"error：col name '{col_name}' does no exist。colnames are below: {available_cols}")
        
    target_data = df[col_name]
    total_sum = target_data.sum()
    
    if target_data.dtype == object:
        target_data = target_data.str.replace('%', '').str.replace(',', '').astype(float)
    
    stem_vector = (target_data / total_sum).values
    return label, stem_vector


def calculate_mse(r, n, stem_data):
    """
    核心计算逻辑：利用合流超几何函数替代数值积分
    避免了 comb 导致的数值溢出，并将 O(N) 的积分计算压缩为向量化 O(1)
    """
    x_vals = np.arange(1, n)
    
    if abs(r) < 1e-9:
        # r 趋近于 0 时的解析解
        b = 1.0 / x_vals
    else:
        # 核心解析解：完全替代原有 quad 积分与 comb 组合数
        hyp_term = hyp1f1(n - x_vals, n, -2 * r)
        b = (n / (x_vals * (n - x_vals))) * (1 - hyp_term) / (1 - np.exp(-2 * r))

    # 分箱逻辑：将计算出的数组平均分为 10 组并计算均值
    d = n / 10.0
    c = np.array([np.mean(b[int(i*d):int((i+1)*d)]) for i in range(10)])
    m = np.mean(c)
    y = 0.1 * c / m
    
    return np.mean((y - stem_data)**2)


def grid_search(r_range, step, n, stem_data):
    """
    执行两阶段网格搜索寻找最小 MSE
    """
    best_r, min_mse = r_range[0], float('inf')
    for r in np.arange(r_range[0], r_range[1] + step, step):
        current_r = 1e-9 if abs(r) < 1e-12 else r
        mse = calculate_mse(current_r, n, stem_data)
        if mse < min_mse:
            min_mse, best_r = mse, r
    return best_r, min_mse


def main():
    # 参数处理保持不变
    parser = argparse.ArgumentParser(description="calculate gama value with snp proporation.")
    parser.add_argument("-i", "--input", type=str, required=True, help="input CSV file (eg. PARS_cds_...csv)")
    parser.add_argument("-o", "--output", type=str, required=True, help="output txt file")
    parser.add_argument("-n", type=int, default=790, help="model config n")
    parser.add_argument("-l", "--label", type=str, default=None, help="manual label (optional)")
    parser.add_argument("-c", "--colname", type=str, default="SNPs", help="Column name to calculate (default: SNPs)")
    
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"error: can't find file {args.input}")
        return

    print(f"--- preprocess: {args.input} ---")
    label, stem_data = preprocess_data(args.input, args.colname, args.label)
    print(f"extract label: {label}")

    print(f"--- calculating: (n={args.n}) ---")
    
    # 网格搜索逻辑保持不变
    r_coarse, mse_coarse = grid_search((1, 100), 1, args.n, stem_data)
    r_best, mse_best = grid_search((r_coarse - 0.99, r_coarse + 0.99), 0.01, args.n, stem_data)
    
    if abs(r_best - 0.01) < 1e-5:
        r_neg_c, mse_neg_c = grid_search((-100, -1), 1, args.n, stem_data)
        r_neg_f, mse_neg_f = grid_search((r_neg_c - 0.99, r_neg_c + 0.99), 0.01, args.n, stem_data)
        if mse_neg_f < mse_best:
            r_best = r_neg_f

    print(f"complete. gama r: {r_best:.4f}")
    with open(args.output, "a") as f:
        f.write(f"{label}\t{r_best:.4f}\n")

if __name__ == "__main__":
    main()
