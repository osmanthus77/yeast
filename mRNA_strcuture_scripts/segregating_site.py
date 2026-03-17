import pandas as pd
import argparse
import sys

def calculate_segregating_sites(input_file, output_file):
    try:
        # 1. 读取输入文件 (TSV格式)
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
        
        # 检查必要的列
        required_cols = ['gene', 'structure', 'fold_dot_length', 'fold_left_length', 'fold_right_length']
        for col in required_cols:
            if col not in df.columns:
                print(f"error: the below col not in input file： '{col}'")
                sys.exit(1)

        # 2. 定义 Stem 和 Loop 标记
        df['is_loop'] = df['structure'].get_label_name = df['structure'].str.contains('loop', case=False, na=False)
        df['is_stem'] = df['structure'].get_label_name = df['structure'].str.contains('stem', case=False, na=False)

        # 3. 按基因分组统计
        # - count(): 统计行数即为 SNP 数量
        # - first(): 长度信息取第一行即可
        results = []
        for gene, group in df.groupby('gene'):
            # 统计 SNP 行数
            loop_snps = group['is_loop'].sum()
            stem_snps = group['is_stem'].sum()
            
            # 获取长度信息 (取该基因内的第一个值)
            loop_len = group['fold_dot_length'].iloc[0]
            stem_len = group['fold_left_length'].iloc[0] + group['fold_right_length'].iloc[0]
            
            # 计算 Segregating Sites (SNP / Length)
            # 避免除以 0
            loop_ss = loop_snps / loop_len if loop_len > 0 else 0
            stem_ss = stem_snps / stem_len if stem_len > 0 else 0
            
            results.append({
                'gene': gene,
                'stem_length': stem_len,
                'stem_snps': stem_snps,
                'stem_segregating_sites': stem_ss,
                'loop_length': loop_len,
                'loop_snps': loop_snps,
                'loop_segregating_sites': loop_ss
            })

        # 4. 转换为 DataFrame 并输出
        output_df = pd.DataFrame(results)
        
        # 输出格式
        final_cols = [
            'gene', 'stem_length', 'stem_snps', 'stem_segregating_sites',
            'loop_length', 'loop_snps', 'loop_segregating_sites'
        ]
        output_df[final_cols].to_csv(output_file, sep='\t', index=False)
        
        print(f"complete！")

    except Exception as e:
        print(f"error in procee: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="calculating Stem and Loop Segregating Sites")
    parser.add_argument("-i", "--input", required=True, help="input tsv file")
    parser.add_argument("-o", "--output", required=True, help="output tsv file")
    
    args = parser.parse_args()
    calculate_segregating_sites(args.input, args.output)

if __name__ == "__main__":
    main()
