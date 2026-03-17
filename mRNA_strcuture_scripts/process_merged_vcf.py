import sys
import csv
import os
from bisect import bisect_right
from collections import defaultdict

## Usage: python process_vcf.py <input_vcf> <root_dir> <output_vcf>


def load_all_indices(root_dir, sample_names_from_vcf):
    """
    read all sample intervals data, and make a dic
    structure: { sample: { chrom: [(start1, end1), (start2, end2)] } }
    """
    master_index = {}
    
    for name in sample_names_from_vcf:
        folder_name = "S288cvs{}".format(name)
        tsv_path = os.path.join(root_dir, folder_name, "intervals.tsv")
        chrom_map = defaultdict(list)
        
        if os.path.exists(tsv_path):
            with open(tsv_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    chrom_map[row['Chromosome']].append((int(row['Start']), int(row['End'])))
            
            # sort for binary search
            for chrom in chrom_map:
                chrom_map[chrom].sort()
        
        master_index[name] = chrom_map
    return master_index

def is_in_range(pos, sorted_intervals):
    """Binary Search"""
    if not sorted_intervals:
        return False
    idx = bisect_right(sorted_intervals, (pos, float('inf'))) - 1
    if idx >= 0:
        start, end = sorted_intervals[idx]
        if start <= pos <= end:
            return True
    return False

def main():
    if len(sys.argv) < 4:
        print("Usage: python process_vcf.py <input_vcf> <root_dir> <output_vcf>")
        sys.exit(1)

    vcf_in = sys.argv[1]
    root_dir = sys.argv[2]
    vcf_out = sys.argv[3]

    # whitelist of meta lines to keep
    keep_meta = set([
        "##fileformat=VCFv4.1",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "##contig=<ID=I>",
        "##contig=<ID=II>",
        "##contig=<ID=III>",
        "##contig=<ID=IV>",
        "##contig=<ID=V>",
        "##contig=<ID=VI>",
        "##contig=<ID=VII>",
        "##contig=<ID=VIII>",
        "##contig=<ID=IX>",
        "##contig=<ID=X>",
        "##contig=<ID=XI>",
        "##contig=<ID=XII>",
        "##contig=<ID=XIII>",
        "##contig=<ID=XIV>",
        "##contig=<ID=XV>",
        "##contig=<ID=XVI>",
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">',
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">'
    ])

    with open(vcf_in, 'r') as fin, open(vcf_out, 'w') as fout:
        
        csv.field_size_limit(sys.maxsize)
        
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')
        
        col_configs = []
        # save (col index, sample, sample's dic)

        for row in reader:
            # 1. delete ## meta
            if row[0].startswith('##'):
                line = '\t'.join(row)
                if line in keep_meta:
                    fout.write(line + "\n")
                continue
            
            # 2. process header #CHROM
            if row[0].startswith('#CHROM'):
                try:
                    fmt_idx = row.index('FORMAT')
                    sample_names_vcf = row[fmt_idx + 1:]
                    
                    # read all sample intervasl
                    master_index = load_all_indices(root_dir, sample_names_vcf)
                    
                    # record each col and its data
                    for i, name in enumerate(sample_names_vcf):
                        col_idx = fmt_idx + 1 + i
                        col_configs.append((col_idx, name, master_index[name]))
                    
                    writer.writerow(row)
                except ValueError:
                    print("Error: no FORMAT col in VCF file")
                    sys.exit(1)
                continue

            # 3. process data row
            chrom = row[0]
            pos = int(row[1])

            # 4. process by col
            # each col and determine
            for col_idx, s_name, s_lookup in col_configs:
                if row[col_idx] == '.':
                    # get intervals of this sample
                    target_intervals = s_lookup.get(chrom, [])
                    
                    # binary search
                    if is_in_range(pos, target_intervals):
                        row[col_idx] = '0'
                    else:
                        row[col_idx] = '.'
            
            writer.writerow(row)
    print("complete")

if __name__ == "__main__":
    main()
