# -*- coding: utf-8 -*-
import yaml
import sys

filter_yml = sys.argv[1]
genes_bed = sys.argv[2]

roman_map = {
    "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6,
    "VII": 7, "VIII": 8, "IX": 9, "X": 10, "XI": 11, "XII": 12,
    "XIII": 13, "XIV": 14, "XV": 15, "XVI": 16
}

records = []

## yml
with open(filter_yml, 'r') as f:
    content = f.read().replace('，', ',')
    gene_data = yaml.safe_load(content)

## change into bed
for gene_id, loc in gene_data.items():
    for chrom, pos_range in loc.items():
        start_str, end_str = pos_range.split('-')
        start = int(start_str) - 1  # 1-based -> 0-based
        end = int(end_str)
        records.append((chrom, start, end, gene_id))

records.sort(key=lambda x: (roman_map.get(x[0], x[0]), x[1]))

with open(genes_bed, 'w') as bed_out:
    for chrom, start, end, gene_id in records:
        bed_out.write("{}\t{}\t{}\t{}\n".format(chrom, start, end, gene_id))
