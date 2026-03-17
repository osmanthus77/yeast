import pandas as pd
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Process and Merge TBLASTN result:Handling adjacent fragmentation and redundant overlap")
    parser.add_argument("-i", "--input", required=True, help="input TSV file")
    parser.add_argument("-o", "--output", required=True, help="merged output file")
    parser.add_argument("-g", "--gap_pct", type=float, default=20.0, help="two adjacent hit Gap threshold (default 20%)")
    parser.add_argument("-v", "--overlap_pct", type=float, default=80.0, help="overlapping threshold (default 80%)")
    return parser.parse_args()

## --- calculate Query sequence union length ---
def calculate_q_union_len(intervals):
    if not intervals: return 0
    # sort by start
    sorted_intervals = sorted([(min(a, b), max(a, b)) for a, b in intervals])
    merged = [sorted_intervals[0]]
    for curr in sorted_intervals[1:]:
        prev = merged[-1]
        if curr[0] <= prev[1]:
            merged[-1] = (prev[0], max(prev[1], curr[1]))
        else:
            merged.append(curr)
    return sum(end - start + 1 for start, end in merged)


## --- merge hits of the same query in the same chrom ---
def process_group(group, gap_limit, overlap_threshold_pct):
    """
    group: DataFrame contains the same group hits
    gap_limit: gap threshold (qlen * 3 * gap_pct)
    overlap_threshold_pct: hit oveplapping threshold
    """
    # sort by start (sstart)
    #  sstart <= send
    group = group.copy()
    group['real_sstart'] = group[['sstart', 'send']].min(axis=1)
    group['real_send'] = group[['sstart', 'send']].max(axis=1)
    hits = group.sort_values(by='real_sstart').to_dict('records')
    
    if not hits: return []

    merged_clusters = []
    current_cluster = [hits[0]]

    for i in range(1, len(hits)):
        curr = hits[i]
        prev = current_cluster[-1]
        
        # --- case 1: overlapping ratio calculation---
        q_overlap = max(0, min(prev['real_send'], curr['real_send']) - max(prev['real_sstart'], curr['real_sstart']))
        # shorter one to calculate ratio
        min_len = min(prev['len_aa'], curr['len_aa'])
        overlap_ratio = (q_overlap / min_len * 3) * 100 if min_len > 0 else 0

        if overlap_ratio >= overlap_threshold_pct:
            # overlapping ratio high(>=80), only Bitscore highest hit
            if curr['bitscore'] > prev['bitscore']:
                current_cluster[-1] = curr
            continue # next hit

        # --- case2: calculate gap(when no overlap)---
        gap = curr['real_sstart'] - prev['real_send']
        
        if gap <= gap_limit:
            # <= gap threshold, add to this cluster
            current_cluster.append(curr)
        else:
            # > gap threshold, end and start a new cluster
            merged_clusters.append(summarize_cluster(current_cluster))
            current_cluster = [curr]
    
    # the last cluster
    merged_clusters.append(summarize_cluster(current_cluster))
    return merged_clusters

def summarize_cluster(cluster):
    """merge Hits of the same cluster into one hit"""
    first = cluster[0]
    
    # pos: the left and right
    s_min = min(h['real_sstart'] for h in cluster)
    s_max = max(h['real_send'] for h in cluster)
    
    # calculate query union length coverage
    q_intervals = [(h['qstart'], h['qend']) for h in cluster]
    union_len = calculate_q_union_len(q_intervals)
    
    # cal len_aa and bitscor
    total_len_aa = sum(h['len_aa'] for h in cluster)
    total_bitscore = sum(h['bitscore'] for h in cluster)
    # pident
    weighted_pident = sum(h['pident'] * h['len_aa'] for h in cluster) / total_len_aa if total_len_aa > 0 else 0
    
    len_ratio = (total_len_aa / first['qlen']) * 100
    qcovs = (union_len / first['qlen']) * 100
    min_evalue = min(h['evalue'] for h in cluster)
    
    return {
        'qseqid': first['qseqid'],
        'qlen': first['qlen'],
        'sstrand': first['sstrand'],
        'sseqid': first['sseqid'],
        'slen': first['slen'],
        'chrom': first['chrom'],
        'sstart': s_min,
        'send': s_max,
        'len_aa': int(total_len_aa),
        'len_ratio': len_ratio,
        'bitscore': total_bitscore,
        'pident': weighted_pident,
        'evalue': min_evalue,
        'qcovs': qcovs
    }
    

def main():
    args = parse_args()
    
    # input header
    input_cols = [
        'qseqid', 'qstart', 'qend', 'qlen', 'sstrand',
        'sseqid', 'slen', 'chrom', 'sstart', 'send',
        'len_aa', 'len_ratio', 'bitscore', 'pident', 'evalue', 'qcovs'
    ]
    
    # read input
    try:
        df = pd.read_csv(args.input, sep='\t', names=input_cols)
    except Exception as e:
        print("input file reading error: {}".format(e))
        sys.exit(1)

    final_output = []

    # group by query and chrom
    for (q_id, chrom), group in df.groupby(['qseqid', 'chrom']):
        qlen = group['qlen'].iloc[0]
        # cal gap threshold: qlen * 3 * gap_pct%
        gap_limit = qlen * 3 * (args.gap_pct / 100.0)
        
        # merge
        clusters = process_group(group, gap_limit, args.overlap_pct)
        
        for c in clusters:
                c['len_ratio'] = "{:.1f}".format(c['len_ratio'])
                c['bitscore'] = int(round(c['bitscore']))
                c['pident'] = "{:.3f}".format(c['pident'])
                c['evalue'] = "{:.2e}".format(c['evalue'])
                c['qcovs'] = int(round(c['qcovs']))
                final_output.append(c)
    
    # output header
    output_cols = [
        'qseqid', 'qlen', 'sstrand', 'sseqid', 'slen', 'chrom',
        'sstart', 'send', 'len_aa', 'len_ratio', 'bitscore',
        'pident', 'evalue', 'qcovs'
    ]
    
    out_df = pd.DataFrame(final_output)
    if not out_df.empty:
        out_df[output_cols].to_csv(args.output, sep='\t', index=False, header=False)

if __name__ == "__main__":
    main()
