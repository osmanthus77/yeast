import os
import sys
import csv

## usage:python extract_intervals.py <root_path>

def main():
    if len(sys.argv) < 2:
        print("usage: python extract_intervals.py <root_path>")
        sys.exit(1)

    root_path = sys.argv[1]

    for dirpath, dirnames, filenames in os.walk(root_path):
        intervals = []
        fas_files = [f for f in filenames if f.endswith('.fas')]
        
        if not fas_files:
            continue

        for file_name in fas_files:
            try:
                # format: S288c.I.+.24176-24603.fas
                parts = file_name.split('.')
                chrom = parts[1]
                coord_part = parts[-2]
                start, end = map(int, coord_part.split('-'))
                intervals.append([chrom, start, end])
            except Exception:
                continue
        if intervals:
            # sort by chrom and start pos
            intervals.sort(key=lambda x: (x[0], x[1]))
            output_file = os.path.join(dirpath, "intervals.tsv")
            with open(output_file, 'w') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(['Chromosome', 'Start', 'End'])
                writer.writerows(intervals)
            print("complete: {}".format(output_file))

if __name__ == "__main__":
    main()
