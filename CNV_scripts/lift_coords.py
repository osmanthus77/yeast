import sys
import re

def lift_coordinates(input_file, output_file):
    try:
        with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
            # header
            # header = "qseqid\tqlen\tsstrand\tsseqid\tslen\tchrom\tsstart\tsend\tlen_aa\tlen_ratio\tbitscore\tpident\tevalue\tqcovs\n"
            # fout.write(header)
            
            rows = []
            
            for line in fin:
                if not line.strip():
                    continue
                cols = line.strip().split()
                
                if len(cols) < 11:
                    continue

                qseqid, qstart, qend, qlen, sstrand, sseqid, slen, sstart, send, bitscore, pident, evalue, qcovs = cols
                
                #  sseqid (format like V:334809-334888)
                match = re.match(r'(.+):(\d+)-(\d+)', sseqid)
                if not match:
                    print("Warning: Could not parse sseqid formatting: {}".format(sseqid))
                    continue
                
                chrom = match.group(1)
                start0 = int(match.group(2))
                # end0 = int(match.group(3)) # no need end0

                
                sstart_val = int(sstart)
                send_val = int(send)

                # lift coords
                if sstrand == '+':
                    sstart_lift = start0 + sstart_val - 1
                    send_lift = start0 + send_val - 1
                elif sstrand == '-':
                    sstart_lift = start0 + send_val - 1
                    send_lift = start0 + sstart_val - 1
                else:
                    print("Warning: Unknown strand '{}' in line: {}".format(sstrand, line.strip()))
                    continue
                len_aa = int((send_lift - sstart_lift + 1)/3)
                len_ratio = round(len_aa / int(qlen) * 100, 1)
                
                # output
                output_row = [
                    qseqid, qstart, qend, qlen, sstrand, sseqid, slen,
                    chrom, str(sstart_lift), str(send_lift), len_aa, len_ratio,
                    bitscore, pident, evalue, qcovs
                ]
                
                rows.append(output_row)
            
            rows.sort(key=lambda x: (x[0], x[7], int(x[8])))
                
            for r in rows:
                fout.write("\t".join(map(str, r)) + "\n")

        #print(f"Success! Processed coordinates saved to: {output_file}")

    except FileNotFoundError:
        print("Error: Input file not found.")
    except Exception as e:
        print("An error occurred: {}".format(e))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python lift_coords.py <input_file> <output_file>")
    else:
        lift_coordinates(sys.argv[1], sys.argv[2])
