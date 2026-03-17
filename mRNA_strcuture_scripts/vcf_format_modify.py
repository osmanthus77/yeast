#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import re

if len(sys.argv) != 4:
    sys.stderr.write(
        "Usage: python vcf_format_modify_with_gap.py ref_alignment.fasta input.vcf output.tsv\n"
    )
    sys.exit(1)

aln_fasta = sys.argv[1]
input_vcf = sys.argv[2]
output_tsv = sys.argv[3]

# --------------------------------------------------
# ref in alignment fas
# --------------------------------------------------
ref_seq = []
with open(aln_fasta) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if ref_seq:
                break
        else:
            ref_seq.append(line)

ref_seq = "".join(ref_seq)

# alignment column -> ref block coordinate
aln2ref = {}
ref_pos = 0
for aln_pos, base in enumerate(ref_seq, start=1):
    if base != "-":
        ref_pos += 1
        aln2ref[aln_pos] = ref_pos
    else:
        aln2ref[aln_pos] = None

# --------------------------------------------------
# modify VCF
# --------------------------------------------------
chrom = None
block_start = None
contig_written = False

with open(input_vcf) as fin, open(output_tsv, "w") as fout:
    for line in fin:
        line = line.rstrip("\n")
            
        # -----------------------------
        # meta header
        # -----------------------------
        if line.startswith("##"):
            if line.startswith("##contig"):
                continue
            fout.write(line + "\n")
            continue
        
        # -----------------------------
        # column header
        # -----------------------------
        if line.startswith("#CHROM"):
            header = line.split("\t")
            ref_field = header[9]
            sample = header[10]

            # S288c.IV(+):536918-538149
            m = re.search(r"\.(\w+)\([^)]*\):(\d+)-(\d+)", ref_field)

            chrom = m.group(1)
            block_start = int(m.group(2))
            
            if not contig_written:
                fout.write("##contig=<ID={}>\n".format(chrom))
                contig_written = True
            
            header[0] = "#CHROM"
            header[9] = "S288c"
            header[10] = sample.split(".")[0]
            fout.write("\t".join(header[:11]) + "\n")
            continue

        # -----------------------------
        # data line
        # -----------------------------
        fields = line.split("\t")
        aln_pos = int(fields[1])

        ref_block_pos = aln2ref.get(aln_pos)

        # gap of ref in alignment fas, ignore
        if ref_block_pos is None:
            continue

        genome_pos = block_start + ref_block_pos - 1

        fields[0] = chrom
        fields[1] = str(genome_pos)
        fout.write("\t".join(fields) + "\n")
