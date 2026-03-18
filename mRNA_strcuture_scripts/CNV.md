# Gene Copy Number Analysis in yeast

- [Gene Copy Number Analysis in yeast](#gene-copy-number-analysis-in-yeast)
  - [Get duplicated regions in genomes](#get-duplicated-regions-in-genomes)
  - [BLAST](#blast)
  - [Process blast result](#process-blast-result)
    - [Lift coordinates](#lift-coordinates)
    - [Merge hit](#merge-hit)
    - [Filter](#filter)
  - [Statistics on copies](#statistics-on-copies)


## Get duplicated regions in genomes

```shell
mkdir -p ~/project/yeast/cnv/genomes
cd ~/project/yeast/cnv

echo "==> get duplicated regions..."
for f in genomes/*.Final.Fasta; do
    basename "$f" ".Final.Fasta"
done > genome.lst

parallel -j 24 '
    mkdir -p repeat/{}
    bash scripts/pgr-repeat.sh ${genome_dir}/{}.Final.fasta repeat/{}
    perl scripts/duplicate_merge.pl -i repeat/{}/duplicated_regions.rg -o repeat/{}/duplicated_regions_merge.rg
    pgr fa range genomes/{}.Final.fasta -r repeat/{}/duplicated_regions_merge.rg -o repeat/{}/duplicated_regions.fa
' :::: genome.lst

```

## BLAST

```shell
cd ~/project/yeast/cnv

echo "==> BLAST..."

parallel -j 12 '
    mkdir -p db/{}
    makeblastdb \
        -in repeat/{}/duplicated_regions.fa \
        -dbtype nucl \
        -out db/{}/{}
' :::: genome.lst

parallel -j 24 --colsep '\t' '
    out_file="repeat/{}/{}.tsv"

    tblastn \
        -query S288c.pep.all.fa \
        -db db/{}/{} \
        -evalue 1e-5 \
        -outfmt "6 qseqid qstart qend qlen sstrand sseqid slen sstart send bitscore pident evalue qcovs" \
    | sed "s/plus/+/g;s/minus/-/g" \
    > "$out_file"
    
    [ ! -s "$out_file" ] && rm "$out_file"
' :::: genome.lst 

```

## Process blast result

### Lift coordinates

```shell
cd ~/project/yeast/cnv


# header:qseqid  qlen  sstrand  sseqid  slen  chrom  sstart  send  len_aa  len_ratio bitscore  pident  evalue  qcovs
echo "==> lift coordinates..."

parallel -j 24 --linebuffer --bar '
    python scripts/lift_coords.py repeat/{}/{}.tsv repeat/{}/{}_lift.tsv
' :::: genome.lst
```

### Merge hit

```shell
cd ~/project/yeast/cnv

echo "==> merging hits..."

parallel -j 24 --linebuffer --bar '
    python scripts/hit_merge.py -i repeat/{}/{}_lift.tsv -o repeat/{}/{}_merge.tsv
' :::: genome.lst
```


### Filter

```shell
# filter by pident、qcovs
cd ~/project/yeast/cnv

echo "==> filter BLAST results..."

parallel -j 24 --linebuffer --bar '
    cat repeat/{}/{}_merge.tsv |
        tsv-filter --le 13:1.0e-20 --ge 12:90  --ge 14:90 |
        sort -k1,1 -k6,6 -k7,7n |
        uniq > repeat/{}/{}_filter.tsv
' :::: genome.lst
```

## Statistics on copies

```shell
cd ~/project/yeast/cnv

cat genome.lst | xargs -I {} mkdir -p copy_num/{}

# all genes
parallel -j 24 --linebuffer --bar "
    cat repeat/{}/{}_filter.tsv | 
        cut -f 1 | sort | uniq -c |
        sort -k2,2n |
        sed -E 's/^ *([0-9]+) +/\1\t/' \
        > copy_num/{}/{}_copy_num.tsv
" :::: genome.lst

# gene copy = 1
parallel -j 24 --linebuffer --bar "
    cat copy_num/{}/{}_copy_num.tsv | 
        grep '^ *1\t' \
        > copy_num/{}/{}_copy_1.tsv
" :::: genome.lst

# gene copies >= 2
parallel -j 24 --linebuffer --bar "
    cat copy_num/{}/{}_copy_num.tsv | 
        grep -v '^ *1\t' \
        > copy_num/{}/{}_copies.tsv
" :::: genome.lst

# genome and its number of multiple copy genes
parallel -j 24 --linebuffer --bar "
    wc -l copy_num/{}/{}_copies.tsv |
        sed -E 's/^ *([0-9]+) +/\1\t/' |
        tr '/' '\t' | cut -f 1,3 |
        sed 's/.tsv//' \
        >> genome_gene.tsv
" :::: genome.lst
```

