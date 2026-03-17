## Gene Copy Number Analysis in yeast

# working directory:
# all genome fasta in one directory(genomes)
# pgr-repeat.sh in 'scripts' directory
# S288c pep_all fasta in working directory
# notice: in genome fasta file, there no '.' in sequence id. It should be modified to '_'


## usage: bash cnv.sh 

## options
while getopts "g:s:" opt; do
    case $opt in
        g) genome_dir=$OPTARG ;;
        s) genome_suffix=$OPTARG ;;
    esac
done


## --- get duplicated regions ---
echo "==> get duplicated regions..."
for f in ${genome_dir}/*"$genome_suffix"; do
    basename "$f" "$genome_suffix"
done > genome.lst

parallel -j 24 '
    mkdir -p repeat/{}
    bash scripts/pgr-repeat.sh ${genome_dir}/{}.Final.fasta repeat/{}
    perl scripts/duplicate_merge.pl -i repeat/{}/duplicated_regions.rg -o repeat/{}/duplicated_regions_merge.rg
    pgr fa range genomes/{}.Final.fasta -r repeat/{}/duplicated_regions_merge.rg -o repeat/{}/duplicated_regions.fa
' :::: genome.lst

## --- BLAST ---
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


# - lift coordinates -
# header:qseqid  qlen  sstrand  sseqid  slen  chrom  sstart  send  len_aa  len_ratio bitscore  pident  evalue  qcovs
echo "==> lift coordinates..."

parallel -j 24 --linebuffer --bar '
    python scripts/lift_coords.py repeat/{}/{}.tsv repeat/{}/{}_lift.tsv
' :::: genome.lst


## - merge hit -
echo "==> merging hits..."

parallel -j 24 --linebuffer --bar '
    python scripts/hit_merge.py -i repeat/{}/{}_lift.tsv -o repeat/{}/{}_merge.tsv
' :::: genome.lst



## - filter -
# pident、qcovs
echo "==> filter BLAST results..."

parallel -j 24 --linebuffer --bar '
    cat repeat/{}/{}_merge.tsv |
        tsv-filter --le 13:1.0e-20 --ge 12:90  --ge 14:90 |
        sort -k1,1 -k6,6 -k7,7n |
        uniq > repeat/{}/{}_filter.tsv
' :::: genome.lst


## --- Statistics on copies --
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


