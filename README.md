# mRNA structure in yeast

reference from https://github.com/wang-q/pars/blob/master/README.md#cut-mrna-alignments-and-extract-snp-list

## 0 Installation

```shell
# homebrew
brew install parallel pigz wget aria2 pv
brew install bcftools blast samtools mafft

brew tap brewsci/bio
brew install raxml

brew tap wang-q/tap
brew install faops lastz multiz sparsemem intspan

curl -fsSL https://raw.githubusercontent.com/wang-q/App-Egaz/master/share/check_dep.sh | bash

# perl modules
cpanm App::Fasops App::Rangeops App::Egaz
cpanm Statistics::ChisqIndep

# R packages
parallel -j 1 -k --line-buffer '
    Rscript -e '\'' if (!requireNamespace("{}", quietly = TRUE)) { install.packages("{}", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN") } '\''
    ' ::: \
        getopt gsubfn RSQLite sqldf sm remotes \
        extrafont ggplot2 scales gridExtra pander \
        readr plyr dplyr proto reshape ape

```

## 1 data download

genome download from article [*From genotype to phenotype with 1,086 near telomere-to-telomere yeast genomes*](https://www.nature.com/articles/s41586-025-09637-0)






## 2 prepare sequences

```shell
cd /scratch/wangq/wsn/T2T_pars/


# prep assembly
egaz template \
    Assemblies \
    --prep -o GENOMES \
    --min 1000 --about 1_000_000 \
    -v --repeatmasker "--species Fungi --parallel 12"

bash GENOMES/0_prep.sh


```

## 3 Align

```shell
cd /scratch/wangq/wsn/T2T_pars/alignment

egaz template \
    GENOMES/S288c \
    $(
        cat ../NewID_Xunhua.list |
            tr -s '[:space:]' '\n' |
            sed 's/^/GENOMES\//'
    ) \
    GENOMES/spar GENOMES/seub \
    --multi -o n791 \
    -v --parallel 12

bash n128/1_pair.sh


# clean
find . -mindepth 1 -maxdepth 3 -type d -name "*_raw"   | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr
```


## blast

```shell

cat DMS_map_seq_prediction_summary_with_sequences.csv | sed '1d' | cut -d "," -f 1,16 | sed -E 's/^/>/g' | sed -E 's/,/\n/g' >> sce_genes_dms.fasta

blastn -task blastn -evalue 1e-3 -num_threads 4 -num_descriptions 10 -num_alignments 10 -outfmt 0 \
    -dust yes -soft_masking true \
    -db S288c.fa -query sce_genes_dms.fasta -out sce_genes_dms.blast

cat DMS_map_seq_prediction_summary_with_sequences.csv | sed '1d' | awk -F ',' '{print $1"\t"$16"\t"$15}'  > DMS_structure/sce_genes_dms_folded.tab

```

```perl
# read gene score
{
    my $file_score = path( $dir_pars, "sce_Score.tab" )->absolute->stringify;
    print "Read $file_score\n";

    my (@non_exist);

    open my $fh, '<', $file_score;
    while ( my $line = <$fh> ) {
        chomp $line;
        my ( $gene, undef, $score ) = split /\t/, $line;

        if ( exists $gene_info_of->{$gene} ) {
            my @scores = split /\;/, $score;
            $gene_info_of->{$gene}{mF_score}    = App::Fasops::Common::mean(@scores);
            $gene_info_of->{$gene}{pars_scores} = [@scores];
        }
        else {
            push @non_exist, $gene;
        }
    }
    close $fh;

    if (@non_exist) {
        print scalar @non_exist, " genes don't exist in $file_gene\n";
        print join( " ", @non_exist, "\n" );
    }
    print "\n";
}



```



## gene filter

### create protein coding gene list

```shell
cd /scratch/wangq/wsn/T2T_pars/gene_filter_dms

# sgd/saccharomyces_cerevisiae.gff → gene list
cat ../sgd/saccharomyces_cerevisiae.gff |
    perl -nla -e '
        next if /^#/;
        next unless $F[2] eq q{gene};
        my $annotation = $F[8];
        $annotation =~ /ID=(.*);Name=/ or next;
        my $ID = $1;
        my $chr = $F[0];
        $chr =~ s/^chr//i;
        next if $chr eq q{mt}; # Skip genes on mitochondria
        print join qq{,}, $ID, qq{$chr($F[6]):$F[3]-$F[4]};
    ' \
    > gene_list.csv

mkdir -p genes
cat gene_list.csv |
    parallel --colsep ',' --no-run-if-empty --linebuffer -k -j 12 '
        >&2 echo {1}
        echo {2} | spanr cover stdin -o genes/{1}.yml
    '
spanr merge genes/*.yml -o genes.merge.yml
rm -fr genes

# overlapped regions
cut -d, -f 2 gene_list.csv |
    spanr coverage -m 2 stdin -o overlapped.yml

spanr statop \
    ../blast/S288c.sizes \
    genes.merge.yml overlapped.yml \
    --op intersect --all -o stdout |
    grep -v "^key" |
    perl -nla -F, -e '
        $F[4] == 0 and print $F[0];
    ' \
    > non-overlapped.lst

# DMS gene
cat non-overlapped.lst |
    grep -Fx -f <(cut -f 1 ../blast_dms/sce_genes_dms.blast.tsv) \
    > DMS-non-overlapped.lst

cat ../blast_dms/sce_genes_dms.blast.tsv |
    perl -nla -e '
        next if /^#/;
        my $ID = $F[0];
        my $chr = $F[2];
        next if $chr eq q{mt}; # Skip genes on mitochondria
        print join qq{,}, $ID, qq{$chr($F[5]):$F[3]-$F[4]};
    ' \
    > DMS_gene_list.csv

mkdir DMS
cat DMS_gene_list.csv |
    parallel --colsep ',' --no-run-if-empty --linebuffer -k -j 12 '
        echo {1}
        echo {2} | spanr cover stdin -o DMS/{1}.yml
    '
spanr merge DMS/*.yml -o DMS.merge.yml
rm -fr DMS

spanr some genes.merge.yml DMS-non-overlapped.lst -o genes.non-overlapped.yml
#spanr split mRNAs.non-overlapped.yml -o mRNAs

spanr some DMS.merge.yml DMS-non-overlapped.lst -o DMS.non-overlapped.yml
spanr split DMS.non-overlapped.yml -o DMS_new -s .yml


```

### intact mRNAs

```shell
cd /scratch/wangq/wsn/T2T_pars/gene_filter_dms

gzip -dcf ../alignment/n791_2/Scer_n790_Seub_refined/*.gz |
    pigz > Scer_n790_Seub.fas.gz

fasops covers -n S288c Scer_n790_Seub.fas.gz -o Scer_n790_Seub.yml

spanr statop \
        ../blast_dms/S288c.sizes \
        genes.non-overlapped.yml Scer_n790_Seub.json \
        --op intersect --all -o stdout |
        grep -v "^key" |
        perl -nla -F, -e '
            $F[2] == $F[4] and print $F[0];
        ' \
        > Scer_n790_Seub.intact.lst

spanr statop \
        ../blast_dms/S288c.sizes \
        genes.non-overlapped.yml Scer_n790_Seub.json \
        --op intersect --all -o Scer_n790_Seub.csv


wc -l *.lst ../blast_dms/*.tsv* |
    grep -v "total$" |
    datamash reverse -W |
    (echo -e "File\tCount" && cat) |
    mlr --itsv --omd cat


# outgroup Spar

gzip -dcf ../alignment/n791/Scer_n790_Spar_refined/*.gz | pigz > Scer_n790_Spar.fas.gz
fasops covers -n S288c Scer_n790_Spar.fas.gz -o Scer_n790_Spar.yml

spanr statop \
        ../blast_dms/S288c.sizes \
        genes.non-overlapped.yml Scer_n790_Spar.json \
        --op intersect --all -o stdout |
        grep -v "^key" |
        perl -nla -F, -e '
            $F[2] == $F[4] and print $F[0];
        ' \
        > Scer_n790_Spar.intact.lst

```

in mac
```shell

for i in $(cat list.txt); do
    yq -P -o=json eval ${i}.yml > ${i}.json
    rm ${i}.yml
done
rm list.txt

spanr statop ../blast/S288c.sizes genes.non-overlapped.yml Scer_n790_Seub.yml \
        --op intersect --all -o stdout |
        grep -v "^key" |
        perl -nla -F, -e '
            $F[2] == $F[4] and print $F[0];
        ' \
        > Scer_n790_Seub.intact.lst


yq -P -o=json eval Scer_n128_Seub.yml > Scer_n128_Seub.json



```


### Cut mRNA alignments and extract SNP list

```shell
cd /scratch/wangq/wsn/T2T_pars/gene_filter_dms

# slice fas.gz
mkdir DMS_Scer_n790_Spar
cat Scer_n790_Spar.intact.lst | parallel --no-run-if-empty --linebuffer -k -j 12 "
    fasops slice Scer_n790_Spar.fas.gz DMS_yml/{}.yml -n S288c -o DMS_Scer_n790_Spar/{}.fas
"

# SNP
mkdir SNP_Scer_n790_Spar
cat Scer_n790_Spar.intact.lst |
    parallel --no-run-if-empty --linebuffer -k -j 12 "
        fasops vars --outgroup --nocomplex DMS_Scer_n790_Spar/{}.fas -o stdout |
             sed 's/\$/\t{}/' \
            > SNP_Scer_n790_Spar/{}.tsv
        "

#loccation,REF,ALT,mutant_to,freq,occured,gene
cat SNP_Scer_n790_Spar/*.tsv |
    tsv-select -f 5,6,7,9,10,8,14 \
    > Scer_n790_Spar.SNPs.tsv


```