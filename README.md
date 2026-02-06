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


## 4 blast

```shell
cat DMS_map_seq_prediction_summary_with_sequences.csv | sed '1d' | cut -d "," -f 1,16 | sed -E 's/^/>/g' | sed -E 's/,/\n/g' >> sce_genes_dms.fasta

cd /scratch/wangq/wsn/T2T_pars/blast_dms

cat ../GENOMES/S288c/{I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI,Mito}.fa > S288c.fa
perl -nl -i -e '/^>/ or $_ = uc $_; print'  S288c.fa
faops size S288c.fa > S288c.sizes
makeblastdb -dbtype nucl -in S288c.fa -parse_seqids

# blast DMS transcripts
blastn -task blastn -evalue 1e-3 -num_threads 4 -num_descriptions 10 -num_alignments 10 -outfmt 0 \
    -dust yes -soft_masking true \
    -db S288c.fa -query sce_genes_dms.fasta -out sce_genes_dms.blast
perl ~/Scripts/pars/blastn_transcript.pl -f sce_genes_dms.blast -m 0

# blast PARS transcripts
cd /scratch/wangq/wsn/T2T_pars/blast_pars
cp ../blast_dms/S288c* .
blastn -task blastn -evalue 1e-3 -num_threads 4 -num_descriptions 10 -num_alignments 10 -outfmt 0 \
    -dust yes -soft_masking true \
    -db S288c.fa -query ../PARS10/sce_genes.fasta -out sce_genes_pars.blast
perl ../scripts/blastn_transcript.pl -f sce_genes_pars.blast -m 0
```



## 5 gene filter

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
    ../blast_dms/S288c.sizes \
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
        next if $chr eq q{mito}; # Skip genes on mitochondria
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
spanr split DMS.non-overlapped.yml -o DMS

## -- run the same steps with PARS data --

# non-overlapped region

# PARS gene
cd /scratch/wangq/wsn/T2T_pars/gene_filter_pars
cp ../gene_filter_dms/gene_list.csv ../gene_filter_dms/genes.merge.yml ../gene_filter_dms/overlapped.yml ../gene_filter_dms/non-overlapped.lst .

cat non-overlapped.lst |
    grep -Fx -f <(cut -f 1 ../blast_pars/sce_genes_pars.blast.tsv) \
    > PARS-non-overlapped.lst

cat ../blast_pars/sce_genes_pars.blast.tsv |
    perl -nla -e '
        next if /^#/;
        my $ID = $F[0];
        my $chr = $F[2];
        next if $chr eq q{mito}; # Skip genes on mitochondria
        print join qq{,}, $ID, qq{$chr($F[5]):$F[3]-$F[4]};
    ' \
    > PARS_gene_list.csv

spanr some genes.merge.yml DMS-non-overlapped.lst -o genes.non-overlapped.yml

```

### extract SNP

```shell
cd /scratch/wangq/wsn/T2T_pars/gene_filter_dms
mkdir -p Scer_n790_pair

ls ../alignment/n791/Pairwise > n791.lst
echo -e "I\nII\nIII\nIV\nV\nVI\nVII\nVIII\nIX\nX\nXI\nXII\nXIII\nXIV\nXV\nXVI" > chromosome.lst

for strain in $(cat n791.lst); do
    for i in $(cat chromosome.lst); do
        gzip -dcf ../alignment/n791/Pairwise/${strain}/mafSynNet/${i}.synNet.maf.gz |
            pigz >> Scer_n790_pair/${strain}.maf.gz
    done
done

cat n791.lst | parallel --no-run-if-empty --linebuffer -k -j 12 "
    fasops maf2fas Scer_n790_pair/{}.maf.gz -o Scer_n790_pair/{}.fas
    gzip Scer_n790_pair/{}.fas
    rm -f Scer_n790_pair/{}.maf.gz
"

#mkdir -p Scer_n790_yml Scer_n790_json
#cat n791.lst | parallel --no-run-if-empty --linebuffer -k -j 12 "
#   fasops covers -n S288c Scer_n790_pair/{}.fas.gz -o Scer_n790_yml/{}.yml
#    yq -P -o=json eval Scer_n790_yml/{}.yml > Scer_n790_json/{}.json
#"

mkdir -p Scer_n790_pair_refine
bsub -q mpi -n 48 bash Scer_n790_pair_refine/refine.sh
#rm -fr Scer_n790_pair

# fasta → block split
mkdir -p Scer_n790_pair_split
cat n791.lst | parallel --no-run-if-empty --linebuffer -k -j 12 "
    mkdir -p Scer_n790_pair_split/{}
    fasops split -r Scer_n790_pair_refine/{}.fas.gz -o Scer_n790_pair_split/{}
"

# extract SNPs by snp-sites
mkdir -p Scer_n790_vcf
cat n791.lst | parallel --no-run-if-empty --linebuffer -k -j 12 "
    ls Scer_n790_pair_split/{} > Scer_n790_pair_split/{}/fasta.lst
    sed -i '/fasta.lst/d' Scer_n790_pair_split/{}/fasta.lst
    sed -i 's/\.fas//g' Scer_n790_pair_split/{}/fasta.lst
"

for pair in $(cat n791.lst); do
    mkdir -p Scer_n790_vcf/${pair}
    for fas in $(cat Scer_n790_pair_split/${pair}/fasta.lst);do
        snp-sites -v Scer_n790_pair_split/${pair}/${fas}.fas -o Scer_n790_vcf/${pair}/${fas}.vcf
    done
done

# modify the CHROM and POS in vcf file
cat n791.lst | parallel --no-run-if-empty --linebuffer -k -j 12 "
    ls Scer_n790_vcf/{} > Scer_n790_vcf/{}/vcf.lst
    sed -i '/vcf.lst/d' Scer_n790_vcf/{}/vcf.lst
    sed -i 's/\.vcf//g' Scer_n790_vcf/{}/vcf.lst
"

for pair in $(cat n791.lst); do
    mkdir -p Scer_n790_vcf_modify/${pair}
    for vcf in $(cat Scer_n790_vcf/${pair}/vcf.lst);do
        python vcf_format_modify.py Scer_n790_pair_split/${pair}/${vcf}.fas Scer_n790_vcf/${pair}/${vcf}.vcf Scer_n790_vcf_modify/${pair}/${vcf}.vcf
    done
done

# vcf concat and sort
mkdir Scer_n790_vcf_concat
bsub -q mpi -n 48 bash Scer_n790_vcf_concat/concat.sh
find Scer_n790_vcf_concat -mindepth 2 -maxdepth 2 -name '*_sample.sort.vcf.gz' | sort | grep -v "S288cvsspar" > Scer_seub.lst
find Scer_n790_vcf_concat -mindepth 2 -maxdepth 2 -name '*_sample.sort.vcf.gz' | sort | grep -v "S288cvsseub" > Scer_spar.lst

mkdir Scer_n790_vcf_merge
bcftools merge -l Scer_seub.lst --output-type z --output Scer_n790_vcf_merge/Scer_seub.vcf.gz
bcftools merge -l Scer_spar.lst --output-type z --output Scer_n790_vcf_merge/Scer_spar.vcf.gz


# -- filter by DMS data --
python yml2bed.py genes.non-overlapped.yml genes.non-overlapped.bed
bgzip genes.non-overlapped.bed
tabix -p bed genes.non-overlapped.bed.gz
# outgroup seub
tabix -p vcf Scer_n790_vcf_merge/Scer_seub.vcf.gz
bcftools annotate \
  -a genes.non-overlapped.bed.gz \
  -c CHROM,FROM,TO,ID \
  -R genes.non-overlapped.bed.gz \
  -Oz -o Scer_n790_vcf_merge/Scer_seub_DMS.vcf.gz \
  Scer_n790_vcf_merge/Scer_seub.vcf.gz 

# outgroup spar
tabix -p vcf Scer_n790_vcf_merge/Scer_spar.vcf.gz
bcftools annotate \
  -a genes.non-overlapped.bed.gz \
  -c CHROM,FROM,TO,ID \
  -R genes.non-overlapped.bed.gz \
  -Oz -o Scer_n790_vcf_merge/Scer_spar_DMS.vcf.gz \
  Scer_n790_vcf_merge/Scer_spar.vcf.gz 

# process vcf 
for i in $(cat n791.lst);do 
    python Scer_n790_pair_split/extract_intervals.py Scer_n790_pair_split/${i}
done

bgzip -d Scer_n790_vcf_merge/Scer_spar_DMS.vcf.gz
python Scer_n790_vcf_merge/process_merged_vcf.py Scer_n790_vcf_merge/Scer_spar_DMS.vcf Scer_n790_pair_split Scer_n790_vcf_merge/Scer_spar_DMS_processed.vcf
#bgzip Scer_n790_vcf_merge/Scer_spar_DMS_processed.vcf

bgzip -d Scer_n790_vcf_merge/Scer_seub_DMS.vcf.gz
python Scer_n790_vcf_merge/process_merged_vcf.py Scer_n790_vcf_merge/Scer_seub_DMS.vcf Scer_n790_pair_split Scer_n790_vcf_merge/Scer_seub_DMS_processed.vcf
#bgzip Scer_n790_vcf_merge/Scer_seub_DMS_processed.vcf

for i in spar seub; do
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%TGT]\n' Scer_n790_vcf_merge/Scer_${i}_DMS_processed.vcf > Scer_n790_vcf_merge/Scer_${i}_DMS.snp
    perl Scer_n790_vcf_merge/polarize_snp_vcf.pl Scer_n790_vcf_merge/Scer_${i}_DMS.snp > Scer_n790_vcf_merge/Scer_${i}_DMS.SNPs.tsv
done

wc -l Scer_n790_vcf_merge/*.SNPs.tsv|
    grep -v "total$" |
    datamash reverse -W |
    (echo -e "File\tCount" && cat) |
    mlr --itsv --omd cat
```
| File | Count |
| --- | --- |
| Scer_n790_vcf_merge/Scer_seub_DMS.SNPs.tsv | 2162768 |
| Scer_n790_vcf_merge/Scer_spar_DMS.SNPs.tsv | 1996511 |


use PARS data to filter and extract SNP

```shell
mkdir -p /scratch/wangq/wsn/T2T_pars/gene_filter_pars
cd /scratch/wangq/wsn/T2T_pars/gene_filter_pars
ln -s ../gene_filter_dms/Scer_n790_vcf_merge Scer_n790_vcf_merge
ln -s ../gene_filter_dms/Scer_n790_pair_split Scer_n790_pair_split
cp ../gene_filter_dms/*py .

# -- filter by PARS data --
python yml2bed.py genes.non-overlapped.yml genes.non-overlapped.bed
bgzip genes.non-overlapped.bed
tabix -p bed genes.non-overlapped.bed.gz

mkdir -p SNP_PARS
# outgroup seub
bcftools annotate \
  -a genes.non-overlapped.bed.gz \
  -c CHROM,FROM,TO,ID \
  -R genes.non-overlapped.bed.gz \
  -Ov -o SNP_PARS/Scer_seub_PARS.vcf \
  Scer_n790_vcf_merge/Scer_seub.vcf.gz 

# outgroup spar
bcftools annotate \
  -a genes.non-overlapped.bed.gz \
  -c CHROM,FROM,TO,ID \
  -R genes.non-overlapped.bed.gz \
  -Ov -o SNP_PARS/Scer_spar_PARS.vcf \
  Scer_n790_vcf_merge/Scer_spar.vcf.gz 

# process vcf 
for i in spar seub; do
    python SNP_PARS/process_merged_vcf.py SNP_PARS/Scer_${i}_PARS.vcf Scer_n790_pair_split SNP_PARS/Scer_${i}_PARS_processed.vcf
done

for i in spar seub; do
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%TGT]\n' SNP_PARS/Scer_${i}_PARS_processed.vcf > SNP_PARS/Scer_${i}_PARS.snp
    perl SNP_PARS/polarize_snp_vcf.pl SNP_PARS/Scer_${i}_PARS.snp > SNP_PARS/Scer_${i}_PARS.SNPs.tsv
done

## TO RUN!!
for i in spar seub; do
    perl SNP_PARS/polarize_snp_vcf.pl SNP_PARS/Scer_${i}_PARS.snp > SNP_PARS/Scer_${i}_PARS.SNPs.tsv
done

wc -l SNP_PARS/*.SNPs.tsv|
    grep -v "total$" |
    datamash reverse -W |
    (echo -e "File\tCount" && cat) |
    mlr --itsv --omd cat

```
| File | Count |
| --- | --- |
| SNP_PARS/Scer_seub_PARS.SNPs.tsv | 1069926 |
| SNP_PARS/Scer_spar_PARS.SNPs.tsv | 960961 |


## 6 VEP

```shell
mkdir -p /scratch/wangq/wsn/T2T_pars/vep
cd /scratch/wangq/wsn/T2T_pars/vep

## -- DMS --
cp ../gene_filter_dms/Scer_n790_vcf_merge/Scer_spar_DMS.SNPs.tsv ../gene_filter_dms/Scer_n790_vcf_merge/Scer_seub_DMS.SNPs.tsv /scratch/wangq/wsn/T2T_pars/vep

for i in spar seub; do
    cat Scer_${i}_DMS.SNPs.tsv | 
        perl -nla -F"\t" -e '
            my $loc = $F[0];
            $loc =~ /^(.*):(.*)/;
            my $chr = $1;
            my $pos = $2;
            print qq{$chr\t$pos\t$pos\t$F[1]\t$F[2]};
        ' \
        > Scer_${i}_DMS.upload.tsv
done

## -- PARS --
cp ../gene_filter_pars/SNP_PARS/Scer_spar_PARS.SNPs.tsv ../gene_filter_pars/SNP_PARS/Scer_seub_PARS.SNPs.tsv .

for i in spar seub; do
    cat Scer_${i}_PARS.SNPs.tsv | 
        perl -nla -F"\t" -e '
            my $loc = $F[0];
            $loc =~ /^(.*):(.*)/;
            my $chr = $1;
            my $pos = $2;
            print qq{$chr\t$pos\t$pos\t$F[1]\t$F[2]};
        ' \
        > Scer_${i}_PARS.upload.tsv
done

wc -l *.upload.tsv |
    grep -v "total$" |
    datamash reverse -W |
    (echo -e "File\tCount" && cat) |
    mlr --itsv --omd cat
```
| File | Count |
| --- | --- |
| Scer_seub_DMS.upload.tsv | 2162768 |
| Scer_spar_DMS.upload.tsv | 1996511 |
| Scer_seub_PARS.upload.tsv | 1069926 |
| Scer_spar_PARS.upload.tsv | 960961 |

upload `upload.tsv` to https://asia.ensembl.org/Tools/VEP


```shell
cd /scratch/wangq/wsn/T2T_pars/vep

for i in spar seub; do
    cat Scer_${i}_DMS.vep.txt | 
        perl -nla -F"\t" -e '
            next if /^#/;
            my $loca = $F[1];
            $loca =~ /^(.*)-[0-9]+/;
            my $ID = $1;
            #location,allele,gene,consequence,CDS_position,amino_acids,codons,existing_variation
            print qq{$ID\t$F[2]\t$F[6]\t$F[3]\t$F[15]\t$F[17]\t$F[18]\t$F[19]};
        ' \
    > Scer_${i}_DMS.vep.tsv
done

for i in spar seub; do
    cat Scer_${i}_PARS.vep.txt | 
        perl -nla -F"\t" -e '
            next if /^#/;
            my $loca = $F[1];
            $loca =~ /^(.*)-[0-9]+/;
            my $ID = $1;
            #location,allele,gene,consequence,CDS_position,amino_acids,codons,existing_variation
            print qq{$ID\t$F[2]\t$F[6]\t$F[3]\t$F[15]\t$F[17]\t$F[18]\t$F[19]};
        ' \
    > Scer_${i}_PARS.vep.tsv
done

wc -l *.vep.tsv |
     grep -v "total$" |
     datamash reverse -W |
     (echo -e "File\tCount" && cat) |
     mlr --itsv --omd cat

```
| File | Count |
| --- | --- |
| Scer_seub_DMS.vep.tsv | 2163988 |
| Scer_spar_DMS.vep.tsv | 1997641 |
| Scer_seub_PARS.vep.tsv | 1070134 |
| Scer_spar_PARS.vep.tsv | 961142 |

## 7 process DMS data

```shell
cd /scratch/wangq/wsn/T2T_pars/DMS_process

perl ../scripts/blastn_transcript.pl -f ../blast_dms/sce_genes_dms.blast -m 0

cp ../../DMS_map_seq_prediction_summary_with_sequences.csv DMS/DMS_summary.csv
cat DMS/DMS_summary.csv | 
    perl -nla -F"," -e '
        next if $. == 1;
        print qq{$F[0]\t$F[15]\t$F[14]};
    ' \
> DMS/sce_genes_folded.tab


for i in spar seub; do
    cat ../vep/Scer_${i}_DMS.vep.tsv | 
        tsv-select -f 1 |
        sort -u \
        > Scer_${i}_DMS.snp.rg
    
    perl ../scripts/read_fold.pl \
        --dms DMS \
        --gene sce_genes_dms.blast.tsv \
        --pos Scer_${i}_DMS.snp.rg \
        > Scer_${i}_DMS_fail_pos.txt

    perl ../scripts/process_vars_in_fold.pl --file Scer_${i}_DMS.gene_variation.yml
done

## -- gene --
cat sce_genes_dms.blast.tsv |
    perl -nla -e '
        $F[2] eq "Mito" and next;
        print qq{$F[2]:$F[3]-$F[4]};
    ' |
    sort |
    spanr cover stdin -o sce_genes.yml

## -- intron --
cat ../sgd/orf_coding_all.fasta |
    perl -n -MAlignDB::IntSpan -e '
        />/ or next;
        /Chr\s+(\w+)\s+from\s+([\d,-]+)/ or next;
        $1 eq "Mito" and next;

        my $chr = $1;
        my $range = $2;
        my @ranges = sort { $a <=> $b } grep {/^\d+$/} split /,|\-/, $range;
        my $intspan = AlignDB::IntSpan->new()->add_range(@ranges);
        my $hole = $intspan->holes;

        printf qq{%s:%s\n}, $chr, $hole->as_string if $hole->is_not_empty;
    ' |
    spanr cover stdin -o sce_intron.yml

# produce orf_genomic set
cat ../sgd/orf_genomic_all.fasta |
    perl -n -e '
        />/ or next;
        /Chr\s+(\w+)\s+from\s+(\d+)\-(\d+)/ or next;
        $1 eq "Mito" and next;

        if ($2 == $3) {
            print qq{$1:$2\n};
        }
        elsif ($2 < $3) {
            print qq{$1:$2-$3\n};
        }
        else {
            print qq{$1:$3-$2\n};
        }
    ' |
    spanr cover stdin -o sce_orf_genomic.yml

## -- mRNA、utr、CDS --
spanr compare --op diff sce_genes.yml sce_orf_genomic.yml -o sce_utr.yml
spanr compare --op diff sce_genes.yml sce_intron.yml -o sce_mRNA.yml
spanr compare --op diff sce_mRNA.yml sce_utr.yml -o sce_cds.yml

for NAME in genes intron orf_genomic utr mRNA cds; do
    spanr stat ../blast_dms/S288c.sizes "sce_${NAME}.yml" --all |
        sed '1 s/^/Name,/' |
        sed "2 s/^/${NAME},/"
done |
    tsv-uniq |
    mlr --icsv --omd cat

```
| Name | chrLength | size | coverage |
| --- | --- | --- | --- |
| genes | 12071326 | 7534730 | 0.6242 |
| intron | 12071326 | 65519 | 0.0054 |
| orf_genomic | 12071326 | 8897088 | 0.7370 |
| utr | 12071326 | 8341 | 0.0007 |
| mRNA | 12071326 | 7474103 | 0.6192 |
| cds | 12071326 | 7465762 | 0.6185 |


```shell
cd /scratch/wangq/wsn/T2T_pars/PARS_process

perl ../scripts/blastn_transcript.pl -f ../blast_pars/sce_genes_pars.blast -m 0

for i in spar; do
    cat ../vep/Scer_${i}_PARS.vep.tsv | 
        tsv-select -f 1 |
        sort -u \
        > Scer_${i}_PARS.snp.rg
    
    perl ../scripts/read_fold.pl \
        --dms ../PARS10 \
        --gene sce_genes_pars.blast.tsv \
        --pos Scer_${i}_PARS.snp.rg \
        > Scer_${i}_PARS_fail_pos.txt

    perl ../scripts/process_vars_in_fold.pl --file Scer_${i}_PARS.gene_variation.yml
done

## -- gene --
cat sce_genes_pars.blast.tsv |
    perl -nla -e '
        $F[2] eq "Mito" and next;
        print qq{$F[2]:$F[3]-$F[4]};
    ' |
    sort |
    spanr cover stdin -o sce_genes.yml

## -- intron --
cat ../sgd/orf_coding_all.fasta |
    perl -n -MAlignDB::IntSpan -e '
        />/ or next;
        /Chr\s+(\w+)\s+from\s+([\d,-]+)/ or next;
        $1 eq "Mito" and next;

        my $chr = $1;
        my $range = $2;
        my @ranges = sort { $a <=> $b } grep {/^\d+$/} split /,|\-/, $range;
        my $intspan = AlignDB::IntSpan->new()->add_range(@ranges);
        my $hole = $intspan->holes;

        printf qq{%s:%s\n}, $chr, $hole->as_string if $hole->is_not_empty;
    ' |
    spanr cover stdin -o sce_intron.yml

# produce orf_genomic set
cat ../sgd/orf_genomic_all.fasta |
    perl -n -e '
        />/ or next;
        /Chr\s+(\w+)\s+from\s+(\d+)\-(\d+)/ or next;
        $1 eq "Mito" and next;

        if ($2 == $3) {
            print qq{$1:$2\n};
        }
        elsif ($2 < $3) {
            print qq{$1:$2-$3\n};
        }
        else {
            print qq{$1:$3-$2\n};
        }
    ' |
    spanr cover stdin -o sce_orf_genomic.yml

## -- mRNA、utr、CDS --
spanr compare --op diff sce_genes.yml sce_orf_genomic.yml -o sce_utr.yml
spanr compare --op diff sce_genes.yml sce_intron.yml -o sce_mRNA.yml
spanr compare --op diff sce_mRNA.yml sce_utr.yml -o sce_cds.yml

for NAME in genes intron orf_genomic utr mRNA cds; do
    spanr stat ../blast_dms/S288c.sizes "sce_${NAME}.yml" --all |
        sed '1 s/^/Name,/' |
        sed "2 s/^/${NAME},/"
done |
    tsv-uniq |
    mlr --icsv --omd cat

```


## 8 SNP analysis

### count per gene GC content 

```shell
cd /scratch/wangq/wsn/T2T_pars/

for i in spar seub; do
    mkdir -p result/Scer_${i}

    perl scripts/count_ACGT_percent.pl \
        --file DMS_process/Scer_${i}_DMS.gene_variation.process.yml \
        --varfold DMS_process/Scer_${i}_DMS.gene_variation.fold_class.tsv \
        --output result/Scer_${i}/fold_class.tsv
    
    datamash check < result/Scer_${i}/fold_class.tsv
done

# 5505 lines, 43 fields

```

### count SNPs and gene

```shell
cd /scratch/wangq/wsn/T2T_pars/

for i in spar seub; do
    mkdir -p result/Scer_${i}

    tsv-join -z \
        vep/Scer_${i}_DMS.SNPs.tsv \
        -f vep/Scer_${i}_DMS.vep.tsv \
        --key-fields 1 \
        --append-fields 2-8 |
        perl -nla -F"\t" -e '
            if ($F[8] eq "-" || $F[6] eq $F[8]){
                splice @F, 8, 1;
                my $line = join ("\t", @F);
                print qq{$line};
            }
            BEGIN{
                print qq{location\tREF\tALT\tmutant_to\tfreq\toccured\tgene\tallele\tconsequence\tCDS_position\tamino_acids\tcodons\texisting_variation};
            }
        ' \
        > result/Scer_${i}/SNPs.vep.tsv
    datamash check < result/Scer_${i}/SNPs.vep.tsv
done
# 1996096 lines, 13 fields
# 2162295 lines, 13 fields

for i in spar seub; do
    cat result/Scer_${i}/SNPs.vep.tsv |
        tsv-join \
            -f result/Scer_${i}/fold_class.tsv \
            -H --key-fields gene \
            --append-fields 2-43 \
        > result/Scer_${i}/SNPs.fold_class.tsv
    datamash check < result/Scer_${i}/SNPs.fold_class.tsv
done
# 1996096 lines, 55 fields
# 2162295 lines, 55 fields

for i in spar seub; do
    cat DMS_process/Scer_${i}_DMS.gene_variation.var_pars.tsv |
        tsv-select -H -e gene |
        sed '1 s/^name/location/' |
        tsv-join \
            -f result/Scer_${i}/SNPs.fold_class.tsv\
            -H --key-fields location \
            --append-fields 2-55 \
        > result/Scer_${i}/data_SNPs_DMS_mRNA.tsv
    datamash check < result/Scer_${i}/data_SNPs_DMS_mRNA.tsv
done
# 1995731 lines, 61 fields
# 2162289 lines, 61 fields

cat result/Scer_spar/data_SNPs_DMS_mRNA.tsv |
    tsv-summarize -H --count --group-by consequence
#consequence     count
#stop_lost       2326
#missense_variant        1136629
#synonymous_variant      777211
#stop_gained     53956
#splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant        408
#splice_polypyrimidine_tract_variant,intron_variant      605
#intron_variant  16063
#splice_region_variant,intron_variant    152
#splice_donor_variant    116
#start_lost      2925
#missense_variant,splice_region_variant  177
#stop_retained_variant   1683
#splice_region_variant,synonymous_variant        113
#splice_acceptor_variant 102
#splice_donor_region_variant,intron_variant      160
#intergenic_variant      3016
#splice_donor_5th_base_variant,intron_variant    51
#stop_gained,splice_region_variant       17
#start_lost,splice_region_variant        8
#upstream_gene_variant   6
#splice_acceptor_variant,synonymous_variant      1
#start_lost,synonymous_variant   1
#splice_donor_region_variant,synonymous_variant  2
#5_prime_UTR_variant     1
#missense_variant,splice_donor_region_variant    1


for i in spar seub; do
    cat result/Scer_${i}/data_SNPs_DMS_mRNA.tsv |
        tsv-filter -H --str-ne "CDS_position:-" \
        > result/Scer_${i}/data_SNPs_DMS_cds.tsv

    cat result/Scer_${i}/data_SNPs_DMS_mRNA.tsv |
        tsv-filter -H --str-eq "CDS_position:-" \
        > result/Scer_${i}/data_SNPs_DMS_utr.tsv

    cat result/Scer_${i}/data_SNPs_DMS_mRNA.tsv |
        tsv-filter -H --or \
            --str-eq "consequence:stop_retained_variant" \
            --str-eq "consequence:synonymous_variant" \
        > result/Scer_${i}/data_SNPs_DMS_syn.tsv

    cat result/Scer_${i}/data_SNPs_DMS_mRNA.tsv |
        tsv-filter -H --or \
            --str-eq "consequence:missense_variant" \
            --str-eq "consequence:start_lost" \
            --str-eq "consequence:stop_gained" \
            --str-eq "consequence:stop_lost" \
        > result/Scer_${i}/data_SNPs_DMS_nsy.tsv
done


for i in spar seub; do
    printf "Area\tSNPs\tGenes\n" > result/Scer_${i}/count.tsv
    for AREA in mRNA cds utr syn nsy; do
        echo ${AREA}
        cat result/Scer_${i}/data_SNPs_DMS_${AREA}.tsv |
            tsv-summarize -H --count |
            sed '1d'
        cat result/Scer_${i}/data_SNPs_DMS_${AREA}.tsv |
            tsv-summarize -H --unique-count gene |
            sed '1d'
    done |
    paste - - - \
    >> result/Scer_${i}/count.tsv
done

```


### count A/T <-> G/C

```shell
cd /scratch/wangq/wsn/T2T_pars/result

for i in spar seub; do
    Rscript ../scripts/count_AT_GC_chi_square.R \
    -i /scratch/wangq/wsn/T2T_pars/result \
    -s DMS \
    -n Scer_${i}
done

```

### count stem length selection

```shell
cd /scratch/wangq/wsn/T2T_pars/result/Scer_spar
mkdir -p freq_10/stem_length

perl ../../scripts/count_position_gene.pl \
    --file ../../DMS_process/Scer_spar_DMS.gene_variation.process.yml \
    --origin data_SNPs_DMS_mRNA.tsv \
    --output data_SNPs_DMS_mRNA_pos.tsv

Rscript ../../scripts/count_AT_GC_gene_trait.R \
    -s DMS \
    -i /scratch/wangq/wsn/T2T_pars/result \
    -n Scer_spar


for i in syn nsy; do
    perl ../../scripts/count_position_gene.pl \
        --file ../../DMS_process/Scer_spar_DMS.gene_variation.process.yml \
        --origin data_SNPs_DMS_${i}.tsv \
        --output data_SNPs_DMS_${i}_pos.tsv

    Rscript ../../scripts/count_AT_GC_gene_trait.R \
        -s DMS \
        -i /scratch/wangq/wsn/T2T_pars/result \
        -n Scer_spar \
        -r ${i}
done

# calculate stem/loop length without SNP
cat data_SNPs_DMS_mRNA.tsv |
    perl -nl -a -F"\t" -e 'print qq{$F[12]};' |
    sort -u |
    perl -nl -a -F"\t" -e 'next if /gene/; print qq{$F[0]}; BEGIN{print qq{gene};}' \
    > mRNA.gene.list.tsv

for i in stem loop; do
    perl ../../scripts/count_structure_length_gene.pl \
        --file ../../DMS_process/Scer_spar_DMS.gene_variation.process.yml \
        --name mRNA.gene.list.tsv \
        --structure ${i} \
        --output ${i}_length_mRNA.tsv
done

```

### count codon gene

```shell
cd /scratch/wangq/wsn/T2T_pars/result/Scer_spar

perl ../../scripts/count_codon_gene.pl \
    --origin data_SNPs_DMS_syn.tsv \
    --output data_SNPs_DMS_syn_codon.tsv

Rscript ../../scripts/count_AT_GC_codon_chi.R \
    -s DMS \
    -i /scratch/wangq/wsn/T2T_pars/result \
    -n Scer_spar

```


### summary

```shell
cd /scratch/wangq/wsn/T2T_pars/result/Scer_spar




# calculate gama value
for i in mRNA nsy syn tRNA 4D; do
    python gama.py -i gama/DMS_${i}_stat_SNPs_freq_10.csv -n 790 -o gama/DMS_${i}_gama.txt
    cat DMS_${i}_gama.txt >> gama/figBC.txt
done

for i in mRNA nsy syn; do
    for type in stem loop;do 
        for GC in AT_GC GC_AT; do
            python gama.py -i gamaE/DMS_${i}_stat_${type}_${GC}_freq_10.csv \
                -n 790 -l ${type}_${GC} \
                -o gamaE/DMS_${i}_stat_${type}_${GC}_gama.txt
        done
    done

    cat gamaE/DMS_${i}_stat_${type}_${GC}_gama.txt >> gamaE/DMS_${i}_gama.txt
done

for i in mRNA nsy syn; do
    for type in stem loop;do 
        for GC in AT_GC GC_AT; do
            cat gamaE/DMS_${i}_stat_${type}_${GC}_gama.txt >> gamaE/DMS_${i}_gama.txt
        done
    done
done

for i in mRNA nsy syn; do
    for type in stem loop;do 
        for GC in AT_GC GC_AT; do
            for length in 1 2 3 4 5 6 7 8 9 10;do
                python gama.py -i figF/DMS_${i}_stat_${type}_${GC}_freq_10.csv \
                    -n 790 -l ${type}_${GC} \
                    -c ${i}_${length} \
                    -o gamaF/DMS_${i}_stat_${type}_${GC}_gama_${length}.txt
            done
        done
    done
done

for i in mRNA ; do
    for type in stem loop;do 
        for GC in AT_GC GC_AT; do
                tr < figF/DMS_${i}_stat_${type}_${GC}_freq_10.tsv "\t" ","  >> DMS_${i}_stat_${type}_${GC}_freq_10.csv
        done
    done
done



```