mkdir -p /scratch/wangq/wsn/n128/gene_filter/Scer_n128_pair
cd /scratch/wangq/wsn/n128/gene_filter
ln -s /share/home/wangq/data/mrna-structure/sgd ../sgd

ls /share/home/wangq/data/mrna-structure/alignment/n128/Pairwise > n128.lst
echo -e "I\nII\nIII\nIV\nV\nVI\nVII\nVIII\nIX\nX\nXI\nXII\nXIII\nXIV\nXV\nXVI" > chromosome.lst
ln -s /share/home/wangq/data/mrna-structure/alignment/n128/Pairwise Pairwise

for strain in $(cat n128.lst); do
    for i in $(cat chromosome.lst); do
        gzip -dcf Pairwise/${strain}/mafSynNet/${i}.synNet.maf.gz |
            pigz >> Scer_n128_pair/${strain}.maf.gz
    done
done

cat n128.lst | parallel --no-run-if-empty --linebuffer -k -j 12 "
    fasops maf2fas Scer_n128_pair/{}.maf.gz -o Scer_n128_pair/{}.fas
    gzip Scer_n128_pair/{}.fas
    rm -f Scer_n128_pair/{}.maf.gz
"

mkdir -p Scer_n128_pair_refine
bsub -q mpi -n 48 bash Scer_n128_pair_refine/refine.sh

mkdir -p Scer_n128_pair_split
cat n128.lst | parallel --no-run-if-empty --linebuffer -k -j 12 "
    mkdir -p Scer_n128_pair_split/{}
    fasops split -r Scer_n128_pair_refine/{}.fas.gz -o Scer_n128_pair_split/{}
"

mkdir -p Scer_n128_vcf
cat n128.lst | parallel --no-run-if-empty --linebuffer -k -j 12 "
    ls Scer_n128_pair_split/{} > Scer_n128_pair_split/{}/fasta.lst
    sed -i '/fasta.lst/d' Scer_n128_pair_split/{}/fasta.lst
    sed -i 's/\.fas//g' Scer_n128_pair_split/{}/fasta.lst
"

for pair in $(cat n128.lst); do
    mkdir -p Scer_n128_vcf/${pair}
    for fas in $(cat Scer_n128_pair_split/${pair}/fasta.lst);do
        snp-sites -v Scer_n128_pair_split/${pair}/${fas}.fas -o Scer_n128_vcf/${pair}/${fas}.vcf
    done
done

cat n128.lst | parallel --no-run-if-empty --linebuffer -k -j 12 "
    ls Scer_n128_vcf/{} > Scer_n128_vcf/{}/vcf.lst
    sed -i '/vcf.lst/d' Scer_n128_vcf/{}/vcf.lst
    sed -i 's/\.vcf//g' Scer_n128_vcf/{}/vcf.lst
"

for pair in $(cat n128.lst); do
    mkdir -p Scer_n128_vcf_modify/${pair}
    for vcf in $(cat Scer_n128_vcf/${pair}/vcf.lst);do
        python vcf_format_modify.py Scer_n128_pair_split/${pair}/${vcf}.fas Scer_n128_vcf/${pair}/${vcf}.vcf Scer_n128_vcf_modify/${pair}/${vcf}.vcf
    done
done

mkdir Scer_n128_vcf_concat
bsub -q mpi -n 48 bash Scer_n128_vcf_concat/concat.sh
find Scer_n128_vcf_concat -mindepth 2 -maxdepth 2 -name '*_sample.sort.vcf.gz' | sort | grep -v "S288cvsspar" > Scer_seub.lst
find Scer_n128_vcf_concat -mindepth 2 -maxdepth 2 -name '*_sample.sort.vcf.gz' | sort | grep -v "S288cvsseub" > Scer_spar.lst

mkdir Scer_n128_vcf_merge
bcftools merge -l Scer_seub.lst --output-type z --output Scer_n128_vcf_merge/Scer_seub.vcf.gz
bcftools merge -l Scer_spar.lst --output-type z --output Scer_n128_vcf_merge/Scer_spar.vcf.gz

# outgroup seub
tabix -p vcf Scer_n128_vcf_merge/Scer_seub.vcf.gz
bcftools annotate \
    -a genes.non-overlapped.bed.gz \
    -c CHROM,FROM,TO,ID \
    -R genes.non-overlapped.bed.gz \
    -Ov -o Scer_n128_vcf_merge/Scer_seub_DMS.vcf \
    Scer_n128_vcf_merge/Scer_seub.vcf.gz 

# outgroup spar
tabix -p vcf Scer_n128_vcf_merge/Scer_spar.vcf.gz
bcftools annotate \
    -a genes.non-overlapped.bed.gz \
    -c CHROM,FROM,TO,ID \
    -R genes.non-overlapped.bed.gz \
    -Ov -o Scer_n128_vcf_merge/Scer_spar_DMS.vcf \
    Scer_n128_vcf_merge/Scer_spar.vcf.gz 

for i in $(cat n128.lst);do 
    python Scer_n128_pair_split/extract_intervals.py Scer_n128_pair_split/${i}
done

python Scer_n128_vcf_merge/process_merged_vcf.py Scer_n128_vcf_merge/Scer_spar_DMS.vcf Scer_n128_pair_split Scer_n128_vcf_merge/Scer_spar_DMS_processed.vcf

python Scer_n128_vcf_merge/process_merged_vcf.py Scer_n128_vcf_merge/Scer_seub_DMS.vcf Scer_n128_pair_split Scer_n128_vcf_merge/Scer_seub_DMS_processed.vcf


for i in spar seub; do
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%TGT]\n' Scer_n128_vcf_merge/Scer_${i}_DMS_processed.vcf > Scer_n128_vcf_merge/Scer_${i}_DMS.snp
    perl Scer_n128_vcf_merge/polarize_snp_vcf.pl Scer_n128_vcf_merge/Scer_${i}_DMS.snp > Scer_n128_vcf_merge/Scer_${i}_DMS.SNPs.tsv
done


mkdir -p /scratch/wangq/wsn/n128/vep
cd /scratch/wangq/wsn/n128/vep

cp ../gene_filter/Scer_n128_vcf_merge/Scer_spar_DMS.SNPs.tsv ../gene_filter/Scer_n128_vcf_merge/Scer_seub_DMS.SNPs.tsv .

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


