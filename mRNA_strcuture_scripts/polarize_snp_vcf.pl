#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(uniq);

## usage: perl Scer_n790_vcf_merge/polarize_snp_vcf.pl Scer_n790_vcf_merge/test_snp.tsv > Scer_n790_vcf_merge/test.SNPs.tsv

my $input_file = $ARGV[0];

open my $fh, '<', $input_file or die "Could not open '$input_file': $!\n";

while (my $line = <$fh>) {
    chomp $line;
    next if $line =~ /^\s*$/; # ignore space

    my @cols = split(/\s+/, $line);
    
    # 0:CHR, 1:POS, 2:ID, 3:REF, 4:ALT, 5:Outgroup, 6...:Samples
    my ($chr, $pos, $gene_id, $ref, $alt, $outgroup) = @cols[0..5];

    # 1 if outgroup base is . ignore this row
    next if $outgroup eq '.';
    
    # 2 get all alleles
    my @samples_raw = @cols[6..$#cols];
    
    ## filter snps with sample missing rate >= 1%
    my $missing_count = 0;
    for my $s (@samples_raw) {
        $missing_count++ if $s eq '.';
    }
    my $total_samples = scalar(@samples_raw);
    my $missing_rate = $missing_count / $total_samples;

    next if $missing_rate >= 0.01;

    # remove missing samples
    my @valid_samples = grep { $_ ne '.' } @samples_raw;
    
    my @all_observed_alleles = grep { $_ ne '.' } ($ref, $alt, $outgroup, @valid_samples);
    my @unique_alleles = uniq(@all_observed_alleles);
    
    my @in_group = map { uc $_ } ($ref, @valid_samples);
    my @internal_alleles = uniq(@in_group);
    
    # base type of ref and sample just 1, ignore(not snp)
    next if scalar(@internal_alleles) < 2;
    
    # base type >= 3, ignore this row
    next if scalar(@unique_alleles) >= 3;
    
    # 3 detemine mutant_to
    my $derived_base = "";
    for my $a (@unique_alleles) {
        if ($a ne $outgroup) {
            $derived_base = $a;
            last;
        }
    }
    my $mutant_to = "$outgroup->$derived_base";
    
    # 4 get all output base(ref + valid_sample)
    my @output_bases_list = ($ref, @valid_samples);
    my $base_string = join("", @output_bases_list);
    
    # 5 calculate freq
    my $freq = 0;
    for my $b (@output_bases_list) {
        $freq++ if $b eq $derived_base;
    }
    
    next if $freq == 0;
    
    # 6 output
    my $location = "$chr:$pos";

    # chr:pos/ref/alt/mutant_to/freq/occured/gene_id
    printf "%s\t%s\t%s\t%s\t%d\t%s\t%s\n",
           $location, $ref, $alt, $mutant_to, $freq, $base_string, $gene_id;
}

close $fh;
