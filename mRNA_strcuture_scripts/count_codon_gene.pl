#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw();
use FindBin;
use YAML::Syck qw();

use Path::Tiny;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use List::Util qw/min/;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS

    perl ~/Scripts/pars/program/count_codon_gene.pl
    --origin data_SNPs_PARS_cds.update.csv
    --output data_SNPs_PARS_cds.update_codon.csv

=cut

Getopt::Long::GetOptions(
    'help|?'      => sub { Getopt::Long::HelpMessage(0) },
    'origin|or=s' => \my $origin,
    'output|o=s'  => \my $output,
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Process codon info...");

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#

open my $csv_fh, '<', $origin;

open my $out_fh, '>', $output;

while (<$csv_fh>) {
    chomp;
    s/"//g;
    my @snp = split /\t/, $_;
    
    if ( $snp[12] eq "gene" ) {

        splice @snp, 18, 0, "Codon_to";
        splice @snp, 18, 0, "Codon_pos";
        my $snp = join "\t", @snp;

        my $header = join "\t", @snp;
        print {$out_fh} $header, "\n";

    }
    else {

        my @codon;
        @codon = split "/",  $snp[17] if defined( $snp[17] );
        my @mutant;
        @mutant = split "->", $snp[9] if defined( $snp[9] );
        my $codon_to;
        my $codon_pos;

        if ( $snp[6] eq '+' ) {
            if ( $codon[0] =~ m/$mutant[0]/ ) {
                $codon_to  = "$codon[0]->$codon[1]";
                $codon_pos = index( $codon[0], $mutant[0] );
            }
            else {
                $codon_to  = "$codon[1]->$codon[0]";
                $codon_pos = index( $codon[1], $mutant[0] );
            }
        }
        else {
            $mutant[0] =~ tr/ATGC/TACG/;
            $mutant[1] =~ tr/ATGC/TACG/;
            if ( $codon[0] =~ m/$mutant[0]/ ) {
                $codon_to  = "$codon[0]->$codon[1]";
                $codon_pos = index( $codon[0], $mutant[0] );
            }
            else {
                $codon_to  = "$codon[1]->$codon[0]";
                $codon_pos = index( $codon[1], $mutant[0] );
            }
        }
        $codon_to = uc($codon_to);

        splice @snp, 18, 0, $codon_pos if defined( $snp[17] );
        splice @snp, 19, 0, $codon_to;
        my $snp = join "\t", @snp;

        print {$out_fh} $snp, "\n";

    }
}

close $csv_fh;
close $out_fh;

__END__
