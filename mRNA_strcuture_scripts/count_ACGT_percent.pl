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
use App::Rangeops;

use Statistics::ChisqIndep;
use POSIX;

#use List::Util;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS

    perl ~/Scripts/pars/program/count_ACGT_percent.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --varfold ~/data/mrna-structure/process/$NAME.gene_variation.fold_class.tsv --output $NAME.gene_variation.fold_class.csv

=cut

Getopt::Long::GetOptions(
    'help|?'       => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'     => \my $file,
    'varfold|vf=s' => \my $varfold,
    'output|o=s'   => \my $output,
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Process folding info...");

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
print "Load $file\n";
my $gene_info_of = YAML::Syck::LoadFile($file);

open my $tsv_fh, '<', $varfold;

open my $out_fh, '>', $output;

while (<$tsv_fh>) {
    chomp;
    s/	/,/g;
    my @snp = split /,/, $_;

    if ( $snp[0] eq "gene" ) {
        $snp[8]  = "stem_A_num";
        $snp[9] = "stem_A_per";
        $snp[10] = "stem_C_num";
        $snp[11] = "stem_C_per";
        $snp[12] = "stem_G_num";
        $snp[13] = "stem_G_per";
        $snp[14] = "stem_U_num";
        $snp[15] = "stem_U_per";
        $snp[16] = "loop_A_num";
        $snp[17] = "loop_A_per";
        $snp[18] = "loop_C_num";
        $snp[19] = "loop_C_per";
        $snp[20] = "loop_G_num";
        $snp[21] = "loop_G_per";
        $snp[22] = "loop_U_num";
        $snp[23] = "loop_U_per";
        $snp[24] = "A_num";
        $snp[25] = "A_per";
        $snp[26] = "C_num";
        $snp[27] = "C_per";
        $snp[28] = "G_num";
        $snp[29] = "G_per";
        $snp[30] = "U_num";
        $snp[31] = "U_per";
        $snp[32] = "stem_AU_num";
        $snp[33] = "stem_CG_num";
        $snp[34] = "loop_AU_num";
        $snp[35] = "loop_CG_num";
        $snp[36] = "AU_num";
        $snp[37] = "CG_num";
        $snp[38] = "stem_CG_content";
        $snp[39] = "loop_CG_content";
        $snp[40] = "CG_content";
        $snp[41] = "X2";
        $snp[42] = "P";

        my $snp = join "\t", @snp;
        print {$out_fh} $snp, "\n";

    }
    else {

        if ( grep ( $snp[0], keys %{$gene_info_of} ) )
        {    # $snp[0] represents snp in which gene

            my $info = $gene_info_of->{ $snp[0] };

            my $seq = $info->{seq};

            my $stem = AlignDB::IntSpan->new;

            $stem->add( $info->{fold_left} );
            $stem->add( $info->{fold_right} );

            my $loop = AlignDB::IntSpan->new( $info->{fold_dot} );

            my @stem_position = $stem->sets();
            my @loop_position = $loop->sets();

            my @stem_seq = ();
            foreach my $ranges1 (@stem_position) {
                my $max1          = $ranges1->max();
                my $min1          = $ranges1->min();
                my $length1       = $max1 - $min1 + 1;
                my $stem_seq_part = substr( $seq, $min1 - 1, $length1 );
                push( @stem_seq, $stem_seq_part );
            }
            my @loop_seq = ();
            foreach my $ranges2 (@loop_position) {
                my $max2          = $ranges2->max();
                my $min2          = $ranges2->min();
                my $length2       = $max2 - $min2 + 1;
                my $loop_seq_part = substr( $seq, $min2 - 1, $length2 );
                push( @loop_seq, $loop_seq_part );
            }

            my $stem_seq = join "", @stem_seq;
            my $loop_seq = join "", @loop_seq;

            my $stem_A_num      = ( $stem_seq =~ s/A/A/g );
            my $stem_A_pre      = $stem_A_num / length($stem_seq);
            my $stem_C_num      = ( $stem_seq =~ s/C/C/g );
            my $stem_C_pre      = $stem_C_num / length($stem_seq);
            my $stem_G_num      = ( $stem_seq =~ s/G/G/g );
            my $stem_G_pre      = $stem_G_num / length($stem_seq);
            my $stem_U_num      = ( $stem_seq =~ s/U/U/g );
            my $stem_U_pre      = $stem_U_num / length($stem_seq);
            my $loop_A_num      = ( $loop_seq =~ s/A/A/g );
            my $loop_A_pre      = $loop_A_num / length($loop_seq);
            my $loop_C_num      = ( $loop_seq =~ s/C/C/g );
            my $loop_C_pre      = $loop_C_num / length($loop_seq);
            my $loop_G_num      = ( $loop_seq =~ s/G/G/g );
            my $loop_G_pre      = $loop_G_num / length($loop_seq);
            my $loop_U_num      = ( $loop_seq =~ s/U/U/g );
            my $loop_U_pre      = $loop_U_num / length($loop_seq);
            my $A_num           = ( $seq =~ s/A/A/g );
            my $A_pre           = $A_num / length($seq);
            my $C_num           = ( $seq =~ s/C/C/g );
            my $C_pre           = $C_num / length($seq);
            my $G_num           = ( $seq =~ s/G/G/g );
            my $G_pre           = $G_num / length($seq);
            my $U_num           = ( $seq =~ s/U/U/g );
            my $U_pre           = $U_num / length($seq);
            my $stem_AU_num     = $stem_A_num + $stem_U_num;
            my $stem_CG_num     = $stem_C_num + $stem_G_num;
            my $loop_AU_num     = $loop_A_num + $loop_U_num;
            my $loop_CG_num     = $loop_C_num + $loop_G_num;
            my $AU_num          = $A_num + $U_num;
            my $CG_num          = $C_num + $G_num;
            my $stem_CG_content = $stem_CG_num / length($stem_seq);
            my $loop_CG_content = $loop_CG_num / length($loop_seq);
            my $CG_content      = $CG_num / length($seq);

            my $stem_AU_num_th =
              ( $stem_AU_num + $loop_AU_num ) * ( $stem_AU_num + $stem_CG_num )
              / ( $stem_AU_num + $stem_CG_num + $loop_AU_num + $loop_CG_num );
            my $loop_AU_num_th =
              ( $stem_AU_num + $loop_AU_num ) * ( $loop_AU_num + $loop_CG_num )
              / ( $stem_AU_num + $stem_CG_num + $loop_AU_num + $loop_CG_num );
            my $stem_CG_num_th =
              ( $stem_CG_num + $loop_CG_num ) * ( $stem_AU_num + $stem_CG_num )
              / ( $stem_AU_num + $stem_CG_num + $loop_AU_num + $loop_CG_num );
            my $loop_CG_num_th =
              ( $stem_CG_num + $loop_CG_num ) * ( $loop_AU_num + $loop_CG_num )
              / ( $stem_AU_num + $stem_CG_num + $loop_AU_num + $loop_CG_num );
            my $X2 =
              ( ( $stem_AU_num - $stem_AU_num_th ) *
                  ( $stem_AU_num - $stem_AU_num_th ) /
                  $stem_AU_num_th ) +
              ( ( $stem_CG_num - $stem_CG_num_th ) *
                  ( $stem_CG_num - $stem_CG_num_th ) /
                  $stem_CG_num_th ) +
              ( ( $loop_AU_num - $loop_AU_num_th ) *
                  ( $loop_AU_num - $loop_AU_num_th ) /
                  $loop_AU_num_th ) +
              ( ( $loop_CG_num - $loop_CG_num_th ) *
                  ( $loop_CG_num - $loop_CG_num_th ) /
                  $loop_CG_num_th );

            my $obs = [ [ $stem_AU_num, $stem_CG_num ],
                [ $loop_AU_num, $loop_CG_num ] ];
            my $chi = Statistics::ChisqIndep->new;
            $chi->load_data($obs);
            $chi->print_summary();
            my $P = ${$chi}{'p_value'};

            my @sum = (
                $stem_A_num,      $stem_A_pre,      $stem_C_num,
                $stem_C_pre,      $stem_G_num,      $stem_G_pre,
                $stem_U_num,      $stem_U_pre,      $loop_A_num,
                $loop_A_pre,      $loop_C_num,      $loop_C_pre,
                $loop_G_num,      $loop_G_pre,      $loop_U_num,
                $loop_U_pre,      $A_num,           $A_pre,
                $C_num,           $C_pre,           $G_num,
                $G_pre,           $U_num,           $U_pre,
                $stem_AU_num,     $stem_CG_num,     $loop_AU_num,
                $loop_CG_num,     $AU_num,          $CG_num,
                $stem_CG_content, $loop_CG_content, $CG_content,
                $X2,              $P
            );
            push( @snp, @sum );

        }

        my $snp_processed = join "\t", @snp;
        print {$out_fh} $snp_processed, "\n";
    }

}

close $tsv_fh;
close $out_fh;

__END__
