#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my $input_file;
my $output_file;

GetOptions(
    'i|input=s'  => \$input_file,
    'o|output=s' => \$output_file,
) or die "usage: $0 [-i input.txt] [-o output.txt] \n";

my $in_fh;
if ($input_file) {
    open($in_fh, '<', $input_file) or die "can't open input file: '$input_file': $!";
} else {
    $in_fh = \*STDIN;
}

my $out_fh;
if ($output_file) {
    open($out_fh, '>', $output_file) or die "can't create output file: '$output_file': $!";
} else {
    $out_fh = \*STDOUT;
}


## ---- main ----
my ($prev_chr, $prev_start, $prev_end);

sub print_range {
    my ($chr, $start, $end) = @_;
    print $out_fh "$chr:$start-$end\n";
}

while (<$in_fh>) {
    chomp;
    next if /^\s*$/;

    my ($chr, $range) = split /:/, $_;

    # single pos
    if ($range !~ /-/) {
        # directly output
        if (defined $prev_chr) {
            print_range($prev_chr, $prev_start, $prev_end);
            undef $prev_chr;
        }
        print $out_fh "$chr:$range\n";
        next;
    }

    # range
    my ($start, $end) = split /-/, $range;

    if (!defined $prev_chr) {
        ($prev_chr, $prev_start, $prev_end) = ($chr, $start, $end);
        next;
    }

    if ($chr eq $prev_chr && ($start - $prev_end) <= 30) {
        $prev_end = $end if $end > $prev_end;
    } else {
        print_range($prev_chr, $prev_start, $prev_end);
        ($prev_chr, $prev_start, $prev_end) = ($chr, $start, $end);
    }
}

print_range($prev_chr, $prev_start, $prev_end) if defined $prev_chr;

close $in_fh;
close $out_fh;
