#!/usr/bin/perl -w

use strict;
use File::Basename;


my $file = shift @ARGV;
print STDERR "Processing $file...\n";
my $name = basename($file);
my $sample = (split /_/, $name)[0];

open (VCF, "<", "$file") or die "couldn't open vcf $file: $?\n";
open (OUT, ">", "${sample}_SBpositions.txt") or die "couldn't open output file: $?\n";

my %counts;

while (my $line = <VCF>) {
    if ($line =~ /^#/) {
	next;
    }
    my ($pos, $field) = (split /\t/, $line)[1,9];
    my $val = (split /:/, $field)[2];
    if ($val >= 13) {
	print OUT "$pos\t$val\n";
    }
}
