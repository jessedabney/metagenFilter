#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

my $output;
GetOptions ( "output=s" => \$output);

my @files = @ARGV; 
my %data;
my $sampnum = 0;

foreach my $file (@files) {
	open READ, $file or die "couldn't open file: $file, $?!\n";
	print STDERR "processing $file...\n";
	while (my $line = <READ>) {
		my ($gene, $snps, $cov, $value) = split /\t/, $line;
		if ($value =~ /na/){
			$value = -1;
		}
		$data{$gene}{"snps"} += $snps;
		$data{$gene}{"cov"} += $cov;
		$data{$gene}{"value"} += $value;
	}
	$sampnum++;
}
close READ;

open AVERAGE, ">${output}" or die "couldn't open file for output, $?\n";
foreach my $gene (sort keys %data) {
	print AVERAGE join ("\t",
		$gene, 
		$data{$gene}{"snps"} / $sampnum, 
		$data{$gene}{"cov"} / $sampnum);
	print AVERAGE "\t";
	printf(AVERAGE "%.9f", $data{$gene}{"value"} / $sampnum);
	print AVERAGE "\n";
}
close AVERAGE;

