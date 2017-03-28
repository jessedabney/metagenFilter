#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

my $output;
GetOptions ( "output=s" => \$output);

my @files = @ARGV; 
my %data;
my $refName;
my $sampnum = 0;

foreach my $file (@files) {
	open READ, $file or die "couldn't open file: $file, $?!\n";
	print STDERR "processing $file...\n";
	while (my $line = <READ>) {
		my ($refID, $position, $snps, $cov, $value) = split /\t/, $line;
		if ($value =~ /na/){
			$value = -1;
		}
		$refName = $refID;
		$data{$position}{"snps"} += $snps;
		$data{$position}{"cov"} += $cov;
		$data{$position}{"value"} += $value;
	}
	$sampnum++;
}
close READ;

open AVERAGE, ">>${output}" or die "couldn't open file for output, $?\n";
foreach my $position (sort {$a <=> $b} keys %data) {
	print AVERAGE join ("\t",
		$refName,
		$position, 
		$data{$position}{"snps"} / $sampnum, 
		$data{$position}{"cov"} / $sampnum);
	print AVERAGE "\t";
	printf(AVERAGE "%.9f", $data{$position}{"value"} / $sampnum);
	print AVERAGE "\n";
}
close AVERAGE;

