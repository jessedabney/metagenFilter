#!/usr/bin/perl -w 
#make_window_table_wiscdata.pl
#Jesse Dabney
#06.26.17

#usage: ./make_window_table_wiscdata.pl
#ALERT: this script was designed for specific samples, so should be modified before running


use strict;

my %overview;
my @stats = ("pi", "theta", "td");
my $outFile = "window_table_wiscdata.txt";
open(OUT, ">", "/home/jdabney/projects/sputum/analysis/wisc_data/$outFile") or die "Couldn't open $outFile for printing: $!\n";

#read in overview file
my $infoFile = "/home/jdabney/projects/sputum/data/wisc_data_overview.txt";
open(INFO, "<", "$infoFile") or die "Couldn't open $infoFile for reading: $!\n";

my $header = <INFO>; #strip off header
while (my $line = <INFO>) {
	chomp $line;
	my($patient, $source, $capture, $sample) = split /\t/, $line;
	$overview{$sample}{"source"} = $source;
	$overview{$sample}{"capture"} = $capture;
	$overview{$sample}{"patient"} = $patient;
}

#find, open and read in each subsample average file.
my @windows;
my $n = 0;
open (MISS, ">>", "/home/jdabney/projects/sputum/analysis/wisc_data/missing_wisc_samples.txt") or die "Couldn't open missing data file for output: $!\n";
foreach my $sample (sort keys %overview) {
	my $dir = "/home/jdabney/projects/sputum/data/metgen_filtered/$sample";
	unless (-d $dir) {
		print MISS "$sample\n";
		next;
	} else {
		chdir "$dir" or die "Couldn't change into $sample directory: $!\n";
		my $file;
		foreach my $stat (@stats) {
			if ($sample =~ /^E/) {
				my $file = "${sample}_subAvg.${stat}";
			} elsif ($sample =~ /^sput/) {
				my $file = "${sample}_subAvg_windows_withNA.${stat}";
			} else { 
				die "weird name for $sample\n"; 
			}
			my %data;
			open (FILE, "<", "$file") or die "Couldn't open $file for reading: $!\n";
			print STDERR "working on $sample\n";
			while (my $line = <FILE>) {
				chomp $line;
				my ($win, $val) = (split /\t/, $line)[1,4];
				$data{$win} = $val;
				if (scalar @windows < 426) {
					push @windows, $win;
				}	
			}
			if ($n < 1) {
				my $header = "sample\tpatient\tsource\tcapture\tstat\t";
				print OUT "$header", join("\t", @windows), "\n";
				$n++;
			}
			print STDERR "printing $sample\n";
			print OUT "$sample\t$overview{$sample}{'patient'}\t$overview{$sample}{'source'}\t$overview{$sample}{'capture'}\t$stat";
			foreach my $win (sort {$a <=> $b} keys %data) {
				print OUT "\t$data{$win}";
			}
			print OUT "\n";
		}
	}
}
