#!/usr/bin/perl -w 
#make table for gene data
#Jesse Dabney
#01.09.16

#ALERT: this script is set up for specific samples, so should be modified before running


use strict;
use Sort::Naturally;

my %overview;
my @genenames;
my @bovisgenes;
my @stats = ("pi", "theta");
#my $outFile = "gene_table_brownsamples_withNA.txt";
my $outFile2 = "gene_table_wisc_data.txt";

#read in overview file
my $infoFile = "/home/jdabney/projects/sputum/data/wisc_data_overview.txt";
open(INFO, "<", "$infoFile") or die "Couldn't open $infoFile for reading: $!\n";
print STDERR "getting overview data...\n";
#my $header = <INFO>; #strip off header
while (my $line = <INFO>) {
	chomp $line;
	#my($patient, $score, $source, $capture, $sample) = (split /\t/, $line)[0,1,3,4,5];
	my ($patient, $source, $capture, $sample) = split /\t/, $line;
	$overview{$sample}{"source"} = $source;
	#$overview{$sample}{"score"} = $score;
	$overview{$sample}{"capture"} = $capture;
	$overview{$sample}{"patient"} = $patient;
}
print STDERR "getting gene names...\n";

#get list of Mtb genes
# open(TBGENES, "<", "/home/jdabney/projects/sputum/data/mtb_gene_names.txt") or die "couldn't open gene name list: $!\n";
# while (my $line = <TBGENES>) {
# 	chomp $line;
# 	push @genenames, $line;
# }
# my @mtb_genes = nsort @genenames; 

#get list of mbovis genes
open(MBGENES, "<", "/home/jdabney/projects/sputum/data/mbovis_genenames.txt") or die "couldn't open mbovis gene list: $!\n";
while (my $line = <MBGENES>) {
	chomp $line;
	push @bovisgenes, $line;
}
my @mbovis_genes = nsort @bovisgenes; 


#open(OUT, ">>", "/home/jdabney/projects/sputum/analysis/genes/$outFile") or die "Couldn't open $outFile for printing: $!\n";
open (OUT2, ">", "/home/jdabney/projects/sputum/analysis/wisc_data/$outFile2");

#find, open and read in each subsample average file.
print STDERR "processing subsample average files...\n";
my $n = 0;
open (MISS, ">>", "/home/jdabney/projects/sputum/analysis/genes/missing_wisc_samples_gene.txt") or die "Couldn't open missing data file for output: $!\n";
foreach my $sample (sort keys %overview) {
	my $dir = "/home/jdabney/projects/sputum/data/metgen_filtered/$sample";
	unless (-d $dir) {
		print MISS "$sample\n";
		next;
	} else {
		print STDERR "processing ${sample}...\n";
		chdir "$dir" or die "Couldn't change into $sample directory: $!\n";
		foreach my $stat (@stats) {
			my $file;
			if ($sample =~ /^E/) {
				$file = "${sample}_subAvg_genes.${stat}";
			} elsif ($sample =~ /^sput/) {
				$file = "${sample}_subAvg_genes_withNA.${stat}";
			} else {
				die "weird name in $sample\n";
			}
			my %data;
			open (FILE, "<", "$file") or die "Couldn't open $file for reading: $!\n";
			while (my $line = <FILE>) {
				chomp $line;
				my ($gene, $val) = (split /\t/, $line)[0,3];
				$data{$gene} = $val;
			}
			if ($n < 1) {
				my $header = "sample\tpatient\tsource\tcapture\tstat\t";
				#print OUT "$header", join("\t", @mtb_genes), "\n";
				print OUT2 "$header", join("\t", @mbovis_genes), "\n";	
				$n++;
			}
			# if ($sample =~ /^E/) {
			# 	print OUT join ("\t", 
			# 				$sample,
			# 				$overview{$sample}{'score'},
			# 				$overview{$sample}{'patient'},
			# 				$overview{$sample}{'source'},
			# 				$overview{$sample}{'capture'},
			# 				$stat);
			# 	foreach my $gene (nsort keys %data) {
			# 		print OUT "\t$data{$gene}";
			# 	}
			# 	print OUT "\n";
			# } elsif ($sample =~ /^sput/) {
			print OUT2 join ("\t", 
						$sample,
						#$overview{$sample}{'score'},
						$overview{$sample}{'patient'},
						$overview{$sample}{'source'},
						$overview{$sample}{'capture'},
						$stat);
			foreach my $gene (nsort keys %data) {
				print OUT2 "\t$data{$gene}";
			}
			print OUT2 "\n";
		}
	}
}
