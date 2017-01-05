#!/usr/bin/perl -w
#mpileup2pop_nosub.pl 
#Jesse Dabney
#01.04.16

#This script runs the final filtered mpileup from the metagenFilter script through popoolation without creating subsamples.

use strict;
#use Cwd;

#my $samplefile = "home/jdabney/projects/sputum/samplenames.txt";
my $sample = shift @ARGV;
#my $pwd = cwd();

#open (SAMPLES, "<", "$samplefile") or die "Couldn't open $samplefile: $?\n";

#while (my $sample = <SAMPLES>) {
#	chomp $sample;
	chdir "$sample" or die "Couldn't change dir to $sample: $?\n";
	if ($sample =~ /^E/) {
		print STDERR "Processing ${sample} as M.tb...\n";
		system "~/projects/sputum/scripts/metagenFilter/automatePoolSeqDiversityStats_noSubsample.py -p ${sample} -g ~/projects/sputum/scripts/metagenFilter/tbdb_gichrom.gtf";
	} elsif ($sample =~ /^s/) {
		print STDERR "Prcessing ${sample} as M.bovis...\n";
		system "~/projects/sputum/scripts/metagenFilter/automatePoolSeqDiversityStats_noSubsample.py -p ${sample} -g /mnt/PepPop_export/data/mbovis/mbovis_AF2122_transcripts.gtf";
	} else {
		die "Unexpected sample name: $sample...\n";
	}
	#chdir "$pwd";
#}