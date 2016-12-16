#!/usr/bin/perl -w

#metagenomic filtering stats


use strict;
use File::Basename;

my $progDir = "/mnt/PepPop_export/PepPrograms";
my $krakenDir = "/opt/PepPrograms";
my $krakenDB = "/home/yudong/DB";
my $dataDir = "/mnt/PepPop_export/data";
my $sam = "${progDir}/samtools-1.2/samtools";

my @files = @ARGV; #list of bam files
open (STATS, ">>", "kraken_stats.txt") or die "couldn't open file for stats output: $?\n";
print STATS "sample\ttotalSeqs\tnumberClassified\ttotalMTBC\tfinalMTBC\n";

foreach my $file (@files) {

	#start working on bam file, convert to fastq
	my $name = basename($file); 
	#print STDERR "$name\n";
	my ($sample, $temp1, $temp2) = split /\./, $name;
	my $fqOut = "${sample}.realn.fastq";
	my $filteredBam = "${sample}_classSeqs.bam";
	my $outVCF = "${sample}.vcf";
	my $mpileup1 = "${sample}_classSeqs.mpileup";
	my $indelGTF = "${sample}_indelRegs.gtf";
	my $mpileup2 = "${sample}_classSeqs_noIndel.mpileup";
	my $mpileup3 = "${sample}_classSeqs_noIndel_remRegs.mpileup";
	my $mpileup4 = "${sample}_classSeqs_noIndel_remRegs_noSB.mpileup";

	print STDERR "Processing $sample [bam2fq]...\n";
	system "$sam bam2fq $file > $fqOut";
	#system "gzip $fqOut"; #can compress here, use other kraken option

#call kraken on fastqs
	
	#use below if inputting compressed fq
	#system "${progDir}/kraken/kraken -preload --threads 8 --fastq-input --gzip-compressed --db ${progDir}/kraken/DB/ ${fqOut}.gz > ${sample}.kraken";
	print STDERR "Processing $sample [kraken]...\n";
	system "${krakenDir}/kraken/kraken --threads 8 --fastq-input --db $krakenDB ${fqOut} > ${sample}.kraken";
#translate kraken file
	print STDERR "Processing $sample [kraken translate]...\n";
	system "${krakenDir}/kraken/kraken-translate --db $krakenDB ${sample}.kraken > ${sample}.kraken.labels";

#work on .kraken files
	open (IN, "<", "${sample}.kraken") or die "couldn't open $file: $?\n";
	print STDERR "Processing $sample [kraken output]...\n";
	my $numClass = 0;
	my $numUnclass = 0;
	# my $lengthsClass = 0;
	# my $lengthsUnclass = 0;
	# my $totalLength = 0;

	while (my $line = <IN>) {
		chomp $line;
		my ($class, $name, $tax, $length, $taxChain) = split "\t", $line;
		if ($class eq "C") {
			$numClass++;
			#$lengthsClass += $length;
		} elsif ($class eq "U") {
			$numUnclass++;
			#$lengthsUnclass += $length;
		} else {
			warn "found an unexpected value in 1st column: ${sample}.kraken";
			next;
		}
		#$totalLength += $length;
	}
	my $totalSeqs = $numClass + $numUnclass;
	# my $fracClass = $numClass / $totalSeqs;
	# my $avgLength = $totalLength / $totalSeqs;
	# my $avgLengthC = $lengthsClass / $numClass;
	# my $avgLengthU = $lengthsUnclass / $numUnclass || "no seqs";
	#print STDOUT "Total sequences: $totalSeqs\nFraction classified: $fracClass\nAverage Length total: $avgLength\nAverage Length C: $avgLengthC\nAverage Length U: $avgLengthU\n";
	close IN;
	
#now work on output from kraken-translate
	print STDERR "Processing $sample [kraken labels]...\n";
	open (LABELS, "<", "${sample}.kraken.labels") or die "couldn't open labels file: $?\n";
	#open (ALLNAMES, ">>", "${sample}_MTBC_allSeqs.txt") or die "couldn't open file for all the seqs: $?\n";
	my %seqNames;
	my $totalMTBC = 0;
	while (my $line = <LABELS>) {
		chomp $line;
		my ($seq, $taxList) = split /\t/, $line;
	#now filter out non MTBC hits
		my @taxa = split /;/, $taxList;
		if (scalar @taxa < 9) {
			next;
		} elsif ($taxa[7] ne "Mycobacterium") {
			next;
		} else {
		#	print ALLNAMES "$seqName\n";
		 	$totalMTBC++;
		# }
	#now count the number if times each sequence appears, removing the trailing /[12]
		my $seqName = (split /\//, $seq)[0];
		$seqNames{$seqName}++;
		}
	}
#now filter against sequences where only 1 of the pair was assigned
	open (NAMES, ">>", "${sample}_MTBC_seqNames.txt") or die "couldn't open output for sequence names: $?\n";
	my $numMTBC = 0;
	foreach my $seqName (keys %seqNames) {
		if ($seqNames{$seqName} != 2) {
			next;
		} else {
			print NAMES "$seqName\n";
			$numMTBC++;
		}
	}

#now print some stats
	my $finalMTBC = 2 * $numMTBC;
	#print STDERR "Single end seqs: $diff\n";
	print STATS "$sample\t$totalSeqs\t$numClass\t$totalMTBC\t$finalMTBC\n";

#now filter the original bam file by sequence names
	print STDERR "Processing $sample [Picard FilterSamReads]...\n";
	system "java -Xmx8g -jar ${progDir}/picard-tools-1.138/picard.jar FilterSamReads INPUT=${sample}.realn.bam OUTPUT=$filteredBam READ_LIST_FILE=${sample}_MTBC_seqNames.txt FILTER=includeReadList";

#now convert filtered bam to mpileup and do some further filtering on the mpileup
	print STDERR "Processing $sample [mpileup generation and filtering]...\n";
	#open (INBAM, "<", "${sample}_MTBCseqs.bam") or die "couldn't open BAM to generate mpileup: $?\n";
	system "$sam mpileup -B -Q 20 -f ${dataDir}/mtuberculosis/MtbNCBIH37Rv.fa $filteredBam > $mpileup1";
	system "perl ${progDir}/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --input $mpileup1 --output $indelGTF";
	system "perl ${progDir}/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input $mpileup1 --gtf $indelGTF --output $mpileup2";
	system "perl ${progDir}/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input $mpileup2 --gtf ${dataDir}/mtuberculosis/150423_removeRegions.gtf --output $mpileup3";
#remove all but the final mpileup (.noIndel.remRegs.mpileup)
	unlink "$mpileup1";
	unlink "$mpileup2";
	unlink "$fqOut";
	unlink "${sample}.kraken.label";

#strand bias filter
	print STDERR "Processing $sample [generateing VCF]...\n";
	system "$sam mpileup -B -Q 20 -f ${dataDir}/mtuberculosis/MtbNCBIH37Rv.fa -uv -t DP,DP4,SP $filteredBam > $outVCF";
	print STDERR "Processing $sample [strand bias filter]...\n";
	open (VCF, "<", "$outVCF") or die "couldn't open $outVCF: $?\n";
	#open (SBOUT, ">", "${sample}_SB_positions.txt") or die "couldn't open output file for strand bias positions: $?\n";
	my %SBpos;
	while (my $line = <VCF>) {
    	if ($line =~ /^#/) {
			next;
    	}
    	my ($pos, $field) = (split /\t/, $line)[1,9];
    	my $val = (split /:/, $field)[2];
    	if ($val >= 13) {
			$SBpos{$pos}=$val;
    	} else {
    		next;
    	}
    }
    open (MPILEUP, "<", "$mpileup3") or die "couldn't open mpileup for strand bias filtering: $?\n";
    open (MPOUT, ">>", "$mpileup4") or die "couldn't open file for strand bias output: $?\n";
    while (my $line = <MPILEUP>) {
    	my $position = (split /\t/, $line)[1];
    	if (-e $SBpos{$position}) {
    		next;
    	} else {
    		print MPOUT "$line"
    	}
    }
#now call automatePoolSeqDiversityStats.py script
	system "~/projects/sputum/scripts/metagenFilter/automatePoolSeqDiversityStats.py -p ${sample} -g ~/projects/sputum/scripts/metagenFilter/tbdb_gichrom.gtf"

}



