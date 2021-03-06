#Metagenomic filtering


###To run:

1. RGAPepPipe [https://github.com/tracysmith/RGAPepPipe]. Requires fastq files, reference genome.

  * Produces *.realn.bam files that are used in the next script.

2. metagenFilter.pl. Requires bam file, reference, gtf of CDS regions, and gtf of regions to remove.

  * converts bam to fastq
  * runs fastq files through Kraken, keeping seqs that are assigned to the Mycobacterium genus (higher up than MTBC).
  * creates mpileup
  * remove indels, conserved and repetitive regions, and stran biased positions
  * calls script to run PoPoolation:
    * subsampling
    * pi, theta and Tajima's D statistic on 100kb windows, with 10kb steps
    * pi, theta stats on genes
  * ouputs:
    * bunch of subsamples and associated files
    * average values across subsamples of all values.


3. make_windows_table.pl Requirements: need to point it at sample overview file and subsample averages.

  * make a table with data for all samples, including metadata and stats based on windows analysis

4. make_gene_table.pl Requirements: need to point at sample overview, subsample averages and list of gene names

  * make a table with data for all samples, including metadata and stats for each gene. 