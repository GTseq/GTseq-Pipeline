#!/usr/bin/perl
# GTseq_HashSeqs.pl adopted from Mike Miller's RAD sequencing pipeline for use in GTseq
# Usage: $ GTseq_HashSeqs.pl something.fastq > something.hash
# collapses unique reads into a single fasta file and counts the occurrences of each.

$file = $ARGV[0];
$lib = $ARGV[1];

open(FILE, "<$file")
	or die;

while (<FILE>) {
		
	$seq_line = <FILE>;
	chomp($seq_line);
	$tags{$seq_line}++;

	$_ = <FILE>;
	$_ = <FILE>;
}
close FILE;


$z = 1;
foreach $key (sort { $tags{$b} <=> $tags{$a} } keys %tags) {

	print ">$lib;$z;$tags{$key}\n";
	print "$key\n";

	$z++;
}


