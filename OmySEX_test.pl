#!/usr/bin/perl
# OmySEX_test.pl
# by Nate 
# Sex fish by GITseq

use strict; use warnings;

my @Files = `ls *.fastq`;
chomp ( @Files );

foreach my $samples (@Files){
	open (FILE, "<$samples") or die;
	my $reads = `wc -l $samples`;
	$reads =~ s/\s.*fastq//;
	$reads = ($reads/4) * 0.6;
	my $thresh = ($reads * 0.0012)/5;
	$samples =~ s/.fastq//;
	print "$samples,";
	my $counts = 0;
	my $sex = "NA";
	while (<FILE>) {
		chomp;
		if(/ATGTGTTCATATGCCAG/){$counts++}
	}
if(($counts > $thresh) && ($reads > 5000)) {$sex = "M"}
elsif (($counts < $thresh) && ($reads > 5000)) {$sex = "F"}
close FILE;
print "$sex\n";
}


