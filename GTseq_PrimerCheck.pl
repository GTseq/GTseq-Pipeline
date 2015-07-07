#!/usr/bin/perl
# GTseq_PrimerCheck.pl
# Checks a list of forward and reverse primers for possible problems with primer-dimer artifacts in multiplex PCR.
# Usage: provide a list of primers in the format Name\tFWD-Primer\tREV-Primer.

use strict; use warnings;

die "Provide a list of Primers\n" unless @ARGV == 1;

my %FWD_Search = ();
my %REV_Search = ();

open (FILE, "<$ARGV[0]") or die "Error opening $ARGV[0]\n";

while (<FILE>) {
	chomp;
	my @info = split "\t", $_;
	$FWD_Search{$info[0]} = substr $info[1], -9;
	$FWD_Search{$info[0]} =~ tr/ACGT/TGCA/;
	$FWD_Search{$info[0]} = reverse $FWD_Search{$info[0]};
	$REV_Search{$info[0]} = substr $info[2], -9;
	$REV_Search{$info[0]} =~ tr/ACGT/TGCA/;
	$REV_Search{$info[0]} = reverse $REV_Search{$info[0]};
	#print "$info[0],$FWD_Search{$info[0]},$REV_Search{$info[0]}\n";  #for testing...
	}
close FILE;

open (FILE2, "<$ARGV[0]") or die;

while (<FILE2>) {
	chomp;
	my @info2 = split "\t", $_;
	$info2[2] = substr $info2[2], 14;
	$info2[1] = substr $info2[1], 14;

	foreach my $loci (sort keys %FWD_Search) {
	if ($info2[2] =~ m/$FWD_Search{$loci}/) {print "FWD\t$loci\t$FWD_Search{$loci}\tREV\t$info2[0]\t$info2[2]\n";}
	#else {print "No Match\n"}
	if ($info2[1] =~ m/$REV_Search{$loci}/) {print "REV\t$loci\t$REV_Search{$loci}\tFWD\t$info2[0]\t$info2[1]\n";}
	#else {print "No Match\n"}
	}
}
close FILE2;
