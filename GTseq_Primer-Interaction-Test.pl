#!/usr/bin/perl
#GTseq_Primer-Interaction-Test.pl
#Test .hash file to count Primer Pair sequences in readthroughs
# Provide a tab delimited file of LocusName\tFWD-Primer sequence\tREV-Primer sequence\n and a .hash file

use strict; use warnings;

die "usage: provide <tab delimited txt file of locus\tfwd\trev seqs> and <*.hash>\n" unless @ARGV == 2;

my @assays = ();
my %fwd_seq = ();
my %fwd_RC = ();
my %rev_seq = ();
my %rev_RC = ();
my %PrimerComb = ();
my %PrimerCombFW = ();

#read in assay and allele information and push to arrays...

open(SEQ, "<$ARGV[0]") or die "error reading $ARGV[0]\n";

while (<SEQ>) {
	chomp;
	my @info = split(/\t/, $_);
	push @assays, $info[0];
	$fwd_seq{$info[0]} = $info[1];
	$fwd_RC{$info[0]} = reverse $fwd_seq{$info[0]};
	$fwd_RC{$info[0]} =~ tr/ACGT/TGCA/;
	$fwd_RC{$info[0]} = substr $fwd_RC{$info[0]}, 0, 10;
	$rev_seq{$info[0]} = $info[2];
	$rev_RC{$info[0]} = reverse $rev_seq{$info[0]};
	$rev_RC{$info[0]} =~ tr/ACGT/TGCA/;
	$rev_RC{$info[0]} = substr $rev_RC{$info[0]}, 0, 10;
	#print "$fwd_seq{$info[0]};$rev_seq{$info[0]};$fwd_RC{$info[0]};$rev_RC{$info[0]}\n";  #Testing...
		}
close SEQ;

foreach my $Locus (@assays) {

open (HASH, "<$ARGV[1]") or die "error reading $ARGV[1]\n";

	while (<HASH>) {
		my $hash = $_;
		my $R1seq = <HASH>;
		chomp ($hash);
		my @info = split(/;/, $hash);
		my $count = $info[2];
		if ($R1seq =~ m/$fwd_seq{$Locus}/) {
			foreach my $Revs (@assays) {
				if ($R1seq =~ m/$rev_RC{$Revs}/) {
					my $pair = "$Locus $Revs";
					if (exists $PrimerComb{$pair}) {$PrimerComb{$pair} = $PrimerComb{$pair} + $count;}
					else {$PrimerComb{$pair} = $count;}
					}
				}
			foreach my $FWDRC (@assays) {
				if ($R1seq =~ m/$fwd_RC{$FWDRC}/) {
					my $pair = "$Locus $FWDRC";
					if (exists $PrimerCombFW{$pair}) {$PrimerCombFW{$pair} = $PrimerCombFW{$pair} + $count;}
					else {$PrimerCombFW{$pair} = $count;}
					}
				}
		}
	}
}
close HASH;

#print headers:
print "Combination\tCounts\n";

foreach my $Combos (sort keys %PrimerComb) {
	print "$Combos\t$PrimerComb{$Combos}\n";
	}
print "******FWD Primer Combinations********\n";
foreach my $Combos2 (sort keys %PrimerCombFW) {
	print "$Combos2\t$PrimerCombFW{$Combos2}\n";
	}
