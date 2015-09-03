#!/usr/bin/perl
#GTseq_Primer-Interaction-Test.pl
#Test .hash file to count Primer Pair sequences in readthroughs
# Provide a tab delimited file of LocusName\tFWD-Primer sequence\tREV-Primer sequence\n and a .hash file

use strict; use warnings;

die "usage: provide <tab delimited txt file of locus\tfwd\trev seqs> and <*.hash>\n" unless @ARGV == 2;

my @assays = ();
my %fwd_seq = ();
my %fwd_RC = ();
my %fwd_start = ();
my %rev_seq = ();
my %rev_RC = ();
my %rev_start = ();
my %rev_RC_Full = ();
my %OT_PrimerComb = ();
my %DC_Artifact = ();
my %PrimerComb = ();
my %PrimerCombFW = ();
my %PrimerCombRV = ();
my %PrimerCombRV_RV = ();

#read in assay and allele information and push to arrays...

open(SEQ, "<$ARGV[0]") or die "error reading $ARGV[0]\n";

while (<SEQ>) {
	chomp;
	my @info = split(/\t/, $_);
	push @assays, $info[0];
	$fwd_seq{$info[0]} = $info[1];
	$fwd_start{$info[0]} = substr $fwd_seq{$info[0]}, 0, 10;
	$fwd_RC{$info[0]} = reverse $fwd_seq{$info[0]};
	$fwd_RC{$info[0]} =~ tr/ACGT/TGCA/;
	$fwd_RC{$info[0]} = substr $fwd_RC{$info[0]}, -10;
	$rev_seq{$info[0]} = $info[2];
	$rev_start{$info[0]} = substr $rev_seq{$info[0]}, 0, 10;
	$rev_RC{$info[0]} = reverse $rev_seq{$info[0]};
	$rev_RC{$info[0]} =~ tr/ACGT/TGCA/;
	$rev_RC_Full{$info[0]} = $rev_RC{$info[0]};
	$rev_RC{$info[0]} = substr $rev_RC{$info[0]}, -10;
	#print "$info[0];$fwd_seq{$info[0]};$rev_seq{$info[0]};$fwd_RC{$info[0]};$rev_RC{$info[0]};$fwd_start{$info[0]};$rev_start{$info[0]};$rev_RC_Full{$info[0]}\n";  #Testing...
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
				if (($R1seq =~ m/$rev_RC_Full{$Revs}/) && ($Locus !~ m/$Revs/)) {
					my $pair = "$Locus $Revs";
					if (exists $DC_Artifact{$pair}) {$DC_Artifact{$pair} = $DC_Artifact{$pair} + $count;}
					else {$DC_Artifact{$pair} = $count;}
					}
				elsif (($R1seq =~ m/$rev_RC{$Revs}/) && ($Locus !~ m/$Revs/)) {
					my $pair = "$Locus $Revs";
					if (exists $PrimerComb{$pair}) {$PrimerComb{$pair} = $PrimerComb{$pair} + $count;}
					else {$PrimerComb{$pair} = $count;}
					}
				elsif (($R1seq =~ m/$rev_RC{$Revs}/) && ($Locus =~ m/$Revs/)) {
					my $pair = "$Locus $Revs";
					if (exists $OT_PrimerComb{$pair}) {$OT_PrimerComb{$pair} = $OT_PrimerComb{$pair} + $count;}
					else {$OT_PrimerComb{$pair} = $count;}
					}
				}
			foreach my $FWDRC (@assays) {
				if (($R1seq =~ m/$fwd_RC{$FWDRC}/) && ($Locus !~ m/$FWDRC/)) {
					my $pair = "$Locus $FWDRC";
					if (exists $PrimerCombFW{$pair}) {$PrimerCombFW{$pair} = $PrimerCombFW{$pair} + $count;}
					else {$PrimerCombFW{$pair} = $count;}
					}
				}
			}
		elsif ($R1seq =~ m/$rev_RC_Full{$Locus}/) {
			foreach my $FWD_starts (@assays) {
				if (($R1seq =~ m/$fwd_start{$FWD_starts}/) && ($Locus !~ m/$FWD_starts/)) {
					my $pair = "$Locus $FWD_starts";
					if (exists $PrimerCombRV{$pair}) {$PrimerCombRV{$pair} = $PrimerCombRV{$pair} + $count;}
					else {$PrimerCombRV{$pair} = $count;}
					}
				}
			foreach my $revs2 (@assays) {
				if (($R1seq =~ m/$rev_start{$revs2}/) && ($Locus !~ m/$revs2/)) {
					my $pair = "$Locus $revs2";
					if (exists $PrimerCombRV_RV{$pair}) {$PrimerCombRV_RV{$pair} = $PrimerCombRV_RV{$pair} + $count;}
					else {$PrimerCombRV_RV{$pair} = $count;}
					}
				}
			}
	}
}

#print headers:
print "***Proper On-Target Primer combinations\nCombination\tCounts\n";
foreach my $Combos1 (sort keys %OT_PrimerComb) {
	print "$Combos1\t$OT_PrimerComb{$Combos1}\n";
	}
print "***Double-Complement Primer-Artifacts***\n";
foreach my $CombosDC (sort keys %DC_Artifact) {
	print "$CombosDC\t$DC_Artifact{$CombosDC}\n";
	}
print "****FWD-REV Primer mis-primes****\n";
foreach my $Combos (sort keys %PrimerComb) {
	print "$Combos\t$PrimerComb{$Combos}\n";
	}
print "******FWD-FWD Primer mis-primes********\n";
foreach my $Combos2 (sort keys %PrimerCombFW) {
	print "$Combos2\t$PrimerCombFW{$Combos2}\n";
	}
print "****REV-FWD Primer mis-primes****\n";
foreach my $Combos3 (sort keys %PrimerCombRV) {
	print "$Combos3\t$PrimerCombRV{$Combos3}\n";
	}
print "******REV-REV Primer mis-primes********\n";
foreach my $Combos4 (sort keys %PrimerCombRV_RV) {
	print "$Combos4\t$PrimerCombRV_RV{$Combos4}\n";
	}
