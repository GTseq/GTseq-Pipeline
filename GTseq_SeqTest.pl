#!/usr/bin/perl
#GTseq_SeqTest.pl by Nate Campbell
#Test .hash file to see how many times loci fwd primers and probe seqs occur.

use strict; use warnings;

die "usage: provide <tab delimited txt file of locus, fwd, probe1, probe2 seqs> and <*.hash>\n" unless @ARGV == 2;

my @assays = ();
my @fwd_seq = ();
my @probe1 = ();
my @probe2 = ();
my @probe1RC = ();
my @probe2RC = ();

my @fwd_count = ();
my @probe_count = ();
my @both_count = ();

#read in assay and allele information and push to arrays...

open(SEQ, "<$ARGV[0]") or die "error reading $ARGV[0]\n";

while (<SEQ>) {
	chomp;
	my @info = split(/\t/, $_);
	push @assays, $info[0];
	push @fwd_seq, $info[1];
	push @probe1, $info[2];
	push @probe2, $info[3];
	my $p1RC = reverse $info[2];
	my $p2RC = reverse $info[3];
	$p1RC =~ tr/ACGT/TGCA/;
	$p1RC =~ tr/\]\[/\[\]/;
	$p2RC =~ tr/ACGT/TGCA/;
	$p2RC =~ tr/\]\[/\[\]/;
	push @probe1RC, $p1RC;
	push @probe2RC, $p2RC;
		}
close SEQ;

my $Targets = @assays;

for (my $i = 0; $i < $Targets; $i++){
$fwd_count[$i] = 0;
$probe_count[$i] = 0;
$both_count[$i] = 0;

open(HASH, "<$ARGV[1]") or die "error reading $ARGV[1]\n";

	while (<HASH>) {
		my $hash = $_;
		my $R1_seq = <HASH>;
		chomp ($hash);
		my @info = split(/;/, $hash);
		my $count = $info[2];

			if ($R1_seq =~ m/$fwd_seq[$i]/){
			$count = $fwd_count[$i] + $count;
			$fwd_count[$i] = $count;
			}
			$count = $info[2];
			if ($R1_seq =~ m/$probe1[$i]|$probe2[$i]|$probe1RC[$i]|$probe2RC[$i]/){
			$count = $probe_count[$i] + $count;
			$probe_count[$i] = $count;
			}
			$count = $info[2];
			if (($R1_seq =~ m/$fwd_seq[$i]/) && ($R1_seq =~ m/$probe1[$i]|$probe2[$i]|$probe1RC[$i]|$probe2RC[$i]/)) {
			$count = $both_count[$i] + $count;
			$both_count[$i] = $count;
			}
		}
	}
close HASH;

# print "$Allele1_Count[0]\n"; #testing...

for (my $j = 0; $j < $Targets; $j++){
	print "$assays[$j],$fwd_count[$j],$probe_count[$j],$both_count[$j]\n";
	}
