#!/usr/bin/perl
# GTseq_ErrorReport_v3.pl
# by Nate Campbell
# Execute this script in a directory containing .genos files generated using the GTseq_Genotyper_v3.pl script.
# This script summarizes possible null alleles and their sequences within the genotyped samples.  In some instances these sequences are the 
# result of co-amplification of a similar locus and should be ignored.  In other cases they may be the result of an undetected allele (null allele)
# within some of the genotyped samples and a modification of the input file is warranted.  It is up to the user to make these determinations and modify
# the input file accordingly.
# Usage: Include genotyping input file as $ARGV[0] on command line;

use strict; use warnings;

my @File_List = `ls *genos`;

chomp ( @File_List );

print "Locus,A1_probe,A2_probe,Possible_Null,Samples_Flagged\n";

my %Nulls = (); #store loci as values and sequences as keys...
my %NullCT = ();  #store counts as values and sequences as keys...

foreach my $file (@File_List) {
	open (FILE, "<$file") or die "Error opening $file\n";
	while (<FILE>) {
		if ($. > 1) {
			my @info = split ",", $_;
			if ($info[11] =~ m/Fail/) {
				$Nulls{$info[13]} = $info[0];
				if (exists $NullCT{$info[13]}) {$NullCT{$info[13]}++}
				else {$NullCT{$info[13]}++}
				}
			}
		}
	close FILE;
	}

my $input_file = $ARGV[0];

foreach my $seqs (sort {$Nulls{$a} cmp $Nulls{$b} } keys %Nulls) {
my $locus = $Nulls{$seqs};
my $counts = $NullCT{$seqs};
my $A1 = 0;
my $A2 = 0;

open (INPUT, "<$input_file") or die "Error opening input file\n";
	while (<INPUT>) {
		my @info2 = split ",", $_;
		if ($info2[0] =~ m/$locus/) {$A1 = $info2[3]; $A2 = $info2[4];}
		}
	close INPUT;
print "$locus,$A1,$A2,$seqs,$counts\n";
}
