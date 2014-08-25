#!/usr/bin/perl
#GTseq_GenoCompile.pl by Nate Campbell
#Compile genotypes and provide summary data from individual genotype files from GTseq analysis...

use strict; use warnings;

my @Files = `ls *genos`;
chomp ( @Files );
print "Sample,Raw Reads,On-Target Reads,%On-Target,%GT";

open (FILE1, "<$Files[0]") or die;
	while (<FILE1>) {
		my @info1 = split ",", $_;
		my $assay1 = $info1[0];
		print ",$assay1";
		}
close FILE1;

print "\n";

foreach my $samples (@Files) {
	my $on_target = 0;
	my $fastq = $samples;
	$fastq =~ s/.genos/.fastq/;
	my $raw = `wc -l $fastq`;
	$raw =~ s/\s.*fastq//;
	$raw = $raw/4;
	my $GT_pct = 0;
	my $num_targets = 0;
	my $sample_name = $samples;
	$sample_name =~ s/.genos//;
	print "$sample_name,$raw,";
	open (FILE, "<$samples") or die;
	while (<FILE>) {
		$num_targets++;
		chomp;
		my @info = split ",", $_;
		my $assay = $info[0];
		my $geno = $info[4];
		if ($geno =~ m/NA/) {$GT_pct++}
		my $count1 = $info[1];
		$count1 =~ s/.=//;
		my $count2 = $info[2];
		$count2 =~ s/.=//;
		$on_target = $on_target + $count1 + $count2;
			}
		close FILE;
		$GT_pct = 100-($GT_pct/$num_targets*100);
		my $OT_pct = $on_target/$raw*100;
		print "$on_target,$OT_pct,$GT_pct,";
		
	open (FILE, "<$samples") or die;
	while (<FILE>) {
		chomp;
		my @info = split ",", $_;
		my $assay = $info[0];
		my $geno = $info[4];
		print "$geno,";
			}
		print "\n"; 
		close FILE;
		}
