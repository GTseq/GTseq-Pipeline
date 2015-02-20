#!/usr/bin/perl
#GenoCompile_v2.pl by Nate Campbell
#Compile genotypes from individual genotype files from GTseq_Genotyper_v2 output ".genos" files...
#This version utilizes the expanded output from the GTseq_Genotyper_v2 script to gather summary data and does not
#require the individual fastq files.
#For optional output formats use arguments. N for numeric genotypes or C for allele counts.  Defaults to SNP genotypes.
#Optional filtered output requires 2 argument values.  Genotype output type and a genotyping threshold [S,N, or C] [90]
# example: $ GTseq_GenoCompile_v2.pl S 90 (outputs SNP genotypes for individual .genos files with 90% or higher genotyping percentage)
# genotypes for individuals with less than the threshold genotyping percentage are converted to "00".

use strict; use warnings;

my $flag = "S";
my $geno_thresh = 0;
if (@ARGV == 1) {$flag = $ARGV[0];}
if (@ARGV == 2) {$flag = $ARGV[0]; $geno_thresh = $ARGV[1];}


my @Files = `ls *genos`;
chomp ( @Files );
print "Sample,Raw Reads,On-Target Reads,%On-Target,%GT";

open (FILE1, "<$Files[0]") or die;
	while (<FILE1>) {
		if ($. > 1){
			my @info1 = split ",", $_;
			my $assay1 = $info1[0];
			print ",$assay1";
			}
		}
close FILE1;

print "\n";

foreach my $samples (@Files) {
	my $raw = 0;
	my $on_target = 0;
	my $GT_pct = 0;
	my $num_targets = 0;
	my $sample_name = $samples;
	$sample_name =~ s/.genos//;
	print "$sample_name,";
	open (FILE, "<$samples") or die;
	while (<FILE>) {
	if ($. == 1) {my @summary = split ",", $_; $summary[1] =~ s/Raw-Reads\://; print "$summary[1],"; $raw = $summary[1];}
	elsif ($. > 1) {
		$num_targets++;
		chomp;
		my @info = split ",", $_;
		my $assay = $info[0];
		my $geno = $info[4];
		if ($geno =~ m/NA|00/) {$GT_pct++}
		my $count1 = $info[1];
		$count1 =~ s/.*=//;
		my $count2 = $info[2];
		$count2 =~ s/.*=//;
		$on_target = $on_target + $count1 + $count2;
				}
			}
		close FILE;
		$GT_pct = 100-($GT_pct/$num_targets*100);
		my $OT_pct = $on_target/$raw*100;
		my $Print_GT_pct = sprintf("%.2f", $GT_pct);
		$OT_pct = sprintf("%.2f", $OT_pct);
		print "$on_target,$OT_pct,$Print_GT_pct,";
		
	open (FILE, "<$samples") or die;
	while (<FILE>) {
		if ($. > 1) {
		chomp;
		my @info1 = split ",", $_;
		my $geno = $info1[4];
		my $L_count = 0;
		$info1[1] =~ s/.*=//;
		$info1[2] =~ s/.*=//;
		$L_count = $info1[1] + $info1[2];
		my $NumGT = "00";
		if($info1[5] =~ m/A1HOM/) {$NumGT = "11";}
		elsif($info1[5] =~ m/HET/) {$NumGT = "12";}
		elsif($info1[5] =~ m/A2HOM/) {$NumGT = "22";}
		elsif($info1[5] =~ m/NA/) {$NumGT = "00";}
		
		if(($flag =~ m/S/) && ($GT_pct >= $geno_thresh)) {print "$geno,";}
		elsif(($flag =~ m/C/) && ($GT_pct >= $geno_thresh)) {print "$L_count,";}
		elsif(($flag =~ m/N/) && ($GT_pct >= $geno_thresh)) {print "$NumGT,";}
		else {print "00,";}
				}
			}
		print "\n"; 
		close FILE;
	}
