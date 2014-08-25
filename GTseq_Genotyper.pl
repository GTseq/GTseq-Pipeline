#!/usr/bin/perl
#GTseq_Genotyper.pl by Nate Campbell
#Genotype GTseq libraries from individual .fastq files
#Requires .csv file of [Locus Name,Allele1,Allele2,ProbeSeq1,ProbeSeq2] for all loci as $ARGV[0] and an individual fastq file.

use strict; use warnings;

die "usage: provide csv file of locus, alleles, and probe seqs plus individual fastq\n" unless @ARGV == 2;

my @assays = ();
my @allele1name = ();
my @allele2name = ();
my @probeA1 = ();
my @probeA2 = ();
my @probeA1_RC = ();
my @probeA2_RC = ();
my @Allele1_Count = ();
my @Allele2_Count = ();

#read in assay and allele information and push to arrays...

open(PROBES, "<$ARGV[0]") or die "error reading $ARGV[0]\n";

while (<PROBES>) {
	chomp;
	my @info = split(/,/, $_);
	push @assays, $info[0];
	push @allele1name, $info[1];
	push @allele2name, $info[2];
	push @probeA1, $info[3];
	push @probeA2, $info[4];
	my $revcomp1 = reverse $info[3];
	$revcomp1 =~ tr/[ACGT]/[TGCA]/;
	$revcomp1 =~ tr/\]\[/\[\]/;
	push @probeA1_RC, $revcomp1;
	my $revcomp2 = reverse $info[4];
	$revcomp2 =~ tr/[ACGT]/[TGCA]/;
	$revcomp2 =~ tr/\]\[/\[\]/;
	push @probeA2_RC, $revcomp2;
		}
close PROBES;

my $ArraySize = @assays;

#Initialize Allele counts at zero...
for (my $h = 0; $h < $ArraySize; $h++) {
	$Allele1_Count[$h] = 0;
	$Allele2_Count[$h] = 0;
	}

#Count alleles...
	open(FASTQ, "<$ARGV[1]") or die "error reading $ARGV[1]\n";

	while (<FASTQ>) {
		my $R1_ID1 = $_;
		my $R1_seq = <FASTQ>;
		my $R1_ID2 = <FASTQ>;
		my $R1_qual = <FASTQ>;
		chomp ($R1_seq);

		for (my $i = 0; $i < $ArraySize; $i++){

			if ($R1_seq =~ m/$probeA1[$i]|$probeA1_RC[$i]/){
			$Allele1_Count[$i]++;
			}
			elsif ($R1_seq =~ m/$probeA2[$i]|$probeA2_RC[$i]/){
			$Allele2_Count[$i]++;
			}
		}
	}
close FASTQ;

# print "$Allele1_Count[0]\n"; #testing...

for (my $j = 0; $j < $ArraySize; $j++){
	
	my $A1fix = 0;
	my $A2fix = 0;
	my $geno = "NA";  #Initialize genotype variable at "NA"...
	my $genoclass = "NA"; #Initialize genotype classification at "NA"...
	
#Fix allele counts to non-zero number for division ratio calculation...
	if ($Allele1_Count[$j] == 0) {$A1fix = 0.1}
	else {$A1fix = $Allele1_Count[$j]}
	if ($Allele2_Count[$j] == 0) {$A2fix = 0.1}
	else {$A2fix = $Allele2_Count[$j]}
	my $ratio = $A1fix/$A2fix;

	if ($Allele1_Count[$j] + $Allele2_Count[$j] < 10) {$geno = "NA"; $genoclass = "NA";} #Set genotypes of low allele count loci to "NA"
	elsif ($ratio >= 10) {$geno = "$allele1name[$j]$allele1name[$j]"; $genoclass = "A1HOM";} #Allele1 Homozygotes
	elsif ($ratio <= 0.1) {$geno = "$allele2name[$j]$allele2name[$j]"; $genoclass = "A2HOM";} #Allele2 Homozygotes
	elsif ($ratio <= 0.2) {$geno = "NA"; $genoclass = "NA";} #In-betweeners
	elsif ($ratio <= 5) {$geno = "$allele1name[$j]$allele2name[$j]"; $genoclass = "HET";} #Heterozygotes

	print "$assays[$j],$allele1name[$j]=$Allele1_Count[$j],$allele2name[$j]=$Allele2_Count[$j],$ratio,$geno,$genoclass\n";
		}
