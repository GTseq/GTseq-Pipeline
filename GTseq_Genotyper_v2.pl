#!/usr/bin/perl
#GTseq_Genotyper_v2.pl by Nate Campbell
#Genotype GTseq libraries from individual .fastq files
#100x faster and more accurate than the original version of the GTseq_Genotyper.pl script.  This version makes use of the forward primer sequences for each locus.
#Requires .csv file of [Locus Name,Allele1,Allele2,ProbeSeq1,ProbeSeq2,FWD_Primer,A1_correction,A2_correction] for all loci as $ARGV[0] and an individual fastq file.
#Correction values for input file are optional.  If correction values are not used, then allele correction values will be zero at all loci.
#Output files contain a header line with summary information followed by locus specific data.  Locus specific output fields are 
# LocusName, Allele1_counts, Allele2_counts, A1/A2-ratio, Genotype, Genotype_Class, A1_correction value, A2_correction value, On_target reads, Locus OT_percentage, Locus On-target reads as percentage of total on-target reads
# The penultimate field is defined as: #reads beginning with forward primer sequence and containing a probe sequence / #All reads beginning with forward primer sequence *100
# The final field is defined as: total number of on-target reads for locus / Overall on-target reads in panel *100

use strict; use warnings;

die "usage: provide csv file of locus, alleles, and probe seqs plus individual fastq\n" unless @ARGV == 2;

my %On_Target = ();
my %Off_Target = ();
my %F_Primer = ();
my %F_PrimerKey = ();
my %allele1name = ();
my %allele2name = ();
my %probeA1 = ();
my %probeA2 = ();
my %probeA1_RC = ();
my %probeA2_RC = ();
my %Allele1_Count = ();
my %Allele2_Count = ();
my %A1_corr = ();
my %A2_corr = ();
my $ArraySize = 0;
my $OT_Reads = 0;
my $Raw_Reads = 0;

#read in assay and allele information and push to hashes...

open(PROBES, "<$ARGV[0]") or die "error reading $ARGV[0]\n";

while (<PROBES>) {
	chomp;
	my @info = split(/,/, $_);
	$F_Primer{$info[0]} = substr $info[5], 0, 14;
	$F_PrimerKey{$F_Primer{$info[0]}} = $info[0];
	$allele1name{$info[0]} = $info[1];
	$allele2name{$info[0]} = $info[2];
	$probeA1{$info[0]} = $info[3];
	$probeA2{$info[0]} = $info[4];
	if (exists $info[6]) {$A1_corr{$info[0]} = $info[6]}
	else {$A1_corr{$info[0]} = 0}
	if (exists $info[7]) {$A2_corr{$info[0]} = $info[7]}
	else {$A2_corr{$info[0]} = 0}
	my $revcomp1 = reverse $info[3];
	$revcomp1 =~ tr/[ACGT]/[TGCA]/;
	$revcomp1 =~ tr/\]\[/\[\]/;
	$probeA1_RC{$info[0]} = $revcomp1;
	my $revcomp2 = reverse $info[4];
	$revcomp2 =~ tr/[ACGT]/[TGCA]/;
	$revcomp2 =~ tr/\]\[/\[\]/;
	$probeA2_RC{$info[0]} = $revcomp2;
	#Initialize Allele counts at zero...
	$Allele1_Count{$info[0]} = 0;
	$Allele2_Count{$info[0]} = 0;
	$On_Target{$info[0]} = 0;
	$Off_Target{$info[0]} = 0;
	$ArraySize++;
		}
close PROBES;


#Count alleles...
	open(FASTQ, "<$ARGV[1]") or die "error reading $ARGV[1]\n";

	while (<FASTQ>) {
		my $R1_ID1 = $_;
		my $R1_seq = <FASTQ>;
		my $R1_ID2 = <FASTQ>;
		my $R1_qual = <FASTQ>;
		chomp ($R1_seq);
		
		my $FP_seq = substr $R1_seq, 0, 14;
		$Raw_Reads++;
		if(exists $F_PrimerKey{$FP_seq}) {
			my $target = $F_PrimerKey{$FP_seq};
			if ($R1_seq =~ m/$probeA1{$target}|$probeA1_RC{$target}/){
			$Allele1_Count{$target}++; $On_Target{$target}++; $OT_Reads++;
			}
			elsif ($R1_seq =~ m/$probeA2{$target}|$probeA2_RC{$target}/){
			$Allele2_Count{$target}++; $On_Target{$target}++; $OT_Reads++;
			}
			else {$Off_Target{$target}++}
		}
	}
close FASTQ;

my $OT_Percentage = $OT_Reads/$Raw_Reads * 100;
$OT_Percentage = sprintf("%.1f", $OT_Percentage);
# print header line...
print "$ARGV[1],Raw-Reads:$Raw_Reads,On-Target reads:$OT_Reads,%On-Target:$OT_Percentage\n";

foreach my $loci (sort keys %F_Primer){
	my $A1fix = 0;
	my $A2fix = 0;
	my $sum_xy = $Allele1_Count{$loci} + $Allele2_Count{$loci};
	$Allele1_Count{$loci} = $Allele1_Count{$loci} - ($sum_xy / 4 * $A1_corr{$loci});
	if ($Allele1_Count{$loci} < 0) {$Allele1_Count{$loci} = 0}
	$Allele2_Count{$loci} = $Allele2_Count{$loci} - ($sum_xy / 4 * $A2_corr{$loci});
	if ($Allele2_Count{$loci} < 0) {$Allele2_Count{$loci} = 0}
	$Allele1_Count{$loci} = int ( $Allele1_Count{$loci} );
	$Allele2_Count{$loci} = int ( $Allele2_Count{$loci} );
	my $geno = "00";  #Initialize genotype variable at "00"...
	my $genoclass = "NA"; #Initialize genotype classification at "NA"...
	
#Fix allele counts to non-zero number for division ratio calculation...
	if ($Allele1_Count{$loci} == 0) {$A1fix = 0.1}
	else {$A1fix = $Allele1_Count{$loci}}
	if ($Allele2_Count{$loci} == 0) {$A2fix = 0.1}
	else {$A2fix = $Allele2_Count{$loci}}
	my $ratio = $A1fix/$A2fix;
	$ratio = sprintf("%.3f", $ratio);

	if ($Allele1_Count{$loci} + $Allele2_Count{$loci} < 10) {$geno = "00"; $genoclass = "NA";} #Set genotypes of low allele count loci to "00"
	elsif ($ratio >= 10) {$geno = "$allele1name{$loci}$allele1name{$loci}"; $genoclass = "A1HOM";} #Allele1 Homozygotes
	elsif ($ratio <= 0.1) {$geno = "$allele2name{$loci}$allele2name{$loci}"; $genoclass = "A2HOM";} #Allele2 Homozygotes
	elsif ($ratio <= 0.2) {$geno = "00"; $genoclass = "NA";} #In-betweeners
	elsif ($ratio <= 5) {$geno = "$allele1name{$loci}$allele2name{$loci}"; $genoclass = "HET";} #Heterozygotes
	
	if($Off_Target{$loci} == 0) {$Off_Target{$loci} = 0.01}
	my $On_Target_Per = ($On_Target{$loci}/($Off_Target{$loci} + $On_Target{$loci})) * 100; 
	my $Per_of_AllOTreads = $On_Target{$loci}/$OT_Reads * 100;
	$On_Target_Per = sprintf("%.1f", $On_Target_Per);
	$Per_of_AllOTreads = sprintf("%.3f", $Per_of_AllOTreads);

	print "$loci,$allele1name{$loci}=$Allele1_Count{$loci},$allele2name{$loci}=$Allele2_Count{$loci},$ratio,$geno,$genoclass,$A1_corr{$loci},$A2_corr{$loci},$On_Target{$loci},$On_Target_Per,$Per_of_AllOTreads\n";
		}
