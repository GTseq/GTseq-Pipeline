#!/usr/bin/perl
#GTseq_Genotyper_v3.pl by Nate Campbell
#Genotype GTseq libraries from individual .fastq files
#100x faster and more accurate than the original version of the GTseq_Genotyper.pl script.  This version makes use of the forward primer sequences for each locus.
#Requires .csv file of [Locus Name,Allele1,Allele2,ProbeSeq1,ProbeSeq2,FWD_Primer,A1_correction,A2_correction] for all loci as $ARGV[0] and an individual fastq file.
#Correction values for input file are optional.  If correction values are not used, then allele correction values will be zero at all loci.
#Output files contain a header line with summary information followed by locus specific data.  Locus specific output fields are 
# LocusName, Allele1_counts, Allele2_counts, A1/A2-ratio, Genotype, Genotype_Class, A1_correction value, A2_correction value, On_target reads, Locus OT_percentage, Locus On-target reads as percentage of total on-target reads
# The penultimate field is defined as: #reads beginning with forward primer sequence and containing a probe sequence / #All reads beginning with forward primer sequence *100
# The final field is defined as: total number of on-target reads for locus / Overall on-target reads in panel *100
# This version also outputs the IFI score (Individual fuzziness index) for each individual sample.  This is a measure of DNA cross contamination and is calculated
# using read counts from background signal at homozygous and No-Call loci.  Low scores are better than high scores.
# version 3 update: Includes fuzzy matching for detection of possible null alleles.  Requires the String::Approx perl extension.
# Possible null alleles can be summarized using the GTseq_ErrorReport_v3.pl script.  Keep in mind that some loci with flagged null alleles
# could be the result of co-amplification of a paralogous locus.
# In some cases the reported fuzzy match sequence doesn't match either probe sequence at all.  This happens when the location of the fuzzy match doesn't match
# the position of the exact probe match on the amplicon.  This probably means that there is amplification of a PSV with these primers.

use strict; use warnings;

use String::Approx 'amatch';

die "usage: provide csv file of locus, alleles, and probe seqs plus individual fastq\n" unless @ARGV == 2;

my %On_Target = ();
my %TargetStart = ();
my %ProbeLength = ();
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
my %print_line = ();
my %fuzzy_match = ();
my %fuzzy_seqs = ();
my %fuzzyCT = ();
my $ArraySize = 0;
my $OT_Reads = 0;
my $Raw_Reads = 0;
my %Null_Error = ();

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
	$fuzzy_match{$info[0]} = 0;
	$Null_Error{$info[0]} = "Pass";
	$TargetStart{$info[0]} = 20;
	$ProbeLength{$info[0]} = 20;
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
			$TargetStart{$target} = $-[0];
			$Allele1_Count{$target}++; $On_Target{$target}++; $OT_Reads++;
			}
			elsif ($R1_seq =~ m/$probeA2{$target}|$probeA2_RC{$target}/){
			$TargetStart{$target} = $-[0];
			$Allele2_Count{$target}++; $On_Target{$target}++; $OT_Reads++;
			}
			else {
				$Off_Target{$target}++;
				my $fuzz = 0;
				my $fuzzy_A1 = $probeA1{$target};
				my $subnumber = 0;
				my $fuzzy_set = "[ \"I1\", \"S1\", \"D1\" ]";
				if ($fuzzy_A1 =~ m/\[|\.|\?|\+/) {$fuzzy_A1 =~ s/\[..\]/N/g; $fuzzy_A1 =~ tr/[.?+]/[NNN]/; $subnumber = $fuzzy_A1 =~ tr/N/N/ + 1; $fuzzy_set = "[ \"I1\", \"S$subnumber\", \"D1\" ]";}
				my $fuzzy_A1RC = $probeA1_RC{$target};
				if ($fuzzy_A1RC =~ m/\[|\.|\?|\+/) {$fuzzy_A1RC =~ s/\[..\]/N/g; $fuzzy_A1RC =~ tr/[.?+]/[NNN]/; $subnumber = $fuzzy_A1RC =~ tr/N/N/ + 1; $fuzzy_set = "[ \"I1\", \"S$subnumber\", \"D1\" ]";}
				my $fuzzy_A2 = $probeA2{$target};
				if ($fuzzy_A2 =~ m/\[|\.|\?|\+/) {$fuzzy_A2 =~ s/\[..\]/N/g; $fuzzy_A2 =~ tr/[.?+]/[NNN]/; $subnumber = $fuzzy_A2 =~ tr/N/N/ + 1; $fuzzy_set = "[ \"I1\", \"S$subnumber\", \"D1\" ]";}
				my $fuzzy_A2RC = $probeA2_RC{$target};
				if ($fuzzy_A2RC =~ m/\[|\.|\?|\+/) {$fuzzy_A2RC =~ s/\[..\]/N/g; $fuzzy_A2RC =~ tr/[.?+]/[NNN]/; $subnumber = $fuzzy_A2RC =~ tr/N/N/ + 1; $fuzzy_set = "[ \"I1\", \"S$subnumber\", \"D1\" ]";}
				if (length $fuzzy_A1 >= length $fuzzy_A2) {$ProbeLength{$target} = length $fuzzy_A1}
				else {$ProbeLength{$target} = length $fuzzy_A2;}
				#print "$fuzzy_set\t$fuzzy_A1\t$fuzzy_A1RC\t$fuzzy_A2\t$fuzzy_A2RC\n";  #...Testing...
				if (amatch($fuzzy_A1, $fuzzy_set, $R1_seq)) {
					my $fuzzy_str = substr $R1_seq, $TargetStart{$target}, $ProbeLength{$target};
					if (length $fuzzy_str > 10) {
						$fuzz++;
						if (exists $fuzzy_seqs{$fuzzy_str}) {$fuzzyCT{$fuzzy_str}++}
						else {$fuzzy_seqs{$fuzzy_str} = $target; $fuzzyCT{$fuzzy_str}++}
						}
				}
				elsif (amatch($fuzzy_A1RC, $fuzzy_set, $R1_seq)) {
					my $fuzzy_str = substr $R1_seq, $TargetStart{$target}, $ProbeLength{$target};
					$fuzzy_str = reverse $fuzzy_str;
					$fuzzy_str =~ tr/ACGT/TGCA/;
					if (length $fuzzy_str > 10) {
						$fuzz++;
						if (exists $fuzzy_seqs{$fuzzy_str}) {$fuzzyCT{$fuzzy_str}++}
						else {$fuzzy_seqs{$fuzzy_str} = $target; $fuzzyCT{$fuzzy_str}++}
						}
				}
				elsif (amatch($fuzzy_A2, $fuzzy_set, $R1_seq)) {
					my $fuzzy_str = substr $R1_seq, $TargetStart{$target}, $ProbeLength{$target};
					if (length $fuzzy_str > 10) {
						$fuzz++;
						if (exists $fuzzy_seqs{$fuzzy_str}) {$fuzzyCT{$fuzzy_str}++}
						else {$fuzzy_seqs{$fuzzy_str} = $target; $fuzzyCT{$fuzzy_str}++}
						}
				}
				elsif (amatch($fuzzy_A2RC, $fuzzy_set, $R1_seq)) {
					my $fuzzy_str = substr $R1_seq, $TargetStart{$target}, $ProbeLength{$target};
					$fuzzy_str = reverse $fuzzy_str;
					$fuzzy_str =~ tr/ACGT/TGCA/;
					if (length $fuzzy_str > 10) {
						$fuzz++;
						if (exists $fuzzy_seqs{$fuzzy_str}) {$fuzzyCT{$fuzzy_str}++}
						else {$fuzzy_seqs{$fuzzy_str} = $target; $fuzzyCT{$fuzzy_str}++}
						}
				}
				if ($fuzz > 0) {$fuzzy_match{$target}++}
				}
		}
	}
close FASTQ;

if ($OT_Reads == 0) {$OT_Reads = 1}
my $OT_Percentage = $OT_Reads/$Raw_Reads * 100 unless $Raw_Reads == 0;
if ($Raw_Reads == 0) {$OT_Percentage = 0.0}
$OT_Percentage = sprintf("%.1f", $OT_Percentage);
# print header line...
print "$ARGV[1],Raw-Reads:$Raw_Reads,On-Target reads:$OT_Reads,%On-Target:$OT_Percentage,";

my $HOM_CT = 0;
my $BKGRD_CT = 0;
my $IFI = 0;

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
	elsif ($ratio >= 10) {$geno = "$allele1name{$loci}$allele1name{$loci}"; $genoclass = "A1HOM"; #Allele1 Homozygotes
		$HOM_CT = $HOM_CT + $Allele1_Count{$loci}; $BKGRD_CT = $BKGRD_CT + $Allele2_Count{$loci};} 
	elsif (($ratio < 10) && ($ratio > 5)) {$geno = "00"; $genoclass = "NA"; #In-betweeners
		$HOM_CT = $HOM_CT + $Allele1_Count{$loci}; $BKGRD_CT = $BKGRD_CT + $Allele2_Count{$loci};}
	elsif ($ratio <= 0.1) {$geno = "$allele2name{$loci}$allele2name{$loci}"; $genoclass = "A2HOM"; #Allele2 Homozygotes
		$HOM_CT = $HOM_CT + $Allele2_Count{$loci}; $BKGRD_CT = $BKGRD_CT + $Allele1_Count{$loci};} 
	elsif ($ratio <= 0.5) {$geno = "00"; $genoclass = "NA"; #In-betweeners
		$HOM_CT = $HOM_CT + $Allele2_Count{$loci}; $BKGRD_CT = $BKGRD_CT + $Allele1_Count{$loci};}
	elsif ($ratio <= 2) {$geno = "$allele1name{$loci}$allele2name{$loci}"; $genoclass = "HET";} #Heterozygotes
	
#Calculate the presence of possible null alleles using fuzzy match ratio...
	
	if ($sum_xy == 0) {$sum_xy = 0.1};
	my $fuzz_ratio = $fuzzy_match{$loci}/$sum_xy;
	if (($genoclass !~ "HET") && ($fuzz_ratio > 0.8) && ($sum_xy > 10)) {$Null_Error{$loci} = "Fail"}
	
	if($Off_Target{$loci} == 0) {$Off_Target{$loci} = 0.01}
	my $On_Target_Per = ($On_Target{$loci}/($Off_Target{$loci} + $On_Target{$loci})) * 100; 
	my $Per_of_AllOTreads = $On_Target{$loci}/$OT_Reads * 100;
	$On_Target_Per = sprintf("%.1f", $On_Target_Per);
	$Per_of_AllOTreads = sprintf("%.3f", $Per_of_AllOTreads);
	$print_line{$loci} = "$loci,$allele1name{$loci}=$Allele1_Count{$loci},$allele2name{$loci}=$Allele2_Count{$loci},$ratio,$geno,$genoclass,$A1_corr{$loci},$A2_corr{$loci},$On_Target{$loci},$On_Target_Per,$Per_of_AllOTreads,$Null_Error{$loci},$fuzzy_match{$loci}";
		}
		
if ($HOM_CT == 0) {$HOM_CT = 1}
$IFI = $BKGRD_CT/$HOM_CT * 100;
$IFI = sprintf("%.2f", $IFI);
print "IFI_score:$IFI\n";

#Sort common fuzzy match sequences...
my %fuzz_common = ();
my %fuzz_commonCT = ();

foreach my $sorted_fuzzies (sort { $fuzzyCT{$a} <=> $fuzzyCT{$b} } keys %fuzzyCT) {
	my $fuzzlocus = $fuzzy_seqs{$sorted_fuzzies};
	if ($Null_Error{$fuzzlocus} =~ m/Fail/) {
		if ((exists $fuzz_common{$fuzzlocus}) && ($fuzz_commonCT{$fuzzlocus} < $fuzzyCT{$sorted_fuzzies})) {
			$fuzz_common{$fuzzlocus} = $sorted_fuzzies;
			$fuzz_commonCT{$fuzzlocus} = $fuzzyCT{$sorted_fuzzies};
		}
		else {
			$fuzz_common{$fuzzlocus} = $sorted_fuzzies;
			$fuzz_commonCT{$fuzzlocus} = $fuzzyCT{$sorted_fuzzies};
			}
		}
	}

foreach my $loci2 (sort keys %F_Primer) {
	if ($Null_Error{$loci2} =~ m/Fail/) {$print_line{$loci2} = "$print_line{$loci2},$fuzz_common{$loci2},$fuzz_commonCT{$loci2}";}
	print "$print_line{$loci2}\n";
	}
