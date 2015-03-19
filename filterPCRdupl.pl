#!/usr/bin/perl

# filterPCRdupl.pl
# November 2010	
# Author: Linnéa Smeds (linnea.smeds@ebc.uu.se)

use strict;
use warnings;
use Getopt::Long;
use List::Util qw[min max];


my $usage = "
# filterPCRdupl.pl
# November 2010	
# Author: Linnéa Smeds (linnea.smeds\@ebc.uu.se)
# -----------------------------------------------------------------------------
# Description: Extract all non-redundant read pairs in a fastq file pair, by 
# comparing the first N bases in both reads between all pairs. If there are 
# several copies of the same pair, only the pair with the highest quality 
# score is kept. Also saves a histogram of the copy distribution.
# Usage: perl filterPCRdupl.pl -fastq1=file1 -fastq2=file2 [-prefix=s -cmp=N]

-fastq1=file \t Fastq files with the read pairs, (pair1 in first file and pair2 
-fastq2=file \t in the second). Reads must have the same order in both files.
-prefix=string \t Prefix for the output files. The filtered fastq files will be
 \t\t named prefix_uniq1.fastq and prefix_uniq2.fastq, and the 
 \t\t histogram with copy distribution will be named prefix_copy.hist.
 \t\t [same prefix as file1].
-cmp=number \t The number of bases from the start of the reads in each pair
 \t\t that will be used for comparison [50].
-h \t\t Print this help message.\n";
				

# Starting time
my $time = time;

# Input parameters
my ($read1, $read2, $prefix, $baseComp, $help);
GetOptions(
 	"fastq1=s" => \$read1,
 	"fastq2=s" => \$read2,
   	"cmp=i" => \$baseComp,
   	"prefix=s" => \$prefix,
	"h" => \$help);

#Checking input, set default if not given
unless($read1 && $read2) {
	die $usage . "\n";
}
if($help) {
	die $usage . "\n";
}
unless($baseComp) {
	$baseComp=50;
}
unless($prefix) {
	$prefix = $read1;
	$prefix =~ s/\.\w+//;
}
unless(-e $read1) {
	die "Error: File $read1 doesn't exist!\n";
}
unless(-e $read2) {
	die "Error: File $read2 doesn't exist!\n";
}

my %reads =();

open(RD1, $read1);
open(RD2, $read2);

print "uniqueFastqPairs started " . localtime() . "\n";

$| = 1;
print "Going through all read pairs...\n";
my($hd1, $hd2, $seq1, $seq2, $plus1, $plus2, $qual1, $qual2);
my ($totCnt, $duplCnt)=(0,0);
while(<RD1>) {
	$hd1=$_;
	chomp($seq1=<RD1>);
	$plus1=<RD1>;
	$qual1=<RD1>;

	$hd2 =<RD2>;
	chomp($seq2=<RD2>);
	$plus2=<RD2>;
	$qual2=<RD2>;
	
	if(min(length($seq1), length($seq2))<$baseComp) {
		my $smallest =min(length($seq1), length($seq2));
		die "Error: Try to compare $baseComp first bases, but the read length is $smallest.\n".
			"Try decreasing -cmp below $smallest\n";
	}
	my $key = substr($seq1, 0, $baseComp) . substr($seq2, 0, $baseComp);
	my @t1 = split(//, $qual1);
	my @t2 = split(//, $qual2);

	my $score=0;
	for(my $i=0; $i<scalar(@t1); $i++) {
		my $temp1 = ord($t1[$i])-64;
		$score+=$temp1;	
	}
	for(my $i=0; $i<scalar(@t2); $i++) {
		my $temp2 = ord($t2[$i])-64;
		$score+=$temp2;	
	}

	if(defined $reads{$key}) {
		if($score > $reads{$key}{'scr'}) {
			 $reads{$key}{'id'}=$hd1;
			 $reads{$key}{'scr'}=$score;
		}
		$reads{$key}{'cpy'}++;
		$duplCnt++;
	}
	else {
		$reads{$key}{'scr'}=$score;
		$reads{$key}{'id'}=$hd1;
		$reads{$key}{'cpy'}=1;
	}
	$totCnt++;
}
close(RD1);
close(RD2);

my $curr_time = time-$time;
if($curr_time<60) {
	print "\tdone in $curr_time seconds.\n";
}
else {
	$curr_time=int($curr_time/60 + 0.5);
	print "\tdone in $curr_time minutes.\n";
}
$curr_time = time;

print "Writing unique read pairs to output...\n";

my $out1 = $prefix . "_uniq1.fastq";
my $out2 = $prefix . "_uniq2.fastq";
my %copies = ();

open(OUT1, ">$out1");
open(OUT2, ">$out2");
open(RD1, $read1);
open(RD2, $read2);

while(<RD1>) {
	$hd1=$_;
	chomp($seq1=<RD1>);
	$plus1=<RD1>;
	$qual1=<RD1>;

	$hd2 =<RD2>;
	chomp($seq2=<RD2>);
	$plus2=<RD2>;
	$qual2=<RD2>;

	my $key = substr($seq1, 0, $baseComp) . substr($seq2, 0, $baseComp);


	if(defined $reads{$key} && $reads{$key}{'id'} eq $hd1) {
		print OUT1 $hd1 . $seq1 . "\n" . $plus1 . $qual1;
		print OUT2 $hd2 . $seq2 . "\n" . $plus2 . $qual2;
		if(defined $copies{$reads{$key}{'cpy'}}) {
			$copies{$reads{$key}{'cpy'}}++;
		}
		else {
			$copies{$reads{$key}{'cpy'}}=1;
		}
		delete $reads{$key};	
	}
}
close(RD1);
close(RD2);
close(OUT1);
close(OUT2);

$curr_time = time-$curr_time;
if($curr_time<60) {
	print "\tdone in $curr_time seconds.\n";
}
else {
	$curr_time=int($curr_time/60 + 0.5);
	print "\tdone in $curr_time minutes.\n";
}
$curr_time = time;


print "Writing copy distribution histogram...\n";
my $distOut = $prefix ."_copy.hist";
open(OUT, ">$distOut");
foreach my $key (sort {$a<=>$b} keys %copies) {
	print OUT $key ."\t". $copies{$key} ."\n";
}
close(OUT);
print "\tdone!\n";

my $frac = 100*($duplCnt/$totCnt);
$frac = sprintf "%.2f", $frac;
print "$totCnt pairs in total, $duplCnt ($frac%) duplicated pairs removed.\n";

$time = (time-$time);
if($time<60) {
	print "Total time elapsed: $time seconds.\n";
}
else {
	$time=int($time/60 + 0.5);
	print "Total time elapsed: $time minutes.\n";
}


