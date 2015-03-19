#!/usr/bin/perl

# condetri.pl
# September 2010, mod Oct 2011, Mar 2012
# Author: Linnéa Smeds (linnea.smeds@ebc.uu.se)

use strict;
use warnings;
use Getopt::Long;


my $usage = "# condetri.pl
# September 2010
# Version 2.1, Mar 2012
# Author: Linnéa Smeds (linnea.smeds\@ebc.uu.se)
# ---------------------------------------------------------------------------------
# Description: Trim reads from the 3'-end and extract reads (or read pairs) of
# good quality. If the reads are paired, the filtering is done pairwise, and 
# if one read in a pair has low quality, the remaning read is saved as single end.
# Usage: perl condetri.pl -fastq1=file1 [-fastq2=file2 -prefix=s -cutfirst=i 
# -rmN -hq=i -lq=i -frac=i -lfrac=i -minlen=i -mh=i -ml=i -sc=i -pb=s]

-fastq1=file \t Fastq file. If a second file is given, the files are trimmed
-fastq2=file \t as a pair. The reads must have the same order in both files.
-prefix=string \t Prefix for the output file(s). The filtered fastq file(s) will
 \t\t be named prefix_trim1.fastq (and prefix_trim2.fastq if present). For pairs,
 \t\t a third file will be given with unpaired reads (reads from pairs where one 
 \t\t low quality read has been removed).
-cutfirst=i \t Remove i first bases from the 5'end [0].
-rmN \t\t Remove non-ATCG bases from 5'end before any trimming [no].
-hq=i \t\t Hiqh quality threshold [25].
-lq=i \t\t Low quality threshold [10].
-frac=[0,1]\t Fraction of read that must exceed hq [0.8].
-lfrac=[0,1]\t Maximum fraction of reads with qual<lq [0].
-minlen=i \t Min allowed read length [50].
-mh=i \t\t When this no of sequential hq bases is reached, the trimming stops [5].
-ml=i \t\t Max no of lq bases allowed after a stretch of hq bases from 3'-end [1].
-sc=i\t\t Scoring table, Score=ASCII-sc, 64 for Solexa and Illumina 1.3+ and
 \t\t Illumina 1.5+, 33 for Sanger and newer Illumina data (from 1.8+). [64].
-pb=fa|fq\t Print the removed low quality reads to a fasta file or a fastq
\t\t file (for investigation/evaluation purposes) [no].
-q \t\t Print Solexa/Illumina scoring table (64 offset).
-q33 \t\t Print Sanger/new Illumina scoring table (33 offset).
-h \t\t Print this help message.\n";
		
# Starting time
my $time = time;

# Input parameters
my ($read1,$read2,$prefix,$cutfirst,$removeN,$HQlim,$lowlim,$minfrac,$lfrac,$minReadLen,
	$maxNoHQ,$maxNoLQ,$scoring,$printBad,$tbl,$tbl33,$help);
GetOptions(
  	"fastq1=s" => \$read1,
   	"fastq2=s" => \$read2,
  	"prefix=s" => \$prefix,
	"cutfirst=i" => \$cutfirst,
	"rmN" => \$removeN,
  	"hq=i" => \$HQlim,
	"lq=i" => \$lowlim,
	"frac=s" => \$minfrac,
	"lfrac=s" => \$lfrac,
	"minlen=i" => \$minReadLen,
	"mh=i" => \$maxNoHQ,
   	"ml=i" => \$maxNoLQ,
	"sc=i" => \$scoring,
	"pb=s" => \$printBad,
	"q" => \$tbl,
	"q33" => \$tbl33,
	"h" => \$help);
	

#--------------------------------------------------------------------------------
#Checking input, set default if not given
if($tbl) {
	&table64();
	exit;
}
if($tbl33) {
	&table33();
	exit;
}
unless($read1) {
	die $usage . "\n";
}
if($help) {
	die $usage . "\n";
}
unless($prefix) {
	$prefix = $read1;
	$prefix =~ s/\.\w+//;
}
unless($cutfirst) {
	$cutfirst=0;
}
unless($HQlim) {
	$HQlim=25;
}
unless($lowlim) {
	$lowlim=10;
}
if($minfrac) {
	if($minfrac<0 || $minfrac>1) {
		die "Error: frac must be between 0 and 1.\n";
	}
}
else{
	$minfrac=0.8;
}
if($lfrac) {
	if($lfrac<0 || $lfrac>1) {
		die "Error: lfrac must be between 0 and 1.\n";
	}
}
else{
	$lfrac=0;
}
unless($maxNoHQ) {
	$maxNoHQ=5;
}
unless($maxNoLQ) {
	$maxNoLQ=1;
}
unless($minReadLen) {
	$minReadLen=50;
}
unless($scoring) {
	$scoring=64;
}
if($printBad) {
	unless($printBad eq "fa" || $printBad eq "fq") {
		die "Error: pb must be either \"fa\" or \"fq\".\n";
	}
}
else {
	$printBad="";
}
unless(-e $read1) {
	die "Error: File $read1 doesn't exist!\n";
}
print "\ncondetri.pl started " . localtime() . "\n";
print "------------------------------------------------------------------\n";
print "Settings: ";
if($cutfirst>0) {
	print "Cutfirst=$cutfirst "; 
}
print "HQ=$HQlim LQ=$lowlim Frac=$minfrac Lfrac=$lfrac MH=$maxNoHQ ML=$maxNoLQ Minlen=$minReadLen Scoring=$scoring\n";
if($removeN) {
	print "Remove non-ATCG characters from 5'-end before trimming\n";
}
if($printBad ne "") {
	print "Print removed reads to a $printBad file.\n";
}


my ($totNoReads, $pairReads, $unpairedReads, $badReads) = (0,0,0,0);
my ($totNoBases, $pairBases, $unpairedBases, $firstbases) = (0,0,0,0);
#--------------------------------------------------------------------------------
# Trimming paired files
if($read2) {
	unless(-e $read2) {
		die "Error: File $read2 doesn't exist!\n";
	}
	my $out1 = $prefix . "_trim1.fastq";
	my $out2 = $prefix . "_trim2.fastq";
	my $out3 = $prefix . "_trim_unpaired.fastq";
	my $out4 = $prefix . "_badreads.".$printBad;
	open(IN1, $read1);
	open(IN2, $read2);
	open(OUT1, ">$out1");
	open(OUT2, ">$out2");
	open(OUT3, ">$out3");
	if($printBad ne "") {
		open(OUT4, ">$out4");
	}
	my($head1,$seq1,$plus1,$qual1,$head2,$seq2,$plus2,$qual2);
	print "Processing...\n";

	while(my $line = <IN1>) {
		$head1 = $line;
		chomp($seq1 = <IN1>);
		$plus1 = <IN1>;
		chomp($qual1 = <IN1>);
		
		$head2 = <IN2>;
		chomp($seq2 = <IN2>);
		$plus2 = <IN2>;
		chomp($qual2 = <IN2>);
		
		my($newseq1, $newscore1, $newseq2, $newscore2)=($seq1,$qual1,$seq2,$qual2);	

		# Cut from 5'-end
		if($cutfirst>0) {
			$newseq1 = substr($seq1, $cutfirst, length($seq1)-$cutfirst);
			$newscore1 = substr($qual1, $cutfirst, length($qual1)-$cutfirst);
			$newseq2 = substr($seq2, $cutfirst, length($seq2)-$cutfirst);
			$newscore2 = substr($qual2, $cutfirst, length($qual2)-$cutfirst);
		}
		
		# Remove non-ATCG characters from 5'end
		if($removeN) {
			my ($tmp1, $tmp2);
			($newseq1, $newscore1, $tmp1) = &removeNs($seq1, $qual1);
			($newseq2, $newscore2, $tmp2) = &removeNs($seq2, $qual2);
			$firstbases=$firstbases+$tmp1+$tmp2;
		}
		
		# Trim both reads
		($newseq1, $newscore1) = &trimEnd($newseq1,$newscore1, $HQlim, $lowlim, $maxNoHQ, $maxNoLQ, $minReadLen, $scoring);
		($newseq2, $newscore2) = &trimEnd($newseq2,$newscore2, $HQlim, $lowlim, $maxNoHQ, $maxNoLQ, $minReadLen, $scoring);

		# Check if reads are ok, print good reads
		if(readOK($newscore1, $HQlim, $lowlim, $minfrac, $lfrac, $scoring, $minReadLen)) {
			if(&readOK($newscore2, $HQlim, $lowlim, $minfrac, $lfrac, $scoring, $minReadLen)) {
				print OUT1 $head1 . $newseq1 ."\n" . $plus1 . $newscore1 . "\n";
				print OUT2 $head2 . $newseq2 ."\n" . $plus2 . $newscore2 . "\n";
				$pairReads+=2;
				$pairBases+=length($newseq1)+length($newseq2);
			}
			else {
				print OUT3 $head1 . $newseq1 ."\n" . $plus1 . $newscore1 . "\n";
				$unpairedReads++;
				$unpairedBases+=length($newseq1);
				if($printBad eq "fa") {
					print OUT4 ">".$head2.$seq2."\n";
					$badReads++;
				}
				elsif($printBad eq "fq") {
					print OUT4 $head2.$seq2."\n".$plus2.$qual2."\n";
					$badReads++;
				}
			}
		}
		else {
			if(&readOK($newscore2, $HQlim, $lowlim, $minfrac, $lfrac, $scoring, $minReadLen)) {
				print OUT3 $head2 . $newseq2 ."\n" . $plus2 . $newscore2 . "\n";
				$unpairedReads++;
				$unpairedBases+=length($newseq2);
				if($printBad eq "fa") {
					print OUT4 ">".$head1.$seq1."\n";
					$badReads++;
				}
				elsif($printBad eq "fq") {
					print OUT4 $head1.$seq1."\n".$plus1.$qual1."\n";
					$badReads++;
				}
			}
			else {
				if($printBad eq "fa" ) {
					print OUT4 ">".$head1.$seq1."\n>".$head2.$seq2."\n";
					$badReads+=2;
				}
				elsif($printBad eq "fq") {
					print OUT4 $head1.$seq1."\n".$plus1.$qual1."\n".
								$head2.$seq2."\n".$plus2.$qual2."\n";
					$badReads+=2;
				}
			}
		}
		$totNoBases+=length($seq1)+length($seq2);
		$totNoReads+=2;
		if ($totNoReads%100000==0) {
			print "$totNoReads reads processed\r";
		}
	}
}
# Trimming single end files
else {
	my $out1 = $prefix . "_trim.fastq";
	my $out4 = $prefix . "_badreads.".$printBad;
	open(IN1, $read1);
	open(OUT1, ">$out1");
	if($printBad ne "") {
		open(OUT4, ">$out4");
	}	

	print "Processing...\n";
	my($head1,$seq1,$plus1,$qual1);

	while(my $line = <IN1>) {
		$head1 = $line;
		chomp($seq1 = <IN1>);
		$plus1 = <IN1>;
		chomp($qual1 = <IN1>);

		my($newseq1, $newscore1) = ($seq1,$qual1);

		# Cut from 5'-end
		if($cutfirst>0) {
			$newseq1 = substr($seq1, $cutfirst, length($seq1)-$cutfirst);
			$newscore1 = substr($qual1, $cutfirst, length($qual1)-$cutfirst);
		}
		
		# Remove non-ATCG characters from 5'end
		if($removeN) {
			my $tmp1;
			($newseq1, $newscore1, $tmp1) = &removeNs($newseq1, $newscore1);
			$firstbases+=$tmp1;
		}
		
		($newseq1, $newscore1) = &trimEnd($newseq1, $newscore1, $HQlim, $lowlim, $maxNoHQ, $maxNoLQ, $minReadLen, $scoring);
		if(readOK($newscore1, $HQlim, $lowlim, $minfrac, $lfrac, $scoring, $minReadLen)) {
			print OUT1 $head1 . $newseq1 ."\n" . $plus1 . $newscore1 . "\n";
			$unpairedReads++;
			$unpairedBases+=length($newseq1);
		}
		else {
			if($printBad eq "fa") {
				print OUT4 ">".$head1.$seq1."\n";
				$badReads++;
			}
			elsif($printBad eq "fq") {
				print OUT4 $head1.$seq1."\n".$plus1.$qual1."\n";
				$badReads++;
			}
		}
		$totNoBases+=length($seq1);
		$totNoReads+=1;
		if ($totNoReads%100000==0) {
			print "$totNoReads reads processed\r";
		}
	}
	
}

$totNoBases+=$cutfirst*$totNoReads;
$firstbases+=$cutfirst*$totNoReads;
#--------------------------------------------------------------------------------
# Print statistics to table
open(STATS, ">$prefix".".stats");
print STATS $prefix."\t".$totNoReads."\t".$totNoBases."\t".$pairReads."\t".$pairBases.
		"\t".$unpairedReads."\t".$unpairedBases."\n";

print "\nDone!\n";
print "------------------------------------------------------------------\n";
print "$totNoReads reads with $totNoBases bases in input files\n";
if($firstbases>0) {
	print "$firstbases bases was removed from the 5'-end before trimming\n";
}
if($read2) {
	my $percent = 100*($pairReads/$totNoReads);
	$percent = sprintf "%.2f", $percent;
	print "$pairReads ($percent%) reads with $pairBases bases saved in pair files\n";
	$percent = 100*($unpairedReads/$totNoReads);
	$percent = sprintf "%.2f", $percent;
	print "$unpairedReads ($percent%) reads with $unpairedBases bases saved in unpaired file\n";
	print "  due to low quality of the other read in the pair\n";
}
else {
	my $percent = 100*($unpairedReads/$totNoReads);
	$percent = sprintf "%.2f", $percent;
	print "$unpairedReads ($percent%) reads with $unpairedBases bases saved\n";
}
if($printBad ne "") {
	print "$badReads removed reads was printed to a $printBad file\n";  
}
print "------------------------------------------------------------------\n";
$time = time-$time;
if($time<60) {
	print "Total time elapsed: $time sec.\n";
}
else {
	$time = int($time/60 + 0.5);
	print "Total time elapsed: $time min.\n";
}

#--------------------------------------------------------------------------------
# Subroutines

sub removeNs {
	my $seq = shift;
	my $qual = shift;
	
#	print "seq är $seq\n";

	my @s = split("", $seq);
	my @t = split("", $qual);
	my $cnt = 0;
	
	while(scalar(@s)>0 && $s[0] !~ m/[ATCG]/i) {
		shift(@s);
		shift(@t);
		$cnt++;	
	}
	$seq = join("", @s);
	$qual = join("", @t);	
	return ($seq, $qual, $cnt);
}

sub trimEnd { 

	my $seq = shift;
	my $qual = shift;
	my $HQ = shift;
	my $LQ = shift;
	my $maxHQ = shift;
	my $maxLQ = shift;
	my $len = shift;
	my $sc = shift;

	if($seq ne "") {
		my $LQ_flag = 0;
		my $HQ_warn = 0;
		my $LQinHQ_flag = "no";
		my ($qual_end, $seq_end) = ("","");

		my @t = split("", $qual);
		my @s = split("", $seq);

		while(scalar(@t)>$len && $HQ_warn<=$maxHQ) {

			if (ord($t[scalar(@t)-1])-$sc < $HQ) {
			
				if ($HQ_warn > 0 && $LQ_flag <= $maxLQ && ord($t[scalar(@t)-1])-$sc > $LQ) {
					$qual_end = pop(@t).$qual_end;
					$seq_end = pop(@s).$seq_end;
					$LQinHQ_flag = "yes";
					$LQ_flag++;

				}
				else {
					pop(@t);
					pop(@s);
					($qual_end, $seq_end) = ("","");
					$HQ_warn = 0;
				}
				
			}
			else {
				$qual_end = pop(@t).$qual_end;
				$seq_end = pop(@s).$seq_end;
				if($LQinHQ_flag eq "yes") {
					$HQ_warn = 1;
					$LQinHQ_flag = "no";
				}
				else {
					$HQ_warn++;
				}
			}
		}
		$seq = join("", @s);
		$qual = join("", @t);
		$seq .= $seq_end;
		$qual .= $qual_end;
	
	}
	return ($seq, $qual);
}
sub readOK {

	my $qual = shift;
	my $HQ = shift;
	my $LQ = shift;
	my $frac = shift;
	my $lfrac = shift;
	my $sc = shift;
	my $minlen = shift;
		
	if(length($qual)<$minlen) {
		return 0;
	}

	my @t = split("", $qual);

	my $score_cnt=0;
	my $bad_cnt=0;

	for (my $i=0; $i<scalar(@t)-1; $i++) {
		if(ord($t[$i])-$sc >= $HQ) {
			$score_cnt++;
		}
		if(ord($t[$i])-$sc < $lowlim) {
			if($lfrac==0) {
				return 0;
			}
			$bad_cnt++;
		}
	}
	if($score_cnt/scalar(@t) >= $frac && $bad_cnt/scalar(@t)<=$lfrac) {
		return 1;
	}
	else {
		return 0;
	}
}
sub table64 {
print "Char	ASCII	Char-64	P(error)
;	59	-5	0.7597
<	60	-4	0.7153
=	61	-3	0.6661
>	62	-2	0.6131
?	63	-1	0.5573
@	64	0	0.5000
A	65	1	0.4427
B	66	2	0.3869
C	67	3	0.3339
D	68	4	0.2847
E	69	5	0.2403
F	70	6	0.2008
G	71	7	0.1663
H	72	8	0.1368
I	73	9	0.1118
J	74	10	0.0909
K	75	11	0.0736
L	76	12	0.0594
M	77	13	0.0477
N	78	14	0.0383
O	79	15	0.0307
P	80	16	0.0245
Q	81	17	0.0196
R	82	18	0.0156
S	83	19	0.0124
T	84	20	0.0099
U	85	21	0.0079
V	86	22	0.0063
W	87	23	0.0050
X	88	24	0.0040
Y	89	25	0.0032
Z	90	26	0.0025
[	91	27	0.0020
\\	92	28	0.0016
]	93	29	0.0013
^	94	30	0.0010
_	95	31	0.0008
`	96	32	0.0006
a	97	33	0.0005
b	98	34	0.0004
c	99	35	0.0003
d	100	36	0.0003
e	101	37	0.0002
f	102	38	0.0002
g	103	39	0.0001
h	104	40	0.0001
";
}
sub table33 {
print "Char	ASCII	Char-33	P(error)
!	33	0	0.5000
\"	34	1	0.4427
#	35	2	0.3869
\$	36	3	0.3339
%	37	4	0.2847
&	38	5	0.2403
'	39	6	0.2008
(	40	7	0.1663
)	41	8	0.1368
*	42	9	0.1118
+	43	10	0.0909
,	44	11	0.0736
-	45	12	0.0594
.	46	13	0.0477
/	47	14	0.0383
0	48	15	0.0307
1	49	16	0.0245
2	50	17	0.0196
3	51	18	0.0156
4	52	19	0.0124
5	53	20	0.0099
6	54	21	0.0079
7	55	22	0.0063
8	56	23	0.0050
9	57	24	0.0040
:	58	25	0.0032
;	59	26	0.0025
<	60	27	0.0020
=	61	28	0.0016
>	62	29	0.0013
?	63	30	0.0010
@	64	31	0.0008
A	65	32	0.0006
B	66	33	0.0005
C	67	34	0.0004
D	68	35	0.0003
E	69	36	0.0003
F	70	37	0.0002
G	71	38	0.0002
H	72	39	0.0001
I	73	40	0.0001
J	74	41	0.0001
";
}
