#!/usr/bin/env perl
use strict;
use Getopt::Long;

my $sep = "\t";
my $IP = '[INFO]';
my @args = @ARGV;
my $type;
my $minIdent;
my $minALen;
my $minALenFrac;
my $maxLenDiff;
my $sameFrame;
my $maxDist;
my $outFile;

GetOptions(
	"type:s"      => \$type,
	"identity:f"  => \$minIdent,
	"alen:i"      => \$minALen,
	"alen-frac:f" => \$minALenFrac,
	"len-diff:f"  => \$maxLenDiff,
	"same-frame!" => \$sameFrame,
	"bound-dist:i"=> \$maxDist,
	"out:s"       => \$outFile
);

my $inFile = shift or &usage();

#set default values
$type = lc($type) if($type);
$minIdent = 0 unless($minIdent);
$minALen  = 0 unless($minALen);
$minALenFrac = 0 unless($minALenFrac);
$maxLenDiff = 1 unless(defined $maxLenDiff);
$maxDist = 10 unless(defined $maxDist);
$outFile = '-' unless($outFile);

open(I,"< $inFile") or die "cannot open $inFile:$!";
open(O,"> $outFile") or die "cannot open $outFile:$!";
print O "#Produced by $0 @args\n";
print O join($sep, qw/#qExon tExon reciprocal pair_id/),"\n";
my $preFile = '';
my $preAln = '';
my %bestExonPairs;
my $counter = 0;
while(<I>)
{
	next if /^#/ or /^pair_id$sep/;
	chomp;
	my @fields = split $sep;
	# check if a new alignment
	if($fields[15] ne $preFile or $fields[16] ne $preAln)
	{
		output() if(%bestExonPairs);
		$preFile = $fields[15];
		$preAln  = $fields[16];
		%bestExonPairs = ();
	}

	# filter the pairs
	next if $type and $fields[5] ne $type;
	next if $fields[6] < $minIdent;
	next if $fields[7] < $minALen;
	next if $fields[8] < $minALenFrac;
	next if $fields[9] > $maxLenDiff;
	next if $fields[14] > $maxDist;
	next if $sameFrame and not_same_frame($fields[10]);
	# store/update the pair
	my $tmp = $fields[7]*$fields[6]; # number of identical bases
	my $qT = $bestExonPairs{$fields[1]}->{$fields[3]} || [];
	my $tQ = $bestExonPairs{$fields[3]}->{$fields[1]} || [];
	if(!exists($qT->[$fields[2]])
		or $qT->[$fields[2]]->[1] < $tmp)
	{
		$bestExonPairs{$fields[1]}->{$fields[3]}->[$fields[2]] =
		[$fields[4], $tmp, $fields[0]];
	}
	if(!exists($tQ->[$fields[4]]) or $tQ->[$fields[4]]->[1] < $tmp )
	{
		$bestExonPairs{$fields[3]}->{$fields[1]}->[$fields[4]] =
		[$fields[2], $tmp, $fields[0]];
	}

	warn "$IP $counter exon pairs have been recorded\n"
	if(++$counter % 10000 == 0); 
}
output() if(%bestExonPairs);
close O;
close I;

warn "Work done!\nIn total $counter exon pairs have been recorded\n";

exit 0;

# output the selected best exons
sub output
{
	my @seqNames = sort {$a cmp $b} keys %bestExonPairs;
	# output results for each pair
	for(my $i = 0; $i < $#seqNames; $i++)
	{
		my $q = $seqNames[$i];
		for(my $j = $i+1; $j <= $#seqNames; $j++)
		{
			my $t = $seqNames[$j];
			my $qT=$bestExonPairs{$q}->{$t};
			my $tQ=$bestExonPairs{$t}->{$q};
			# loop over each exon in $q and $t
			# $q first
			my $targetExon;
			for(my $i = 1; $i <= $#$qT; $i++)
			{
				next unless($qT->[$i]); # no homolog for this exon
				$targetExon = $qT->[$i]->[0];
				my $reciprocal = 0;
				if(exists($tQ->[$targetExon]->[0]) and
					$tQ->[$targetExon]->[0] == $i)
				{
					$reciprocal = 1;
					# also delete the target exon from the tQ relation
					$tQ->[$targetExon] = undef;
				}

				print O join($sep, "$q:$i", "$t:$targetExon",
					$reciprocal, $qT->[$i]->[2]),"\n";
			}
			# $t now, reciprocal ones have been deleted, so only
			# non-reciprocal ones left
			for(my $i = 1; $i <= $#$tQ; $i++)
			{
				next unless($tQ->[$i]); # no exon or it is reciprocal
				$targetExon = $tQ->[$i]->[0];
				print O join($sep, "$t:$i", "$q:$targetExon",
					0, $tQ->[$i]->[2]),"\n";
			}
		}
	}

	return 1;
}

sub not_same_frame
{
	my $frames = shift;

	my ($f1, $f2) = split ',', $frames;

	return 1 if $f1 ne $f2;
	return 0;
}

sub usage
{
	print <<USAGE;
Usage: [options] <homologous-exons-file>

This program selects the best homologous exon for each exon. The file
<homologous-exons-file> should be the output from
homologous_exons_by_gene_alignment.pl, which has been grouped by
alignment (i.e., the sequences from the same alignment are in
contiguous lines).

Options:

Here are criteria for a pair of exons to regard homologous, pairs
failing any of these will be disgarded in output:

--type: the alignment type, 5 for aligned at 5' ends, 3 for 3' ends,
and 'both' for both ends. Default is any type of alignment

--identity: the minimal identity for the aligned region, in
percentage, i.e., range from 0 to 100. Default is 0.

--alen: the minimal length of aligned region. Default is 0.

--alen-frac: the minimal fraction of aligned region, calculated as
aligned-region-len/length-of-longer-exon. Deault is 0.

--len-diff: the maximal difference of lengths between the two exons in
fraction (ranging from 0 to 1). Default is 1.

--same-frame: a switch option, requiring the two exons have the same 
translation frames if provided.

--bound-dist: the maximal distance (in bps) between the two exons at the
unaligned end (for type is 5 or 3 only) in the alignment. Default is
10.

Output options:

--out: the file to store/display output. Default is screen.

Output format:

qExon   tExon    reciprocal  pair_id
seq1:1  seq2:1     1          1
seq1:2  seq2:2     1          2
seq2:3  seq1:2     0          4
seq1:3  seq2:4     1          7
... ...

The output contains four fields separated by <tab>, 'qExon' and
'tExon' are the homologous exon pair, each having the format
'sequence-name:exon-number'. The third field 'reciprocal', when = 1,
denotes that the best homologous exon of the tExon in the homologous
sequence is the qExon, i.e., reciprocally best. When this is true,
note that we do not list the best exon for the tExon in another line
in order to save space. The last field is the unique exon pair id
direcly copied from the input file.

Created:
Thu Jan 19 22:07:11 EST 2017

USAGE

	exit 1;
}
