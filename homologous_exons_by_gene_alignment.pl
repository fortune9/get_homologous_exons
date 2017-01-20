#!/usr/bin/env perl

use warnings;
use strict;
#use MyPackage::Functions qw/test_if_file_exist/;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Getopt::Long;
use File::Temp qw/tempdir tempfile/;

my $IP = '[INFO]'; # info prefix
my $WP = '[WARN]'; # warning prefix
my $EP = '[ERROR]'; # error prefix
my $minimalBlockNum = 1; # the minimal number of blocks in a sequence to analyze
my $alnFile;
my @structureFiles;
my $outFile;
my $format;
my $minusStop;
my $ignoreAln;
my $cdsStartCol;
my $noEnds;

GetOptions(
		"aln-file=s"       => \$alnFile,
		"format:s"         => \$format,
		"structure-file=s@" => \@structureFiles,
		"ignore!"          => \$ignoreAln,
		"out-file:s"       => \$outFile,
		"minus-stop!"      => \$minusStop,
		"no-ends!" => \$noEnds,
		"cds-offset-col:i" => \$cdsStartCol
);

&usage() unless($alnFile and @structureFiles);

#test_if_file_exist($outFile) if($outFile);

$outFile ||= '-';
$format  ||= 'clustalw';
#$ignoreAln = 0 unless(defined $ignoreAln);
#$minusStop = 1 unless(defined $minusStop);
my $stopLen = 3; # the length to be substracted from the sum of block lengths

# parse the structure file
my $structureParser = _parse_structure_file(@structureFiles);

# create a temporary directory for usage when necessary
#my $tmpDir = tempdir("homologous_exons.$$.XXXXXXXXX", DIR => '/tmp', CLEANUP => 1);
my $tmpDir = tempdir("homologous_exons.$$.XXXXXXXXX", CLEANUP => 1);
# got the alignment file names
my @alnFiles;
if(-d $alnFile) # a directory
{
	@alnFiles = <$alnFile/*>;
}else
{
	if($alnFile =~ /\.tar\.gz$/i)
	{
		system("tar -xzf $alnFile -C $tmpDir") and die "extracttion failed:$!";
		@alnFiles = <$tmpDir/*>;
	}elsif($alnFile =~ /\.tar$/i)
	{
		system("tar -xf $alnFile -C $tmpDir") and die "extracttion failed:$!";
		@alnFiles = <$tmpDir/*>;
	}else # just a simple file
	{
		push @alnFiles, $alnFile;
	}
}

warn "$IP ".scalar(@alnFiles)." alignment files are found\n";

my $counter = 0;
my $alignmentCounter = 0; # record the total alignment analyzed
my $successAlnCounter = 0;
my $pairId = 0; # the unique id for each homologous exon pair

open(OUT,"> $outFile") or die "can not open $outFile:$!";
print OUT "#produced by $0\n";
print OUT join("\t",qw/pair_id seq_5 exon_5 seq_3 exon_3 align_type 
	identity aligned_length aligned_length_frac length_diff frames
    aligned_bound bound_5 bound_3 exon_dist file aln_counter/),"\n";
warn "$IP finding homologous exons ...\n";
foreach my $file (@alnFiles)
{
	my $alnIn = Bio::AlignIO->new(-file => $file, -format => $format);
	warn "Can not open alignemnt $file\n" and next unless($alnIn);
	
	my $baseFile = $file;
	$baseFile =~ s/.*\///;
	
	# find the homologous splice junctions
	my $alnCounter = 0;
	while(my $aln = $alnIn->next_aln)
	{
		$alnCounter++;
		# get the homologous exon pairs with aligned criterion
		my ($homoExonPairArrayRef) = &find_homologous_exon_pairs($aln);

		unless($homoExonPairArrayRef)
		{
			warn "$WP no homolog exons found in $file [$alnCounter]\n";
			next;
		}
		
		# output the result
		foreach my $pairNode (@$homoExonPairArrayRef)
		{
			print OUT join("\t",@$pairNode,$baseFile,$alnCounter),"\n";
		}

		warn "$IP $successAlnCounter alignments have been processed\n" 
		if(++$successAlnCounter % 100 == 0);
	} # end of the alignment file
	
	$alignmentCounter += $alnCounter;
	$counter++;
}


close OUT;

# remove temporary directory
warn "$IP Removing the temporaty directory ...\n";
unlink <$tmpDir/*>;
rmdir($tmpDir) or warn "$WP can not remove $tmpDir\n";

warn '*' x 60, "\nThe whole work is successfully finished\n";
warn '=' x 60, "\nTotally $counter files($alignmentCounter alignments [$successAlnCounter successes]) are processed\n";

exit 0;

sub find_homologous_exon_pairs
{
	my $aln          = shift;

	# find the exons aligned at same positions
	my ($alignedPosRef,$exonHash) = &aligned_positions($aln);
	return undef unless($alignedPosRef);

	# got the exons aligned at one end
	my $alignedExons5 = $alignedPosRef->{'5'}; # exons aligned at 5' end
	my $alignedExons3 = $alignedPosRef->{'3'}; # exons aligned at 3' end

	# merge the two results to get final result, particularly for
	# exons aligned at both ends
	my %homoExonPairs;
	while( my ($alnPos,$exonsRef) = each %$alignedExons5 )
	{
		next unless($#$exonsRef > 0);# more than one exons exist
		for(my $i = 0; $i < $#$exonsRef; $i++)
		{
			my $qExon = $exonsRef->[$i];
			for(my $j = $i+1; $j <= $#$exonsRef; $j++)
			{
				my $tExon = $exonsRef->[$j];
				$homoExonPairs{$qExon}->{$tExon} = '5';
			}
		}# end of a given position
	}# end of 5' end alignment
	
	# update this by 3' end positions
	while( my ($alnPos,$exonsRef) = each %$alignedExons3 )
	{
		next unless($#$exonsRef > 0);# need more than one exons for following analysis
		for(my $i = 0; $i < $#$exonsRef; $i++)
		{
			my $qExon = $exonsRef->[$i];
			for(my $j = $i+1; $j <= $#$exonsRef; $j++)
			{
				my $tExon = $exonsRef->[$j];
				if(exists $homoExonPairs{$qExon}->{$tExon}) # the 5' end also align
				{
					$homoExonPairs{$qExon}->{$tExon} = 'both';
				}elsif(exists $homoExonPairs{$tExon}->{$qExon}) # the same, by another name order
				{
					$homoExonPairs{$tExon}->{$qExon} = 'both';
				}else # only aligned at 3' end
				{
					$homoExonPairs{$qExon}->{$tExon} = '3';
				}
			}
		}# end of a given position
	}# end of the 3' end alignment

	my @homoExonPairInfo;
	# get the final outputed result
	while( my ($qExon,$tExonRef) = each %homoExonPairs )
	{
		my ($qName,$qOrder,$qFrame,$qLength,$qAlnPos5,$qAlnPos3) =
		@{$exonHash->{'by_id'}->[$qExon]};
		while( my ($tExon,$alignType) = each %$tExonRef )
		{
			my ($tName,$tOrder,$tFrame,$tLength,$tAlnPos5,$tAlnPos3) =
			@{$exonHash->{'by_id'}->[$tExon]};
			# create a new alignment containing only these two sequences
			my $newAln = _get_sub_aln($aln,$qName,$tName);
			my $lengthDiff = _get_length_diff($qLength,$tLength);
			my $frameString;
			my $identity;
			my $alignedLen;
			my $alignedPos; # alignment position where exons are aligned
			my $aligned5Pos; # the position of the other unaligned end near 5'
			my $aligned3Pos; # the position of the other unaligned end near 3'
			my $distance; # the distance for the boundaries at the unaligned ends, 
			# gap-only columns are excluded before calculation
			my $firstExon; # the exon in homologous pair that is
			               # closer 5' end
			my $secondExon; # the exon in homologous pair that is
			                # farther to 5' end
			if($alignType =~ /both/i)
			{
				$frameString = "$qFrame,$tFrame";
				warn "$EP The 3' ends between $qName:$qOrder and"
				." $tName:$tOrder don't align\n" if($qAlnPos3 != $tAlnPos3);
				$alignedPos  = $qAlnPos5; # using the left aligned boundary
				$aligned5Pos = $qAlnPos3; # this and next should be equal
				$aligned3Pos = $tAlnPos3;
				$firstExon  = "$qName\t$qOrder"; # the proximal exon first
				$secondExon = "$tName\t$tOrder";
				my $subAln;
				eval{
				   $subAln = $newAln->slice($alignedPos,$aligned5Pos);
				};
				if($@)
				{
					warn $@;
					next;
				}
				$identity = $subAln->average_percentage_identity;
				$distance = 0;
				$alignedLen = $subAln->remove_gaps->length;
			}elsif($alignType =~ /5/i) # aligned at 5' ends
			{
				$alignedPos  = $qAlnPos5;
				if($qAlnPos3 <= $tAlnPos3) # q is before t
				{
					$aligned5Pos = $qAlnPos3;
					$aligned3Pos = $tAlnPos3;
					$firstExon  = "$qName\t$qOrder";
					$secondExon = "$tName\t$tOrder";
					$frameString = "$qFrame,$tFrame";
				}else # q is after t
				{
					$aligned5Pos = $tAlnPos3;
					$aligned3Pos = $qAlnPos3;
					$firstExon  = "$tName\t$tOrder";
					$secondExon = "$qName\t$qOrder";
					$frameString = "$tFrame,$qFrame";
				}
				my $subAln;
				eval {
				      $subAln = $newAln->slice($alignedPos,$aligned5Pos);
				};
				if($@)
				{
					warn $@;
					next;
				}
				$identity = $subAln->average_percentage_identity;
				my $distanceAln = $newAln->slice($aligned5Pos,$aligned3Pos);
				$distance = $distanceAln->length - 1;
				$alignedLen = $subAln->remove_gaps->length;
			}elsif($alignType =~ /3/i) # aligned on the 3'ends
			{
				$alignedPos = $qAlnPos3;
				if($qAlnPos5 <= $tAlnPos5) # q is before t
				{
					$aligned5Pos = $qAlnPos5;
					$aligned3Pos = $tAlnPos5;
					$firstExon  = "$qName\t$qOrder";
					$secondExon = "$tName\t$tOrder";
					$frameString = "$qFrame,$tFrame";
				}else # t is before q
				{
					$aligned5Pos = $tAlnPos5;
					$aligned3Pos = $qAlnPos5;
					$firstExon  = "$tName\t$tOrder";
					$secondExon = "$qName\t$qOrder";
					$frameString = "$tFrame,$qFrame";
				}
				my $subAln;
				eval {
					$subAln = $newAln->slice($aligned3Pos,$alignedPos);
				};
				if($@)
				{
					warn $@;
					next;
				}
				$identity = $subAln->average_percentage_identity;
				my $distanceAln = $newAln->slice($aligned5Pos,$aligned3Pos);
				$distance = $distanceAln->length - 1;
				$alignedLen = $subAln->remove_gaps->length;
			}else # this should not happen
			{
				warn "$EP Unknow homologous exon pair type:$alignType".
				" for $qName:$qOrder vs $tName:$tOrder\n";
				next;
			}

			# store the result
			my $alignedFrac = sprintf("%.5f",
				$alignedLen/_max($qLength,$tLength));
			push
			@homoExonPairInfo,[++$pairId,$firstExon,$secondExon,$alignType,
			                   sprintf("%.5f", $identity),$alignedLen,
							   $alignedFrac,$lengthDiff,$frameString,$alignedPos,
							   $aligned5Pos,$aligned3Pos,$distance];

		}# end of inner while
	}# end of outer while

	return \@homoExonPairInfo;
}

sub _get_sub_aln
{
	my $aln = shift;

	my $newAln = Bio::SimpleAlign->new();
	
	foreach my $id (@_)
	{
		my ($seq) = $aln->each_seq_with_id($id);
		$newAln->add_seq($seq);
	}
	
	return $newAln;
}

sub aligned_positions
{
	my $aln = shift;
	
	my $exonHash = {};
	my %alignedEndPos;

	my $leftSeqNum = $aln->num_sequences;
	my $uniqExonId = 0;
	my $success = 0; # record the number of sequence successfully analyzed
	foreach my $seq ($aln->each_seq)
	{
		last if($leftSeqNum < 2); # less than two sequences available for this alignment
		my $id = $seq->id();
		my $spliceInfo = $structureParser->{lc($id)};
		unless($spliceInfo) # no structure information
		{
			last if $ignoreAln;
			$leftSeqNum--;
			next;
		}

		my ($exonOffset,$blockString,$cdsOffSet) = @$spliceInfo;
		my @blocks = split(',',$blockString);
		if(scalar(@blocks) < $minimalBlockNum) # only one block
		{
			last if($ignoreAln);
			$leftSeqNum--;
			next;
		}
		
		# check the length for this sequence
		my $blockLength = _calculate_block_length(@blocks);
		my $seqLen = $seq->end - $seq->start + 1;
		if($seqLen != $blockLength)
		{
			warn("$EP The sum of block lengths($blockLength) in $id is not equal that in alignment ($seqLen)\n");
			last if($ignoreAln);
			$leftSeqNum--;
			next;
		}
		
		# get the aligned positions for each exon
		my $cumulativeLen = 0;
		my $blockNum = scalar(@blocks);
		for(my $i = 0; $i < $blockNum; $i++)
		{
			$cumulativeLen += $blocks[$i];
			next if($i == $blockNum-1 and $minusStop and $blocks[$i]
				<= $stopLen); # ignore last block if it contains
			# (part) stop codon only
			# ignore two blocks/exons at two ends if necessary
			next if($noEnds and ($i == 0 or $i == $blockNum-1 ));
			# which boundary to check, the 5' or 3'
			my $end3 = $cumulativeLen; # 3' end
			my $end5 = $cumulativeLen - $blocks[$i] + 1; # 5'end
			if($end3 > $blockLength) # this occurs when the last
			# exon/block only contains part of stop codon, i.e., shorter
			# than 3 nts
			{
				$end3 = $blockLength;
			}
			my $alnPos5; # 5' end of the exon
			my $alnPos3; # 3' end
			eval{
			 $alnPos5 = $aln->column_from_residue_number($id,$end5);
			 $alnPos3 = $aln->column_from_residue_number($id,$end3);
			};
			if($@)
			{
				warn $@;
				next;
			}
	
			my $exonNum = $exonOffset + $i;
			$uniqExonId++;
			my $frame;
			if($cdsOffSet > $cumulativeLen) # no coding region in this exon
			{
				$frame = -1;
			}elsif($cdsOffSet > $cumulativeLen - $blocks[$i]) #CDS starts from this exon
			{
				$frame = 0;
			}else # start coding in preceeding exon
			{
				$frame = ($cumulativeLen - $blocks[$i] - $cdsOffSet +
					1)%3; # number of codon bases in preceding exon
			}
			
			# seqname, exonid, coding-frame, exon-length, 5ss pos, 3ss pos	
			$exonHash->{'by_id'}->[$uniqExonId] =
			[$id,$exonNum,$frame,$blocks[$i],$alnPos5,$alnPos3];
			# $exonHash->{'by_name'}->{$id}->{$exonNum} = $uniqExonId;
			
			# record the aligned exon
			push @{$alignedEndPos{'5'}->{$alnPos5}}, $uniqExonId;
			push @{$alignedEndPos{'3'}->{$alnPos3}}, $uniqExonId;
			#$alignedEndPos->{'3ss'}->{$alnPos3ss}->{$id} = $uniqExonId;
			
		} # end for this sequence

		$success++;
	} # end for this alignment

	return undef unless($success >= 2);
	return \%alignedEndPos,$exonHash;
}

sub _get_length_diff
{
	my $len1 = shift;
	my $len2 = shift;

	return sprintf("%.5f", abs($len1 - $len2)/_max($len1,$len2));
}

sub _max
{
	my @els = @_;

	my $max;

	foreach my $el (@els)
	{
		$max = $el unless(defined $max);
		$max = $el if($max < $el);
	}

	return $max;
}

sub usage
{
	print <<USAGE;

	Usage: $0 <--option=value> ...

	This program finds the homologous exons based on sequence
	alignments and exon-intron structure information.

	The criteria for determining homologous exons are:
	1) homologous exon pairs must be aligned on at least one end (5'
	   or 3' or both)

	The output for each homologous exon pair will include the
	following parameters (field name is shown in []):
	a) matching type in the alignment, at both ends (both), just 5' (5) 
	   or 3' end (3) [align_type]
	b) the sequence identity between the two exons after excluding
	   gaps in alignment [identity]
	c) the length of aligned region for the two exons after excluding gaps 
	   [aligned_length]
	d) the length difference fraction between the exon pair, calculated as 
	   abs(diff)/longer_exon_length [length_diff]
	e) the translation frames (if any) for the two exons. The frame
	   denotes the number bases included in precedding exon, so
	   frame=0, =1, =2 means that the exon starts with a complete
	   codon, a codon with 2 bases, and a codon with 1 base,
	   respectively [frames]
	f) the position in the alignment where the two exons are aligned, for 
	   the 'both' type, the left (5') end is reported. [aligned_bound]
	g) the position of the other (un)aligned exon ends in the
	   alignment which is near 5' end. For 'align_type' is 'both',
	   this end is aligned and the value in this field will equal that
	   in the field 'bound_3'. [bound_5]
	h) similar to the field 'bound_5', but mark one of the two exon
	   ends distant from 5' end. [bound_3]
	i) the distance between the two (un)aligned exon ends in alignment. If
	   'align_type' is 'both', the distance is 0. Note the gap-only columns
	   are excluded before calculation. [exon_dist]

	Options:
	--aln-file: the file/directory containing sequence alignments.
	Valid values for this option include a directory, a tared file or
	a single file containing many alignments. Note the sequence names
	in alignments must match those in the 'structure-file' (below).
	
	--format:  sequence alignment format, default is clustalw.

	--structure-file: the file containing exon-intron structure
	information for each sequence. Sequences
	without structure information will be skipped. See below for the
	requirement on the format of the structure-file. This option can
	be specified multiple times to input the structure information
	distributed in multiple files.

	--cds-offset-col: at which column in the structure-file it
	specifies the cds start in each sequence. If not provided, the
	first base in each sequence is assumed as the start of CDS. Note
	the first column should be specified as 0, second as 1, etc. The
	information is used to calculate translation frames for each exon.

	--ignore: a switch option. If provided, sequence alignments
	containing any sequence that has only one exon or has no structure
	information will be wholy skipped. In default, only the relevant 
	sequences will be ignored and the remaining sequences in the
	alignment will be used.
	
	--out-file: the file to store the result. In default, it is screen.
	
	--minus-stop: a switch option. If provided, it indicates that the
	last exon for a sequence should minus 3 bases to trim the stop
	codon. This is usually the case for CDS alignment, so don't forget
	this option when necessary.

	--no-ends: a switch option. If provided, the most 5' and 3' exons
	for each sequence will be ignored.
	
	Format of structure-file:
	#seq_name  number_of_first_exon  exon_sizes  cds-start
	seq1         2                   100,200,150    50
	seq2         1                   300,100,150    1
	...		...

	The structure-file should contain at least 3 fields separated by
	<tab>. The fields after the first 3 will be ignored. Optionally,
	another field specified by the option --cds-offset-col may be
	used. The first three fields are: sequence name, the exon number
	of first exon included in the alignment (e.g., 1 if the first exon
	of the sequence is in the alignment), and the sizes of each
	exons.

	Notes:
	1) in the result, one exon of sequence A may match more than one 
	   exon in another sequence B. This can occur when the two exons
	   in sequence A lost the intron in-between during evolution.
	
	Author:  Zhenguo Zhang
	Contact: zhangz.sci\@gmail.com

	Last-modified: 2017-01-17

USAGE

	exit 1;
}

# a function to read into the exon-intron structure information
sub _parse_structure_file
{
	my $sep = "\t";
	my %hash;
	foreach my $file (@_)
	{
		open(IN,"< $file") or die "can not open $file:$!";
		while(<IN>)
		{
			next if(/^#/);
			my @fields = split($sep,$_);
			next unless(scalar(@fields) > 3);
			# get the CDS start offset in the sequence
			my $cdsOffset = 1; # default is the first position of sequence
			$cdsOffset = $fields[$cdsStartCol] if(defined $cdsStartCol);
			# mRNA acc as key, the exon number of first exon in the
			# alignment, a string of exon block sizes, and the start
			# of cds in the included sequence part as value
			$hash{lc($fields[0])} = [@fields[1,2],$cdsOffset];
		}
		close IN;
	}
	return \%hash;
}

sub _calculate_block_length
{
	my @array = @_;

	my $length = 0;

	foreach my $el (@array)
	{
		$length += $el;
	}

	$length -= $stopLen if($minusStop);

	return $length;
}
