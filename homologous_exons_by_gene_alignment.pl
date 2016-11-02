#!/usr/bin/env perl

use warnings;
use strict;
use MyPackage::Functions qw/test_if_file_exist/;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Getopt::Long;
use File::Temp qw/tempdir tempfile/;

my $minimalBlockNum = 3; # the minimal number of blocks in a sequence to analyze
my $alnFile;
my $structureFile;
my $outFile;
my $format;
my $minusStop;
my $ignoreAln;
my $cdsStartCol;

GetOptions(
		"aln-file=s"       => \$alnFile,
		"structure-file=s" => \$structureFile,
		"cds-start-column:i" => \$cdsStartCol, 
		"out-file:s"       => \$outFile,
		"format:s"         => \$format,
		"minus-stop!"      => \$minusStop,
		"ignore!"          => \$ignoreAln
);

&usage() unless($alnFile and $structureFile);

test_if_file_exist($outFile) if($outFile);

$outFile ||= '-';
$format  ||= 'clustalw';
$ignoreAln = 0 unless(defined $ignoreAln);
$minusStop = 1 unless(defined $minusStop);
my $stopLen = 3; # the length to be substracted from the sum of block lengths

# parse the structure file
my $structureParser = _parse_structure_file($structureFile);

# create a temporary directory for usage when necessary
my $tmpDir = tempdir("homologous_exons.$$.XXXXXXXXX", DIR => '/tmp', CLEANUP => 1);
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

warn scalar(@alnFiles)." alignment files are found\n";

my $counter = 0;
my $alignmentCounter = 0; # record the total alignment analyzed
my $successAlnCounter = 0;
my $pairId = 0; # the unique id for each homologous exon pair

open(OUT,"> $outFile") or die "can not open $outFile:$!";
print OUT "#produced by $0\n";
print OUT join("\t",qw/pair_id seq_5 exon_5 seq_3 exon_3 align_type identity aligned_length length_diff frames
                    aligned_bound bound_5 bound_3 exon_dist file aln_counter/),"\n";

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

		next unless($homoExonPairArrayRef); # no homolog exons found in this alignment
		
		# output the result
		foreach my $pairNode (@$homoExonPairArrayRef)
		{
			print OUT join("\t",@$pairNode,$baseFile,$alnCounter),"\n";
		}
		$successAlnCounter++;

	} # end of the alignment file
	
	$alignmentCounter += $alnCounter;

	warn "$counter files has been parsed\n" if(++$counter%100 == 0);
}


close OUT;

# remove temporary directory
warn "Removing the temporaty directory ...\n";
unlink <$tmpDir/*>;
rmdir($tmpDir) or warn "can not remove $tmpDir\n";

warn '*' x 60, "\nThe whole work is successfully finished\n";
warn '=' x 60, "\nTotally $counter files($alignmentCounter alignments [$successAlnCounter successes]) are analyzed\n";


exit 0;


sub find_homologous_exon_pairs
{
	my $aln          = shift;

	# find the exons aligned at same positions
	my ($alignedPosRef,$exonHash) = &aligned_positions($aln);
	return undef unless($alignedPosRef);

	# got the homologous exon pairs
	my $alignedPos5ss = $alignedPosRef->{'5ss'};
	my $alignedPos3ss = $alignedPosRef->{'3ss'};

	# merge the two results and get the final result
	my %homoExonPairs;
	while( my ($alnPos,$exonsRef) = each %$alignedPos5ss )
	{
		next unless($#$exonsRef > 0);# more than one exons exist
		for(my $i = 0; $i <= $#$exonsRef; $i++)
		{
			my $qExon = $exonsRef->[$i];
			for(my $j = $i+1; $j <= $#$exonsRef; $j++)
			{
				my $tExon = $exonsRef->[$j];
				$homoExonPairs{$qExon}->{$tExon} = '5ss';
			}
		}# end of a given position
	}# end of 5'ss mapping
	
	# update this by 3ss positions
	while( my ($alnPos,$exonsRef) = each %$alignedPos3ss )
	{
		next unless($#$exonsRef > 0);# need more than one exons for following analysis
		for(my $i = 0; $i <= $#$exonsRef; $i++)
		{
			my $qExon = $exonsRef->[$i];
			for(my $j = $i+1; $j <= $#$exonsRef; $j++)
			{
				my $tExon = $exonsRef->[$j];
				if(exists $homoExonPairs{$qExon}->{$tExon}) # the 5'ss exon boundary also align
				{
					$homoExonPairs{$qExon}->{$tExon} = 'both';
				}elsif(exists $homoExonPairs{$tExon}->{$qExon}) # the same, by another name
				{
					$homoExonPairs{$tExon}->{$qExon} = 'both';
				}else
				{
					$homoExonPairs{$qExon}->{$tExon} = '3ss';
				}
			}
		}# end of a given position
	}# end of the 3'ss mapping

	my @homoExonPairInfo;
	# get the final outputed result
	while( my ($qExon,$tExonRef) = each %homoExonPairs )
	{
		my ($qName,$qOrder,$qFrame,$qLength,$qAln5ssPos,$qAln3ssPos) =
		@{$exonHash->{'by_id'}->[$qExon]};
		while( my ($tExon,$alignType) = each %$tExonRef )
		{
			my ($tName,$tOrder,$tFrame,$tLength,$tAln5ssPos,$tAln3ssPos) =
			@{$exonHash->{'by_id'}->[$tExon]};
			# create a new alignment containing only these two sequences
			my $newAln = _get_sub_aln($aln,$qName,$tName);
			my $lengthDiff = _get_length_diff($qLength,$tLength);
			my $frameString;
			my $identity;
			my $alignedLen;
			my $alignedPos;
			my $aligned5Pos; # the other unaligned pos near 5'
			my $aligned3Pos; # the other unaligned pos near 3'
			my $distance; # the distance for the boundaries at the unaligned side, gap-only columns
			my $firstExon; # the exon with the other boundary close to the aligned boundary
			my $secondExon; # the exon with the other boundary distant to the aligned boundary
			# are excluded before calculation
			if($alignType =~ /both/i)
			{
				$frameString = "$qFrame,$tFrame";
				warn "Th----------------------------------\n" if($qAln3ssPos != $tAln3ssPos);
				$alignedPos  = $qAln3ssPos; # using the left aligned boundary
				$aligned5Pos = $qAln5ssPos; # this two should be equal
				$aligned3Pos = $tAln5ssPos;
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
			}elsif($alignType =~ /3ss/i)
			{
				$alignedPos  = $qAln3ssPos;
				if($qAln5ssPos <= $tAln5ssPos)
				{
					$aligned5Pos = $qAln5ssPos;
					$aligned3Pos  = $tAln5ssPos;
					$firstExon  = "$qName\t$qOrder";
					$secondExon = "$tName\t$tOrder";
					$frameString = "$qFrame,$tFrame";
				}else
				{
					$aligned5Pos = $tAln5ssPos;
					$aligned3Pos = $qAln5ssPos;
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
			}elsif($alignType =~ /5ss/i)
			{
				$alignedPos = $qAln5ssPos;
				if($qAln3ssPos <= $tAln3ssPos) # equal is actually impossible
				{
					$aligned5Pos = $qAln3ssPos;
					$aligned3Pos = $tAln3ssPos;
					$firstExon  = "$qName\t$qOrder";
					$secondExon = "$tName\t$tOrder";
					$frameString = "$qFrame,$tFrame";
				}else
				{
					$aligned5Pos = $tAln3ssPos;
					$aligned3Pos = $qAln3ssPos;
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
				warn "Unknow homologous exon pair type:$alignType\n";
				next;
			}

			# store the result
			push
			@homoExonPairInfo,[++$pairId,$firstExon,$secondExon,$alignType,
			                   $identity,$alignedLen,$lengthDiff,$frameString,$alignedPos,
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
	my %alnSplicePos;

	my $leftSeqNum = $aln->no_sequences;
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
			warn("The sum of block lengths($blockLength) in $id is not equal that in alignment ($seqLen)\n");
			last if($ignoreAln);
			$leftSeqNum--;
			next;
		}
		
		my $cumulativeLen = 0;
		my $blockNum = scalar(@blocks);
		for(my $i = 0; $i < $blockNum - 1; $i++)
		{
			$cumulativeLen += $blocks[$i];
			next unless($i > 0); # ignore the first block
			# which boundary to check, the 5' or 3'
			my $ss5Pos = $cumulativeLen;
			my $ss3Pos = $cumulativeLen - $blocks[$i] + 1;
			if($ss5Pos > $blockLength) # this occurs when the last exon/block only contain partial stop codon
			{
				$ss5Pos = $blockLength;
			}
			my $alnPos5ss;
			my $alnPos3ss;
			eval{
			 $alnPos5ss = $aln->column_from_residue_number($id,$ss5Pos);
			 $alnPos3ss = $aln->column_from_residue_number($id,$ss3Pos);
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
			}elsif($cdsOffSet > $cumulativeLen - $blocks[$i]) # start from this exon
			{
				$frame = 0;
			}else # start coding in preceeding exon
			{
				$frame = ($cumulativeLen - $blocks[$i] - $cdsOffSet + 1)%3;
			}
			
			# seqname, exonid, coding-frame, exon-length, 5ss pos, 3ss pos	
			$exonHash->{'by_id'}->[$uniqExonId] =
			[$id,$exonNum,$frame,$blocks[$i],$alnPos5ss,$alnPos3ss];
			$exonHash->{'by_name'}->{$id}->{$exonNum} = $uniqExonId;
			
			# record the aligned exon
			push @{$alnSplicePos{'5ss'}->{$alnPos5ss}}, $uniqExonId;
			push @{$alnSplicePos{'3ss'}->{$alnPos3ss}}, $uniqExonId;
			#$alnSplicePos->{'3ss'}->{$alnPos3ss}->{$id} = $uniqExonId;
			
		} # end for this sequence

		$success++;
	} # end for this alignment

	return undef unless($success >= 2);
	return \%alnSplicePos,$exonHash;
}

sub _get_length_diff
{
	my $len1 = shift;
	my $len2 = shift;

	return abs($len1 - $len2)/_max($len1,$len2);
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

	This script can determine the homologous exons from the multiple sequence alignment. The
	criteria for a pair of exons to be homologous:
	1) the exon pairs must be aligned at least in one boundary

	Meanwhile, these parameters will be calculated to easily filter the results[field]:
	a) the type of alignment, match at both boundaries(both), or 5'ss(5ss) or 3'ss(3ss)[align_type]
	b) the identity between the two exons, note the gaps are excluded in calculation
	c) the length of aligned region for the two exons without gaps[aligned_length]
	d) the difference of exon lengths, abs(diff)/longgest_exon_length[length_diff]
	e) the frames for the two exons[frames]
	f) the alignment position where the exons are aligned, for the 'both' type, the left(5') bound
	   chosen.[aligned_bound]
	g) another unaligned position in the alignment near 5' end, in fact this should be aligned
	   in the 'both' type, and the same as bound_3 field.[bound_5]
	h) another unaligned position in the alignment near 3' end, in fact this should be aligned   
	   in the 'both' type, and the same as bound_5 field.[bound_3]
	i) the distance between the two unaligned boundaries in the alignemnt. Note the gap-only columns
	   are excluded before counting.[exon_dist]

	Options:
	--aln-file: the file containing alignment files. This can be a directory, a tared file or a
	normal file.

	--structure-file: the file containing the structure information for each sequence in the
	alignment. The alignment without structure information will be ignored. The length of the
	sequence in alignments will be checked if it equals the sum of lengths in structure.
	A example file:
	/database/ftp/NCBI/genomes/Human/Build36.2/human_36.2_RefSeq_translated_protein.info.slim

	--cds-start-colunmn: which column in the structure above is the CDS start position for each
	transcript. If not provided, it will assume the first position in each sequence is the CDS
	start, this is usually the case for the CDS sequence. Original columin is 0.

	--ignore: a switch option indicative to ignore the alignment in which if there is at least one
	sequence with only one exon block or this sequence does not have structure information. If this
	option is false, then only this sequence is ignored, and the left sequences are analyzed. The
	default is false.
	
	--format:  the multiple sequence alignment format, default is clustalw.
	
	--out-file: the file to store the result, the default is to print to stdout.
	
	--minus-stop: a switch option indicative whether it is necessary fo substracting stop codon
	length from the sum of block lengths. This is usually true for the CDS sequence, because they
	include the stop sequence in it. The alignment in which the length of a sequence does not match
	the sum of blocks in the structure file will be ignored. The default value is true.
	
	Notifications:
	1) the two marginal exons/blocks will not be considered in analysis.
	2) in the result, there may be a certain exon matched with another two exons of the same
	   transcript. There are at least two causes for this, one is that the two exons are merged
	   (loss of the intron between them), and the other is loss or splicing out of one exon in the
	   former transcript. One method to distinguish the two is compare the exon lengths among them.
	   If the former case is true, we should find that the sum of two exons (nearly) equal to the
	   exon matched. If the latter case is true, then only one of the two exons in the same
	   transcript equal to the matched exon in another transcript.
	
	Author:  Zhenguo Zhang
	Contact: fortunezzg\@gmail.com

	Last-modified: 2008-10-15

USAGE

	exit 1;
}

sub _parse_structure_file
{
	my $file = shift;
	
	my $sep = "\t";
	my %hash;
	
	open(IN,"< $file") or die "can not open $file:$!";
	while(<IN>)
	{
		next if(/^#/);
		my @fields = split($sep,$_);
		next unless(scalar(@fields) > 2);
		# get the CDS start offset in the sequence
		my $cdsOffset = 1; # default is the first position of sequence
		$cdsOffset = $fields[$cdsStartCol] if(defined $cdsStartCol);
		# mRNA acc as key, the start exon and block strings as value
		$hash{lc($fields[0])} = [@fields[1,2],$cdsOffset];
	}
	close IN;

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
