#! /usr/bin/perl

#Tristan LEFEBURE

#1june2022
# - add option to only select a section of the reference sequence
# - to do so modified hoxw the revcom seq are handle: they are first, instead than at the end
# - add the blastx

#9mai2022
# - clean up code

#23fev2016
#- allow to switch between query or hit extraction

#25nov2016
#-change colonm selection
#-

use strict;
use warnings;
use Bio::Index::Fasta;
use Getopt::Long;
use YAML;
# use List::Util qw[min max];

#options
my $min_length = 100;
my $help;
my $extract = 'query';
my @region;
my $blastx;
my $yaml;

GetOptions(
	"length=i" => \$min_length,
	"help|?" => \$help,
	"extract=s" => \$extract,
	"region=s{2,2}" => \@region,
	"yaml=s" => \$yaml,
	"blastx" => \$blastx
);

my $usage = 
"
A script to mass extract sequences given BLAST alignments to a 
reference sequence

Usage: $0 <blast report file> <fasta index> <output>

with:
  -blast report file: comes from a regular blast with -outfmt 6
  -fasta index: produced by index_fasta.pl

Options:
  -help 

  -length: minimum length for a hit to be included, default=$min_length.

  -extract <query|target> : what do you want to extract ? [$extract]

  -region 100 900 : extract the hits in that given region, it should
	            be in nucleotide coordinates if you used blastx. 
		    Be carefull: if there is indels up or downstream
		    you will extract more/less than expected. You may
		    want to enlarge the region to avoid missing
		    sequences

  -blastx, if using blastx the query is in nucleotide space,
           the target ans the length in aa space. Using this option
           puts everybody in nucleotide space

  -yaml <FILE> : print the sequence to ba extracted coordinates
           into a YAML file
\n";


#declaration
my %coord; #AN: start->dd, end->dd
my %seq_per_family; #family->nbr de seq

if ($#ARGV<1) {
	print $usage;
	exit;
}

###set the column numbers

my $pname; my $pstart; my $pstop; 
my $pRstart; my $pRstop;

if($extract eq 'query') {
	$pname = 0;
	$pstart = 6;
	$pstop = 7;
	$pRstart = 8;
	$pRstop = 9;
	
}
elsif($extract eq 'target') {
	$pname = 1;
	$pstart = 8;
	$pstop = 9;
	$pRstart = 6;
	$pRstop = 7;


} else {
	die "I don't understand the option extract : $extract\n";
}


#the index
my $inx = Bio::Index::Fasta->new($ARGV[1]) or die "Can't open the index file\n";
my $out = Bio::SeqIO->new(-file => ">$ARGV[2]", -format => 'fasta');

# Parse the pairwise alignment


# VOCABULARY:
#
# ref:     ----------------------------
#          |  |    |        |    |
#          1  m0   R0       R1   m1
   
# region:          |--------|
# target:     --------------------
#             |    |        |    |
#             t0   t0p      t1p  t1
#
# m: match position
# R: region
# t: target
#
# If the region is included in the target sequence:
# t0p = t0 + (R0 - m0)
# t1p = t1 - (m1 - R1)

open BLAST, "$ARGV[0]" or die "Can't open $ARGV[1]\n";
while (<BLAST>) {
	chomp;
	my @cols = split /\t/;

	# Change the scale if this is a blastx reaport
	if($blastx) {
		#move the length and the target
		$cols[3] = $cols[3] * 3;
		$cols[8] = $cols[8] * 3 - 2;
		$cols[9] = $cols[9] * 3;
	}
	

	my $length = $cols[3];
	next if($length < $min_length);

	my $t0 = $cols[$pstart];
	my $t1 = $cols[$pstop];
	my $id = $cols[$pname];	 
	my $m0 = $cols[$pRstart];
	my $m1 = $cols[$pRstop];

	# First, handle the revcom
	if($t0 > $t1) {
		# reverse the coordinates
		# I need the length... and need to get the sequence to get it...
		my $seq;
        	eval{ $seq = $inx->fetch($id) };
        	if ($@) {
                	print "Problem with $id in the index, not included\n$@\n";
                	next;
       		 }
		my $length = $seq->length;
		$t0 = $length - $t0 +1;
		$t1 = $length - $t1 +1;
		# and store the strand
		$coord{$id}{strand} = '-';
	} else {
		$coord{$id}{strand} = '+';
	}


	# If no sub-regions are specified,
	# just store the hit coordinates,
	# and the strand 

	if($#region < 1) {
		$coord{$id}{start} = $t0; 
		$coord{$id}{end} = $t1;

	} else {

		# A sub-region is declared
		my $R0 = $region[0];
		my $R1 = $region[1];
	
		# skip if the hit does not overlap the region
		next if $m0 > $R1;

		my $t0p; my $t1p;
		# calculate the left coordinate:	
		# if there if something to remove on the left
		if($R0 > $m0) {	
			$t0p = $t0 + ($R0 - $m0);	
		} else { 
			$t0p = $t0;
		}

		# same for the right coordinates
		if($R1 < $m1) {
			$t1p = $t1 - ($m1 - $R1);
		} else {
			$t1p = $t1;
		}
		
		# The sequence may now be too small
		my $newlength = $t1p - $t0p;
		next if $newlength < $min_length;

		# Store the coordinates:
		$coord{$id}{start} = $t0p; 
		$coord{$id}{end} = $t1p;
	}
}

# Dump the hash to a YAML file
if($yaml) {
	open OUT, ">$yaml";
	print OUT Dump(%coord);
}



#for each hit, cut and export
print "Exporting the sequences:\n";
my $n = 0;
foreach my $id (keys %coord) {
	
	# Some hit have a strand but no coordinates
	# because they did not contain the region,
	# skip them

	next unless( exists $coord{$id}{start} ); 

	printf "  %s            \r", $id;	
	my $seq;
        eval{ $seq = $inx->fetch($id) };
        if ($@) {
                print "Problem with $id in the index, not included\n$@\n";
                next;
        }

        # Extract the alignment
        if($coord{$id}{strand} eq '+') {
            my $seq_to_print = $seq->trunc($coord{$id}{start},$coord{$id}{end});
            $out->write_seq($seq_to_print);
	    ++$n;	
        }
        elsif($coord{$id}{strand} eq '-') {
             my $seq_to_print = $seq->revcom->trunc($coord{$id}{start},$coord{$id}{end});
                $out->write_seq($seq_to_print);
		++$n;
        }
        else {
            print "Problem with $id\n";
        }
}
print "\n";
print "Exported $n sequences to $ARGV[2]\n";

