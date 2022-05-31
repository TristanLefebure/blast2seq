#! /usr/bin/perl

#9mai2022
# - clean up code

#23fev2016
#- allow to switch between query or hit extraction

#25nov2016
#-changge colonm selection
#-

use strict;
use warnings;
use Bio::Index::Fasta;
use Getopt::Long;
# use List::Util qw[min max];

#options
# my $min_seq_per_family = 2;
my $min_length = 100;
my $help;
my $extract = 'query';
# my $block_start;
# my $block_stop;

GetOptions(
# 	"start=s" => \$block_start;
# 	"end=s" => \$block_end;
	"length=i" => \$min_length,
	"help|?" => \$help,
	"extract=s" => \$extract
);

my $usage = "Usage: $0 <blast report file> <fasta index> <output>
with:
-blast report file: comes from a regular blast with -outfmt 6
-fasta index: produced by index_fasta.pl

Options:
  -help
  -length: minimum length for a hit to be included, default=$min_length.
   Carefull: if parsing a blastx alignment, this length is in aa units.
  -extract <query|target> : what do you want to extract ? [$extract]\n";

#2015-04-28: based n blast2fam_al,.pl

#declaration
# my %families; #AN->family
my %coord; #AN: start->dd, end->dd
my %seq_per_family; #family->nbr de seq
# my %fam_an; #family->array of AN

if ($#ARGV<1) {
	print $usage;
	exit;
}

###set the column numbers

my $pname; my $pstart; my $pstop;

if($extract eq 'query') {
	$pname = 0;
	$pstart = 6;
	$pstop = 7;
}
elsif($extract eq 'target') {
	$pname = 1;
	$pstart = 8;
	$pstop = 9;
} else {
	die "I don't understand the option extract : $extract\n";
}


#the index
my $inx = Bio::Index::Fasta->new($ARGV[1]) or die "Can't open the index file\n";
my $out = Bio::SeqIO->new(-file => ">$ARGV[2]", -format => 'fasta');

#the pairwise alignment
my $n = 0;
open BLAST, "$ARGV[0]" or die "Can't open $ARGV[1]\n";
while (<BLAST>) {
# 	++$n;
# 	if ($n == 1) { next }
	chomp;
	my @cols = split /\t/;
	my $length = $cols[3];
	if($length >= $min_length) {
		$coord{$cols[$pname]}{start} = $cols[$pstart]; #nov2016 change 1 -> 0
		$coord{$cols[$pname]}{end} = $cols[$pstop];
		if($cols[$pstart] > $cols[$pstop]) {
                    $coord{$cols[$pname]}{strand} = '-';
		}
		else {
                    $coord{$cols[$pname]}{strand} = '+';
		}
	}
	else { print "Seq $cols[$pname] as a hit too short: not included\n" }
}


#for each family
foreach my $id (keys %coord) {
	print "$id\n";
	my $seq;
        eval{ $seq = $inx->fetch($id) };
        if ($@) {
                print "Problem with $id in the index, not included\n$@\n";
                next;
        }
        #extract the alignment
        if($coord{$id}{strand} eq '+') {
            my $seq_to_print = $seq->trunc($coord{$id}{start},$coord{$id}{end});
            $out->write_seq($seq_to_print);
        }
        elsif($coord{$id}{strand} eq '-') {
             my $seq_to_print = $seq->trunc($coord{$id}{end},$coord{$id}{start});
                $out->write_seq($seq_to_print->revcom);
        }
        else {
            print "Problem with $id\n";
        }
}
