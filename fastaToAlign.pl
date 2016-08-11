#!/usr/bin/perl -w

#Purpose: Format fasta gene sequence file for DNAsp as alignment file with
#symbols

use strict;
use warnings;
use Getopt::Std;
use List::MoreUtils qw{ any };
use Data::Dumper qw(Dumper);
use Bio::SeqIO;


# Get global options

our($opt_i,$opt_o,$opt_r);
getopts('i:o:r:');
my $in = $opt_i;
my $out = $opt_o || "output.csv";
my $ref=$opt_r || "B73";

my $outdir="output/";

#Open INPUT (a fasta file. First seq ">" is reference)

my $seqio = Bio::SeqIO->new(-file => $in, '-format' => 'Fasta');

my %hash;

open(FASTA, $in);


my $currentHead=<FASTA>;
chomp($currentHead);
#print $currentHead;
my $currentSeq="";

while(<FASTA>) {
    chomp($_);
    if ($_ =~  m/^>/ ) {
		#When new header, store last one in hash
		$hash{$currentHead}=$currentSeq;
		#Reassign
		$currentHead=$_;
		$currentSeq="";
		}
	else{
		#Read sequence into string $currentSeq
		$currentSeq=$currentSeq.$_;
		}	
    }

close FASTA;
my $change=".";
print "Looking for sequences to edit\n";

foreach my $head (keys %hash){
	print $head,"\n";
	if (index($head, $ref) == -1) {
		my $s=$hash{$head};
		#empty hash
		delete $hash{$head};
		my @seq=split("",$s);
		#place to store edited sequence
		my @newseq=();
		###########Now check each nucleotide
		
		for (my $i=0;$i<scalar @seq;$i++){
		print $seq[$i];
			#There is a deleted seq and we need to put "-"
			if ($seq[$i] eq "{"){
				$i++; #skip the opening curly bracket
				while ($seq[$i] ne "}"){
					if ($seq[$i] eq ("/"||"-")){
						$i++;
					}
					else{
						push (@newseq,"-");
						$i++;
					}
				}
				#$i++; #skip curly closing bracket
			}

			#There is a SNP and we need to keep the nucleotide letters
			elsif ($seq[$i] eq "["){
				$i++; #skip opening square bracket
				while($seq[$i] ne "]"){
					push (@newseq,$seq[$i]);
				}
				#$i++; #skip closing square bracket
				}
			else{
				push (@newseq,".");
			}	
		}
	my $stringSeq=join('',@newseq);
	$hash{$head}=$stringSeq;
	print Dumper(\%hash);
	die;		
	}

}

print Dumper(\%hash);
	


die;


while (my $line=<IN>){
	$line=chomp($line);
	#Check if it is a header line beginning with >
	if (index($line, ">") != -1){
		#When new header, store last one in hash
		$hash{$currentHead}={$currentSeq};
		#Reassign
		$currentHead=$line;
		$currentSeq="";
		}
	else{
		#Read sequence into string $currentSeq
		$currentSeq=$currentSeq.$line;
		}	
	}
close IN;
	
print Dumper(\%hash);
	

	