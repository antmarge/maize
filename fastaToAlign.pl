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
my $out = $opt_o;
my $ref=$opt_r || "B73";

if (!$out){
	my @temp=split(/\./,$in);
	print $in;
	foreach(@temp){
		print $_,"\n";
		}
	my $name=$temp[0];
	$out="${name}_align.fa";
	}

#Open INPUT (a fasta file. First seq ">" is reference)

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


#################################################
print "Looking for sequences to edit\n";

open OUT,">",$out;


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
					if ($seq[$i] eq ("/"||"-")){
						$i++;
					}
					else{
						push (@newseq,$seq[$i]);
						$i++;
						}
				}
				#$i++; #skip closing square bracket
				}
			else{
				push (@newseq,".");
			}	
		}
	my $stringSeq=join('',@newseq);
	$stringSeq =~ s/(.{1,80})/$1\n/gs;
	$hash{$head}=$stringSeq;
	}
	
	else{
		my $refSeq=$hash{$head};
		$refSeq=~ s/(.{1,80})/$1\n/gs;
		
		print OUT $head,"\n";
		print OUT $refSeq,"\n";
		#$hash{$head}=$refSeq;
		}

}

#print Dumper(\%hash);
	
open OUT,">",$out;

foreach my $k(keys %hash){
	print OUT $k,"\n";
	print OUT $hash{$k},"\n";
	}
close OUT;



#while (my $line=<IN>){
#	$line=chomp($line);
#	#Check if it is a header line beginning with >
	#if (index($line, ">") != -1){
#		#When new header, store last one in hash
#		$hash{$currentHead}={$currentSeq};
#		#Reassign
#		$currentHead=$line;
#		$currentSeq="";
#		}
#	else{
#		#Read sequence into string $currentSeq
#		$currentSeq=$currentSeq.$line;
#		}	
#	}
#close IN;
	
#print Dumper(\%hash);
