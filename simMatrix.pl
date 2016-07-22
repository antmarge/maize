#!/usr/bin/perl -w

#Purpose: Create similarity matrix to compare SNPs for different samples

#TO DO:
#Add recursive function to query all possible comparisons
#Allow gene list input

use strict;
use warnings;
use Getopt::Std;
use List::MoreUtils qw{ any };
use Data::Dumper qw(Dumper);



# Get global options

our($opt_i, $opt_c,$opt_o,$opt_s,$opt_g,$opt_h,$opt_m);
getopts('hi:c:o:s:g:m');
my $in = $opt_i;
my $ind = $opt_c;
my $out = $opt_o || "output.csv";
my $s = $opt_s;
my $g = $opt_g;
my $h = $opt_h;
my $m= $opt_m; #if specified then include missing "./." values. Default: ignore.

#usage: perl simMatrix.pl -i master.csv -c 6,7,29,34,39,44,49,54,59,64,69,74 -s 0,1
#usage: perl simMatrix.pl -i master.csv -s 6,7

#my @ind=split(",",$ind);
sub cleaner{
	my $line=$_[0];
	chomp($line);
	$line =~ s/\x0d{0,1}\x0a{0,1}\Z//s; 
	return $line;
	}
	
open (IN, "<", $in) or print "Cannot open file $in !" and exit;

my $headerLine=<IN>;
$headerLine=cleaner($headerLine); #gets rid of carriage returns (^M)
my @head=split("\t",$headerLine);

my $i=0;
if ($h){
	foreach (@head){
		print $i,"-",$_,"\n";
		$i++;
		}
	exit;
	}
	
my @sim=split(",",$s);
my %input;

while (<IN>){
	my $line=cleaner($_); #gets rid of carriage returns (^M)
	my @fields=split("\t",$line);
	my $geneName=$fields[1];
	if ($g and ($g ne $geneName)){
		next;
		}
	my $pos=$fields[3];
	#If the col index is a comparison col then check for N/N format first
	#Then store it into the @mod array.
	my @mod;
	for (my $i=0; $i<scalar @fields;$i++){
		my $f=$fields[$i];
		# if the field entry is a snp then check and modify if needed
		if (any { $_ eq $i} @sim){	
			if ($f!~/\//){
				$f=$f."/".$f;
				}
			}
		push (@mod,$f);
		}
		if (($m and any { $_ eq "./."} @mod)){
			next;
		}		
	$input{$pos}=\@mod;
	#The %input (hash) contains ALL column entries for a given position (key)
	}
	
close IN;

#print Dumper \%input;
#die;


#### Subroutine for comparing two columns (all genes, or indiv genes)
sub sim{
	#Accepts indices of columns to be compared 
	my $c1=shift;
	my $c2=shift;
	
	my $sum=0;
	my $total=0;
	foreach (values %input){
		my @vals=@{$_}; #de-reference the array stored in value of hash
		#Now @snps is an array of the SNPs being compared
		if ($vals[$c1] eq $vals[$c2]){
			$sum++;
			}
		$total++;
		}
	my $score=$sum/$total;
	my $cleanScore=sprintf("%.2f",($score*100));
	return $cleanScore;
}


if (!$g){
	$g="all genes";
	}

#### Loop to make all possible comparisons between "samples"
#Just need column indices, not dealing with data

#my @compNames=@cols[@sim];
my $n=scalar @sim;
my $ri=0;
my @out; #store double array (matrix) of results here
my @colnames=@head[@sim];
unshift (@colnames," "); #need empty space for corner of table

push (@out,\@colnames); #column names for output matrix

foreach my $r(@sim){
	my @row;
	my $rname=$head[$r];
	push (@row,$rname);
	my $result;
	for (my $ci=0;$ci<$n;$ci++){
		#my $result="X";
		#if ($ci>=$ri){
		#foreach my $c(@sim){
			my $c=$sim[$ci];
			$result=sim($r,$c);
			#}
		push (@row,$result);
		}
	push (@out,\@row);
	$ri++;
		
	}

#Now print output array

open OUT,">",$out;
foreach my $row(@out){
	my @rowArray=@{$row};
	my $string=join(",",@rowArray);
	print OUT $string,"\n";
	}
close OUT;
