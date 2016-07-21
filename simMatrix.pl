#!/usr/bin/perl -w

#Purpose: Create similarity matrix to compare SNPs for different samples

#TO DO:
#Add recursive function to query all possible comparisons
#Allow gene list input

use strict;
use Getopt::Std;
use Data::Dumper qw(Dumper);

# Get global options

our($opt_i, $opt_c,$opt_o,$opt_s,$opt_g,$opt_h);
getopts('hi:c:o:s:g:');
my $in = $opt_i;
my $ind = $opt_c;
my $out = $opt_o || "output.csv";
my $s = $opt_s;
my $g = $opt_g;
my $h = $opt_h;

#usage: perl simMatrix.pl -i master.csv -c 6,7,29,34,39,44,49,54,59,64,69,74 -s 0,1
#usage: perl simMatrix.pl -i master.csv -s 6,7

#my @ind=split(",",$ind);

open (IN, "<", $in);

my $headerLine=<IN>;
my @cols=split("\t",$headerLine);

my $i=0;
if ($h){
	foreach (@cols){
		print $i,"-",$_,"\n";
		$i++;
		}
	exit;
	}
	
my @sim=split(",",$s);
my $c1=$sim[0];
my $c2=$sim[1];

my @compNames=@cols[@sim];
my $compNum=scalar @compNames;

my %input;

while (<IN>){
	my @fields=split("\t",$_);
	my $geneName=$fields[1];
	if ($g and ($g ne $geneName)){
		next;
		}
	my $pos=$fields[3];
	my @temp=@fields[@sim];
	my @imp; #after checking for N/N format, store here
	foreach my $snp(@temp){
		if ($snp!~/\//){
			$snp=$snp."/".$snp;
			}
		push (@imp,$snp);
		}		
	$input{$pos}=\@imp;
	}
close IN;

#print Dumper \%input;
#die;

if (!$g){
	$g="all genes";
	}


my $sum=0;
my $total=0;
foreach (values %input){
	my @snps=@{$_}; #de-reference the array stored in value of hash
	#Now @snps is an array of the SNPs being compared
	if ($snps[0] eq $snps[1]){
		$sum++;
		}
	$total++;
	}
my $score=$sum/$total;
my $cleanScore=sprintf("%.2f",($score*100));
print "\nSNP similarity between ",$compNames[0], " and ",$compNames[1],
"\nfor $g: ",$cleanScore,"%\n";

print "\n";
	
#print Dumper \%input