#!/usr/bin/perl -w

#Purpose: Create similarity matrix to compare SNPs for different samples

#TO DO:
#Allow gene list input---hash of hash

use strict;
use warnings;
use Getopt::Std;
use List::MoreUtils qw{ any };
use Data::Dumper qw(Dumper);



# Get global options

our($opt_i, $opt_c,$opt_o,$opt_s,$opt_g,$opt_h,$opt_m,$opt_l);
getopts('hi:c:o:s:g:ml:');
my $in = $opt_i;
my $ind = $opt_c;
my $out = $opt_o || "output.csv";
my $s = $opt_s;
my $g = $opt_g; #single gene, comma separated list, or file
my $h = $opt_h;
my $m= $opt_m; #if specified then include missing "./." values. Default: ignore.
my $l=$opt_l;
#perl ../maize/testdoublehash.pl -i tabmaster.txt -s 6,7,29,34,39 -o doubleHash1.csv 
my $outdir="output/";

#my @ind=split(",",$ind);
sub cleaner{
	my $line=$_[0];
	chomp($line);
	$line =~ s/\x0d{0,1}\x0a{0,1}\Z//s; 
	return $line;
	}

#Open gene file if exists

my @genes; #an array of all the genes to query for snp similarity
my @labels;
if (defined $g){
	if ($g=~/.txt/){ 
	# -g is a line delimited text file of genes
		open (GEN,"<",$g);
		while (<GEN>){
			my $gene=cleaner($_);
			my @gene=split("\t",$gene);
			push (@genes,$gene[1]);
			push (@labels,$gene[0]);
			}
		close GEN;
		}
	else{ 
	# -g is a comma separated list or single gene in the command line
		@genes=split(",",$g);
		@labels=split(",",$l);

		}
	}
	push (@genes,"all");
	push (@labels,"all");
	
	


#Open master file of SNPs and store them in the input hash (%input)

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
#this is going to be a double hash with gene name as key for 
#value that is key2:value2 pos:arrayOfInfo

while (<IN>){
	my $line=cleaner($_); #gets rid of carriage returns (^M)
	my @fields=split("\t",$line);
	my $geneName=$fields[1];
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
	$input{"all"}{$pos}=\@mod;
	$input{$geneName}{$pos}=\@mod;
	#The %input (hash) contains ALL column entries for a given position (key)
	}
	
close IN;

open SNP,">","snpNum.csv";

foreach my $gene (@genes){
	if (exists $input{$gene}){
	my %hash=%{$input{$gene}};
	my $countpergene=scalar (keys %hash);
	print SNP $gene,",",$countpergene,"\n";
	}
	
	}
close SNP;
die;
#print Dumper \%input;
#die;
my $allCount=0;
my $allSum=0;

sub simGene{
	#Accepts indices of columns to be compared and gene to group by
	my $c1=shift;
	my $c2=shift;
	my $geneKeep=shift;
	my $geneSum=0;
	my $geneCount=0;
	my %geneVals=%{$input{$geneKeep}};
	foreach (values %geneVals){
		my @vals=@{$_}; #de-reference the array stored in value of hash
		#Now @snps is an array of the SNPs being compared
		my $snp1=$vals[$c1];
		my $snp2=$vals[$c2];
		if ($m and (($snp1 eq "./.") || ($snp2 eq "./."))){
			next;
		}
		if ($snp1 eq $snp2){
			$geneSum++;
			}
		$geneCount++;
		}
	#my $score=$sum/$total;
	#print $c1, " and " ,$c2,"\t",$total,"\n";
	#my $cleanScore=sprintf("%.2f",($score*100));
	return ($geneSum,$geneCount);
}


#### Loop to make all possible comparisons between "samples"
#Just need column indices, not dealing with data

	my $gsi=0;

	foreach my $gene(@genes){
		my $label=$labels[$gsi];
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
					my $c=$sim[$ci];
					my ($sum,$count)=simGene($r,$c,$gene);

				push (@row,$result);
				}
			push (@out,\@row);
			$ri++;
			}

		#Now print output array
		my $outfile=join("",($label,$out));
		open (OUT,">",$outfile) or print "Cannot open $outfile \n" and exit;
		foreach my $row(@out){
			my @rowArray=@{$row};
			my $string=join(",",@rowArray);
			print OUT $string,"\n";
			}
		close OUT;
		$gsi++;
		}
