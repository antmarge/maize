#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
use List::Util qw(sum);
use Bio::SeqIO;
use Data::Random qw(:all);
use List::Util qw(sum);
use List::BinarySearch qw( :all );
use List::BinarySearch::XS;
use List::MoreUtils qw(uniq);
use File::Path;
use File::Basename;
use feature qw/say/;
use autodie;
use Data::Dumper qw(Dumper);


#ASSIGN INPUTS TO VARIABLES
our ($main,$ret,$out,$col,$key);
GetOptions(
'm:s' => \$main,
'r:s' => \$ret,
'o:s' =>\$out,
'c:s'=>\$col,
'k:s'=> \$key,
);


#usage: perl getCols.pl -m path/to/main -r path/to/retrival -k 1,0 -c 1,2


if (!$out) {$out="comb".$main}

sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}

my @keys=split(",",$key);
my $mainKey=$keys[0];
my $retKey=$keys[1];
my @cols=split(",",$col);

print "Main file: $main\n";
print "\tIndex for main key: $mainKey\n";
print "Retrieval file: $ret\n";
print "\tIndex for retrieval key: $retKey\n";
print "Output file: $out\n";

#GET DATA OUT OF MAIN FILE AND STORE IT INTO A HASH



sub mean {
	return sum(@_)/@_;
}

sub cleaner{
	my $line=$_[0];
	chomp($line);
	$line =~ s/\x0d{0,1}\x0a{0,1}\Z//s; 
	return $line;
	}

my %mainHash;
my %retHash;
my %outHash;

####Store main file in hash	
open (MAIN, '<', $main) or die "Could not open '$main' file \n";
	
my $line=<MAIN>;
my $mainline=cleaner($line); #gets rid of carriage returns (^M)
my @mainheader=split(",",$mainline);
#print $line,"\n",$mainline,"\n";

print "key column for main file: ", $mainheader[$mainKey],"\n";

while (my $entry = <MAIN>) {
	$entry=cleaner($entry);
	my @fields=split(',',$entry);
	unless(defined($fields[$mainKey])){
    	next;
		}
	$mainHash{$fields[$mainKey]}=\@fields;    
}

close MAIN;

####Store retrieval file in hash
open (RET, '<', $ret) or die "Could not open '$ret' file \n";
	
$line=<RET>;
$line=cleaner($line); #gets rid of carriage returns (^M)
my @retheader=split(",",$line);
print "key column for ret file: ", $retheader[$retKey],"\n";

@retheader=@retheader[@cols];

while (my $entry = <RET>) {
	$entry=cleaner($entry);
	my @retfields=split(',',$entry);
	my $key=$retfields[$retKey];
	if ( ! defined($key)){
    	next;
		}
	my @retcols=@retfields[@cols];
	$retHash{$retfields[$retKey]}=\@retcols;    
}

close RET;

#print Dumper \%retHash;
#die;
#Combine retrieval columns with main file columns

push (@mainheader,@retheader);
my $missCount=0;

foreach my $mk (keys %mainHash){
	my @maincols=@{$mainHash{$mk}};
	if (exists $retHash{$mk}){
		my @retcols=@{$retHash{$mk}};
		push (@maincols,@retcols);
		}
	else{
		my @missing=("NA","NA");
		push (@maincols,@missing);
		$missCount++;
		}
	$outHash{$mk}=\@maincols;
	}

open (OUT,">",$out);
print OUT join(",",@mainheader),"\n";
foreach my $k (keys %outHash){
	my @vals=@{$outHash{$k}};
	print OUT join(",",@vals);
	print OUT "\n";
}
close OUT;

